import hashlib
import json
import logging
import os
import re
import threading
import time
import uuid
from pathlib import Path

from flask import Flask, g, has_request_context, jsonify, render_template, request, send_file
from werkzeug.exceptions import RequestEntityTooLarge
from werkzeug.utils import secure_filename

from cleanup_uploads import cleanup_uploads
from utils import (
    RamachandranManager,
    fetch_alphafold_model,
    fetch_structure_file,
    generate_csv,
    generate_pdf_report,
    get_phi_psi,
    is_valid_uniprot_accession,
    parse_structure,
)

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = "uploads"
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16MB limit
app.config["PDB_FETCH_TIMEOUT_SECONDS"] = 15
app.config["ALPHAFOLD_FETCH_TIMEOUT_SECONDS"] = 15

os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)

# Keep compatibility with classic 4-char PDB IDs and upcoming extended IDs.
LEGACY_PDB_ID_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")
EXTENDED_PDB_ID_RE = re.compile(r"^pdb_[A-Za-z0-9]{8}$", flags=re.IGNORECASE)
ALLOWED_UPLOAD_SUFFIXES = {".pdb", ".cif", ".mmcif"}
CLEANUP_INTERVAL_SECONDS = 3600
CLEANUP_MAX_AGE_HOURS = 24
REFERENCE_CACHE_MAX_AGE_SECONDS = 86400

_last_cleanup_ts = 0.0
_cleanup_lock = threading.Lock()

# Initialize the Ramachandran data manager
# Assuming data is in 'data' directory relative to the project root
rama_manager = RamachandranManager(data_directory="data")

REFERENCE_DATA = {
    key: {
        "phi": value["grid"]["phi"],
        "psi": value["grid"]["psi"],
        "z": value["grid"]["z"],
        "levels": value["levels"],
    }
    for key, value in rama_manager.rama_data.items()
}
REFERENCE_JSON = json.dumps(REFERENCE_DATA, separators=(",", ":"), sort_keys=True)
REFERENCE_ETAG = hashlib.sha256(REFERENCE_JSON.encode("utf-8")).hexdigest()


class RequestFormatter(logging.Formatter):
    def format(self, record):
        if has_request_context():
            record.request_id = getattr(g, "request_id", "-")
        else:
            record.request_id = "-"
        return super().format(record)


def configure_logging():
    formatter = RequestFormatter("%(asctime)s %(levelname)s request_id=%(request_id)s %(name)s: %(message)s")
    root_logger = logging.getLogger()

    if root_logger.handlers:
        for handler in root_logger.handlers:
            handler.setFormatter(formatter)
    else:
        root_handler = logging.StreamHandler()
        root_handler.setFormatter(formatter)
        root_logger.addHandler(root_handler)

    root_logger.setLevel(logging.INFO)

    if not app.logger.handlers:
        app_handler = logging.StreamHandler()
        app_handler.setFormatter(formatter)
        app.logger.addHandler(app_handler)

    for handler in app.logger.handlers:
        handler.setFormatter(formatter)
    app.logger.setLevel(logging.INFO)
    app.logger.propagate = False


def maybe_cleanup_upload_folder():
    global _last_cleanup_ts

    now = time.time()
    # Fast-path check avoids lock contention when cleanup is not due.
    if now - _last_cleanup_ts < CLEANUP_INTERVAL_SECONDS:
        return

    if not _cleanup_lock.acquire(blocking=False):
        return

    try:
        now = time.time()
        # Re-check under lock: another worker may have run cleanup already.
        if now - _last_cleanup_ts < CLEANUP_INTERVAL_SECONDS:
            return

        app.logger.info("Running periodic upload cleanup.")
        cleanup_uploads(
            uploads_dir=app.config["UPLOAD_FOLDER"],
            max_age_hours=CLEANUP_MAX_AGE_HOURS,
            logger=app.logger,
        )
        _last_cleanup_ts = now
    finally:
        _cleanup_lock.release()


configure_logging()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/healthz")
def healthz():
    """Health check endpoint for monitoring and orchestration tools."""
    return jsonify({"status": "healthy"}), 200


@app.route("/info")
def info():
    """Returns application metadata."""
    return jsonify(
        {
            "name": "ramachandran",
            "version": "0.1.0",
            "description": "A web-based Ramachandran plot analysis and visualization tool for protein structures.",
        }
    ), 200


@app.route("/reference")
def reference():
    response = app.response_class(REFERENCE_JSON, mimetype="application/json")
    response.set_etag(REFERENCE_ETAG)
    response.cache_control.public = True
    response.cache_control.max_age = REFERENCE_CACHE_MAX_AGE_SECONDS
    response.make_conditional(request)
    return response


@app.route("/robots.txt")
def robots_txt():
    """Standard instructions for web crawlers."""
    return (
        "User-agent: *\nDisallow: /process\nDisallow: /download/",
        200,
        {"Content-Type": "text/plain"},
    )


@app.route("/favicon.ico")
def favicon():
    """Serve the favicon to prevent 404 logs in browsers."""
    return app.send_static_file("favicon.ico")


@app.errorhandler(RequestEntityTooLarge)
def handle_file_too_large(_error):
    return jsonify({"error": "Uploaded file is too large. Maximum allowed size is 16 MB."}), 413


def is_valid_pdb_id(value):
    if not value:
        return False
    return bool(LEGACY_PDB_ID_RE.fullmatch(value) or EXTENDED_PDB_ID_RE.fullmatch(value))


def is_allowed_upload(filename):
    suffix = Path(filename).suffix.lower()
    return suffix in ALLOWED_UPLOAD_SUFFIXES


@app.before_request
def attach_request_id():
    external_request_id = (request.headers.get("X-Request-ID") or "").strip()
    g.request_id = external_request_id[:64] if external_request_id else uuid.uuid4().hex


@app.after_request
def add_request_id_header(response):
    if getattr(g, "request_id", None):
        response.headers["X-Request-ID"] = g.request_id
    return response


@app.route("/process", methods=["POST"])
def process():
    maybe_cleanup_upload_folder()

    structure_id = (request.form.get("pdb_id") or "").strip()
    file = request.files.get("pdb_file")

    filepath = None

    if structure_id:
        if is_valid_pdb_id(structure_id):
            filepath, fetch_error = fetch_structure_file(
                structure_id,
                output_dir=app.config["UPLOAD_FOLDER"],
                timeout=app.config["PDB_FETCH_TIMEOUT_SECONDS"],
            )
        elif is_valid_uniprot_accession(structure_id):
            filepath, fetch_error = fetch_alphafold_model(
                structure_id,
                output_dir=app.config["UPLOAD_FOLDER"],
                timeout=app.config["ALPHAFOLD_FETCH_TIMEOUT_SECONDS"],
            )
        else:
            error_message = (
                "Invalid structure identifier format. Use a PDB ID (e.g., '1UBQ' or 'pdb_00001abc') "
                "or a UniProt accession (e.g., 'P69905')."
            )
            return jsonify({"error": error_message}), 400

        if not filepath:
            return jsonify({"error": fetch_error or "Failed to fetch structure from remote source."}), 400
    elif file and file.filename != "":
        filename = secure_filename(file.filename)
        if not filename:
            return jsonify({"error": "Invalid upload filename."}), 400
        if not is_allowed_upload(filename):
            allowed = ", ".join(sorted(ALLOWED_UPLOAD_SUFFIXES))
            return jsonify({"error": f"Unsupported file format. Allowed: {allowed}."}), 400
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], f"{uuid.uuid4()}_{filename}")
        file.save(filepath)
    else:
        return jsonify({"error": "No PDB ID or file provided"}), 400

    try:
        structure = parse_structure(filepath)
        phi_psi_data = get_phi_psi(structure)
        if not phi_psi_data:
            raise ValueError("No protein residues were found in the provided structure.")
        if not any(item["phi"] is not None and item["psi"] is not None for item in phi_psi_data):
            raise ValueError("No valid phi/psi angles could be computed for the provided structure.")

        # Classify each residue
        for item in phi_psi_data:
            if item["phi"] is not None and item["psi"] is not None:
                score, category = rama_manager.classify_phipsi(item["rama_type"], item["phi"], item["psi"])
                item["score"] = score
                item["classification"] = category
            else:
                item["score"] = None
                item["classification"] = None

        # Save results to a JSON file for on-demand downloads
        result_id = str(uuid.uuid4())
        results_path = os.path.join(app.config["UPLOAD_FOLDER"], f"{result_id}_results.json")
        with open(results_path, "w") as f:
            json.dump(
                {
                    "pdb_id": structure_id if structure_id else file.filename,
                    "phi_psi": phi_psi_data,
                },
                f,
            )

        response = {
            "result_id": result_id,
            "pdb_id": structure_id if structure_id else file.filename,
            "phi_psi": phi_psi_data,
        }

        return jsonify(response)

    except ValueError as err:
        app.logger.warning(
            "Rejected structure input: %s",
            err,
            extra={
                "pdb_id": structure_id or None,
                "uploaded_filename": file.filename if file else None,
                "filepath": filepath,
            },
        )
        return jsonify({"error": str(err)}), 400
    except Exception:
        app.logger.exception(
            "Failed to process structure",
            extra={
                "pdb_id": structure_id or None,
                "uploaded_filename": file.filename if file else None,
                "filepath": filepath,
            },
        )
        return jsonify({"error": "Internal error while processing structure. Please try again."}), 500
    finally:
        # Optionally clean up files, but maybe keep them for a bit
        pass


@app.route("/download/csv/<result_id>")
def download_csv(result_id):
    results_path = os.path.join(app.config["UPLOAD_FOLDER"], f"{result_id}_results.json")
    if not os.path.exists(results_path):
        return "Result not found", 404

    with open(results_path) as f:
        data = json.load(f)

    csv_content = generate_csv(data["phi_psi"])
    pdb_id = data["pdb_id"].lower()

    from io import BytesIO

    buf = BytesIO(csv_content.encode("utf-8"))

    return send_file(
        buf,
        mimetype="text/csv",
        as_attachment=True,
        download_name=f"ramachandran_{pdb_id}.csv",
    )


@app.route("/download/pdf/<result_id>")
def download_pdf(result_id):
    results_path = os.path.join(app.config["UPLOAD_FOLDER"], f"{result_id}_results.json")
    if not os.path.exists(results_path):
        return "Result not found", 404

    with open(results_path) as f:
        data = json.load(f)

    pdf_buf = generate_pdf_report(data["phi_psi"], data["pdb_id"], rama_manager)
    pdb_id = data["pdb_id"].lower()

    return send_file(
        pdf_buf,
        mimetype="application/pdf",
        as_attachment=True,
        download_name=f"ramachandran_{pdb_id}.pdf",
    )


if __name__ == "__main__":
    app.run(debug=True, port=5001)
