import os
import uuid
import json
from flask import Flask, render_template, request, jsonify, send_file
from werkzeug.utils import secure_filename
from utils import (
    RamachandranManager,
    parse_structure,
    fetch_structure_file,
    get_phi_psi,
    generate_csv,
    generate_pdf_report,
)

app = Flask(__name__)
app.secret_key = "supersecretkey"
app.config["UPLOAD_FOLDER"] = "uploads"
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # 16MB limit

os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)

# Initialize the Ramachandran data manager
# Assuming data is in 'data' directory relative to the project root
rama_manager = RamachandranManager(data_directory="data")


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
    return jsonify({
        "name": "ramachandran",
        "version": "0.1.0",
        "description": "A web-based Ramachandran plot analysis and visualization tool for protein structures."
    }), 200


@app.route("/robots.txt")
def robots_txt():
    """Standard instructions for web crawlers."""
    return "User-agent: *\nDisallow: /process\nDisallow: /download/", 200, {"Content-Type": "text/plain"}


@app.route("/favicon.ico")
def favicon():
    """Serve the favicon to prevent 404 logs in browsers."""
    return app.send_static_file("favicon.ico")


@app.route("/process", methods=["POST"])
def process():
    pdb_id = request.form.get("pdb_id")
    file = request.files.get("pdb_file")

    filepath = None

    if pdb_id:
        filepath = fetch_structure_file(pdb_id, output_dir=app.config["UPLOAD_FOLDER"])
        if not filepath:
            return jsonify({"error": "Failed to fetch PDB ID from RCSB"}), 400
    elif file and file.filename != "":
        filename = secure_filename(file.filename)
        filepath = os.path.join(
            app.config["UPLOAD_FOLDER"], f"{uuid.uuid4()}_{filename}"
        )
        file.save(filepath)
    else:
        return jsonify({"error": "No PDB ID or file provided"}), 400

    try:
        structure = parse_structure(filepath)
        phi_psi_data = get_phi_psi(structure)

        # Classify each residue
        for item in phi_psi_data:
            if item["phi"] is not None and item["psi"] is not None:
                score, category = rama_manager.classify_phipsi(
                    item["rama_type"], item["phi"], item["psi"]
                )
                item["score"] = score
                item["classification"] = category
            else:
                item["score"] = None
                item["classification"] = None

        # Prepare reference data for the frontend (contours)
        # We only send the grid for the types present in the structure, or all if preferred
        # For simplicity, let's send reference data for the 6 standard types
        reference_data = {}
        for key in rama_manager.rama_data:
            reference_data[key] = rama_manager.rama_data[key]["grid"]
            reference_data[key]["levels"] = rama_manager.rama_data[key]["levels"]

        # Save results to a JSON file for on-demand downloads
        result_id = str(uuid.uuid4())
        results_path = os.path.join(
            app.config["UPLOAD_FOLDER"], f"{result_id}_results.json"
        )
        with open(results_path, "w") as f:
            json.dump(
                {
                    "pdb_id": pdb_id if pdb_id else file.filename,
                    "phi_psi": phi_psi_data,
                },
                f,
            )

        response = {
            "result_id": result_id,
            "pdb_id": pdb_id if pdb_id else file.filename,
            "phi_psi": phi_psi_data,
            "reference": reference_data,
        }

        return jsonify(response)

    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        # Optionally clean up files, but maybe keep them for a bit
        pass


@app.route("/download/csv/<result_id>")
def download_csv(result_id):
    results_path = os.path.join(
        app.config["UPLOAD_FOLDER"], f"{result_id}_results.json"
    )
    if not os.path.exists(results_path):
        return "Result not found", 404

    with open(results_path, "r") as f:
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
    results_path = os.path.join(
        app.config["UPLOAD_FOLDER"], f"{result_id}_results.json"
    )
    if not os.path.exists(results_path):
        return "Result not found", 404

    with open(results_path, "r") as f:
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
