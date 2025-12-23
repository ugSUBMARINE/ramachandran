import os
import uuid
from flask import Flask, render_template, request, jsonify, flash, redirect, url_for
from werkzeug.utils import secure_filename
from utils import RamachandranManager, parse_structure, fetch_pdb_file, get_phi_psi

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

@app.route("/process", methods=["POST"])
def process():
    pdb_id = request.form.get("pdb_id")
    file = request.files.get("pdb_file")
    
    filepath = None
    
    if pdb_id:
        filepath = fetch_pdb_file(pdb_id, output_dir=app.config["UPLOAD_FOLDER"])
        if not filepath:
            return jsonify({"error": "Failed to fetch PDB ID from RCSB"}), 400
    elif file and file.filename != "":
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config["UPLOAD_FOLDER"], f"{uuid.uuid4()}_{filename}")
        file.save(filepath)
    else:
        return jsonify({"error": "No PDB ID or file provided"}), 400
    
    try:
        structure = parse_structure(filepath)
        phi_psi_data = get_phi_psi(structure)
        
        # Classify each residue
        for item in phi_psi_data:
            if item["phi"] is not None and item["psi"] is not None:
                score, category = rama_manager.classify_phipsi(item["rama_type"], item["phi"], item["psi"])
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

        from utils import generate_csv
        csv_data = generate_csv(phi_psi_data)

        response = {
            "pdb_id": pdb_id if pdb_id else file.filename,
            "phi_psi": phi_psi_data,
            "reference": reference_data,
            "csv_data": csv_data
        }
        
        return jsonify(response)
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        # Optionally clean up files, but maybe keep them for a bit
        pass

if __name__ == "__main__":
    app.run(debug=True, port=5001)
