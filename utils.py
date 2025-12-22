import os
import gzip
import urllib.request
import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from scipy import interpolate

# Define the Ramachandran plot types and levels
RAMA_KEYS = ["General", "Ile/Val", "pre-Pro", "Gly", "trans-Pro", "cis-Pro"]
RAMA_TYPES_FILE_MAP = {
    "General": "general-noGPIVpreP",
    "Ile/Val": "ileval-nopreP",
    "pre-Pro": "prepro-noGP",
    "Gly": "gly-sym",
    "trans-Pro": "transpro",
    "cis-Pro": "cispro",
}

LEVELS = {
    "General": [0.0005, 0.02],
    "Ile/Val": [0.001, 0.02],
    "pre-Pro": [0.001, 0.02],
    "Gly": [0.001, 0.02],
    "trans-Pro": [0.001, 0.02],
    "cis-Pro": [0.002, 0.02],
}

class RamachandranManager:
    def __init__(self, data_directory="data"):
        self.data_directory = data_directory
        self.rama_data = self._load_reference_data()

    def _load_reference_data(self):
        rama_data = {}
        dihed_range = np.arange(-181.0, 182.0, 2.0)
        
        for key in RAMA_KEYS:
            filename = os.path.join(self.data_directory, f"rama8000-{RAMA_TYPES_FILE_MAP[key]}.data")
            if not os.path.exists(filename):
                print(f"Warning: Reference data file {filename} not found.")
                continue

            data = []
            with open(filename, "r") as f:
                for line in f:
                    if line and not line.startswith("#"):
                        data.append([float(val) for val in line.split()])
            
            data = np.array(data)
            phi_vals, psi_vals, freq = data.T
            col_idx = ((phi_vals + 180) // 2).astype(int)
            row_idx = ((psi_vals + 180) // 2).astype(int)

            # Build 2D distribution
            rama_dist = np.zeros((180, 180), dtype=float)
            rama_dist[row_idx, col_idx] = freq
            rama_dist = np.pad(rama_dist, 1, mode="wrap")
            
            # create interpolating function
            f_interp = interpolate.RegularGridInterpolator(
                (dihed_range, dihed_range), rama_dist.T
            )

            rama_data[key] = {
                "dist": rama_dist.tolist(), # For JSON serialization if needed
                "levels": LEVELS[key],
                "f": f_interp,
                "grid": {
                    "phi": dihed_range.tolist(),
                    "psi": dihed_range.tolist(),
                    "z": rama_dist.tolist()
                }
            }
        return rama_data

    def classify_phipsi(self, rama_type, phi, psi):
        if rama_type not in self.rama_data:
            return 0.0, "unknown"
            
        f = self.rama_data[rama_type]["f"]
        levels = self.rama_data[rama_type]["levels"]
        
        try:
            value = f((phi, psi))
        except:
            return 0.0, "outlier"

        if value >= levels[1]:
            category = "favoured"
        elif value >= levels[0]:
            category = "allowed"
        else:
            category = "outlier"
            
        return float(value * 100.0), category

def fetch_pdb_file(pdb_id, output_dir="temp_pdb"):
    os.makedirs(output_dir, exist_ok=True)
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz"
    filepath = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    try:
        with urllib.request.urlopen(url) as response:
            if response.status == 200:
                content = gzip.decompress(response.read()).decode("utf-8")
                with open(filepath, "w") as f:
                    f.write(content)
                return filepath
    except Exception as e:
        print(f"Error fetching PDB {pdb_id}: {e}")
    return None

def calc_dihedral(a, b, c, d):
    # Vector implementation similar to the old script
    b1 = b - a
    b2 = c - b
    b3 = d - c
    
    def normalize(v):
        norm = np.linalg.norm(v, axis=-1, keepdims=True)
        return v / norm

    n1 = normalize(np.cross(b1, b2))
    n2 = normalize(np.cross(b2, b3))
    m1 = normalize(np.cross(b2, n1))
    
    x = np.sum(n1 * n2, axis=-1)
    y = np.sum(m1 * n2, axis=-1)
    return np.rad2deg(np.arctan2(y, x))

def get_phi_psi(structure):
    results = []
    
    # Mapping for Ramachandran types
    GENERAL_AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "LYS", "LEU", "MET", "ASN", "GLN", "ARG", "SER", "THR", "TRP", "TYR"]
    
    for model in structure:
        for chain in model:
            residues = [res for res in chain if is_aa(res)]
            
            # Check for chain breaks by CA-CA distance
            for i in range(len(residues)):
                res = residues[i]
                res_name = res.get_resname()
                
                # Default type
                rama_type = "General" if res_name in GENERAL_AA else "unknown"
                if res_name in ["ILE", "VAL"]:
                    rama_type = "Ile/Val"
                elif res_name == "GLY":
                    rama_type = "Gly"
                elif res_name == "PRO":
                    rama_type = "trans-Pro" # Will check for cis later
                
                # Check for pre-Pro
                if i < len(residues) - 1:
                    next_res = residues[i+1]
                    if next_res.get_resname() == "PRO" and res_name not in ["PRO", "GLY"]:
                        rama_type = "pre-Pro"
                
                # Get backbone atoms
                try:
                    # Previous residue C
                    if i > 0:
                        prev_res = residues[i-1]
                        # Check distance for chain break (CA-CA > 4.0)
                        if "CA" in res and "CA" in prev_res:
                            dist = np.linalg.norm(res["CA"].get_coord() - prev_res["CA"].get_coord())
                            if dist > 4.5: # Chain break
                                prev_res_c = None
                            else:
                                prev_res_c = prev_res["C"] if "C" in prev_res else None
                        else:
                            prev_res_c = None
                    else:
                        prev_res_c = None
                        
                    # Current residue atoms
                    n = res["N"]
                    ca = res["CA"]
                    c = res["C"]
                    
                    # Next residue N
                    if i < len(residues) - 1:
                        next_res = residues[i+1]
                        # Check distance for chain break
                        if "CA" in res and "CA" in next_res:
                            dist = np.linalg.norm(res["CA"].get_coord() - next_res["CA"].get_coord())
                            if dist > 4.5:
                                next_res_n = None
                            else:
                                next_res_n = next_res["N"] if "N" in next_res else None
                        else:
                            next_res_n = None
                    else:
                        next_res_n = None
                        
                    phi = None
                    if prev_res_c is not None:
                        phi = calc_dihedral(prev_res_c.get_coord(), n.get_coord(), ca.get_coord(), c.get_coord())
                        
                    psi = None
                    if next_res_n is not None:
                        psi = calc_dihedral(n.get_coord(), ca.get_coord(), c.get_coord(), next_res_n.get_coord())
                        
                    omega = None
                    if prev_res_c is not None:
                        prev_ca = prev_res["CA"]
                        omega = calc_dihedral(prev_ca.get_coord(), prev_res_c.get_coord(), n.get_coord(), ca.get_coord())
                        if res_name == "PRO" and omega is not None and abs(omega) < 90.0:
                            rama_type = "cis-Pro"

                    if phi is not None and psi is not None:
                        results.append({
                            "chain": chain.id,
                            "resSeq": res.id[1],
                            "resName": res_name,
                            "phi": float(phi),
                            "psi": float(psi),
                            "omega": float(omega) if omega is not None else None,
                            "rama_type": rama_type
                        })
                except KeyError:
                    # Missing backbone atoms
                    continue
                    
    return results

def parse_structure(filepath):
    if filepath.endswith(".gz"):
        # Handle compressed files if necessary, but fetch_pdb_file already decompresses
        pass
        
    if filepath.endswith(".cif") or filepath.endswith(".mmcif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)
        
    structure = parser.get_structure("protein", filepath)
    return structure
