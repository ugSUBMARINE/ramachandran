import csv
import gzip
import io
import os
import urllib.request

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from Bio.PDB.Polypeptide import is_aa
from scipy import interpolate

matplotlib.use("Agg")  # Use non-interactive backend

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
            filename = os.path.join(
                self.data_directory, f"rama8000-{RAMA_TYPES_FILE_MAP[key]}.data"
            )
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
                "dist": rama_dist.tolist(),  # For JSON serialization if needed
                "levels": LEVELS[key],
                "f": f_interp,
                "grid": {
                    "phi": dihed_range.tolist(),
                    "psi": dihed_range.tolist(),
                    "z": rama_dist.tolist(),
                },
            }
        return rama_data

    def classify_phipsi(self, rama_type, phi, psi):
        if rama_type not in self.rama_data:
            return 0.0, "unknown"

        f = self.rama_data[rama_type]["f"]
        levels = self.rama_data[rama_type]["levels"]

        try:
            value = f((phi, psi))
        except Exception:
            return 0.0, "outlier"

        if value >= levels[1]:
            category = "favoured"
        elif value >= levels[0]:
            category = "allowed"
        else:
            category = "outlier"

        return float(value * 100.0), category


def fetch_structure_file(pdb_id, output_dir="temp_pdb"):
    os.makedirs(output_dir, exist_ok=True)
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz"
    filepath = os.path.join(output_dir, f"{pdb_id}.cif")

    try:
        with urllib.request.urlopen(url) as response:
            if response.status == 200:
                content = gzip.decompress(response.read()).decode("utf-8")
                with open(filepath, "w") as f:
                    f.write(content)
                return filepath
    except Exception as e:
        print(f"Error fetching structure {pdb_id}: {e}")
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
    GENERAL_AA = [
        "ALA",
        "CYS",
        "ASP",
        "GLU",
        "PHE",
        "HIS",
        "LYS",
        "LEU",
        "MET",
        "MSE",
        "ASN",
        "GLN",
        "ARG",
        "SER",
        "THR",
        "TRP",
        "TYR",
    ]

    for model_nr, model in enumerate(structure):
        # Only process first model
        if model_nr > 0:
            break

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
                    rama_type = "trans-Pro"  # Will check for cis later

                # Check for pre-Pro
                if i < len(residues) - 1:
                    next_res = residues[i + 1]
                    if next_res.get_resname() == "PRO" and res_name not in [
                        "PRO",
                        "GLY",
                    ]:
                        rama_type = "pre-Pro"

                # Get backbone atoms
                try:
                    prev_res_c = None
                    if i > 0:
                        prev_res = residues[i - 1]
                        if "CA" in res and "CA" in prev_res:
                            dist = np.linalg.norm(
                                res["CA"].get_coord() - prev_res["CA"].get_coord()
                            )
                            if dist <= 4.5:
                                prev_res_c = prev_res["C"] if "C" in prev_res else None

                    n = res["N"] if "N" in res else None
                    ca = res["CA"] if "CA" in res else None
                    c = res["C"] if "C" in res else None

                    next_res_n = None
                    if i < len(residues) - 1:
                        next_res = residues[i + 1]
                        if "CA" in res and "CA" in next_res:
                            dist = np.linalg.norm(
                                res["CA"].get_coord() - next_res["CA"].get_coord()
                            )
                            if dist <= 4.5:
                                next_res_n = next_res["N"] if "N" in next_res else None

                    phi = None
                    if (
                        prev_res_c is not None
                        and n is not None
                        and ca is not None
                        and c is not None
                    ):
                        phi = calc_dihedral(
                            prev_res_c.get_coord(),
                            n.get_coord(),
                            ca.get_coord(),
                            c.get_coord(),
                        )

                    psi = None
                    if (
                        n is not None
                        and ca is not None
                        and c is not None
                        and next_res_n is not None
                    ):
                        psi = calc_dihedral(
                            n.get_coord(),
                            ca.get_coord(),
                            c.get_coord(),
                            next_res_n.get_coord(),
                        )

                    omega = None
                    if prev_res_c is not None and n is not None and ca is not None:
                        prev_ca = prev_res["CA"] if "CA" in prev_res else None
                        if prev_ca is not None:
                            omega = calc_dihedral(
                                prev_ca.get_coord(),
                                prev_res_c.get_coord(),
                                n.get_coord(),
                                ca.get_coord(),
                            )
                            if (
                                res_name == "PRO"
                                and omega is not None
                                and abs(omega) < 90.0
                            ):
                                rama_type = "cis-Pro"

                    results.append(
                        {
                            "chain": chain.id,
                            "resSeq": res.id[1],
                            "icode": res.id[2],
                            "resName": res_name,
                            "phi": float(phi) if phi is not None else None,
                            "psi": float(psi) if psi is not None else None,
                            "omega": float(omega) if omega is not None else None,
                            "rama_type": rama_type,
                        }
                    )
                except Exception:
                    # Fallback for unexpected issues with single residues
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


def generate_csv(phi_psi_data):
    output = io.StringIO()
    writer = csv.writer(output)

    # Headers
    writer.writerow(
        [
            "chain",
            "residue number",
            "icode",
            "residue type",
            "phi",
            "psi",
            "omega",
            "interpolated percentage",
            "rama type",
            "category",
        ]
    )

    for p in phi_psi_data:
        # Only include if we have phi and psi
        if p["phi"] is not None and p["psi"] is not None:
            writer.writerow(
                [
                    p["chain"],
                    p["resSeq"],
                    p["icode"].strip(),
                    p["resName"],
                    f"{p['phi']:.3f}",
                    f"{p['psi']:.3f}",
                    f"{p['omega']:.3f}" if p["omega"] is not None else "",
                    f"{p['score']:.4f}" if p.get("score") is not None else "",
                    p["rama_type"],
                    p.get("classification", ""),
                ]
            )

    return output.getvalue()


def generate_pdf_report(phi_psi_data, pdb_id, rama_manager):
    """
    Generates a 6-panel Ramachandran PDF report similar to MolProbity/old_script.
    """
    # Filter out residues without phi/psi
    valid_data = [
        p for p in phi_psi_data if p["phi"] is not None and p["psi"] is not None
    ]
    if not valid_data:
        raise ValueError("No valid phi/psi data to plot.")

    # Create the figure
    fig, axes = plt.subplots(
        3, 2, figsize=(8.0, 12.0), layout="constrained", sharex=True, sharey=True
    )

    # Define color for levels
    # levels_colors = ("#5F5FFF", "#57A1EB") # From old script
    levels_colors = ("#7e48db", "#3036e7")  # Matching the web app's purple/blue

    for i, (ax, rama_type) in enumerate(zip(axes.flat, RAMA_KEYS)):
        # Data for this type
        type_data = [p for p in valid_data if p["rama_type"] == rama_type]

        # 1. Plot reference contours
        if rama_type in rama_manager.rama_data:
            ref = rama_manager.rama_data[rama_type]
            grid = ref["grid"]
            phi_grid, psi_grid = np.meshgrid(grid["phi"], grid["psi"])

            ax.contour(
                phi_grid,
                psi_grid,
                np.array(ref["dist"]),
                levels=ref["levels"],
                colors=levels_colors,
                linewidths=1.0,
                zorder=10,
            )

        # 2. Scatter Points
        if type_data:
            # Non-outliers
            normal_phis = [
                p["phi"] for p in type_data if p.get("classification") != "outlier"
            ]
            normal_psis = [
                p["psi"] for p in type_data if p.get("classification") != "outlier"
            ]

            if normal_phis:
                ax.scatter(
                    normal_phis,
                    normal_psis,
                    marker="o",
                    s=10.0,
                    fc="none",
                    ec="0.4",
                    zorder=20,
                )

            # Outliers
            outliers = [p for p in type_data if p.get("classification") == "outlier"]
            if outliers:
                outlier_phis = [p["phi"] for p in outliers]
                outlier_psis = [p["psi"] for p in outliers]
                ax.scatter(
                    outlier_phis,
                    outlier_psis,
                    marker="o",
                    s=10.0,
                    fc="none",
                    ec="#ef4444",  # Red for outliers
                    zorder=20,
                )

                # Labels for outliers
                for row in outliers:
                    icode_str = f"{row['icode']}".strip()
                    label = f"{row['chain']} {row['resSeq']}{icode_str} {row['resName'].capitalize()}"
                    ax.text(
                        row["phi"] + 5.0,
                        row["psi"],
                        label,
                        fontsize=6,
                        verticalalignment="center",
                        zorder=30,
                    )

        # Aesthetics
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-180, 181, 60))
        ax.set_xlim(-180.0, 180.0)
        ax.set_ylim(-180.0, 180.0)
        ax.set_title(rama_type, fontsize=10, fontweight="bold")

        if i >= 4:
            ax.set_xlabel(r"$\Phi$ (°)")
        if i % 2 == 0:
            ax.set_ylabel(r"$\Psi$ (°)")

        ax.grid(linestyle="dashed", linewidth=0.5, alpha=0.5)

    fig.suptitle(
        f"Ramachandran Report: {pdb_id.upper()}", fontsize=16, fontweight="bold"
    )

    # Save to buffer
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf")
    buf.seek(0)
    plt.close(fig)  # Clean up
    return buf
