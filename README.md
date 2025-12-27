# Ramachandran Plot Analysis Tool

A web-based application for analyzing and visualizing Ramachandran plots of protein structures. This tool calculates backbone dihedral angles (φ and ψ) and classifies residues into **Favored**, **Allowed**, and **Outlier** regions based on high-quality reference data.

## Live Demo

**Try it now**: [ramachandran.onrender.com](https://ramachandran.onrender.com)

*Note: The demo runs on a free tier, so initial startup may take 30-60 seconds.*

## Features

- **PDB/mmCIF Support**: Upload your own structure files or fetch them directly from the RCSB PDB by ID.
- **Accurate Classification**: Uses reference distributions derived from the **Top8000** dataset of high-quality protein structures. These distributions are sourced from the [Richardson Lab's reference data repository](https://github.com/rlabduke/reference_data) and cover 6 specific residue types:
  - General
  - Ile/Val
  - Glycine
  - pre-Proline
  - trans-Proline
  - cis-Proline
- **Interactive Visualization**: Dynamic plotting with contour lines representing favorable and allowed regions.
- **Statistical Summary**: Instant breakdown of residue distributions and outliers.
- **Data Export**: Download analysis results as a CSV file for further study.

## References

This tool uses reference data and validation criteria based on the following publication:

*   **Chen, V. B., et al. (2010).** *MolProbity: all-atom structure validation for macromolecular crystallography.* **Acta Crystallographica Section D**, 66(1), 12–21. [doi:10.1107/S0907444909042073](https://doi.org/10.1107/S0907444909042073)
*   **Williams, C. J., et al. (2018).** *MolProbity: More and better reference data for improved all-atom structure validation.* **Protein Science**, 27(1), 293–301. [doi:10.1002/pro.3330](https://doi.org/10.1002/pro.3330)


The underlying Ramachandran distribution data (rama8000) is provided by the Richardson Lab at Duke University: [github.com/rlabduke/reference_data](https://github.com/rlabduke/reference_data).

## Getting Started

### Prerequisites

- **Python**: 3.13 or higher
- **uv**: Recommended for dependency management (see [uv documentation](https://github.com/astral-sh/uv))

### Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/ramachandran.git
   cd ramachandran
   ```

2. **Install dependencies**:
   Using `uv`:
   ```bash
   uv sync
   ```
   Or using `pip`:
   ```bash
   pip install .
   ```

### Running the Application

Start the Flask development server:

```bash
uv run python app.py
```

The application will be available at `http://localhost:5001`.

## Project Structure

- `app.py`: Flask application and API endpoints.
- `utils.py`: Core logic for PDB parsing, dihedral calculation, and Ramachandran classification.
- `data/`: Reference data files for Ramachandran distributions.
- `static/`: Frontend assets (CSS, JS).
- `templates/`: HTML templates.

## Deployment

This application can be deployed to platforms like Render, Heroku, or similar services.

### Deploying to Render

1. Ensure `gunicorn` is added as a dependency:
   ```bash
   uv add gunicorn
   ```

2. Use the following start command:
   ```bash
   uv run gunicorn app:app
   ```

3. Render will automatically detect and install dependencies from `pyproject.toml` using `uv sync`.

A live demo is available at [ramachandran.onrender.com](https://ramachandran.onrender.com).

## License

This project is open-source and available under the MIT License.
