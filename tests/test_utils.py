import pytest

from app import is_allowed_upload, is_valid_pdb_id
from utils import (
    RamachandranManager,
    fetch_alphafold_model,
    fetch_structure_file,
    generate_csv,
    get_phi_psi,
    is_valid_uniprot_accession,
    parse_structure,
)


def test_is_valid_pdb_id_accepts_legacy_and_extended_formats():
    assert is_valid_pdb_id("1UBQ")
    assert is_valid_pdb_id("pdb_00001abc")
    assert is_valid_pdb_id("PDB_00001ABC")


def test_is_valid_pdb_id_rejects_invalid_formats():
    assert not is_valid_pdb_id("ABCD")
    assert not is_valid_pdb_id("0ABC")
    assert not is_valid_pdb_id("123")


def test_is_valid_uniprot_accession_accepts_known_formats():
    assert is_valid_uniprot_accession("P69905")
    assert is_valid_uniprot_accession("A0A024RBG1")
    assert is_valid_uniprot_accession("q8n158")


def test_is_valid_uniprot_accession_rejects_invalid_formats():
    assert not is_valid_uniprot_accession("1UBQ")
    assert not is_valid_uniprot_accession("INVALID")
    assert not is_valid_uniprot_accession("P1234")


def test_is_allowed_upload_accepts_supported_extensions():
    assert is_allowed_upload("example.pdb")
    assert is_allowed_upload("example.cif")
    assert is_allowed_upload("example.mmcif")
    assert not is_allowed_upload("example.txt")


def test_parse_structure_raises_value_error_for_missing_path(tmp_path):
    missing_path = tmp_path / "missing_file.pdb"
    with pytest.raises(ValueError, match="Unable to parse structure file"):
        parse_structure(str(missing_path))


def test_fetch_structure_file_rejects_empty_pdb_id(tmp_path):
    path, error = fetch_structure_file("   ", output_dir=str(tmp_path))
    assert path is None
    assert error == "PDB ID is required."


def test_fetch_alphafold_model_downloads_latest_prediction(monkeypatch, tmp_path):
    class FakeResponse:
        def __init__(self, body, status=200):
            self._body = body
            self.status = status

        def read(self):
            return self._body

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    def fake_urlopen(url, timeout=15):
        if "api/prediction" in url:
            return FakeResponse(
                b'[{"latestVersion": 2, "cifUrl": "https://alphafold.ebi.ac.uk/files/new.cif"},'
                b' {"latestVersion": 1, "cifUrl": "https://alphafold.ebi.ac.uk/files/old.cif"}]'
            )
        if url.endswith("new.cif"):
            return FakeResponse(b"data_mock_cif")
        raise AssertionError(f"Unexpected URL: {url}")

    monkeypatch.setattr("utils.urllib.request.urlopen", fake_urlopen)

    path, error = fetch_alphafold_model("P69905", output_dir=str(tmp_path))
    assert error is None
    assert path is not None
    assert path.endswith("af-p69905.cif")
    assert (tmp_path / "af-p69905.cif").read_bytes() == b"data_mock_cif"


def test_fetch_alphafold_model_handles_empty_metadata(monkeypatch, tmp_path):
    class FakeResponse:
        def __init__(self, body, status=200):
            self._body = body
            self.status = status

        def read(self):
            return self._body

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr("utils.urllib.request.urlopen", lambda *_args, **_kwargs: FakeResponse(b"[]"))

    path, error = fetch_alphafold_model("P69905", output_dir=str(tmp_path))
    assert path is None
    assert error == "No AlphaFold prediction found for UniProt accession 'P69905'."


def test_get_phi_psi_produces_valid_angles(tripeptide_path):
    structure = parse_structure(str(tripeptide_path))
    data = get_phi_psi(structure)

    assert len(data) == 4
    valid = [item for item in data if item["phi"] is not None and item["psi"] is not None]
    assert len(valid) >= 1
    assert any(item["rama_type"] == "General" for item in data)
    assert any(item["rama_type"] == "Gly" for item in data)


def test_classify_unknown_rama_type_as_outlier():
    manager = RamachandranManager(data_directory="data")
    score, classification = manager.classify_phipsi("unknown", -60.0, -45.0)
    assert score == 0.0
    assert classification == "outlier"


def test_generate_csv_skips_entries_without_complete_angles():
    csv_content = generate_csv(
        [
            {
                "chain": "A",
                "resSeq": 1,
                "icode": " ",
                "resName": "ALA",
                "phi": None,
                "psi": -30.0,
                "omega": None,
                "score": None,
                "rama_type": "General",
                "classification": None,
            },
            {
                "chain": "A",
                "resSeq": 2,
                "icode": " ",
                "resName": "GLY",
                "phi": -60.0,
                "psi": -45.0,
                "omega": 180.0,
                "score": 99.9,
                "rama_type": "Gly",
                "classification": "favoured",
            },
        ]
    )

    lines = csv_content.strip().splitlines()
    assert len(lines) == 2
    assert lines[0].startswith("chain,residue number")
    assert "GLY" in lines[1]
