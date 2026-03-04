import io
import json


def test_reference_endpoint_returns_cacheable_payload(client):
    response = client.get("/reference")

    assert response.status_code == 200
    assert response.headers.get("ETag")
    assert "max-age" in (response.headers.get("Cache-Control") or "")
    assert response.headers.get("X-Request-ID")

    payload = response.get_json()
    assert isinstance(payload, dict)
    assert "General" in payload


def test_reference_endpoint_supports_conditional_get(client):
    first = client.get("/reference")
    etag = first.headers.get("ETag")
    assert etag

    second = client.get("/reference", headers={"If-None-Match": etag})
    assert second.status_code == 304
    assert second.get_data() == b""


def test_process_requires_input(client):
    response = client.post("/process", data={})
    assert response.status_code == 400
    assert response.get_json()["error"] == "No PDB ID or file provided"


def test_process_rejects_invalid_pdb_id(client):
    response = client.post("/process", data={"pdb_id": "ABCD"})
    assert response.status_code == 400
    assert "Invalid structure identifier format" in response.get_json()["error"]


def test_process_rejects_unsupported_upload_extension(client):
    response = client.post(
        "/process",
        data={"pdb_file": (io.BytesIO(b"not-a-structure"), "file.txt")},
        content_type="multipart/form-data",
    )
    assert response.status_code == 400
    assert "Unsupported file format" in response.get_json()["error"]


def test_process_handles_pdb_fetch_errors(client, app_env, monkeypatch):
    app_module, _ = app_env
    monkeypatch.setattr(app_module, "fetch_structure_file", lambda *_, **__: (None, "Fetch failed"))

    response = client.post("/process", data={"pdb_id": "1UBQ"})
    assert response.status_code == 400
    assert response.get_json()["error"] == "Fetch failed"


def test_process_handles_uniprot_fetch_errors(client, app_env, monkeypatch):
    app_module, _ = app_env
    monkeypatch.setattr(app_module, "fetch_alphafold_model", lambda *_, **__: (None, "AF fetch failed"))

    response = client.post("/process", data={"pdb_id": "P69905"})
    assert response.status_code == 400
    assert response.get_json()["error"] == "AF fetch failed"


def test_process_upload_returns_result_id_and_persists_result(client, upload_dir, tripeptide_path):
    with tripeptide_path.open("rb") as handle:
        response = client.post(
            "/process",
            data={"pdb_file": (handle, tripeptide_path.name)},
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    payload = response.get_json()
    assert "result_id" in payload
    assert "phi_psi" in payload
    assert "reference" not in payload
    assert response.headers.get("X-Request-ID")

    result_file = upload_dir / f"{payload['result_id']}_results.json"
    assert result_file.exists()
    persisted = json.loads(result_file.read_text())
    assert persisted["pdb_id"] == tripeptide_path.name
    assert isinstance(persisted["phi_psi"], list)


def test_download_endpoints_serve_csv_and_pdf(client, tripeptide_path):
    with tripeptide_path.open("rb") as handle:
        process_response = client.post(
            "/process",
            data={"pdb_file": (handle, tripeptide_path.name)},
            content_type="multipart/form-data",
        )

    result_id = process_response.get_json()["result_id"]

    csv_response = client.get(f"/download/csv/{result_id}")
    assert csv_response.status_code == 200
    assert csv_response.mimetype == "text/csv"
    assert csv_response.get_data().startswith(b"chain,residue number")

    pdf_response = client.get(f"/download/pdf/{result_id}")
    assert pdf_response.status_code == 200
    assert pdf_response.mimetype == "application/pdf"
    assert pdf_response.get_data().startswith(b"%PDF")


def test_download_endpoints_return_404_for_unknown_result(client):
    csv_response = client.get("/download/csv/does-not-exist")
    pdf_response = client.get("/download/pdf/does-not-exist")

    assert csv_response.status_code == 404
    assert pdf_response.status_code == 404
