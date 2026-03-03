import os
from pathlib import Path

import pytest

# Ensure matplotlib/font caches are writable in sandboxed test environments.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-config")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg-cache")


@pytest.fixture()
def app_env(tmp_path, monkeypatch):
    import app as app_module

    upload_dir = tmp_path / "uploads"
    upload_dir.mkdir(parents=True, exist_ok=True)

    app_module.app.config.update(
        TESTING=True,
        UPLOAD_FOLDER=str(upload_dir),
    )
    monkeypatch.setattr(app_module, "maybe_cleanup_upload_folder", lambda: None)
    return app_module, upload_dir


@pytest.fixture()
def client(app_env):
    app_module, _ = app_env
    with app_module.app.test_client() as test_client:
        yield test_client


@pytest.fixture()
def upload_dir(app_env):
    _, upload_path = app_env
    return upload_path


@pytest.fixture(scope="session")
def tripeptide_path():
    return Path(__file__).parent / "fixtures" / "tripeptide.pdb"
