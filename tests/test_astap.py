"""Finding and fetching the ASTAP solver.

The download is exercised against a mirror built in ``tmp_path`` and served over a
``file://`` URL, so the whole thing runs offline: archive layout, checksum verification,
unpacking, and the executable bit. What is *not* covered is running a real ASTAP binary,
which is what :func:`astap.verify` is for -- so ``verify`` is stubbed where a stub
executable would have to be a real program.

The resolution order is the other half, and matters more than the download: most people
already have ASTAP, and the point of the change is that they no longer have to say where
it is.
"""

import json
import os
import zipfile

import pytest

from seestar_photometry import astap, astrometry


# --- finding an existing install ----------------------------------------------------

@pytest.fixture(autouse=True)
def isolated(monkeypatch, tmp_path):
    """Keep every test off the real cache and the real PATH."""
    monkeypatch.setenv("SEESTAR_ASTAP_DATA", str(tmp_path / "cache"))
    monkeypatch.delenv("ASTAP_EXE", raising=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty-path"))
    return tmp_path


def test_explicit_path_wins(tmp_path):
    exe = tmp_path / "my_astap.exe"
    exe.write_text("")
    assert astap.executable(exe) == exe


def test_env_var_comes_next(tmp_path, monkeypatch):
    exe = tmp_path / "from_env"
    exe.write_text("")
    monkeypatch.setenv("ASTAP_EXE", str(exe))
    assert astap.executable() == exe


def test_path_lookup_finds_a_system_install(tmp_path, monkeypatch):
    """`apt install astap-cli` and the AUR both put `astap_cli` on PATH.

    This is the case that used to need configuring by hand on every non-Windows
    machine, because the default was a hardcoded ``C:\\Program Files`` path.
    """
    bindir = tmp_path / "bin"
    bindir.mkdir()
    name = "astap_cli.exe" if os.name == "nt" else "astap_cli"
    exe = bindir / name
    exe.write_text("")
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir))
    assert astap.executable() == exe


def test_downloaded_copy_is_used_when_nothing_else_exists(tmp_path):
    target = astap.downloaded_executable()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("")
    assert astap.executable() == target
    assert astap.is_installed()


def test_nothing_installed_resolves_to_none():
    # The Windows stock path is the last resort and does not exist on a clean runner;
    # on a developer machine that has ASTAP installed it legitimately does.
    found = astap.executable()
    assert found is None or found.exists()


def test_the_error_names_every_way_out(monkeypatch):
    """A missing solver has three fixes and the message has to give all of them."""
    monkeypatch.setattr(astap, "executable", lambda explicit=None: None)
    with pytest.raises(RuntimeError) as excinfo:
        astrometry.astap_executable()
    message = str(excinfo.value)
    assert "astap-cli" in message          # system package
    assert "astap.download()" in message   # let us fetch it
    assert "solver='local'" in message     # do without


def test_platform_key_is_recognised():
    key = astap.platform_key()
    assert key in astap.EXECUTABLES


# --- the download, against a file:// mirror -------------------------------------------

def build_mirror(directory, with_database=True, corrupt=False):
    """A mirror with a stub binary and a stub star database. Returns its URL."""
    import hashlib

    directory.mkdir(parents=True, exist_ok=True)
    key = astap.platform_key()
    meta = {"version": astap.VERSION, "platforms": {}, "databases": {}}

    binary = directory / f"astap_cli_{key}.zip"
    with zipfile.ZipFile(binary, "w") as zf:
        # A nested directory, because the upstream archives are inconsistent about it.
        zf.writestr(f"astap/{astap.EXECUTABLES[key]}", "#!/bin/sh\nexit 0\n")
    digest = hashlib.sha256(binary.read_bytes()).hexdigest()
    meta["platforms"][key] = {
        "file": binary.name,
        "sha256": "0" * 64 if corrupt else digest,
        "bytes": binary.stat().st_size,
    }

    if with_database:
        db = directory / f"{astap.DEFAULT_DATABASE}.zip"
        with zipfile.ZipFile(db, "w") as zf:
            for i in (1, 2):
                zf.writestr(f"{astap.DEFAULT_DATABASE}_010{i}.1476", "stub")
        meta["databases"][astap.DEFAULT_DATABASE] = {
            "file": db.name,
            "sha256": hashlib.sha256(db.read_bytes()).hexdigest(),
            "bytes": db.stat().st_size,
        }

    (directory / astap.MANIFEST).write_text(json.dumps(meta), encoding="utf-8")
    return directory.as_uri()


@pytest.fixture
def mirror(tmp_path, monkeypatch):
    url = build_mirror(tmp_path / "mirror")
    monkeypatch.setenv("SEESTAR_ASTAP_URL", url)
    # Running the stub would need it to be a real program on every platform.
    monkeypatch.setattr(astap, "verify", lambda path=None: path)
    return url


def test_download_installs_a_runnable_solver(mirror):
    path = astap.download(quiet=True)
    assert path.exists() and path == astap.downloaded_executable()
    assert astap.is_installed()
    assert astap.executable() == path
    if os.name != "nt":
        assert os.access(path, os.X_OK), "the execute bit must be set after unzipping"


def test_download_brings_the_star_database(mirror):
    astap.download(quiet=True)
    assert astap.has_database()
    assert list(astap.database_dir().glob("*.1476"))


def test_download_is_idempotent(mirror):
    first = astap.download(quiet=True)
    stamp = first.stat().st_mtime_ns
    second = astap.download(quiet=True)
    assert first == second
    assert second.stat().st_mtime_ns == stamp, "a second call must not re-fetch"


def test_database_can_be_skipped(mirror):
    astap.download(database=None, quiet=True)
    assert astap.downloaded_executable().exists()
    assert not astap.has_database()


def test_a_corrupt_download_is_refused(tmp_path, monkeypatch):
    """A truncated binary must fail loudly, not be unpacked and run."""
    monkeypatch.setenv("SEESTAR_ASTAP_URL", build_mirror(tmp_path / "bad", corrupt=True))
    monkeypatch.setattr(astap, "verify", lambda path=None: path)
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        astap.download(quiet=True)
    assert not astap.downloaded_executable().exists()


def test_an_unmirrored_database_says_what_there_is(tmp_path, monkeypatch):
    monkeypatch.setenv("SEESTAR_ASTAP_URL",
                       build_mirror(tmp_path / "m", with_database=False))
    monkeypatch.setattr(astap, "verify", lambda path=None: path)
    with pytest.raises(RuntimeError, match="no star database"):
        astap.download(quiet=True)


def test_verify_explains_a_binary_it_cannot_run(tmp_path):
    """The macOS quarantine case, which is the one that actually bites."""
    fake = tmp_path / "not_a_program"
    fake.write_text("this is not an executable")
    fake.chmod(0o644)
    if os.name == "nt":
        pytest.skip("Windows will not refuse to exec a text file the same way")
    with pytest.raises(RuntimeError, match="could not be run"):
        astap.verify(fake)


def test_verify_needs_something_to_verify(monkeypatch):
    # Explicitly nothing installed: a developer machine may well have a real one at the
    # stock Windows path, which is the last resort in the resolution order.
    monkeypatch.setattr(astap, "executable", lambda explicit=None: None)
    with pytest.raises(RuntimeError, match="download"):
        astap.verify()


# --- how solve_astap uses it -----------------------------------------------------------

def test_the_database_flag_is_only_passed_when_we_fetched_one(mirror, monkeypatch):
    """A system ASTAP knows where its own database is; overriding that is presumptuous.

    Checked by capturing the command rather than running it -- what matters is the
    argument list, and a real solve needs a real binary and a real frame.
    """
    seen = {}

    class _Result:
        returncode, stdout, stderr = 1, "", ""

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        return _Result()

    import subprocess

    import conftest

    from seestar_photometry import frames

    monkeypatch.setattr(subprocess, "run", fake_run)

    frames_dir = isolated_dir = astap.cache_dir().parent / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    frame = frames.load_frame(conftest.write_cube(isolated_dir / "f.fit"))

    astap.download(quiet=True)                    # gives us a database
    assert astap.has_database()
    with pytest.raises(RuntimeError, match="ASTAP no solution"):
        astrometry.solve_astap(frame, astap_exe=astap.downloaded_executable())
    assert "-d" in seen["cmd"]
    assert str(astap.database_dir()) in seen["cmd"]

    for stale in astap.database_dir().glob("*.1476"):
        stale.unlink()
    assert not astap.has_database()
    with pytest.raises(RuntimeError, match="ASTAP no solution"):
        astrometry.solve_astap(frame, astap_exe=astap.downloaded_executable(),
                               force=True)
    assert "-d" not in seen["cmd"]
