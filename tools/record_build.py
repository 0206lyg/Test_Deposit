#!/usr/bin/env python3
"""Record or verify the project-local solver, sources and OpenFOAM toolchain."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def command_info(command, args):
    executable = shutil.which(command)
    if not executable:
        raise RuntimeError(f"Required tool is unavailable: {command}")
    result = subprocess.run([executable, *args], text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=20, check=False)
    if result.returncode:
        raise RuntimeError(f"{command} {' '.join(args)} failed ({result.returncode}).")
    return {"path": str(Path(executable).resolve()),
            "version": result.stdout.strip().splitlines()[0] if result.stdout.strip() else "",
            "sha256": digest(executable)}


def source_hashes():
    paths = [p for directory in (ROOT / "src", ROOT / "build/rheoToolSubset")
             for p in directory.rglob("*")
             if p.is_file() and "lnInclude" not in p.parts
             and (p.suffix in {".C", ".H", ".cpp", ".h"}
                  or (p.parent.name == "Make" and p.name in {"files", "options"}))]
    paths += [ROOT / "tools" / name for name in ("prepare_rheotool.py", "fetch_dependencies.py")]
    return {str(p.relative_to(ROOT)): digest(p) for p in sorted(paths) if p.is_file()}


def snapshot():
    for name in ("WM_PROJECT_DIR", "WM_PROJECT_VERSION", "WM_OPTIONS", "FOAM_USER_APPBIN", "FOAM_USER_LIBBIN"):
        if not os.environ.get(name):
            raise RuntimeError(f"{name} is unset. Source tools/load_environment.sh first.")
    if os.environ["WM_PROJECT_VERSION"] != "8":
        raise RuntimeError("The project requires Foundation OpenFOAM-8.")
    app = Path(os.environ["FOAM_USER_APPBIN"]) / "rheoDepositionFoam"
    library = Path(os.environ["FOAM_USER_LIBBIN"]) / "libdepositionConstitutiveEquations.so"
    for path in (app, library):
        if not path.is_file() or ROOT not in path.resolve().parents:
            raise RuntimeError(f"Missing project-local build output: {path}")
    if not os.access(app, os.X_OK):
        raise RuntimeError(f"Solver is not executable: {app}")
    tools = {"compiler": command_info("g++", ["--version"]),
             "mpi_compiler": command_info("mpicxx", ["--showme:version"]),
             "mpi_launcher": command_info(os.environ.get("DEP_MPI_LAUNCHER", "mpirun"), ["--version"])}
    for name in ("wmake", "blockMesh"):
        path = shutil.which(name)
        if not path:
            raise RuntimeError(f"Required OpenFOAM tool is unavailable: {name}")
        tools[name] = {"path": str(Path(path).resolve()), "sha256": digest(path)}
    base_lib = Path(os.environ.get("FOAM_LIBBIN", "")) / "libOpenFOAM.so"
    if not base_lib.is_file():
        raise RuntimeError("Cannot find the active OpenFOAM libOpenFOAM.so.")
    return {"schema_version": 1, "project_root": str(ROOT),
            "openfoam_root": str(Path(os.environ["WM_PROJECT_DIR"]).resolve()),
            "openfoam_version": os.environ["WM_PROJECT_VERSION"],
            "wm_options": os.environ["WM_OPTIONS"], "tools": tools,
            "base_library": {"path": str(base_lib.resolve()), "sha256": digest(base_lib)},
            "compiled_sources": source_hashes(),
            "outputs": {str(p.relative_to(ROOT)): digest(p) for p in (app, library)},
            "package_inputs": {str(p.relative_to(ROOT)): digest(p)
                               for p in (ROOT / "tools/make_case.py", ROOT / "config/environment.sh") if p.is_file()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="reject a missing or stale build")
    args = parser.parse_args()
    manifest = ROOT / "build/build_manifest.json"
    try:
        current = snapshot()
        if args.check:
            if not manifest.is_file():
                raise RuntimeError("No completed build manifest. Run sbatch build_dep.sbatch first.")
            recorded = json.loads(manifest.read_text())
            mismatches = [key for key in ("schema_version", "project_root", "openfoam_root", "openfoam_version",
                          "wm_options", "tools", "base_library", "compiled_sources", "outputs")
                          if recorded.get(key) != current.get(key)]
            if mismatches:
                raise RuntimeError("Build/environment mismatch: " + ", ".join(mismatches)
                                   + ". Rebuild with the current toolchain before running.")
            if recorded.get("package_inputs") != current["package_inputs"]:
                print("Note: case-generator/configuration inputs changed since compilation; current inputs will be used.")
        executable = Path(os.environ["FOAM_USER_APPBIN"]) / "rheoDepositionFoam"
        result = subprocess.run([str(executable), "-help"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=30)
        if result.returncode:
            raise RuntimeError(f"Built solver cannot start (-help exited {result.returncode}). Check shared libraries.")
        if args.check:
            print("Build verified: project solver, compiled sources and OpenFOAM/MPI toolchain match.")
        else:
            current["recorded_at_utc"] = datetime.now(timezone.utc).isoformat()
            manifest.parent.mkdir(parents=True, exist_ok=True)
            temporary = manifest.with_suffix(f".tmp.{os.getpid()}")
            temporary.write_text(json.dumps(current, indent=2) + "\n")
            temporary.replace(manifest)
            print(f"Recorded completed build: {manifest}")
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"Build validation failed: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
