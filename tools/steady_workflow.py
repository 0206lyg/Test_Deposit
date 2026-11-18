#!/usr/bin/env python3
"""Run or submit the planar contraction cases. Every continuation stage must converge."""
import argparse
import contextlib
from datetime import datetime, timezone
from decimal import Decimal
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time

from check_run import atomic_json, checkpoint, complete_time, fields_from_initial, prepare_resume

ROOT = Path(__file__).resolve().parents[1]
BASELINE = ROOT / "config/experiment.json"
MODELS = ("gly75", "cross", "peo")
LEVELS = ("smoke", "coarse", "medium", "fine")
FRACTIONS = (Decimal("0.1"), Decimal("0.3"), Decimal("0.6"), Decimal("1"))
INTERNAL_FIELD = re.compile(r"(?m)^\s*internalField\s+((?:nonuniform|uniform)\s+[^;]*;)")
STAMP = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")


def normalize_flow(text):
    if text not in (".3", "0.3", ".5", "0.5", ".7", "0.7"):
        raise argparse.ArgumentTypeError("Choose 0.3, 0.5 or 0.7 mL/h.")
    return str(Decimal(text))


def positive(text):
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError("Choose a positive integer.")
    return value


def config_hash(config):
    return hashlib.sha256(config.read_bytes()).hexdigest()


def series_path(config):
    path = Path("steady2d")
    if config.resolve() != BASELINE:
        path /= "fit_" + config_hash(config)[:12]
    return path


def cases_for(flow, mesh, config):
    directory = ROOT / "runs" / series_path(config)
    return [directory / f"{key}_Q{flow.replace('.', 'p')}_{mesh}" for key in MODELS]


def stage_plan(final, fluid, flow):
    if fluid != "peo":
        return [(final, fluid, flow)]
    stages = final.with_name(final.name + ".stages")
    rates = [format((Decimal(flow)*fraction).normalize(), "f") for fraction in FRACTIONS]
    plan = [(stages / "01_cross_seed", "cross", rates[0])]
    for i, rate in enumerate(rates):
        path = final if rate == flow else stages / f"{i+2:02d}_peo_Q{rate.replace('.', 'p')}"
        plan.append((path, "peo", rate))
    return plan


def run(command, log=None):
    if log is None:
        subprocess.run(list(map(str, command)), check=True)
    else:
        with log.open("w") as out:
            result = subprocess.run(list(map(str, command)), stdout=out, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{command[0]} exited {result.returncode}. See {log}.")


def transfer(source, dest, iteration):
    """Copy only internal fields on identical meshes. Preserve the destination BCs."""
    if (source/"system/blockMeshDict").read_bytes() != (dest/"system/blockMeshDict").read_bytes():
        raise ValueError("Continuation needs identical source and destination meshes.")
    source_fields = set(fields_from_initial(source))
    target_fields = set(fields_from_initial(dest))
    for name in sorted(source_fields & target_fields):
        src = source / iteration / name
        if src.is_file():
            text = src.read_text()
        else:
            with gzip.open(str(src)+".gz", "rt") as f:
                text = f.read()
        match = INTERNAL_FIELD.search(text)
        target = dest / "0" / name
        original = target.read_text()
        found = INTERNAL_FIELD.search(original)
        if not match or not found:
            raise ValueError(f"Could not find ASCII internalField in {name}.")
        target.write_text(original[:found.start(1)] + match[1] + original[found.end(1):])
    # The pressure-corrected viscous face flux is already mass-conserving at
    # the first elastic stage. At changed Q, OpenFOAM initializes phi from U
    # and the newly prescribed inlet, then applies the pressure correction.
    source_rate = json.loads((source/"metadata.json").read_text())["Q_mL_h"]
    target_rate = json.loads((dest/"metadata.json").read_text())["Q_mL_h"]
    if source_rate == target_rate:
        for name in ("phi", "phi.gz"):
            flux = source/iteration/name
            if flux.is_file():
                shutil.copyfile(flux, dest/"0"/name)
                break


def stage_identity(args, fluid, rate):
    return {"format": "steady2d-v1", "fluid": fluid, "Q_mL_h": rate, "mesh": args.mesh,
            "ranks": args.ranks, "config_sha256": args.config_sha256,
            "stress_scheme": args.stress_scheme, "monitor_interval": args.monitor_interval}


def generate(stage, spec, args, predecessor):
    # Publish only fully initialized meshes. Interrupted setup can be safely retried.
    stage.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=stage.parent, prefix=f".{stage.name}-") as staging:
        temp = Path(staging) / "case"
        config_bytes = args.config.read_bytes()
        if hashlib.sha256(config_bytes).hexdigest() != spec["config_sha256"]:
            raise ValueError("The input configuration changed during this workflow. Resume with its original contents.")
        config_snapshot = Path(staging) / "config.json"
        config_snapshot.write_bytes(config_bytes)
        run([sys.executable, ROOT/"tools/make_case.py", "--fluid", spec["fluid"], "--flow",
             spec["Q_mL_h"], "--mesh", args.mesh, "--case", temp, "--ranks", args.ranks,
             "--max-iters", args.max_iters, "--write-interval", args.write_interval,
             "--monitor-interval", args.monitor_interval, "--config", config_snapshot,
             "--stress-scheme", args.stress_scheme])
        # Keep tool output on failure, including when the unpublished staging tree is removed.
        logs = ROOT / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        for tool in ("blockMesh", "checkMesh"):
            log = logs / f"{stage.name}.{tool}.{STAMP}.log"
            run([tool, "-case", temp], log)
            if tool == "checkMesh" and "Mesh OK." not in log.read_text():
                raise RuntimeError(f"checkMesh did not report Mesh OK. Inspect {log}.")
            shutil.copy(log, temp / f"log.{tool}")
        if predecessor:
            transfer(predecessor[0], temp, predecessor[1])
        run([sys.executable, ROOT/"tools/make_case.py", "--initialize-inlet", temp])
        if args.ranks > 1:
            log = logs / f"{stage.name}.decomposePar.{STAMP}.log"
            run(["decomposePar", "-case", temp], log)
            shutil.copy(log, temp/"log.decomposePar")
        atomic_json(temp/"steady_stage.json", spec)
        temp.rename(stage)


def set_control(stage, key, value):
    run(["foamDictionary", stage/"system/controlDict", "-entry", key, "-set", value],
        stage / f"log.control.{key}")


class Interruption:
    def __init__(self):
        self.requested = False
        self.solving = None

    def signal(self, signum, frame):
        self.requested = True
        if self.solving:
            print("Scheduler stop requested. Writing a restart checkpoint.", flush=True)
            try:
                set_control(self.solving, "stopAt", "writeNow")
            except (OSError, RuntimeError) as exc:
                print(str(exc), file=sys.stderr)

    def check(self):
        if self.requested:
            raise InterruptedError("Job interrupted. Resubmit this workflow with --resume.")


def inspect_run(stage, ranks, code=0, phase="finished", interrupted=False, seconds=None):
    command = [sys.executable, ROOT/"tools/check_run.py", "--case", stage, "--ranks", str(ranks),
               "--exit-code", str(code), "--phase", phase]
    if seconds is not None:
        command += ["--solve-seconds", str(seconds)]
    if interrupted:
        command.append("--interrupted")
    result = subprocess.run(command, check=False, stdout=subprocess.DEVNULL)
    if result.returncode:
        path = stage/"run_status.json"
        status = json.loads(path.read_text())["status"] if path.exists() else "failed"
        raise RuntimeError(f"{stage.name}: {status}. See {path} and {stage}/log.solver.*. "
                           "After an iteration cap, --resume continues from the last complete checkpoint.")
    return json.loads((stage/"run_status.json").read_text())["checkpoint"]


def reconstruct(stage, info, ranks):
    if ranks > 1:
        # Postprocessing reads additional output-only fields as well. Reconstruct all.
        run(["reconstructPar", "-case", stage, "-time", info["time_name"]],
            stage / f"log.reconstructPar.{STAMP}")
        if not complete_time(stage/info["time_name"], info["required_fields"]):
            raise RuntimeError(f"Reconstructed fields incomplete at {stage/info['time_name']}")


def advance(stage, spec, args, predecessor, stop):
    stop.check()
    existing = stage.exists()
    if existing:
        identity = stage / "steady_stage.json"
        if not args.resume or not identity.exists() or json.loads(identity.read_text()) != spec:
            raise ValueError(f"Existing stage {stage} differs in settings or requires --resume.")
    else:
        generate(stage, spec, args, predecessor)
    stop.check()
    if existing and (stage/"postProcessing/deposition/stationary.json").exists():
        try:
            converged = inspect_run(stage, args.ranks)
        except RuntimeError:
            converged = None
        if converged:
            if not complete_time(stage / converged["time_name"], converged["required_fields"]):
                reconstruct(stage, converged, args.ranks)
            print(f"Reusing converged {stage.name} at iteration {converged['time_name']}.")
            return stage, converged["time_name"]

    start = 0
    if existing:
        # Initial fields are complete by construction even if no evolved write exists.
        info = checkpoint(stage, args.ranks)
        if info["time"] > 0:
            prepare_resume(stage, info)
            start = int(info["time"])
        else:
            # Archive diagnostics from a failure before the first periodic write.
            diag = stage/"postProcessing/deposition"
            if diag.exists():
                diag.rename(diag.with_name("deposition.prior."+STAMP))
        for key, value in (("startFrom", "startTime"), ("startTime", start),
                           ("endTime", start + args.max_iters), ("stopAt", "endTime"),
                           ("writeInterval", args.write_interval)):
            set_control(stage, key, value)

    command = [str(Path(os.environ["FOAM_USER_APPBIN"])/"rheoDepositionFoam"), "-case", str(stage)]
    if args.ranks > 1:
        command = [os.environ.get("DEP_MPI_LAUNCHER", "mpirun"), "-np", str(args.ranks)] + command + ["-parallel"]
    log = stage / f"log.solver.{STAMP}"
    print(f"{spec['fluid']} Q={spec['Q_mL_h']} mL/h, {stage}, iteration {start}, log {log}")
    began = time.monotonic()
    with log.open("w") as out:
        stop.solving = stage
        try:
            stop.check()
            result = subprocess.run(command, stdout=out, stderr=subprocess.STDOUT)
        finally:
            stop.solving = None
    elapsed = time.monotonic() - began
    if stop.requested:
        with contextlib.suppress(RuntimeError):
            inspect_run(stage, args.ranks, result.returncode, "solve", True, elapsed)
        stop.check()
    info = inspect_run(stage, args.ranks, result.returncode, seconds=elapsed)
    try:
        reconstruct(stage, info, args.ranks)
    except RuntimeError:
        with contextlib.suppress(RuntimeError):
            inspect_run(stage, args.ranks, 1, "reconstruct")
        raise
    print(f"Converged at iteration {info['time_name']}; solver wall time for this invocation {elapsed:.2f} s.")
    return stage, info["time_name"]


def execute(args):
    config_digest = config_hash(args.config)
    if args.config_sha256 and args.config_sha256 != config_digest:
        raise ValueError("The submitted configuration has changed. Resubmit with its current contents.")
    args.config_sha256 = config_digest
    final = (args.case or cases_for(args.flow, args.mesh, args.config)[MODELS.index(args.fluid)]).resolve()
    final.parent.mkdir(parents=True, exist_ok=True)
    plan = stage_plan(final, args.fluid, args.flow)
    stop = Interruption()
    for sig in (signal.SIGUSR1, signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, stop.signal)
    lock_path = final.with_name("." + final.name + ".lock")
    with lock_path.open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f"Another workflow is running {final}.")
        if not args.resume and any(stage.exists() for stage, _, _ in plan):
            raise ValueError(f"Existing case or continuation stage for {final}. Use --resume or --case.")
        predecessor = None
        for stage, fluid, rate in plan:
            print(f"Stage {stage.name} ({fluid}, {rate} mL/h)", flush=True)
            predecessor = advance(stage, stage_identity(args, fluid, rate), args, predecessor, stop)
        print(f"Target solution: {final/'case.foam'}")


def submit_job(script, args, dependencies=(), flags=()):
    cmd = ["sbatch", "--parsable", *flags]
    if dependencies:
        cmd.extend(["--kill-on-invalid-dep=yes", "--dependency=afterok:" + ":".join(dependencies)])
    result = subprocess.check_output([*cmd, str(ROOT/script), *args], text=True).strip()
    jobid = result.split(";")[0]
    if not jobid.isdigit():
        raise RuntimeError(f"Unexpected sbatch result {result!r}.")
    return jobid


def accepted_or_resume(case, fluid, level, args, digest):
    # Avoid queueing a completed coarse sweep again just to compare meshes.
    plan = stage_plan(case, fluid, args.flow)
    existing = [(p, model, rate) for p, model, rate in plan if p.exists()]
    if not existing:
        return False, False, args.ranks
    records = []
    for path, model, rate in existing:
        identity = path/"steady_stage.json"
        if not identity.is_file():
            raise ValueError(f"Missing stage identity: {identity}")
        spec = json.loads(identity.read_text())
        if (spec.get("format") != "steady2d-v1" or spec.get("fluid") != model or
            spec.get("Q_mL_h") != rate or spec.get("mesh") != level or
            spec.get("config_sha256") != digest or spec.get("stress_scheme") != args.stress_scheme or
            spec.get("monitor_interval") != args.monitor_interval):
            raise ValueError(f"Existing settings differ for {path}; choose a separate case or use the original configuration.")
        records.append(spec)
    ranks = records[-1]["ranks"]
    if (case/"postProcessing/deposition/stationary.json").is_file():
        try:
            info = inspect_run(case, ranks)
            if complete_time(case/info["time_name"], info["required_fields"]):
                return True, False, ranks
        except RuntimeError:
            pass
    if args.ranks is not None and args.ranks != ranks:
        raise ValueError(f"Incomplete {case.name} uses {ranks} MPI ranks. Resume with that count.")
    return False, True, ranks


def submit(args):
    config = args.config.resolve()
    digest = config_hash(config)
    common = ["--config", str(config), "--config-sha256", digest]
    upstream = [args.after_build] if args.after_build else []
    jobs, pending = [], []
    for level in filter(None, [args.mesh, args.compare_mesh]):
        for fluid, case in zip(MODELS, cases_for(args.flow, level, config)):
            accepted, resume, ranks = accepted_or_resume(case, fluid, level, args, digest)
            if accepted:
                print(f"Reusing converged {fluid}, {args.flow} mL/h, {level}: {case}")
            else:
                pending.append((level, fluid, resume, ranks))
    for level, fluid, resume, ranks in pending:
        options = [fluid, args.flow, level, *common, "--max-iters", str(args.max_iters),
                   "--write-interval", str(args.write_interval), "--monitor-interval", str(args.monitor_interval),
                   "--stress-scheme", args.stress_scheme]
        if resume:
            options.append("--resume")
        flags = ["--job-name=dep_"+fluid+"_"+level]
        if ranks is not None:
            flags.append("--ntasks="+str(ranks))
        job = submit_job("run_dep.sbatch", options, upstream, flags)
        jobs.append(job)
        print(f"{fluid}, {args.flow} mL/h, {level}: {job}", flush=True)
    post_args = [args.flow, args.mesh, *common]
    if args.compare_mesh:
        post_args += ["--compare-mesh", args.compare_mesh]
    print("Wall comparison job:", submit_job("postprocess_dep.sbatch", post_args, jobs or upstream))


def post(args):
    if args.config_sha256 and args.config_sha256 != config_hash(args.config):
        raise ValueError("Configuration changed after submission.")
    tag = f"Q{args.flow.replace('.', 'p')}_{args.mesh}"
    out = ROOT/"comparison"/series_path(args.config)/(tag if not args.compare_mesh else tag+"_vs_"+args.compare_mesh)
    cmd = [sys.executable, ROOT/"tools/compare_results.py", *cases_for(args.flow, args.mesh, args.config), "--output", out]
    if args.compare_mesh:
        cmd += ["--refined", *cases_for(args.flow, args.compare_mesh, args.config)]
    run(cmd)
    print(f"Comparison: {out}")


def main():
    action = sys.argv[1] if sys.argv[1:2] in (["submit"], ["post"]) else "solve"
    parser = argparse.ArgumentParser(description=__doc__)
    if action == "solve":
        parser.add_argument("fluid", choices=MODELS)
    parser.add_argument("flow", nargs="?", default="0.5", type=normalize_flow)
    parser.add_argument("mesh", nargs="?", choices=LEVELS, default="coarse")
    parser.add_argument("--config", type=Path, default=Path(os.environ.get("DEP_CONFIG", BASELINE)))
    parser.add_argument("--config-sha256", help=argparse.SUPPRESS)
    if action in ("submit","post"):
        parser.add_argument("--compare-mesh", choices=LEVELS, help="additional, finer mesh level for the same three models")
    if action == "submit":
        parser.add_argument("--after-build", help="wait for this successful SLURM build job ID")
    if action in ("solve","submit"):
        parser.add_argument("--ranks", type=positive, default=(os.getenv("SLURM_NTASKS", os.getenv("DEP_RANKS", "1")) if action=="solve" else None))
        parser.add_argument("--max-iters", "--max-iterations", type=positive, default=12000,
                            help="maximum additional iterations per unfinished stage and invocation")
        parser.add_argument("--write-interval", type=positive, default=200, help="checkpoint interval in steady iterations")
        parser.add_argument("--monitor-interval", type=positive, default=10)
        parser.add_argument("--stress-scheme", choices=("upwind","limitedLinear"), default="upwind")
    if action == "solve":
        parser.add_argument("--case", type=Path, help="override target case directory")
        parser.add_argument("--resume", action="store_true", help="reuse converged stages and restart the first unfinished stage")
    args = parser.parse_args(sys.argv[2:] if action != "solve" else None)
    if getattr(args, "after_build", None) and not args.after_build.isdigit():
        parser.error("--after-build needs a numeric job ID.")
    if getattr(args, "compare_mesh", None) and LEVELS.index(args.compare_mesh) <= LEVELS.index(args.mesh):
        parser.error("--compare-mesh must be finer than the primary mesh.")
    if getattr(args, "ranks", None) and action == "solve" and "SLURM_NTASKS" in os.environ and args.ranks > int(os.environ["SLURM_NTASKS"]):
        parser.error("Requested rank count exceeds the Slurm allocation.")
    try:
        {"submit": submit, "post": post, "solve": execute}[action](args)
    except InterruptedError as exc:
        print(exc, file=sys.stderr)
        return 75
    except (ValueError, RuntimeError, OSError, KeyError, subprocess.CalledProcessError) as exc:
        print(exc, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
