#!/usr/bin/env python3
"""Check saved fields and metric histories; solver completion is not stationarity."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import json
import math
import os
from pathlib import Path
import re
import sys


def atomic_json(path, value):
    temporary = path.with_suffix(f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def time_value(name):
    try:
        value = float(name)
        return value if math.isfinite(value) and value >= 0 else None
    except ValueError:
        return None


def fields_from_initial(case):
    fields = {"U", "p"}
    for path in (case / "0").glob("*"):
        if path.is_file():
            try:
                text = path.read_text(errors="ignore")[:8192]
                if re.search(r"\bclass\s+vol(?:Scalar|Vector|Tensor|SymmTensor)Field\s*;", text):
                    fields.add(path.name)
            except OSError:
                pass
    for name in list(fields):
        if name.startswith("theta_"):
            suffix = name[len("theta"):]
            fields.update({"eigVals" + suffix, "eigVecs" + suffix})
    return sorted(fields)


def complete_time(folder, fields):
    for name in fields:
        path = folder / name
        compressed = folder / (name + ".gz")
        if not path.is_file() and compressed.is_file():
            path = compressed
        if not path.is_file() or path.stat().st_size < 40:
            return False
        try:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as stream:
                head = stream.read(8192)
                if path.suffix == ".gz":
                    tail = head
                    for block in iter(lambda: stream.read(65536), b""):
                        tail = (tail + block)[-8192:]
                else:
                    stream.seek(max(0, path.stat().st_size - 8192))
                    tail = stream.read(8192)
            # OpenFOAM writes this footer only after the complete field. The
            # generated initial dictionaries intentionally have no footer.
            finished = re.search(rb"//\s*\*{5,}\s*//\s*$", tail) if time_value(folder.name) != 0 else tail.rstrip().endswith(b"}")
            if b"FoamFile" not in head or not finished:
                return False
        except (OSError, EOFError):
            return False
    return True


def checkpoint(case, ranks):
    fields = fields_from_initial(case)
    processors = sorted((p for p in case.glob("processor*") if p.is_dir()
                         and re.fullmatch(r"processor\d+", p.name)), key=lambda p: int(p.name[9:]))
    if processors:
        if [p.name for p in processors] != [f"processor{i}" for i in range(ranks)]:
            raise ValueError("Existing processor directories do not match the requested MPI rank count.")
        roots = processors
    else:
        roots = [case]
    times = None
    names = {}
    for root in roots:
        valid = {}
        for path in root.iterdir():
            value = time_value(path.name)
            if path.is_dir() and value is not None and complete_time(path, fields):
                valid[value] = path.name
        times = set(valid) if times is None else times & set(valid)
        names.update(valid)
    if not times:
        raise ValueError("No complete common checkpoint containing all initial fields was found.")
    latest = max(times)
    return {"time": latest, "time_name": names[latest], "parallel": bool(processors),
            "ranks": len(processors) if processors else 1, "required_fields": fields}


def metrics(path):
    rows, invalid, final_valid = {}, 0, True
    with path.open(newline="") as stream:
        reader = csv.DictReader(line for line in stream if line.strip() and not line.startswith("#"))
        if not reader.fieldnames or "time" not in reader.fieldnames:
            raise ValueError("Metrics CSV has no 'time' header.")
        fields = reader.fieldnames
        for raw in reader:
            try:
                if None in raw:
                    raise ValueError("unexpected CSV columns")
                row = {key: float(raw[key]) for key in fields}
                if not all(math.isfinite(value) for value in row.values()):
                    raise ValueError("nonfinite metric")
                rows[row["time"]] = row
                final_valid = True
            except (ValueError, TypeError, KeyError):
                invalid += 1
                final_valid = False
    if not rows:
        raise ValueError("Metrics CSV contains no complete finite rows.")
    return fields, rows, invalid, final_valid


def prepare_resume(case, info):
    if info["time"] <= 0:
        raise ValueError("No evolved checkpoint exists; --resume requires a saved time greater than zero.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    directory = case / "postProcessing/deposition"
    marker = directory / "stationary.json"
    archive_marker = False
    if marker.exists():
        data = json.loads(marker.read_text())
        when = float(data["time"])
        if not math.isfinite(when) or when > info["time"] + 1e-10:
            raise ValueError("Stationarity marker is newer than the complete checkpoint; inspect the case before resuming.")
        archive_marker = True
    history = directory / "metrics.csv"
    history_data = None
    if history.exists():
        fields, rows, invalid, _ = metrics(history)
        if invalid or any(t > info["time"] + 1e-10 for t in rows):
            history_data = fields, rows
    # Finish validation before changing either existing file.
    if archive_marker:
        marker.rename(marker.with_name(f"stationary.prior.{stamp}.json"))
    if history_data:
        fields, rows = history_data
        history.rename(history.with_name(f"metrics.prior.{stamp}.csv"))
        with history.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows[t] for t in sorted(rows) if t <= info["time"] + 1e-10)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--exit-code", type=int, default=0)
    parser.add_argument("--phase", default="finished")
    parser.add_argument("--interrupted", action="store_true")
    parser.add_argument("--checkpoint", action="store_true")
    parser.add_argument("--prepare-resume", action="store_true")
    parser.add_argument("--print-time", action="store_true")
    args = parser.parse_args()
    case = args.case.resolve()
    try:
        if args.ranks < 1 or not case.is_dir():
            raise ValueError("A valid case directory and positive rank count are required.")
        if args.checkpoint or args.prepare_resume:
            info = checkpoint(case, args.ranks)
            if args.prepare_resume:
                prepare_resume(case, info)
            print(info["time_name"] if args.print_time else json.dumps(info))
            return 0
        errors, warnings = [], []
        info, latest_row, marker = None, None, None
        try:
            info = checkpoint(case, args.ranks)
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
        try:
            _, rows, invalid, final_valid = metrics(case / "postProcessing/deposition/metrics.csv")
            latest_row = rows[max(rows)]
            if invalid:
                warnings.append(f"Ignored {invalid} incomplete/nonfinite metric rows.")
            if not final_valid and args.exit_code == 0 and not args.interrupted:
                errors.append("The final metrics row is incomplete or nonfinite.")
        except (OSError, ValueError) as exc:
            errors.append(str(exc))
        marker_path = case / "postProcessing/deposition/stationary.json"
        if marker_path.exists():
            try:
                marker = json.loads(marker_path.read_text())
            except (OSError, ValueError) as exc:
                errors.append(f"Invalid stationarity marker: {exc}")
        stationary = False
        if marker and info and latest_row:
            try:
                t = float(marker["time"])
                stationary = marker.get("status") == "stationary" and math.isfinite(t) \
                    and math.isclose(t, info["time"], rel_tol=1e-8, abs_tol=1e-10) \
                    and math.isclose(t, latest_row["time"], rel_tol=1e-8, abs_tol=1e-10)
                if not stationary:
                    warnings.append("Stationarity marker does not match the latest saved fields and metric row.")
            except (ValueError, TypeError, KeyError):
                errors.append("Stationarity marker is missing a finite time.")
        status = "interrupted" if args.interrupted else "failed" if args.exit_code or errors else "stationary" if stationary else "unconverged"
        payload = {"status": status, "solver_exit_code": args.exit_code, "phase": args.phase,
                   "checked_at_utc": datetime.now(timezone.utc).isoformat(), "case": str(case),
                   "checkpoint": info, "latest_metrics": latest_row,
                   "stationarity_marker": marker, "errors": errors, "warnings": warnings}
        atomic_json(case / "run_status.json", payload)
        print(json.dumps(payload, indent=2, allow_nan=False))
        if status == "unconverged":
            print("Solver completed, but stationarity was NOT established. Inspect histories or resume with a later end time.", file=sys.stderr)
        return 2 if status == "failed" else 0
    except (OSError, ValueError, KeyError) as exc:
        print(f"Run check failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
