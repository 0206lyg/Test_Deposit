#!/usr/bin/env python3
"""Compare two rheoDepositionFoam metric histories without reading mesh data.

Usage:
  python tools/compare_results.py runs/gly75_Q0p5 runs/peo_Q0p5 --output comparison

Each case contains postProcessing/deposition/metrics.csv: a header, `time`,
and numerical SI-valued columns. Restart duplicate times use the last complete
row. By default compare each case's latest saved field with matching metrics;
label STEADY only if both stationary.json markers certify those exact times.
Use --time for a matched-time transient comparison.
Requires numpy and matplotlib (not OpenFOAM or ParaView).
"""
import argparse
import csv
import json
import math
from pathlib import Path
import sys


LABELS = {
    "Qin_m3_s": "Inlet volume flow [m³/s]",
    "Qout_m3_s": "Outlet volume flow [m³/s]",
    "massImbalance": "Relative volume-flow imbalance [-]",
    "deltaP_Pa": "Pressure drop [Pa]",
    "wallShearMean_Pa": "Step-face mean tangential traction [Pa]",
    "wallShearMax_Pa": "Step-face maximum tangential traction [Pa]",
    "downstreamWallShearMean_Pa": "Downstream-sidewall mean traction [Pa]",
    "approachFlux_m3_s": "Positive approach volume flux [m³/s]",
    "approachNetFlux_m3_s": "Net approach volume flux [m³/s]",
    "approachOffset_m": "Approach sampling offset [m]",
    "relativeU": "Relative velocity change [-]",
    "relativeTau": "Relative total extra-stress change [-]",
    "relativeWall": "Relative step-face traction change [-]",
    "relativeApproach": "Relative approach-flux change [-]",
    "relativePressure": "Relative pressure-drop change [-]",
    "maxVelocity_m_s": "Maximum speed [m/s]",
    "maxN1_Pa": "Maximum polymer N1 [Pa]",
    "stableDuration_s": "Continuous stable duration [s]",
}
NOTES = [
    "Approach flux is a finite-offset fluid-transport proxy, not particle capture or deposition flux.",
    "Sharp-corner maximum traction can depend strongly on the mesh; compare averages and resolved spatial profiles.",
    "wallShearMean/Max cover the trimmed contraction step faces; downstreamWallShearMean covers a separate downstream-sidewall region.",
    "Independent latest outputs are solver convergence times, not a matched physical-time transient comparison.",
    "STEADY requires matching stationary.json markers; inspect histories and solver convergence as well.",
]


def read_metrics(path):
    """Ignore incomplete/nonfinite trailing writes; retain complete restart rows."""
    rows, skipped = {}, 0
    with path.open(newline="", encoding="utf-8-sig") as stream:
        data = csv.DictReader(line for line in stream if line.strip() and not line.startswith("#"))
        if not data.fieldnames or "time" not in data.fieldnames:
            raise ValueError(f"{path}: CSV header must include 'time'")
        columns = [key for key in data.fieldnames if key != "time"]
        for raw in data:
            try:
                if None in raw or any(raw.get(key) in (None, "") for key in data.fieldnames):
                    raise ValueError("incomplete row")
                row = {key: float(raw[key]) for key in data.fieldnames}
                if not all(math.isfinite(value) for value in row.values()):
                    raise ValueError("nonfinite row")
            except (TypeError, ValueError):
                skipped += 1
                continue
            rows[row["time"]] = row
    if not rows:
        raise ValueError(f"{path}: no complete finite numerical rows")
    return columns, [rows[t] for t in sorted(rows)], skipped


def common_pairs(left, right):
    """Two-pointer join, with tolerance only for decimal serialization roundoff."""
    pairs, i, j = [], 0, 0
    while i < len(left) and j < len(right):
        a, b = left[i]["time"], right[j]["time"]
        if math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12):
            pairs.append((left[i], right[j]))
            i, j = i + 1, j + 1
        elif a < b:
            i += 1
        else:
            j += 1
    return pairs


def saved_rows(case, rows):
    """Match metrics to serial/reconstructed U and p outputs, if a case is given."""
    if not case.is_dir():
        return rows
    times = []
    for directory in case.iterdir():
        try:
            value = float(directory.name)
        except ValueError:
            continue
        if directory.is_dir() and all(any((directory / (field + suffix)).is_file()
                                           and (directory / (field + suffix)).stat().st_size > 0
                                           for suffix in ("", ".gz")) for field in ("U", "p")):
            times.append(value)
    return [row for row in rows if any(math.isclose(row["time"], t, rel_tol=1e-9, abs_tol=1e-12) for t in times)]


def stationary_marker(case, time_value, relative_path):
    path = case / relative_path if case.is_dir() else case.parent / "stationary.json"
    try:
        marker = json.loads(path.read_text(encoding="utf-8"))
        valid = (case.is_dir() and marker.get("status") == "stationary"
                 and math.isclose(float(marker["time"]), time_value, rel_tol=1e-9, abs_tol=1e-12))
        return {"path": str(path.resolve()), "matches_selected_time": valid, "marker": marker}
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {"path": str(path.resolve()), "matches_selected_time": False, "reason": str(exc)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("newtonian", type=Path, help="Gly75 case directory or metrics CSV")
    parser.add_argument("viscoelastic", type=Path, help="PEO case directory or metrics CSV")
    parser.add_argument("--metrics-file", default="postProcessing/deposition/metrics.csv", help="relative CSV path within each case")
    parser.add_argument("--stationary-file", default="postProcessing/deposition/stationary.json")
    parser.add_argument("--output", type=Path, default=Path("comparison"))
    parser.add_argument("--time", type=float, help="specific common output time; default: each case's latest saved output")
    parser.add_argument("--metrics", help="comma-separated column names; default: all common metrics")
    args = parser.parse_args()
    try:
        import numpy as np
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        parser.error(f"{exc}. Install numpy and matplotlib in the analysis Python environment.")

    paths = [path / args.metrics_file if path.is_dir() else path for path in (args.newtonian, args.viscoelastic)]
    try:
        (cols_a, a, skipped_a), (cols_b, b, skipped_b) = [read_metrics(path) for path in paths]
        saved_a, saved_b = saved_rows(args.newtonian, a), saved_rows(args.viscoelastic, b)
        if not saved_a or not saved_b:
            raise ValueError("No saved serial/reconstructed U/p fields with matching metrics in one or both cases. Reconstruct the final output first.")
        if args.time is not None:
            pairs = common_pairs(saved_a, saved_b)
            pairs = [pair for pair in pairs if math.isclose(pair[0]["time"], args.time, rel_tol=1e-9, abs_tol=1e-12)]
            if not pairs:
                raise ValueError("Requested --time is not saved and complete in both cases; no interpolation was performed.")
            selected_a, selected_b = pairs[-1]
        else:
            selected_a, selected_b = saved_a[-1], saved_b[-1]
        metrics = args.metrics.split(",") if args.metrics else [key for key in cols_a if key in cols_b]
        metrics = [key.strip() for key in metrics]
        missing = [key for key in metrics if key not in cols_a or key not in cols_b or key == "time"]
        if missing or not metrics:
            raise ValueError(f"Missing/invalid requested metrics: {missing}; available: {sorted(set(cols_a) & set(cols_b))}")
    except (OSError, ValueError) as exc:
        parser.error(str(exc))

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    markers = [stationary_marker(case, row["time"], args.stationary_file)
               for case, row in zip((args.newtonian, args.viscoelastic), (selected_a, selected_b))]
    status = "STEADY" if all(marker["matches_selected_time"] for marker in markers) else "UNCONVERGED PREVIEW"
    selection_mode = "matched_time" if args.time is not None else "independent_latest"
    summary = []
    for key in metrics:
        x, y = selected_a[key], selected_b[key]
        summary.append({"metric": key, "Gly75_time_s": selected_a["time"], "PEO_time_s": selected_b["time"], "status": status, "Gly75": x, "PEO": y,
                        "PEO_minus_Gly75": y - x, "PEO_over_Gly75": y / x if x != 0 else None})
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    payload = {"selected_times_s": [selected_a["time"], selected_b["time"]], "selection_mode": selection_mode,
               "status": status, "stationary_markers": markers, "sources": [str(p.resolve()) for p in paths],
               "skipped_rows": [skipped_a, skipped_b], "notes": NOTES, "metrics": summary}
    (output / "comparison.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    # Each metric retains its own physical units; unrelated quantities are not
    # combined on one scale or normalized by a potentially zero reference.
    for page, start in enumerate(range(0, len(metrics), 6), 1):
        subset = metrics[start:start + 6]
        fig, axes = plt.subplots(math.ceil(len(subset) / 2), 2, figsize=(12, 3.2 * math.ceil(len(subset) / 2)), squeeze=False)
        for ax, key in zip(axes.flat, subset):
            for rows, selected, label, color in ((a, selected_a, "Gly75", "#355e9a"), (b, selected_b, "PEO ePTT", "#c36b21")):
                times = np.asarray([row["time"] for row in rows])
                values = np.asarray([row[key] for row in rows])
                ax.plot(times, values, label=label, color=color, linewidth=1.7, marker="." if len(times) < 20 else None)
                ax.axvline(selected["time"], color=color, linestyle=":", linewidth=1)
            ax.set(title=LABELS.get(key, key), xlabel="Simulation time [s]")
            ax.ticklabel_format(axis="both", style="sci", scilimits=(-3, 4), useOffset=False)
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
            # Do not turn roundoff in an imposed constant flow rate into a
            # visually large trend with a 1e-22 offset axis.
            values_all = [row[key] for row in a + b]
            low, high = min(values_all), max(values_all)
            scale = max(abs(low), abs(high))
            if scale > 0 and high - low < 1e-7 * scale:
                ax.set_ylim(low - .025 * scale, high + .025 * scale)
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        for ax in list(axes.flat)[len(subset):]:
            ax.set_visible(False)
        fig.suptitle(f"{status} — Gly75 t={selected_a['time']:.6g} s; PEO t={selected_b['time']:.6g} s", fontsize=13,
                     color="black" if status == "STEADY" else "#b22b27")
        fig.tight_layout(rect=(0, 0, 1, .96))
        fig.savefig(output / f"metric_histories_{page:02d}.png", dpi=180)
        plt.close(fig)
    print(f"{status}: {len(metrics)} metrics, Gly75 t={selected_a['time']:.8g} s, PEO t={selected_b['time']:.8g} s; output: {output}")
    if skipped_a or skipped_b:
        print(f"Ignored incomplete/nonfinite rows: Gly75={skipped_a}, PEO={skipped_b}", file=sys.stderr)
    for note in NOTES:
        print(f"Note: {note}")


if __name__ == "__main__":
    main()
