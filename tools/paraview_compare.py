#!/usr/bin/env pvpython
"""Create paired contraction views and an editable ParaView state (5.10+ API).

Run with ParaView's Python, independently of the OpenFOAM installation:
  pvpython tools/paraview_compare.py runs/gly75_Q0p5 runs/peo_Q0p5 --output comparison/paraview
On a headless machine use a ParaView build with EGL/OSMesa or a working display;
`pvpython --force-offscreen-rendering ...` alone does not supply that backend.

Input: serial/reconstructed OpenFOAM cases, U [m/s], p, and optional N1 [Pa]
and wallShearStressTotal [Pa]. Use reconstructPar before this script for MPI
runs. Defaults to each case's latest field time with matching metrics; only
matching stationary.json markers permit the STEADY label. --time requests a
matched transient snapshot. Each reader is frozen at its selected time.
Coordinates: x along flow, y across width, z through height, contraction x=0.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import sys


def field_exists(directory, name):
    return any((directory / (name + suffix)).is_file() and (directory / (name + suffix)).stat().st_size > 0
               for suffix in ("", ".gz"))


def available_times(case, metrics_file):
    times = {}
    for directory in case.iterdir():
        try:
            value = float(directory.name)
        except ValueError:
            continue
        if directory.is_dir() and math.isfinite(value) and field_exists(directory, "U") and field_exists(directory, "p"):
            times[value] = directory
    # Only metric rows fully written by the postprocessor count as completion
    # markers. Do not select a field directory still being written by a solver.
    metrics = case / metrics_file
    if metrics.is_file():
        complete = []
        with metrics.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(line for line in stream if line.strip() and not line.startswith("#")):
                try:
                    if None in row or any(value in (None, "") for value in row.values()):
                        continue
                    values = [float(value) for value in row.values()]
                    if all(math.isfinite(value) for value in values):
                        complete.append(float(row["time"]))
                except (ValueError, KeyError, TypeError):
                    continue
        times = {t: d for t, d in times.items() if any(math.isclose(t, v, rel_tol=1e-9, abs_tol=1e-12) for v in complete)}
    else:
        print(f"Warning: {metrics} missing; U/p file presence cannot certify a finished write. Use stopped/completed cases.", file=sys.stderr)
    return times


def common_time(a, b, requested):
    pairs = [(x, y) for x in a for y in b if math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-12)]
    if requested is not None:
        pairs = [p for p in pairs if math.isclose(p[0], requested, rel_tol=1e-9, abs_tol=1e-12)]
    if not pairs:
        raise ValueError("No common complete serial/reconstructed field time; check metrics.csv and reconstructPar output.")
    return max(pairs)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("newtonian", type=Path)
    parser.add_argument("viscoelastic", type=Path)
    parser.add_argument("--output", type=Path, default=Path("comparison/paraview"))
    parser.add_argument("--metrics-file", default="postProcessing/deposition/metrics.csv")
    parser.add_argument("--stationary-file", default="postProcessing/deposition/stationary.json")
    parser.add_argument("--time", type=float, help="common transient time; default: each case's latest complete time")
    parser.add_argument("--slice-z", type=float, help="slice height in metres; default: middle of mesh bounds")
    parser.add_argument("--half-window", type=float, default=180e-6, help="half-width of contraction image along x [m]")
    parser.add_argument("--wall-field", default="wallShearStressTotal")
    parser.add_argument("--normal-field", default="N1", help="polymeric tau_xx-tau_yy, in Pa")
    parser.add_argument("--wall-regions", help="comma-separated reader region names; default: patches named wall/step/top/bottom/front/back/side")
    parser.add_argument("--wall-vmax", type=float, help="optional common upper color limit [Pa], recorded in metadata")
    parser.add_argument("--state-only", action="store_true", help="build .pvsm without rendering PNGs")
    parser.add_argument("--no-streamlines", action="store_true")
    args = parser.parse_args()
    if args.half_window <= 0 or (args.wall_vmax is not None and args.wall_vmax <= 0):
        parser.error("--half-window and --wall-vmax must be positive")
    return parser, args


def main():
    parser, args = parse_args()
    try:
        import paraview.simple as pv
    except ImportError:
        parser.error("Run this script with ParaView's pvpython, not ordinary Python. For the GUI, open each case's .foam file directly.")
    cases = [args.newtonian.resolve(), args.viscoelastic.resolve()]
    try:
        time_maps = [available_times(case, args.metrics_file) for case in cases]
        if not all(time_maps):
            raise ValueError("No saved U/p fields with complete metrics in one or both cases; reconstruct final outputs first.")
        times = common_time(*time_maps, args.time) if args.time is not None else tuple(max(values) for values in time_maps)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    pv._DisableFirstRenderCameraReset()
    readers, wall_readers, slices, slice_z = [], [], [], []
    warnings = []
    labels = ["Gly75", "PEO ePTT"]
    markers = []
    for case, time_value in zip(cases, times):
        path = case / args.stationary_file
        try:
            marker = json.loads(path.read_text(encoding="utf-8"))
            valid = (marker.get("status") == "stationary"
                     and math.isclose(float(marker["time"]), time_value, rel_tol=1e-9, abs_tol=1e-12))
            markers.append({"path": str(path), "matches_selected_time": valid, "marker": marker})
        except (OSError, ValueError, TypeError, KeyError) as exc:
            markers.append({"path": str(path), "matches_selected_time": False, "reason": str(exc)})
    status = "STEADY" if all(marker["matches_selected_time"] for marker in markers) else "UNCONVERGED PREVIEW"

    def set_if(proxy, name, value):
        if proxy.GetProperty(name) is not None:
            setattr(proxy, name, value)

    def reader_for(case, title):
        foam = case / "open.foam"
        foam.touch(exist_ok=True)
        reader = pv.OpenFOAMReader(registrationName=title, FileName=str(foam))
        set_if(reader, "SkipZeroTime", 0)
        set_if(reader, "Adddimensionalunitstoarraynames", 0)
        set_if(reader, "Createcelltopointfiltereddata", 1)
        reader.UpdatePipelineInformation()
        return reader

    def freeze(reader, time_value, title):
        frozen = pv.ForceTime(registrationName=title + " fixed time", Input=reader)
        frozen.ForcedTime = time_value
        frozen.IgnorePipelineTime = 1
        frozen.UpdatePipeline()
        return frozen

    wall_regions_by_case = []
    for case, title, time_value in zip(cases, labels, times):
        reader = reader_for(case, title + " internal")
        regions = list(reader.MeshRegions.Available)
        interior = [r for r in regions if r == "internalMesh" or r.endswith("/internalMesh")]
        if not interior:
            raise RuntimeError(f"{case}: no internalMesh; available: {regions}")
        reader.MeshRegions = interior
        fields = list(reader.CellArrays.Available)
        reader.CellArrays = [f for f in ("U", "p", args.normal_field, args.wall_field) if f in fields]
        reader = freeze(reader, time_value, title + " internal")
        if reader.GetDataInformation().GetNumberOfCells() == 0:
            raise RuntimeError(f"{case}: empty internal mesh at {time_value}")
        bounds = reader.GetDataInformation().GetBounds()
        z = args.slice_z if args.slice_z is not None else .5 * (bounds[4] + bounds[5])
        if not bounds[4] < z < bounds[5]:
            raise ValueError(f"slice z={z} must be strictly inside {case}'s z bounds {bounds[4:6]}")
        plane = pv.Slice(registrationName=title + " midheight", Input=reader)
        plane.SliceType = "Plane"
        plane.SliceType.Origin = [0, 0, z]
        plane.SliceType.Normal = [0, 0, 1]
        plane.UpdatePipeline(time_value)
        readers.append(reader)
        slices.append(plane)
        slice_z.append(z)
        if args.wall_regions:
            requested_regions = [r.strip() for r in args.wall_regions.split(",")]
            selected = [r for r in regions if r in requested_regions or r.split("/")[-1] in requested_regions]
        else:
            selected = [r for r in regions if r not in interior and not r.startswith("group/")
                        and re.search(r"wall|step|top|bottom|front|back|side", r, re.I)]
        wall_regions_by_case.append(selected)
        if selected and args.wall_field in fields and field_exists(time_maps[len(readers)-1][time_value], args.wall_field):
            walls = reader_for(case, title + " walls")
            walls.MeshRegions = selected
            walls.CellArrays = [args.wall_field]
            walls = freeze(walls, time_value, title + " walls")
            # Truncate the long inlet/outlet before an oblique wall view: a
            # camera zoom alone would let remote walls obscure the contraction.
            for edge, normal in ((args.half_window, [1, 0, 0]), (-args.half_window, [-1, 0, 0])):
                clipped = pv.Clip(registrationName=title + f" walls x={edge:g}", Input=walls)
                clipped.ClipType = "Plane"
                clipped.ClipType.Origin = [edge, 0, 0]
                clipped.ClipType.Normal = normal
                clipped.Invert = 1
                clipped.UpdatePipeline(time_value)
                walls = clipped
            wall_readers.append(walls)
        else:
            wall_readers.append(None)
            warnings.append(f"{title}: wall field/regions unavailable; wall pair skipped. Available regions: {regions}")

    animation = pv.GetAnimationScene()
    animation.UpdateAnimationUsingDataTimeSteps()
    animation.AnimationTime = times[0]
    pv.GetTimeKeeper().Time = times[0]
    arrays = []

    def information(source, name):
        for association in ("POINTS", "CELLS"):
            collection = source.PointData if association == "POINTS" else source.CellData
            if name in collection.keys():
                return association, collection[name]
        raise ValueError(f"Missing array {name} in {source}")

    def pair_layout(name, sources, field, unit, wall=False):
        layout = pv.CreateLayout(name=name)
        layout.SplitHorizontal(0, .5)
        infos = [information(source, field) for source in sources]
        ranges = [info.GetRange(-1 if info.GetNumberOfComponents() > 1 else 0) for _, info in infos]
        low, high = min(r[0] for r in ranges), max(r[1] for r in ranges)
        if not all(math.isfinite(value) for value in (low, high)):
            raise ValueError(f"Nonfinite color range for {field}")
        if field in ("U", args.wall_field):
            low = 0.0
        elif field == args.normal_field and low < 0 < high:
            low, high = -max(abs(low), abs(high)), max(abs(low), abs(high))
        if wall and args.wall_vmax is not None:
            high = args.wall_vmax
        if high <= low:
            high = low + max(abs(low) * .01, 1e-12)
        lut = pv.GetColorTransferFunction(field)
        lut.RescaleTransferFunction(low, high)
        set_if(lut, "AutomaticRescaleRangeMode", "Never")
        views = []
        for i, (source, title, (association, info)) in enumerate(zip(sources, labels, infos)):
            view = pv.CreateView("RenderView")
            layout.AssignView(i + 1, view)
            view.ViewSize = [900, 700]
            view.ViewTime = times[i]
            set_if(view, "UseColorPaletteForBackground", 0)
            view.Background = [1, 1, 1]
            view.OrientationAxesVisibility = 0
            display = pv.Show(source, view)
            display.Representation = "Surface"
            selection = (association, field, "Magnitude") if info.GetNumberOfComponents() > 1 else (association, field)
            pv.ColorBy(display, selection)
            display.LookupTable = lut
            display.SetScalarBarVisibility(view, True)
            bar = pv.GetScalarBar(lut, view)
            bar.Title = f"{field} [{unit}]"
            bar.ComponentTitle = ""
            bar.TitleColor = [0, 0, 0]
            bar.LabelColor = [0, 0, 0]
            center = [0, 0, slice_z[i]]
            view.CameraFocalPoint = center
            view.CameraParallelProjection = 1
            if wall:
                view.CameraPosition = [2 * args.half_window, -3 * args.half_window, center[2] + 2.5 * args.half_window]
                view.CameraViewUp = [0, 0, 1]
            else:
                view.CameraPosition = [0, 0, center[2] + 1]
                view.CameraViewUp = [0, 1, 0]
            view.CameraParallelScale = args.half_window * 700 / 900
            annotation = pv.Text(registrationName=name + " " + title + " label")
            annotation.Text = status + "\n" + f"{title} | t = {times[i]:.6g} s" + ("" if wall else f" | z = {slice_z[i]*1e6:.3g} µm")
            text_display = pv.Show(annotation, view)
            text_display.Color = [0, 0, 0] if status == "STEADY" else [.7, .1, .1]
            text_display.FontSize = 14
            views.append(view)
        arrays.append({"view": name, "field": field, "common_color_range": [low, high], "units": unit})
        return layout, views

    velocity_layout, velocity_views = pair_layout("Velocity and streamlines", slices, "U", "m/s")
    if not args.no_streamlines:
        for reader, view, title, z, t in zip(readers, velocity_views, labels, slice_z, times):
            points = pv.CellDatatoPointData(registrationName=title + " velocity interpolation", Input=reader)
            points.UpdatePipeline(t)
            # Main-stream seeds plus seeds in both upstream corner regions.
            # These are fluid streamlines at one time, never particle paths.
            for tag, x, lo, hi, count in (("main", -.85 * args.half_window, -85e-6, 85e-6, 30),
                                          ("upper corner", -5e-6, 25e-6, 88e-6, 16),
                                          ("lower corner", -5e-6, -88e-6, -25e-6, 16)):
                seeds = pv.Line(registrationName=title + " " + tag + " seeds")
                seeds.Point1, seeds.Point2, seeds.Resolution = [x, lo, z], [x, hi, z], count
                stream = pv.StreamTracerWithCustomSource(registrationName=title + " " + tag + " streamlines", Input=points, SeedSource=seeds)
                stream.Vectors = ["POINTS", "U"]
                stream.IntegrationDirection = "BOTH"
                stream.MaximumStreamlineLength = 4 * args.half_window
                stream.UpdatePipeline(t)
                representation = pv.Show(stream, view)
                pv.ColorBy(representation, None)
                representation.DiffuseColor = [.12, .12, .12]
                representation.LineWidth = 1.0
    layouts = [("velocity_streamlines.png", velocity_layout)]
    if all(field_exists(time_maps[i][times[i]], args.normal_field) for i in range(2)):
        normal_layout, _ = pair_layout("Polymer normal stress difference", slices, args.normal_field, "Pa")
        layouts.append(("polymer_N1.png", normal_layout))
    else:
        warnings.append("N1 was not written in both selected time folders; N1 comparison skipped (no zero substitution).")
    if all(source is not None for source in wall_readers):
        wall_layout, _ = pair_layout("Total wall tangential traction", wall_readers, args.wall_field, "Pa", wall=True)
        layouts.append(("wall_tangential_traction.png", wall_layout))
    pv.SaveState(str(output / "comparison.pvsm"))
    metadata = {"cases": [str(case) for case in cases], "selected_times_s": list(times), "slice_z_m": slice_z,
                "status": status, "stationary_markers": markers,
                "selection_mode": "matched_time" if args.time is not None else "independent_latest",
                "image_half_window_x_m": args.half_window, "wall_regions": wall_regions_by_case,
                "plots": arrays, "warnings": warnings,
                "notes": ["Streamlines show fluid velocity at the selected time, not finite-particle trajectories.",
                          "N1 is polymeric tau_xx-tau_yy; it is not wall-normal particle force.",
                          "Wall colors use total tangential traction magnitude in Pa, including solvent and polymer contributions.",
                          "Wall surfaces are clipped to the reported x window before rendering; CSV metrics retain their postprocessor-defined domain.",
                          "Independent latest outputs are solver convergence times, not a matched physical-time transient comparison.",
                          "If cases move to another computer, locate their .foam files in the ParaView Load State data-path dialog.",
                          "No particle attachment, detachment or approach flux is inferred from these pictures."]}
    (output / "visualization_metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if not args.state_only:
        for filename, layout in layouts:
            pv.SaveScreenshot(str(output / filename), layout, ImageResolution=[1800, 700])
    for warning in warnings:
        print("Warning: " + warning, file=sys.stderr)
    print(f"{status}: Gly75 t={times[0]:.8g} s; PEO t={times[1]:.8g} s; output: {output}")


if __name__ == "__main__":
    main()
