#!/usr/bin/env python3
"""Regression of the compiled ePTT library against exact homogeneous shear.

Run after loading tools/load_environment.sh and building the project library.
Requires Python's standard library, blockMesh and wmake; no NumPy/SciPy.
Outputs a small audit log and report.json under --workdir.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess


def header(kind, name):
    return f"FoamFile {{version 2.0; format ascii; class {kind}; object {name};}}\n"


def fixture(path, modes, gamma, dt, duration, inner, rotation=None):
    for d in ("0", "constant", "system"):
        (path / d).mkdir(parents=True, exist_ok=True)
    (path / "system/blockMeshDict").write_text(header("dictionary", "blockMeshDict") + """convertToMeters 1;
vertices ((0 0 0)(1 0 0)(1 1 0)(0 1 0)(0 0 1)(1 0 1)(1 1 1)(0 1 1));
blocks (hex (0 1 2 3 4 5 6 7) (3 1 1) simpleGrading (1 1 1));
edges ();
boundary (walls {type wall; faces ((0 4 7 3)(1 2 6 5)(0 1 5 4)(3 7 6 2)(0 3 2 1)(4 5 6 7));});
mergePatchPairs ();
""")
    (path / "system/controlDict").write_text(header("dictionary", "controlDict") + f"""application rheoMaterialPointTest;
startFrom startTime; startTime 0; stopAt endTime; endTime {duration:.17g}; deltaT {dt:.17g};
writeControl timeStep; writeInterval 1000000; purgeWrite 0; writeFormat ascii;
writePrecision 12; writeCompression off; timeFormat general; timePrecision 12;
runTimeModifiable false; shearRate {gamma:.17g}; innerCorrections {inner};
""")
    if rotation is not None:
        with (path / "system/controlDict").open("a") as f:
            f.write("rotation (" + " ".join(f"{v:.17g}" for row in rotation for v in row) + ");\n")
    (path / "system/fvSchemes").write_text(header("dictionary", "fvSchemes") + """ddtSchemes {default Euler;}
gradSchemes {default Gauss linear;}
divSchemes {default Gauss upwind;}
laplacianSchemes {default Gauss linear corrected;}
interpolationSchemes {default linear;}
snGradSchemes {default corrected;}
""")
    # Each component matrix has one diagonal entry. Two linear iterations
    # avoid endless relative-residual iterations at roundoff-level zero stress.
    (path / "system/fvSolution").write_text(header("dictionary", "fvSolution") + """solvers {"theta.*" {solver smoothSolver; smoother symGaussSeidel; tolerance 1e-13; relTol 0; maxIter 2;}}
relaxationFactors {equations {"theta.*" 1;}}
""")
    (path / "0/U").write_text(header("volVectorField", "U") + "dimensions [0 1 -1 0 0 0 0]; internalField uniform (0 0 0); boundaryField {walls {type fixedValue; value uniform (0 0 0);}}\n")
    text = header("dictionary", "constitutiveProperties") + "parameters {type multiMode; models (\n"
    for i, mode in enumerate(modes, 1):
        eta, lam, eps = mode
        text += f"""M{i} {{type PTTLog;
rho rho [1 -3 0 0 0 0 0] 1000;
etaS etaS [1 -1 -1 0 0 0 0] {0.00695 if i == 1 else 0};
etaP etaP [1 -1 -1 0 0 0 0] {eta:.17g};
lambda lambda [0 0 1 0 0 0 0] {lam:.17g};
epsilon epsilon [0 0 0 0 0 0 0] {eps:.17g};
zeta zeta [0 0 0 0 0 0 0] 0;
destructionFunctionType exponential; stabilization coupling;
thermoLambda {{type Constant;}} thermoEta {{type Constant;}}
}}\n"""
        for field, dims in ((f"thetaM{i}", "0 0 0 0 0 0 0"), (f"tauM{i}", "1 -1 -2 0 0 0 0")):
            (path / "0" / field).write_text(header("volSymmTensorField", field) + f"dimensions [{dims}]; internalField uniform (0 0 0 0 0 0); boundaryField {{walls {{type zeroGradient;}}}}\n")
    (path / "constant/constitutiveProperties").write_text(text + ");}\n")


def exact(modes, gamma):
    shear, normal = 0.0, 0.0
    for eta, lam, eps in modes:
        target = 2*eps*(lam*gamma)**2
        lo, hi = 0.0, max(1.0, math.log1p(target))
        for _ in range(90):
            q = .5*(lo + hi)
            if q*math.exp(2*q) < target:
                lo = q
            else:
                hi = q
        f = math.exp(.5*(lo+hi))
        shear += eta*gamma/f
        normal += 2*lam*eta*gamma**2/f**2
    return shear, normal


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    project = Path(__file__).resolve().parents[1]
    ap.add_argument("--workdir", type=Path, default=project / "validation/material_point")
    ap.add_argument("--gamma", type=float, default=1000.)
    ap.add_argument("--dt", type=float, default=1e-5)
    ap.add_argument("--duration", type=float, default=.3)
    ap.add_argument("--inner", type=int, default=4)
    ap.add_argument("--relative-tolerance", type=float, default=.001)
    ap.add_argument("--skip-build", action="store_true")
    args = ap.parse_args()
    if min(args.gamma, args.dt, args.duration, args.relative_tolerance) <= 0 or args.inner < 1:
        ap.error("Rates, times, tolerance and iteration count must be positive")
    if os.environ.get("WM_PROJECT_VERSION") != "8":
        raise SystemExit("Load the project's Foundation OpenFOAM-8 environment first.")
    appbin = os.environ.get("FOAM_USER_APPBIN")
    if not appbin:
        raise SystemExit("FOAM_USER_APPBIN is missing; source tools/load_environment.sh")
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if not args.skip_build:
        with (work / "log.build").open("w") as log:
            subprocess.run(["wmake", str(project / "tests/materialPoint")], stdout=log, stderr=subprocess.STDOUT, check=True)
    executable = Path(appbin) / "rheoMaterialPointTest"
    if not executable.is_file():
        raise SystemExit(f"Executable missing: {executable}")
    config = json.loads((project / "config/experiment.json").read_text())
    baseline = [(m["eta_p_Pa_s"], m["lambda_s"], m["epsilon"]) for m in config["fluids"]["peo"]["modes"]]
    identity = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
    rotation = [[-2/3,2/15,11/15],[2/3,-1/3,2/3],[1/3,14/15,2/15]]
    candidates = {"selected_peo": (baseline, identity),
        "strong_nonlinearity_regression": ([(.13588580220713073, .04, .5323201163935225), (.05716419779286927, .002, .6802672592220357)], identity),
        "selected_peo_rotated_3d": (baseline, rotation)}
    results = {}
    for name, (modes, rotation) in candidates.items():
        case = work / name
        fixture(case, modes, args.gamma, args.dt, args.duration, args.inner, rotation)
        with (case / "log.blockMesh").open("w") as log:
            subprocess.run(["blockMesh", "-case", str(case)], stdout=log, stderr=subprocess.STDOUT, check=True)
        proc = subprocess.Popen([str(executable), "-case", str(case)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        final = full_tensor = None
        algebra_passed = False
        with (case / "log.audit").open("w") as log:
            for line in proc.stdout:
                if line.startswith("AUDIT") or "ERROR" in line or "Warning" in line:
                    log.write(line)
                if line.startswith("AUDIT_FINAL "):
                    final = tuple(map(float, line.split()[1:]))
                if line.startswith("AUDIT_TENSOR "):
                    full_tensor = tuple(map(float, line.split()[1:]))
                if line.startswith("AUDIT_ALGEBRA_PASS "):
                    algebra_passed = True
        if proc.wait() != 0 or final is None or full_tensor is None or not algebra_passed:
            raise SystemExit(f"Material-point executable failed; inspect {case / 'log.audit'}")
        shear, normal = exact(modes, args.gamma)
        reference = [[normal,shear,0.],[shear,0.,0.],[0.,0.,0.]]
        rotated = [[sum(rotation[i][k]*reference[k][l]*rotation[j][l]
                        for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
        expected_tensor = [rotated[i][j] for i,j in ((0,0),(0,1),(0,2),(1,1),(1,2),(2,2))]
        expected = (rotated[0][1],rotated[0][0]-rotated[1][1])
        errors = [abs(a-b)/max(abs(b), 1e-30) for a, b in zip(final[1:], expected)]
        scale=max(abs(v) for v in expected_tensor)
        tensor_errors=[abs(a-b)/max(abs(b),1e-8*scale,1e-30) for a,b in zip(full_tensor,expected_tensor)]
        passed = all(math.isfinite(e) and e <= args.relative_tolerance for e in errors+tensor_errors)
        results[name] = {"passed": passed, "gamma_s-1": args.gamma, "time_s": final[0], "polymer_shear_Pa": final[1], "N1_Pa": final[2], "analytical_polymer_shear_Pa": expected[0], "analytical_N1_Pa": expected[1], "relative_errors": errors, "modes_eta_lambda_epsilon": modes}
        results[name].update({"rotation":rotation,"algebra_checks_passed":algebra_passed,
            "stress_components_order":["xx","xy","xz","yy","yz","zz"],
            "stress_components_Pa":full_tensor,"analytical_stress_components_Pa":expected_tensor,
            "tensor_relative_errors":tensor_errors})
        print(f"{name}: {'PASS' if passed else 'FAIL'}, shear={final[1]:.9g} Pa, N1={final[2]:.9g} Pa, max relative error={max(errors+tensor_errors):.3g}", flush=True)
    (work / "report.json").write_text(json.dumps({"scope": "Actual compiled ePTT constitutive equations in prescribed homogeneous shear; not contraction-flow validation", "dt_s": args.dt, "inner_corrections": args.inner, "relative_tolerance": args.relative_tolerance, "results": results}, indent=2) + "\n")
    if not all(r["passed"] for r in results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
