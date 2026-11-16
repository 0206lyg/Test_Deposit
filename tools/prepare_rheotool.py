#!/usr/bin/env python3
"""Generate the native-only constitutive subset used by rheoDepositionFoam.

No OpenFOAM installation or upstream source is modified. Source model physics
comes from fppimenta/rheoTool commit 81f3e8cf6079e84afd04c5d847389c338404b268.
Upstream files retain their GPL notices; see the copied COPYING file.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

COMMIT = "81f3e8cf6079e84afd04c5d847389c338404b268"
PINNED_MANIFEST_SHA256 = '5888b2a934a4a734ccf04fa42e97a213a65287135c3bd8a6d6b90f50ffc82322'


def replace_function(text, signature, body):
    start = text.index(signature)
    brace = text.index("{", start)
    level = 1
    end = brace + 1
    while level:
        if text[end] == "{":
            level += 1
        elif text[end] == "}":
            level -= 1
        end += 1
    return text[:brace] + "{\n" + body + "\n}" + text[end:]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--upstream", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--export-vendor", type=Path, help="Export the hash-verified minimal upstream source, without Git metadata")
    args = ap.parse_args()
    upstream = args.upstream.resolve()
    if not args.output and not args.export_vendor:
        ap.error("--output or --export-vendor is required")
    if (upstream / ".git").exists():
        revision = subprocess.check_output(["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
        if revision != COMMIT:
            raise SystemExit(f"Expected upstream {COMMIT}, found {revision}")
    root = upstream / "of90/src/libs"
    selected = [
        "constitutiveEquations/constitutiveModel",
        "constitutiveEquations/constitutiveEqs/constitutiveEq",
        "constitutiveEquations/constitutiveEqs/PTT/PTTLog",
        "constitutiveEquations/constitutiveEqs/Giesekus/GiesekusLog",
        "constitutiveEquations/constitutiveEqs/multiMode",
        "constitutiveEquations/constitutiveEqs/GNF/Newtonian",
        "constitutiveEquations/constitutiveEqs/GNF/CarreauYasuda",
        "thermo/thermoFunctions/thermoFunctions",
        "thermo/thermoFunctions/Constant",
    ]
    paths = []
    for directory in selected:
        paths.extend(sorted(p for p in (root / directory).iterdir() if p.suffix in (".C", ".H")))
    paths.extend(root / "constitutiveEquations/constitutiveEqs/utils" / p for p in ("jacobi.H", "boilerLog.H"))
    hashes = {str(source.relative_to(root)): hashlib.sha256(source.read_bytes()).hexdigest() for source in paths}
    audit_hashes = {"of90/src/libs/" + p: h for p, h in hashes.items()}
    audit_hashes["LICENSE"] = hashlib.sha256((upstream / "LICENSE").read_bytes()).hexdigest()
    digest = hashlib.sha256(json.dumps(audit_hashes, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if digest != PINNED_MANIFEST_SHA256:
        raise SystemExit("The selected upstream source bytes differ from the reviewed, pinned version.")
    if args.export_vendor:
        vendor = args.export_vendor.resolve()
        if vendor == upstream or upstream in vendor.parents:
            raise SystemExit("Vendor export must be separate from the upstream checkout.")
        for relative in audit_hashes:
            dest = vendor / relative
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(upstream / relative, dest)
        (vendor / "upstream.source.json").write_text(json.dumps({"upstream": "https://github.com/fppimenta/rheoTool", "commit": COMMIT, "source_sha256": audit_hashes, "manifest_sha256": digest}, indent=2) + "\n")
        print(f"Exported {len(audit_hashes)} verified source files to {vendor}")
    if not args.output:
        return
    output = args.output.resolve()
    if upstream == output or upstream in output.parents:
        raise SystemExit("Output must be separate from the upstream checkout.")
    marker = output / ".deposition-subset"
    if output.exists() and any(output.iterdir()) and not marker.is_file():
        raise SystemExit(f"Refusing to overwrite non-generated directory {output}")
    output.mkdir(parents=True, exist_ok=True)
    marker.write_text(COMMIT + "\n")
    for source in paths:
        relative = source.relative_to(root)
        target = output / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    shutil.copyfile(upstream / "LICENSE", output / "COPYING")
    base = output / "constitutiveEquations/constitutiveEqs/constitutiveEq/constitutiveEq.C"
    text = base.read_text()
    text = text.replace('#include "coupledSolver.H"', '// Native-only subset: block-coupled solver deliberately excluded.')
    text = text.replace('#include "blockOperators.H"', '// Native segregated OpenFOAM linear solvers are used.')
    text = replace_function(text, "void constitutiveEq::divTauImplCoupled", '''    FatalErrorInFunction
        << "Block-coupled stress/momentum solution is not supported by this native-only build."
        << exit(FatalError);''')
    text = replace_function(text, "void constitutiveEq::checkIfCoupledSolver", '''    solveCoupled_ = false;
    const dictionary* coupled = dict.subDictPtr("coupledSolvers");
    if (coupled)
    {
        const dictionary* uptau = coupled->subDictPtr("Uptau");
        if (uptau &&
            (uptau->lookupOrDefault<Switch>("solveCoupledTau", false) ||
             uptau->lookupOrDefault<Switch>("solveCoupledUp", false)))
        {
            FatalErrorInFunction << "This build only supports segregated native OpenFOAM solvers."
                << exit(FatalError);
        }
    }
    (void)tau;''')
    text = replace_function(text, "void constitutiveEq::decomposeGradU", '''    // Fattal--Kupferman decomposition, including its repeated-eigenvalue limit.
    // The upstream denominator regularization alone loses pure-shear startup
    // at A=I. In a degenerate eigenspace B must retain the symmetric block.
    forAll(M, cellI)
    {
        const tensor& values = eigVals[cellI];
        const tensor& m = M[cellI];
        tensor& b = B[cellI];
        tensor& w = omega[cellI];
        b = tensor::zero;
        w = tensor::zero;
        b.xx() = m.xx(); b.yy() = m.yy(); b.zz() = m.zz();
        auto pair = []
        (
            const scalar li, const scalar lj, const scalar mij, const scalar mji,
            scalar& bij, scalar& bji, scalar& wij, scalar& wji
        )
        {
            const scalar scale = max(max(mag(li), mag(lj)), SMALL);
            if (mag(lj - li) <= 1e-8*scale)
            {
                bij = bji = 0.5*(mij + mji);
                wij = 0.5*(mij - mji);
            }
            else
            {
                bij = bji = 0;
                wij = (lj*mij + li*mji)/(lj - li);
            }
            wji = -wij;
        };
        pair(values.xx(), values.yy(), m.xy(), m.yx(), b.xy(), b.yx(), w.xy(), w.yx());
        pair(values.xx(), values.zz(), m.xz(), m.zx(), b.xz(), b.zx(), w.xz(), w.zx());
        pair(values.yy(), values.zz(), m.yz(), m.zy(), b.yz(), b.zy(), w.yz(), w.zy());
    }
    omega = innerP(eigVecs, omega, false);
    B = innerP(eigVecs, B, false);''')
    base.write_text(text)
    # A missing optional thermo dictionary selects Constant upstream but is
    # then dereferenced as nullptr. Pass the parent dictionary in that case.
    thermal = output / "thermo/thermoFunctions/thermoFunctions/newThermoFunction.C"
    text = thermal.read_text()
    old = "cstrIter()(name, mesh, *sDict)"
    if old not in text:
        raise SystemExit("Upstream thermo constructor changed unexpectedly")
    thermal.write_text(text.replace(old, "cstrIter()(name, mesh, sDict ? *sDict : dict)"))
    make = output / "Make"
    make.mkdir(exist_ok=True)
    sources = [str(p.relative_to(root)) for p in paths if p.suffix == ".C"]
    (make / "files").write_text("\n".join(sources) + "\n\nLIB = $(FOAM_USER_LIBBIN)/libdepositionConstitutiveEquations\n")
    (make / "options").write_text('''EXE_INC = \\
    -isystem$(EIGEN_RHEO) \\
    -IconstitutiveEquations/constitutiveEqs/utils \\
    -I$(LIB_SRC)/finiteVolume/lnInclude \\
    -I$(LIB_SRC)/meshTools/lnInclude

LIB_LIBS = \\
    -lfiniteVolume \\
    -lmeshTools
''')
    (output / "PROVENANCE.json").write_text(json.dumps({"upstream": "https://github.com/fppimenta/rheoTool", "commit": COMMIT, "source_subtree": "of90/src/libs", "original_source_sha256": hashes, "patches": ["Remove PETSc/block-coupled base-class dependencies and fail explicitly if requested", "Use the continuous repeated-eigenvalue log-conformation decomposition limit", "Avoid dereferencing a missing optional thermo dictionary"]}, indent=2) + "\n")
    print(f"Prepared {len(sources)} compilation units in {output}; upstream {COMMIT}")


if __name__ == "__main__":
    main()
