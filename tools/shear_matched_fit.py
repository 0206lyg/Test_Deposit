#!/usr/bin/env python3
"""Check the PEO steady-shear fit or generate a shear-equivalent second mode."""
import argparse
import copy
import json
import math
from pathlib import Path

from planar_inlet import shear_factor


def shear_viscosity(fluid, rate):
    return fluid["eta_s_Pa_s"] + sum(
        m["eta_p_Pa_s"] / shear_factor(m, rate) for m in fluid["modes"]
    )


def main():
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", type=Path, default=root / "config/experiment.json")
    ap.add_argument("--lambda2-ms", type=float)
    ap.add_argument("--write-config", type=Path, help="new configuration for the steady launcher")
    args = ap.parse_args()
    config = json.loads(args.base.read_text())
    peo = config["fluids"]["peo"]
    if len(peo["modes"]) != 2 or peo["zeta"] != 0 or peo["destruction_function"] != "exponential":
        ap.error("Requires two-mode exponential PTT with zeta = 0")
    if (args.lambda2_ms is None) != (args.write_config is None):
        ap.error("Provide --lambda2-ms and --write-config together")
    candidate = copy.deepcopy(config)
    fluid = candidate["fluids"]["peo"]
    if args.lambda2_ms is not None:
        lam = args.lambda2_ms / 1000
        if not 0 < lam <= peo["modes"][0]["lambda_s"]:
            ap.error("The second relaxation time must be positive and at most the longest (CaBER) time")
        mode = fluid["modes"][1]
        mode["epsilon"] *= (mode["lambda_s"] / lam) ** 2
        mode["lambda_s"] = lam
        fluid["parameter_provenance"] += (
            f" Sensitivity: lambda2 = {lam:g} s; epsilon2 adjusted to preserve "
            "epsilon2*lambda2^2 and the baseline steady shear viscosity."
        )

    cross = config["fluids"]["cross"]
    print("shear_rate_1_s  Cross_Pa_s  baseline_ePTT_Pa_s  candidate_ePTT_Pa_s")
    for rate in (0, 10, 100, 823, 1372, 1920, 10000):
        eta_cross = cross["eta_inf_Pa_s"] + (cross["eta_zero_Pa_s"] - cross["eta_inf_Pa_s"]) * (
            1 + (cross["k_s"] * rate) ** cross["a"]
        ) ** ((cross["n"] - 1) / cross["a"])
        reference = shear_viscosity(peo, rate)
        eta = shear_viscosity(fluid, rate)
        if not math.isclose(eta, reference, rel_tol=1e-12):
            raise ValueError(f"Steady-shear equivalence failed at {rate}/s")
        print(f"{rate:14g}  {eta_cross:10.7f}  {reference:18.7f}  {eta:19.7f}")
    if args.write_config:
        with args.write_config.open("x") as out:
            json.dump(candidate, out, indent=2, allow_nan=False)
            out.write("\n")
        print(f"Wrote {args.write_config}")


if __name__ == "__main__":
    main()
