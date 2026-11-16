#!/usr/bin/env bash
#SBATCH --job-name=dep_post
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=dep_post_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
dep_root=$(pwd -P)
dep_flow=${1:-0.5}
dep_mesh=${2:-coarse}
case "$dep_flow" in .3|0.3) dep_flow=0.3;; .5|0.5) dep_flow=0.5;; .7|0.7) dep_flow=0.7;; *) echo 'Use Q=0.3, 0.5 or 0.7.' >&2; exit 2;; esac
case "$dep_mesh" in smoke|coarse|medium|fine) ;; *) echo 'Invalid mesh.' >&2; exit 2;; esac
dep_a="$dep_root/runs/gly75_Q${dep_flow/./p}_${dep_mesh}"
dep_b="$dep_root/runs/peo_Q${dep_flow/./p}_${dep_mesh}"
dep_out="$dep_root/comparison/Q${dep_flow/./p}_${dep_mesh}"
dep_py="$dep_root/.venv/bin/python"
[[ -x "$dep_py" ]] || dep_py=${DEP_PYTHON:-python3}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$dep_py" tools/compare_results.py "$dep_a" "$dep_b" --output "$dep_out" \
    --metrics deltaP_Pa,wallShearMean_Pa,downstreamWallShearMean_Pa,approachFlux_m3_s,massImbalance,maxN1_Pa
dep_pv=${DEP_PVPYTHON:-pvpython}
if command -v "$dep_pv" >/dev/null 2>&1; then
    # A state file needs no X server; rendering can be done in desktop ParaView.
    "$dep_pv" --no-mpi tools/paraview_compare.py "$dep_a" "$dep_b" --output "$dep_out/paraview" --state-only
else
    echo 'ParaView pvpython is not loaded. Metric plots are complete; open each case.foam in desktop ParaView, or run tools/paraview_compare.py with your pvpython.'
fi
echo "Comparison outputs: $dep_out"
