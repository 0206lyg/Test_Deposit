#!/usr/bin/env bash
#SBATCH --job-name=dep_gly_steady
#SBATCH --account=lilly-roa3
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --signal=B:USR1@120
#SBATCH --output=dep_gly_steady_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
dep_flow=${1:-0.5}; dep_mesh=${2:-coarse}
if (( $# >= 2 )); then shift 2; else set --; fi
case "$dep_flow" in .3|0.3) dep_flow=0.3;; .5|0.5) dep_flow=0.5;; .7|0.7) dep_flow=0.7;; *) echo 'Use Q=0.3, 0.5 or 0.7.' >&2; exit 2;; esac
dep_case="${DEP_RESULTS_ROOT:-/scratch/gautschi/lee6456/results_deposition}/gly75_Q${dep_flow/./p}_${dep_mesh}_converged"
dep_resume=0
while (( $# )); do
    case "$1" in
        --resume) dep_resume=1; shift;;
        --case) (( $# >= 2 )) || exit 2; dep_case=$2; shift 2;;
        *) echo 'Usage: sbatch run_newtonian_converged.sbatch [Q] [mesh] [--resume] [--case PATH]' >&2; exit 2;;
    esac
done
source tools/load_environment.sh
"$DEP_PYTHON" tools/record_build.py --check
dep_args=()
if (( dep_resume )); then
    "$DEP_PYTHON" tools/configure_newtonian_stop.py --case "$dep_case" --ranks "${SLURM_NTASKS:-4}"
    dep_args+=(--resume)
fi
# exec preserves scheduler signals for the original, tested launcher.
exec bash run_dep.sbatch gly75 "$dep_flow" "$dep_mesh" --case "$dep_case" "${dep_args[@]}"
