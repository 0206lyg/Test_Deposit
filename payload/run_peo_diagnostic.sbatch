#!/usr/bin/env bash
#SBATCH --job-name=dep_peo_diagnose
#SBATCH --account=lilly-roa3
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --signal=B:USR1@120
#SBATCH --output=dep_peo_diagnose_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
dep_flow=${1:-0.5}; dep_mesh=${2:-coarse}
if (( $# >= 2 )); then shift 2; else set --; fi
case "$dep_flow" in .3|0.3) dep_flow=0.3;; .5|0.5) dep_flow=0.5;; .7|0.7) dep_flow=0.7;; *) echo 'Use Q=0.3, 0.5 or 0.7.' >&2; exit 2;; esac
dep_case="${DEP_RESULTS_ROOT:-/scratch/gautschi/lee6456/results_deposition}/peo_Q${dep_flow/./p}_${dep_mesh}_diagnostic"
while (( $# )); do
    case "$1" in
        --case) (( $# >= 2 )) || exit 2; dep_case=$2; shift 2;;
        *) echo 'Usage: sbatch run_peo_diagnostic.sbatch [Q] [mesh] [--case PATH]' >&2; exit 2;;
    esac
done
source tools/load_environment.sh
"$DEP_PYTHON" tools/record_build.py --check
[[ ! -e "$dep_case" ]] || { echo 'Diagnostic case already exists; use a new --case path.' >&2; exit 2; }
dep_child=''
dep_signal() { [[ -z "$dep_child" ]] || kill -USR1 "$dep_child" 2>/dev/null || true; }
dep_collect() {
    dep_code=$?
    trap - EXIT
    "$DEP_PYTHON" "$DEP_ROOT/tools/collect_peo_diagnostics.py" "$dep_case" || true
    exit "$dep_code"
}
trap dep_collect EXIT
trap dep_signal USR1 TERM INT
# Keep the original rheology, mesh, time step and schemes. Only diagnostic
# reporting and a short end-time differ from the failed run. Keep the original
# 0.1 s write interval so output-time alignment does not change the trajectory.
bash run_dep.sbatch peo "$dep_flow" "$dep_mesh" --case "$dep_case" \
    --end-time "${DEP_DIAGNOSTIC_END_TIME:-0.003}" --write-interval 0.1 &
dep_child=$!
set +e
while true; do
    wait "$dep_child"
    dep_code=$?
    if (( dep_code > 128 )) && kill -0 "$dep_child" 2>/dev/null; then continue; fi
    break
done
set -e
dep_child=''
exit "$dep_code"
