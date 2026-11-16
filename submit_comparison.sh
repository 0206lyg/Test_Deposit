#!/usr/bin/env bash
# Submit two independent cases. No job is submitted unless this script is run.
set -Eeuo pipefail
dep_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "$dep_root"
dep_flow=${1:-0.5}; dep_mesh=${2:-coarse}
if (( $# > 2 )); then shift 2; else set --; fi
dep_dependency=""
if (( $# )); then
    [[ $# == 2 && $1 == --after-build && $2 =~ ^[0-9]+$ ]] || {
        echo 'Usage: ./submit_comparison.sh [0.5] [coarse] [--after-build JOB_ID]' >&2; exit 2;
    }
    dep_dependency="--dependency=afterok:$2"
fi
case "$dep_flow" in .3|0.3) dep_flow=0.3;; .5|0.5) dep_flow=0.5;; .7|0.7) dep_flow=0.7;; *) echo 'Flow must be 0.3, 0.5 or 0.7 mL/h.' >&2; exit 2;; esac
case "$dep_mesh" in smoke|coarse|medium|fine) ;; *) echo 'Invalid mesh level.' >&2; exit 2;; esac
command -v sbatch >/dev/null || { echo 'sbatch is unavailable; run this on Gautschi.' >&2; exit 2; }
for dep_fluid in gly75 peo; do
    dep_case="$dep_root/runs/${dep_fluid}_Q${dep_flow/./p}_${dep_mesh}"
    [[ ! -e "$dep_case" ]] || { echo "Existing case: $dep_case. Submit its --resume command explicitly." >&2; exit 2; }
done
dep_args=(--parsable)
[[ -n "$dep_dependency" ]] && dep_args+=("$dep_dependency")
for dep_fluid in gly75 peo; do
    dep_job=$(sbatch "${dep_args[@]}" --job-name="dep_${dep_fluid}" "$dep_root/run_dep.sbatch" "$dep_fluid" "$dep_flow" "$dep_mesh")
    echo "Submitted $dep_fluid: $dep_job"
done
