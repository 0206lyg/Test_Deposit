#!/usr/bin/env bash
#SBATCH --job-name=dep_flow
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=24:00:00
#SBATCH --signal=B:USR1@120
#SBATCH --output=dep_flow_%j.out
set -Eeuo pipefail

usage() {
    echo "Usage: sbatch run_dep.sbatch {gly75|peo|cross} {0.3|0.5|0.7} {smoke|coarse|medium|fine} [--resume] [--end-time SECONDS] [--write-interval SECONDS] [--case PATH]" >&2
    echo 'New-case defaults: smoke end-time 0.0001 s; other meshes 4 s; write interval min(0.1 s, end-time/2). Resume preserves its write interval unless explicitly supplied.' >&2
    exit "${1:-2}"
}
[[ $# == 1 && $1 == --help ]] && usage 0
(( $# >= 3 )) || usage
dep_fluid=$1; dep_flow=$2; dep_mesh=$3; shift 3
case "$dep_fluid" in gly75|peo|cross) ;; *) usage ;; esac
case "$dep_flow" in .3|0.3) dep_flow=0.3;; .5|0.5) dep_flow=0.5;; .7|0.7) dep_flow=0.7;; *) usage;; esac
case "$dep_mesh" in smoke|coarse|medium|fine) ;; *) usage ;; esac
dep_resume=0; dep_end_time=""; dep_write_interval=""; dep_case_arg=""
while (( $# )); do
    case "$1" in
        --resume) dep_resume=1; shift ;;
        --end-time) (( $# >= 2 )) || usage; dep_end_time=$2; shift 2 ;;
        --write-interval) (( $# >= 2 )) || usage; dep_write_interval=$2; shift 2 ;;
        --case) (( $# >= 2 )) || usage; dep_case_arg=$2; shift 2 ;;
        *) usage ;;
    esac
done
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
[[ -f tools/load_environment.sh ]] || { echo 'Submit from /home/lee6456/Deposition.' >&2; exit 2; }
source tools/load_environment.sh
dep_ranks=${SLURM_NTASKS:-${DEP_RANKS:-16}}
[[ "$dep_ranks" =~ ^[1-9][0-9]*$ ]] || { echo 'MPI rank count must be a positive integer.' >&2; exit 2; }
if [[ -n "$dep_case_arg" ]]; then
    DEP_CASE=$("$DEP_PYTHON" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$dep_case_arg")
else
    DEP_CASE="$DEP_ROOT/runs/${dep_fluid}_Q${dep_flow/./p}_${dep_mesh}"
fi
dep_run_id=${SLURM_JOB_ID:-local_$(date -u +%Y%m%dT%H%M%S)}
exec > >(tee "$DEP_ROOT/logs/run_${dep_run_id}.log") 2>&1
echo "Case: $DEP_CASE; fluid: $dep_fluid; Q: $dep_flow mL/h; mesh: $dep_mesh; ranks: $dep_ranks"
"$DEP_PYTHON" tools/record_build.py --check
for dep_command in blockMesh checkMesh decomposePar reconstructPar foamDictionary "$DEP_MPI_LAUNCHER"; do
    command -v "$dep_command" >/dev/null || { echo "Missing command: $dep_command" >&2; exit 2; }
done
if [[ -n "$dep_end_time" ]]; then
    "$DEP_PYTHON" -c 'import math,sys; t=float(sys.argv[1]); assert math.isfinite(t) and t>0, "end time must be positive and finite"' "$dep_end_time"
fi
if [[ -n "$dep_write_interval" ]]; then
    "$DEP_PYTHON" -c 'import math,sys; t=float(sys.argv[1]); assert math.isfinite(t) and t>0, "write interval must be positive and finite"' "$dep_write_interval"
fi

dep_phase=prepare; dep_record_status=0; dep_stop_requested=0; dep_solver_pid=""
dep_finish() {
    local dep_code=$?
    trap - EXIT
    if (( dep_record_status )); then
        local dep_check_args=(--case "$DEP_CASE" --ranks "$dep_ranks" --exit-code "$dep_code" --phase "$dep_phase")
        (( dep_stop_requested )) && dep_check_args+=(--interrupted)
        local dep_check_code=0
        "$DEP_PYTHON" "$DEP_ROOT/tools/check_run.py" "${dep_check_args[@]}" || dep_check_code=$?
        if (( dep_code == 0 && dep_check_code != 0 )); then dep_code=$dep_check_code; fi
    fi
    exit "$dep_code"
}
dep_stop() {
    dep_stop_requested=1
    echo 'Scheduler stop requested: asking OpenFOAM to write the current state and exit.' >&2
    if [[ -n "$dep_solver_pid" ]] && kill -0 "$dep_solver_pid" 2>/dev/null; then
        foamDictionary "$DEP_CASE/system/controlDict" -entry stopAt -set writeNow || true
    else
        exit 75
    fi
}
trap dep_finish EXIT
trap dep_stop USR1 TERM INT

if (( dep_resume )); then
    [[ -d "$DEP_CASE" ]] || { echo "Cannot resume a missing case: $DEP_CASE" >&2; exit 2; }
    "$DEP_PYTHON" - "$DEP_CASE" "$dep_fluid" "$dep_flow" "$dep_mesh" "$dep_ranks" <<'PY'
import json, math, sys
from pathlib import Path
p=Path(sys.argv[1])/'metadata.json'
if not p.is_file(): raise SystemExit('Resume requires the generated metadata.json; case was not changed.')
m=json.loads(p.read_text())
if m.get('fluid_key') != sys.argv[2] or not math.isclose(float(m.get('Q_mL_h',-1)),float(sys.argv[3])) or m.get('mesh_level') != sys.argv[4] or int(m.get('ranks',-1)) != int(sys.argv[5]):
    raise SystemExit('Resume fluid, flow, mesh or rank count differs from metadata.json; case was not changed.')
PY
    dep_checkpoint=$("$DEP_PYTHON" tools/check_run.py --case "$DEP_CASE" --ranks "$dep_ranks" --checkpoint --print-time)
    dep_target_end=${dep_end_time:-$(foamDictionary "$DEP_CASE/system/controlDict" -entry endTime -value)}
    "$DEP_PYTHON" -c 'import sys; c,e=map(float,sys.argv[1:]); assert c>0, "No evolved checkpoint exists"; assert e>c, "Use --end-time greater than the saved checkpoint time"' "$dep_checkpoint" "$dep_target_end"
    # Archive stale stationarity and post-checkpoint history only after validation.
    "$DEP_PYTHON" tools/check_run.py --case "$DEP_CASE" --ranks "$dep_ranks" --prepare-resume >/dev/null
    dep_record_status=1
    foamDictionary "$DEP_CASE/system/controlDict" -entry startFrom -set startTime
    foamDictionary "$DEP_CASE/system/controlDict" -entry startTime -set "$dep_checkpoint"
    foamDictionary "$DEP_CASE/system/controlDict" -entry endTime -set "$dep_target_end"
    foamDictionary "$DEP_CASE/system/controlDict" -entry stopAt -set endTime
    if [[ -n "$dep_write_interval" ]]; then
        foamDictionary "$DEP_CASE/system/controlDict" -entry writeInterval -set "$dep_write_interval"
    fi
    echo "Resuming complete checkpoint $dep_checkpoint s; end-time cap $dep_target_end s."
else
    [[ ! -e "$DEP_CASE" ]] || { echo "Case already exists: $DEP_CASE. Use --resume or a new --case path; nothing was overwritten." >&2; exit 2; }
    if [[ -z "$dep_end_time" ]]; then
        if [[ "$dep_mesh" == smoke ]]; then dep_end_time=0.0001; else dep_end_time=4; fi
    fi
    if [[ -z "$dep_write_interval" ]]; then
        dep_write_interval=$("$DEP_PYTHON" -c 'import sys; print(min(.1,float(sys.argv[1])/2))' "$dep_end_time")
    fi
    "$DEP_PYTHON" tools/make_case.py --fluid "$dep_fluid" --flow "$dep_flow" --mesh "$dep_mesh" --case "$DEP_CASE" --ranks "$dep_ranks" --end-time "$dep_end_time" --write-interval "$dep_write_interval"
    dep_record_status=1
    dep_phase=mesh
    blockMesh -case "$DEP_CASE" > "$DEP_CASE/log.blockMesh.${dep_run_id}" 2>&1
    checkMesh -case "$DEP_CASE" > "$DEP_CASE/log.checkMesh.${dep_run_id}" 2>&1
    if ! grep -Fq 'Mesh OK.' "$DEP_CASE/log.checkMesh.${dep_run_id}"; then
        echo "Mesh validation did not report Mesh OK. Inspect $DEP_CASE/log.checkMesh.${dep_run_id}; the solver was not started." >&2
        exit 3
    fi
fi
dep_phase=decompose
if (( dep_ranks > 1 )); then
    if [[ ! -d "$DEP_CASE/processor0" ]]; then
        dep_decompose_args=(-case "$DEP_CASE")
        if (( dep_resume )); then dep_decompose_args+=(-time "$dep_checkpoint"); fi
        decomposePar "${dep_decompose_args[@]}" > "$DEP_CASE/log.decomposePar.${dep_run_id}" 2>&1
    fi
fi
dep_phase=solve
if (( dep_ranks > 1 )); then
    "$DEP_MPI_LAUNCHER" -np "$dep_ranks" "$FOAM_USER_APPBIN/rheoDepositionFoam" -case "$DEP_CASE" -parallel > "$DEP_CASE/log.solver.${dep_run_id}" 2>&1 &
else
    "$FOAM_USER_APPBIN/rheoDepositionFoam" -case "$DEP_CASE" > "$DEP_CASE/log.solver.${dep_run_id}" 2>&1 &
fi
dep_solver_pid=$!
set +e
while true; do
    wait "$dep_solver_pid"
    dep_solver_code=$?
    if (( dep_solver_code > 128 && dep_stop_requested )) && kill -0 "$dep_solver_pid" 2>/dev/null; then
        continue
    fi
    break
done
set -e
dep_solver_pid=""
(( dep_solver_code == 0 )) || exit "$dep_solver_code"
dep_phase=reconstruct
if (( dep_ranks > 1 )); then
    # Explicitly select the latest complete common checkpoint, never a partial write.
    dep_checkpoint=$("$DEP_PYTHON" tools/check_run.py --case "$DEP_CASE" --ranks "$dep_ranks" --checkpoint --print-time)
    reconstructPar -case "$DEP_CASE" -time "$dep_checkpoint" > "$DEP_CASE/log.reconstructPar.${dep_run_id}" 2>&1
fi
touch "$DEP_CASE/case.foam"
dep_phase=finished
echo "Saved fields: $DEP_CASE/case.foam"
echo "Metrics: $DEP_CASE/postProcessing/deposition/metrics.csv"
if (( dep_stop_requested )); then exit 75; fi
exit 0
