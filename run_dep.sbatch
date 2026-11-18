#!/usr/bin/env bash
#SBATCH --job-name=dep_steady
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --signal=B:USR1@120
#SBATCH --output=dep_steady_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
if [[ ${1:-} == --help ]]; then exec python3 tools/steady_workflow.py --help; fi
source tools/load_environment.sh
"$DEP_PYTHON" tools/record_build.py --check
exec "$DEP_PYTHON" -u tools/steady_workflow.py "$@"
