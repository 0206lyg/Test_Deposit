#!/usr/bin/env bash
#SBATCH --job-name=dep_wall_profiles
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=dep_post_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
dep_py=${DEP_PYTHON:-python3}
if [[ -x .venv/bin/python ]]; then dep_py=.venv/bin/python; fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "$dep_py" tools/steady_workflow.py post "$@"
