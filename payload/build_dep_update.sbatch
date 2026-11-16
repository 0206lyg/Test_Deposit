#!/usr/bin/env bash
#SBATCH --job-name=dep_update_build
#SBATCH --account=lilly-roa3
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=dep_update_build_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
# Updating the add-on does not require reinstalling the plotting environment.
export DEP_BUILD_PLOTS="${DEP_BUILD_PLOTS:-0}"
source tools/load_environment.sh
python3 tests/test_convergence_policy.py
bash build_dep.sbatch
echo 'UPDATE BUILD COMPLETE: convergence stop and PEO diagnosis are installed.'
