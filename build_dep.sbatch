#!/usr/bin/env bash
#SBATCH --job-name=dep_build
#SBATCH --account=lilly-nan6
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=dep_build_%j.out
set -Eeuo pipefail
cd "${DEP_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$(dirname "$0")}}"
[[ -f tools/load_environment.sh ]] || { echo 'Submit from /home/lee6456/Deposition.' >&2; exit 2; }
source tools/load_environment.sh
trap 'echo "Build failed at line $LINENO. See logs/build_${SLURM_JOB_ID:-local}.log; the base OpenFOAM was not modified." >&2' ERR
exec > >(tee "$DEP_ROOT/logs/build_${SLURM_JOB_ID:-local}.log") 2>&1
echo "Base: $WM_PROJECT_DIR; version: $WM_PROJECT_VERSION; options: $WM_OPTIONS"
echo "Project binaries: $FOAM_USER_APPBIN"
command -v wmake g++ "$DEP_PYTHON" "$DEP_MPI_LAUNCHER"
g++ --version | head -1
mpicxx --showme:version
blockMesh -help >/dev/null
"$DEP_PYTHON" tools/fetch_dependencies.py
"$DEP_PYTHON" tools/prepare_rheotool.py --upstream vendor/rheoTool --output build/rheoToolSubset
wmake libso "$DEP_ROOT/build/rheoToolSubset"
wmake "$DEP_ROOT/src/rheoDepositionFoam"
rheoDepositionFoam -help >/dev/null
if ldd "$FOAM_USER_APPBIN/rheoDepositionFoam" | grep -q 'not found'; then
    echo 'Unresolved runtime libraries.' >&2; exit 3
fi
"$DEP_PYTHON" tools/verify_package.py
if [[ "${DEP_RUN_MATERIAL_POINT_TESTS:-1}" == 1 ]]; then
    "$DEP_PYTHON" tests/run_material_point.py
fi
# Plotting is optional for the solver. Keep its environment under this project.
if [[ "${DEP_BUILD_PLOTS:-1}" == 1 ]]; then
    if "$DEP_PYTHON" -m venv "$DEP_ROOT/.venv" && "$DEP_ROOT/.venv/bin/python" -m pip install -r requirements-plot.txt; then
        echo 'Plotting environment ready.'
    else
        echo 'Plotting environment could not be installed; CFD and .foam visualization remain available.' >&2
    fi
fi
"$DEP_PYTHON" tools/record_build.py
echo 'BUILD COMPLETE: rheoDepositionFoam is ready in the project install directory.'
echo 'Next: bash submit_comparison.sh 0.5 coarse'
