# Project-local settings. This file does not modify your shell startup or OpenFOAM.
# These module defaults match the prior Gautschi OpenFOAM-8 workflow.
# Override before sbatch with: export DEP_MODULES='... ...'
DEP_OF_ROOT="${DEP_OF_ROOT:-/home/lee6456/OpenFOAM/OpenFOAM-8}"
DEP_MODULES="${DEP_MODULES-gcc/11.4.1 openmpi/4.1.6}"
DEP_PYTHON="${DEP_PYTHON:-python3}"
DEP_MPI_LAUNCHER="${DEP_MPI_LAUNCHER:-mpirun}"
DEP_BUILD_JOBS="${DEP_BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-8}}"
