#!/usr/bin/env bash
# Source this file; all installation outputs are redirected inside Deposition.
DEP_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "$DEP_ROOT/config/environment.sh"
if [[ -n "$DEP_MODULES" ]]; then
    if ! type module >/dev/null 2>&1; then
        for dep_init in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh; do
            if [[ -f "$dep_init" ]]; then source "$dep_init"; break; fi
        done
    fi
    if ! type module >/dev/null 2>&1; then
        echo 'Module command unavailable. Initialize Gautschi modules, or set DEP_MODULES="" only for an already configured toolchain.' >&2
        return 2
    fi
    read -r -a dep_modules <<< "$DEP_MODULES"
    module purge
    if ! module load "${dep_modules[@]}"; then
        echo 'The configured compiler/MPI modules are unavailable. Use the SAME module pair used to compile your existing OpenFOAM-8; edit config/environment.sh. OpenFOAM was not changed.' >&2
        return 2
    fi
fi
if [[ ! -r "$DEP_OF_ROOT/etc/bashrc" ]]; then
    echo "Cannot read $DEP_OF_ROOT/etc/bashrc; set DEP_OF_ROOT." >&2
    return 2
fi
# OpenFOAM's own shell files access optional unset variables.
dep_was_nounset=0
dep_was_errexit=0
[[ $- == *u* ]] && dep_was_nounset=1
[[ $- == *e* ]] && dep_was_errexit=1
# OF8's startup probes optional files with nonzero statuses; its shell code is
# not written for caller-side errexit. Validate the resulting environment below.
set +eu
source "$DEP_OF_ROOT/etc/bashrc"
dep_source_status=$?
(( dep_was_nounset )) && set -u
(( dep_was_errexit )) && set -e
if (( dep_source_status != 0 )) || [[ ${WM_PROJECT_VERSION:-} != 8 ]]; then
    echo 'This package requires a working Foundation OpenFOAM-8 installation.' >&2
    return 2
fi
export DEP_ROOT DEP_OF_ROOT DEP_PYTHON DEP_MPI_LAUNCHER
export WM_PROJECT_USER_DIR="$DEP_ROOT/install/OpenFOAM-user"
export FOAM_USER_APPBIN="$DEP_ROOT/install/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$DEP_ROOT/install/$WM_OPTIONS/lib"
export FOAM_RUN="$DEP_ROOT/runs"
export EIGEN_RHEO="$DEP_ROOT/deps/eigen-3.4.0"
export WM_NCOMPPROCS="$DEP_BUILD_JOBS"
export PATH="$FOAM_USER_APPBIN:$PATH"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$FOAM_USER_APPBIN" "$FOAM_USER_LIBBIN" "$FOAM_RUN" "$DEP_ROOT/logs"
case "$FOAM_USER_APPBIN:$FOAM_USER_LIBBIN" in
    *"$DEP_OF_ROOT"*) echo 'Refusing to write add-on binaries into the base OpenFOAM tree.' >&2; return 2;;
esac
unset dep_init dep_modules dep_was_nounset dep_was_errexit dep_source_status
