#!/usr/bin/env bash
# Three production solves are independent within each flow rate and mesh.
set -Eeuo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
exec "${DEP_PYTHON:-python3}" tools/steady_workflow.py submit "$@"
