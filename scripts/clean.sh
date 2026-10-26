#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# clean.sh - remove simulation results and (optionally) all generated files
#
# Usage:
#   scripts/clean.sh [--all] [CASE_DIR]           (CASE_DIR defaults to case)
#
#   default  time directories except 0/, processor*/, postProcessing/,
#            solver logs, VTK/, dynamicCode/
#   --all    additionally the mesh (constant/polyMesh), the initialised 0/,
#            all logs, the solver binary (build/) and wmake objects
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

usage() {
    sed -n '3,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

ALL=false
CASE_DIR="${ROOT_DIR}/case"

for arg in "$@"; do
    case "${arg}" in
        --all)     ALL=true ;;
        -h|--help) usage 0 ;;
        -*)        echo "Error: unknown option '${arg}'" >&2; usage 1 ;;
        *)         CASE_DIR="${arg}" ;;
    esac
done

if [[ ! -d "${CASE_DIR}/system" ]]; then
    echo "Error: '${CASE_DIR}' is not an OpenFOAM case (no system/)" >&2
    exit 1
fi

remove() {
    local p
    for p in "$@"; do
        if [[ -e "${p}" || -L "${p}" ]]; then
            echo "  rm ${p}"
            rm -rf -- "${p}"
        fi
    done
}

shopt -s nullglob

echo "Cleaning results in ${CASE_DIR}"

# Time directories: numeric names (0.5, 12, 1e-05, ...), keeping 0/
time_re='^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'
for d in "${CASE_DIR}"/*/; do
    name="$(basename "${d}")"
    if [[ "${name}" =~ ${time_re} && "${name}" != "0" ]]; then
        remove "${d%/}"
    fi
done

remove "${CASE_DIR}"/processor* \
       "${CASE_DIR}/postProcessing" \
       "${CASE_DIR}"/log.interThermalFoam* \
       "${CASE_DIR}/VTK" \
       "${CASE_DIR}/dynamicCode"

if ${ALL}; then
    echo "Cleaning mesh, initial fields and build artefacts"
    remove "${CASE_DIR}/0" \
           "${CASE_DIR}/constant/polyMesh" \
           "${CASE_DIR}"/log.* \
           "${ROOT_DIR}/build"

    # wmake objects: solver/<name>/Make/<WM_OPTIONS>/
    for d in "${ROOT_DIR}"/solver/*/Make/*/; do
        remove "${d%/}"
    done
fi
