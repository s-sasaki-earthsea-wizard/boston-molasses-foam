#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# prepare_anim_case.sh - set up a case directory for an animation run
#
# Usage:
#   scripts/prepare_anim_case.sh SRC_CASE DST_CASE END_TIME WRITE_INTERVAL
#
# Copies the pristine inputs of SRC_CASE (0.orig/, system/, constant/ without
# polyMesh/) into DST_CASE and sets endTime and writeInterval in its
# controlDict. Timestamps are kept, so the mesh and 0/ of DST_CASE are only
# regenerated when the inputs changed. Results in DST_CASE are left alone
# (run-stepN removes them before running).
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
    sed -n '3,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0
[[ $# -eq 4 ]] || usage 1

SRC="$1"
DST="$2"
END_TIME="$3"
WRITE_INTERVAL="$4"

number_re='^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'
for value in "${END_TIME}" "${WRITE_INTERVAL}"; do
    if [[ ! "${value}" =~ ${number_re} ]]; then
        echo "Error: '${value}' is not a positive number" >&2
        exit 1
    fi
done

if [[ ! -d "${SRC}/system" || ! -d "${SRC}/0.orig" ]]; then
    echo "Error: '${SRC}' is not a case with system/ and 0.orig/" >&2
    exit 1
fi
if [[ "$(realpath -m "${SRC}")" == "$(realpath -m "${DST}")" ]]; then
    echo "Error: the animation case must differ from the source case '${SRC}'" >&2
    exit 1
fi
if [[ -d "${DST}" && ! -d "${DST}/system" && -n "$(ls -A "${DST}")" ]]; then
    echo "Error: '${DST}' exists and is not a case directory" >&2
    exit 1
fi

echo "Preparing ${DST} from ${SRC} (endTime ${END_TIME}, writeInterval ${WRITE_INTERVAL})"

mkdir -p "${DST}/constant"
rm -rf "${DST}/0.orig" "${DST}/system"
cp -rp "${SRC}/0.orig" "${SRC}/system" "${DST}/"
for f in "${SRC}"/constant/*; do
    name="$(basename "${f}")"
    [[ "${name}" == "polyMesh" ]] && continue
    rm -rf "${DST}/constant/${name}"
    cp -rp "${f}" "${DST}/constant/"
done

# set_control <key> <value>: replace the top-level "key value;" entry
set_control() {
    local key="$1" value="$2" file="${DST}/system/controlDict"
    local hits
    hits="$(grep -cE "^${key}[[:space:]]+[^;]*;" "${file}" || true)"
    if [[ "${hits}" -ne 1 ]]; then
        echo "Error: entry '${key}' not found exactly once in ${file}" >&2
        exit 1
    fi
    sed -i -E "s/^(${key}[[:space:]]+)[^;]*;/\1${value};/" "${file}"
}

set_control endTime "${END_TIME}"
set_control writeInterval "${WRITE_INTERVAL}"
