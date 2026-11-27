#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# prepare_anim_case.sh - set up a case directory for an animation run
#
# Usage:
#   scripts/prepare_anim_case.sh SRC_CASE DST_CASE
#
# Copies the pristine inputs of SRC_CASE (0.orig/, system/, constant/ without
# polyMesh/) into DST_CASE. Timestamps are kept, so the mesh and 0/ of
# DST_CASE are only regenerated when the inputs changed. Results in DST_CASE
# are left alone (run-stepN removes them before running). endTime and
# writeInterval are set by run-stepN (END_TIME, WRITE_INTERVAL).
# -----------------------------------------------------------------------------
set -euo pipefail

usage() {
    sed -n '3,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0
[[ $# -eq 2 ]] || usage 1

SRC="$1"
DST="$2"

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

echo "Preparing ${DST} from ${SRC}"

mkdir -p "${DST}/constant"
rm -rf "${DST}/0.orig" "${DST}/system"
cp -rp "${SRC}/0.orig" "${SRC}/system" "${DST}/"
for f in "${SRC}"/constant/*; do
    name="$(basename "${f}")"
    [[ "${name}" == "polyMesh" ]] && continue
    rm -rf "${DST}/constant/${name}"
    cp -rp "${f}" "${DST}/constant/"
done
