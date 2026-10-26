#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# set_step.sh - switch constant/transportProperties to a verification step
#
# Usage:
#   scripts/set_step.sh <1|2|3|4> [CASE_DIR]      (CASE_DIR defaults to case)
#
# Steps:
#   1  Isothermal water dam-break         energy off, constant viscosity
#   2  Water + heat exchange with air     energy on,  constant viscosity
#   3  Water + Andrade viscosity mu(T)    energy on,  mu(T) on
#   4  Boston molasses                    energy on,  mu(T) on, molasses data
#
# Existing "key value;" entries are replaced in place (comments, layout and
# other entries are kept). Every replacement must match exactly one entry,
# otherwise the script fails without touching the file.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

usage() {
    sed -n '3,16p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0

STEP="${1:-}"
CASE_DIR="${2:-${ROOT_DIR}/case}"
TP="${CASE_DIR}/constant/transportProperties"

if [[ ! "${STEP}" =~ ^[1-4]$ ]]; then
    echo "Error: step must be 1, 2, 3 or 4 (got '${STEP}')" >&2
    usage 1
fi
if [[ ! -f "${TP}" ]]; then
    echo "Error: ${TP} not found" >&2
    exit 1
fi

# --- Liquid properties --------------------------------------------------------
# Water at 20 degC. B fitted to mu(4 degC)/mu(20 degC) = 1.5673/1.0016 mPa s.
WATER_RHO=998.2       # [kg/m3]
WATER_MU0=1.002e-3    # [Pa s] at T0
WATER_CP=4182         # [J/kg/K]
WATER_KAPPA=0.598     # [W/m/K]
WATER_B=2274          # [K]
WATER_SIGMA=0.0728    # [N/m]

# Cane molasses (Step 4). PLACEHOLDER values of plausible magnitude, to be
# replaced by literature/measured data (mu0 and B dominate the result):
#   B = Ea/R = 1.0e4 K (Ea ~ 83 kJ/mol) -> mu(4 degC)/mu(20 degC) ~ 7
MOLASSES_RHO=1400     # [kg/m3]
MOLASSES_MU0=5.0      # [Pa s] at T0
MOLASSES_CP=2300      # [J/kg/K]
MOLASSES_KAPPA=0.35   # [W/m/K]
MOLASSES_B=1.0e4      # [K]
MOLASSES_SIGMA=0.075  # [N/m]

T_REF=293.15          # [K] reference temperature T0 of mu(T)

case "${STEP}" in
    1)
        LABEL="Isothermal water dam-break"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        SOLVE_ENERGY=no; T_DEPENDENT=no
        ;;
    2)
        LABEL="Water + heat exchange with the atmosphere"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        SOLVE_ENERGY=yes; T_DEPENDENT=no
        ;;
    3)
        LABEL="Water + temperature-dependent viscosity"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        SOLVE_ENERGY=yes; T_DEPENDENT=yes
        ;;
    4)
        LABEL="Boston molasses (high density and viscosity)"
        RHO=${MOLASSES_RHO}; MU0=${MOLASSES_MU0}; CP=${MOLASSES_CP}
        KAPPA=${MOLASSES_KAPPA}; B=${MOLASSES_B}; SIGMA=${MOLASSES_SIGMA}
        SOLVE_ENERGY=yes; T_DEPENDENT=yes
        ;;
esac

# Kinematic viscosity for the Newtonian model of the mixture class (kept
# consistent with mu0; the solver warns if they disagree)
NU="$(awk -v mu="${MU0}" -v rho="${RHO}" 'BEGIN { printf "%.6e", mu/rho }')"

# set_entry <file> <dict/path/key> <value>
#   Replace the value of "key value;" inside the (sub)dictionary path.
#   Supports the layout used in this case: "name" and "{" on separate lines
#   (or "name {"), one entry per line, // comments.
set_entry() {
    local file="$1" path="$2" value="$3"

    if ! awk -v path="${path}" -v value="${value}" '
        function current(   i, s) {
            s = ""
            for (i = 1; i <= depth; i++) s = s stack[i] "/"
            return s
        }
        BEGIN {
            n = split(path, parts, "/")
            key = parts[n]
            prefix = ""
            for (i = 1; i < n; i++) prefix = prefix parts[i] "/"
            depth = 0; hits = 0; pending = ""
        }
        {
            line = $0
            code = line
            sub(/\/\/.*$/, "", code)

            if (code ~ /^[ \t]*[A-Za-z_][A-Za-z0-9_.]*[ \t]*$/) {
                pending = code
                gsub(/[ \t]/, "", pending)
            }
            else if (code ~ /^[ \t]*[A-Za-z_][A-Za-z0-9_.]*[ \t]*\{[ \t]*$/) {
                name = code
                gsub(/[ \t{]/, "", name)
                stack[++depth] = name
            }
            else if (code ~ /^[ \t]*\{[ \t]*$/) {
                stack[++depth] = pending
                pending = ""
            }
            else if (code ~ /^[ \t]*\}[ \t]*$/) {
                depth--
            }
            else if (current() == prefix \
                  && match(code, "^[ \t]*" key "[ \t]+[^;]*;")) {
                # Keep indentation, key alignment and the trailing comment
                # column (when there is room for the new value)
                match(line, "^[ \t]*" key "[ \t]+")
                head = substr(line, 1, RLENGTH)
                rest = substr(line, RLENGTH + 1)
                semi = index(rest, ";")
                tail = substr(rest, semi + 1)
                match(tail, "^[ \t]*")
                width = semi + RLENGTH
                comment = substr(tail, RLENGTH + 1)
                line = head value ";"
                if (comment != "") {
                    pad = width - length(value) - 1
                    if (pad < 1) pad = 1
                    line = line sprintf("%" pad "s", "") comment
                }
                hits++
            }
            print line
        }
        END { exit (hits == 1 ? 0 : 3) }
    ' "${file}" > "${file}.new"; then
        echo "Error: entry '${path}' not found exactly once in ${TP}" >&2
        exit 1
    fi
    mv "${file}.new" "${file}"
}

echo "Step ${STEP}: ${LABEL}"
echo "  -> ${TP}"

# Edit a working copy; the case file is only replaced if every entry matched
WORK="$(mktemp "${TMPDIR:-/tmp}/transportProperties.XXXXXX")"
trap 'rm -f "${WORK}" "${WORK}.new"' EXIT
cp "${TP}" "${WORK}"

set_entry "${WORK}" activeStep                            "${STEP}"
set_entry "${WORK}" water/rho                             "${RHO}"
set_entry "${WORK}" water/nu                              "${NU}"
set_entry "${WORK}" water/Cp                              "${CP}"
set_entry "${WORK}" water/kappa                           "${KAPPA}"
set_entry "${WORK}" sigma                                 "${SIGMA}"
set_entry "${WORK}" energy/solveEnergy                    "${SOLVE_ENERGY}"
set_entry "${WORK}" liquidViscosity/temperatureDependent  "${T_DEPENDENT}"
set_entry "${WORK}" liquidViscosity/mu0                   "${MU0}"
set_entry "${WORK}" liquidViscosity/T0                    "${T_REF}"
set_entry "${WORK}" liquidViscosity/B                     "${B}"

cat "${WORK}" > "${TP}"

printf '  %-22s %s\n' \
    "rho [kg/m3]"          "${RHO}" \
    "mu0 [Pa s]"           "${MU0}" \
    "nu [m2/s]"            "${NU}" \
    "Cp [J/kg/K]"          "${CP}" \
    "kappa [W/m/K]"        "${KAPPA}" \
    "B [K]"                "${B}" \
    "solveEnergy"          "${SOLVE_ENERGY}" \
    "temperatureDependent" "${T_DEPENDENT}"
