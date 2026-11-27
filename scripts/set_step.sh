#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# set_step.sh - switch a case to a verification step
#
# Usage:
#   scripts/set_step.sh <1|2|3|4> [CASE_DIR]      (CASE_DIR defaults to case)
#
# Sets the liquid properties and model switches in constant/transportProperties
# and the initial liquid temperature (liquidT) in system/setFieldsDict.
# Re-run setFields afterwards when liquidT changed (run-stepN does that).
#
# Steps:
#   1  Isothermal water dam-break         energy off, constant viscosity
#   2  Water + heat exchange with air     energy on,  constant viscosity
#   3  Water + Andrade viscosity mu(T)    energy on,  mu(T) on
#   4  Boston molasses                    energy on,  mu(T) on, molasses data
#
# Existing "key value;" entries are replaced in place (comments, layout and
# other entries are kept). Every replacement must match exactly one entry,
# otherwise the script fails without touching either file. A file is only
# written when its content changes, so make sees an unchanged setFieldsDict
# as up to date.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

usage() {
    sed -n '3,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0

STEP="${1:-}"
CASE_DIR="${2:-${ROOT_DIR}/case}"
TP="${CASE_DIR}/constant/transportProperties"
SFD="${CASE_DIR}/system/setFieldsDict"

if [[ ! "${STEP}" =~ ^[1-4]$ ]]; then
    echo "Error: step must be 1, 2, 3 or 4 (got '${STEP}')" >&2
    usage 1
fi
for f in "${TP}" "${SFD}"; do
    if [[ ! -f "${f}" ]]; then
        echo "Error: ${f} not found" >&2
        exit 1
    fi
done

# --- Liquid properties --------------------------------------------------------
# T0 is the reference temperature of mu0 in mu(T) = mu0 exp(B (1/T - 1/T0)),
# TINIT the initial temperature of the liquid block.
#
# Water at 20 degC. B fitted to mu(4 degC)/mu(20 degC) = 1.5673/1.0016 mPa s.
WATER_RHO=998.2       # [kg/m3]
WATER_MU0=1.002e-3    # [Pa s] at T0
WATER_CP=4182         # [J/kg/K]
WATER_KAPPA=0.598     # [W/m/K]
WATER_B=2274          # [K]
WATER_SIGMA=0.0728    # [N/m]
WATER_T0=293.15       # [K]
WATER_TINIT=293.15    # [K]

# Cane molasses (Step 4). PLACEHOLDER values of plausible magnitude, to be
# replaced by literature/measured data (mu0 and B dominate the result):
#   B = Ea/R = 1.0e4 K (Ea ~ 83 kJ/mol) -> mu(4 degC)/mu(20 degC) ~ 7
MOLASSES_RHO=1400     # [kg/m3]
MOLASSES_MU0=5.0      # [Pa s] at T0
MOLASSES_CP=2300      # [J/kg/K]
MOLASSES_KAPPA=0.35   # [W/m/K]
MOLASSES_B=1.0e4      # [K]
MOLASSES_SIGMA=0.075  # [N/m]
MOLASSES_T0=293.15    # [K]
MOLASSES_TINIT=293.15 # [K]

case "${STEP}" in
    1)
        LABEL="Isothermal water dam-break"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        T0=${WATER_T0}; TINIT=${WATER_TINIT}
        SOLVE_ENERGY=no; T_DEPENDENT=no
        ;;
    2)
        LABEL="Water + heat exchange with the atmosphere"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        T0=${WATER_T0}; TINIT=${WATER_TINIT}
        SOLVE_ENERGY=yes; T_DEPENDENT=no
        ;;
    3)
        LABEL="Water + temperature-dependent viscosity"
        RHO=${WATER_RHO}; MU0=${WATER_MU0}; CP=${WATER_CP}
        KAPPA=${WATER_KAPPA}; B=${WATER_B}; SIGMA=${WATER_SIGMA}
        T0=${WATER_T0}; TINIT=${WATER_TINIT}
        SOLVE_ENERGY=yes; T_DEPENDENT=yes
        ;;
    4)
        LABEL="Boston molasses (high density and viscosity)"
        RHO=${MOLASSES_RHO}; MU0=${MOLASSES_MU0}; CP=${MOLASSES_CP}
        KAPPA=${MOLASSES_KAPPA}; B=${MOLASSES_B}; SIGMA=${MOLASSES_SIGMA}
        T0=${MOLASSES_T0}; TINIT=${MOLASSES_TINIT}
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
        echo "Error: entry '${path}' not found exactly once in ${TARGET}" >&2
        exit 1
    fi
    mv "${file}.new" "${file}"
}

# install_if_changed <work copy> <target>
#   Write the target only if the work copy differs (keeps the inode and the
#   timestamp of an unchanged file)
install_if_changed() {
    local work="$1" target="$2"
    if cmp -s "${work}" "${target}"; then
        echo "  -> ${target} (unchanged)"
    else
        cat "${work}" > "${target}"
        echo "  -> ${target}"
    fi
}

echo "Step ${STEP}: ${LABEL}"

# Edit working copies; the case files are only replaced if every entry matched
WORK="$(mktemp "${TMPDIR:-/tmp}/transportProperties.XXXXXX")"
WORK_SFD="$(mktemp "${TMPDIR:-/tmp}/setFieldsDict.XXXXXX")"
trap 'rm -f "${WORK}" "${WORK}.new" "${WORK_SFD}" "${WORK_SFD}.new"' EXIT
cp "${TP}" "${WORK}"
cp "${SFD}" "${WORK_SFD}"

TARGET="${TP}"
set_entry "${WORK}" activeStep                            "${STEP}"
set_entry "${WORK}" water/rho                             "${RHO}"
set_entry "${WORK}" water/nu                              "${NU}"
set_entry "${WORK}" water/Cp                              "${CP}"
set_entry "${WORK}" water/kappa                           "${KAPPA}"
set_entry "${WORK}" sigma                                 "${SIGMA}"
set_entry "${WORK}" energy/solveEnergy                    "${SOLVE_ENERGY}"
set_entry "${WORK}" liquidViscosity/temperatureDependent  "${T_DEPENDENT}"
set_entry "${WORK}" liquidViscosity/mu0                   "${MU0}"
set_entry "${WORK}" liquidViscosity/T0                    "${T0}"
set_entry "${WORK}" liquidViscosity/B                     "${B}"

TARGET="${SFD}"
set_entry "${WORK_SFD}" liquidT                           "${TINIT}"

install_if_changed "${WORK}" "${TP}"
install_if_changed "${WORK_SFD}" "${SFD}"

printf '  %-22s %s\n' \
    "rho [kg/m3]"          "${RHO}" \
    "mu0 [Pa s]"           "${MU0}" \
    "nu [m2/s]"            "${NU}" \
    "Cp [J/kg/K]"          "${CP}" \
    "kappa [W/m/K]"        "${KAPPA}" \
    "B [K]"                "${B}" \
    "T0 [K]"               "${T0}" \
    "initial liquid T [K]" "${TINIT}" \
    "solveEnergy"          "${SOLVE_ENERGY}" \
    "temperatureDependent" "${T_DEPENDENT}"
