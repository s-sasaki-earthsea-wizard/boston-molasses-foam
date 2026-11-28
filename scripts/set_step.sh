#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# set_step.sh - switch a case to a step
#
# Usage:
#   [MESH_REFINE=r] [END_TIME=t] [WRITE_INTERVAL=dt] \
#       scripts/set_step.sh <1|2|3|4|5> [CASE_DIR] [VARIANT]
#
#   CASE_DIR defaults to case. VARIANT selects the molasses viscosity of
#   Steps 4-5: low, central (default) or high; Steps 1-3 accept only central.
#   MESH_REFINE (default 1) multiplies the cell numbers of the step's mesh;
#   END_TIME and WRITE_INTERVAL override the step's endTime / writeInterval.
#
# Sets the liquid properties and model switches in constant/transportProperties,
# the initial liquid temperature (liquidT) in system/setFieldsDict, the domain
# and mesh in system/blockMeshDict, and endTime / writeInterval in
# system/controlDict. Re-run blockMesh and setFields afterwards when those
# changed (run-stepN does that).
#
# Steps:
#   1  Isothermal water dam-break         energy off, constant viscosity
#   2  Water + heat exchange with air     energy on,  constant viscosity
#   3  Water + Andrade viscosity mu(T)    energy on,  mu(T) on
#   4  Boston molasses, wave stage        energy on,  mu(T) on, molasses data
#   5  Boston molasses, spreading         as 4 on a 400 m domain with a
#                                         bottom-refined mesh, 600 s
#   Steps 1-4: 180 m domain, 0.5 m cells, 20 s; Step 5: see blockMeshDict.
#
# Existing "key value;" entries are replaced in place (comments, layout and
# other entries are kept). Every replacement must match exactly one entry,
# otherwise the script fails without touching any file. A file is only
# written when its content changes, so make sees unchanged dictionaries as
# up to date.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

usage() {
    sed -n '3,33p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

[[ "${1:-}" == "-h" || "${1:-}" == "--help" ]] && usage 0

STEP="${1:-}"
CASE_DIR="${2:-${ROOT_DIR}/case}"
VARIANT="${3:-central}"
MESH_REFINE="${MESH_REFINE:-1}"
TP="${CASE_DIR}/constant/transportProperties"
SFD="${CASE_DIR}/system/setFieldsDict"
BMD="${CASE_DIR}/system/blockMeshDict"
CD="${CASE_DIR}/system/controlDict"

if [[ ! "${STEP}" =~ ^[1-5]$ ]]; then
    echo "Error: step must be 1, 2, 3, 4 or 5 (got '${STEP}')" >&2
    usage 1
fi
if [[ ! "${VARIANT}" =~ ^(low|central|high)$ ]]; then
    echo "Error: viscosity variant must be low, central or high (got '${VARIANT}')" >&2
    usage 1
fi
if [[ "${STEP}" -lt 4 && "${VARIANT}" != central ]]; then
    echo "Error: the viscosity variant '${VARIANT}' applies to Steps 4-5 only" >&2
    exit 1
fi
number_re='^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$'
for name in MESH_REFINE END_TIME WRITE_INTERVAL; do
    value="${!name:-1}"
    if [[ ! "${value}" =~ ${number_re} ]] || awk -v v="${value}" 'BEGIN { exit !(v <= 0) }'; then
        echo "Error: ${name} must be a positive number (got '${value}')" >&2
        exit 1
    fi
done
for f in "${TP}" "${SFD}" "${BMD}" "${CD}"; do
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

# Cane final molasses (Step 4). Literature estimates: no source describes the
# 1919 molasses itself. Derivations, ranges and caveats: wiki page
# "Molasses properties"; sources: README, "Molasses property sources".
#   TINIT  about 5 degC above the air (Sharp, Kennedy & Rubinstein 2016);
#          T0 = TINIT, so mu0 is the viscosity at the start
#   mu0    Broadfoot & Steindl (1980) correlation at 40 % purity, shear
#          rate 10 1/s and 10 degC; uncertain by ~x3 through the unknown dry
#          substance, hence three variants:
#            low      20 Pa s   76 % dry substance, 40 % purity
#            central  55 Pa s   78 %, 40 % (reproduces Sharp et al.'s x3)
#            high    200 Pa s   80 %, 35 %
#          Final molasses measured by Shah (2017) near 81 Brix gives
#          48-129 Pa s when moved from 30 to 10 degC with this B.
#          B is the same in all three; it hardly matters until cooling is
#          resolved
#   B      x3 from 10 to 0 degC (Sharp et al. 2016), Ea ~ 71 kJ/mol
#   rho    Rouillard (1985), Brix 80 at 10 degC
#   Cp     Janovskii & Archangelskii (1928) / Rouillard (1985), 78-80 % dry
#          substance at 10 degC
#   kappa  Watson (1986), Brix 80 at 10 degC
#   sigma  no effect at this mesh (capillary length ~2 mm)
MOLASSES_RHO=1430     # [kg/m3]
MOLASSES_MU0_LOW=20       # [Pa s] at T0
MOLASSES_MU0_CENTRAL=55   # [Pa s] at T0
MOLASSES_MU0_HIGH=200     # [Pa s] at T0
MOLASSES_CP=2050      # [J/kg/K]
MOLASSES_KAPPA=0.31   # [W/m/K]
MOLASSES_B=8.5e3      # [K]
MOLASSES_SIGMA=0.075  # [N/m]
MOLASSES_T0=283.15    # [K]
MOLASSES_TINIT=283.15 # [K]

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
    4|5)
        if [[ "${STEP}" == 4 ]]; then
            LABEL="Boston molasses, wave stage, ${VARIANT} viscosity"
        else
            LABEL="Boston molasses, spreading, ${VARIANT} viscosity"
        fi
        case "${VARIANT}" in
            low)     MU0=${MOLASSES_MU0_LOW} ;;
            central) MU0=${MOLASSES_MU0_CENTRAL} ;;
            high)    MU0=${MOLASSES_MU0_HIGH} ;;
        esac
        RHO=${MOLASSES_RHO}; CP=${MOLASSES_CP}
        KAPPA=${MOLASSES_KAPPA}; B=${MOLASSES_B}; SIGMA=${MOLASSES_SIGMA}
        T0=${MOLASSES_T0}; TINIT=${MOLASSES_TINIT}
        SOLVE_ENERGY=yes; T_DEPENDENT=yes
        ;;
esac

# --- Domain, mesh and run time --------------------------------------------------
# Steps 1-4: the wave stage on a 180 m domain with 0.5 m cells; the fronts
#            reach the open end within 9.5-16.5 s, so 20 s.
# Step 5:    the spreading stage, 600 s on a 400 m domain; the front then
#            lies within 10 % of Huppert's viscous current. dx = 0.5 m;
#            dy = 0.05 m up to y = 2.5 m, which resolves a 0.4-1.2 m thick
#            layer with 8-24 cells, then growing geometrically (r = 1.0488,
#            48 cells) to 0.49 m at y = 12 m, the column top, and 0.5 m above.
#            The grading sections are given as (length, cells, last/first).
case "${STEP}" in
    [1-4])
        LENGTH_X=180; BASE_CELLS_X=360; BASE_CELLS_Y=40; GRADING_Y=1
        STEP_END_TIME=20; STEP_WRITE_INTERVAL=0.5
        ;;
    5)
        LENGTH_X=400; BASE_CELLS_X=800; BASE_CELLS_Y=114
        GRADING_Y="((2.5 50 1) (9.5 48 9.3798) (8 16 1))"
        STEP_END_TIME=600; STEP_WRITE_INTERVAL=1
        ;;
esac
END_TIME="${END_TIME:-${STEP_END_TIME}}"
WRITE_INTERVAL="${WRITE_INTERVAL:-${STEP_WRITE_INTERVAL}}"

# scale_cells <cells>: cells * MESH_REFINE; fails unless a whole number >= 1
scale_cells() {
    awk -v n="$1" -v r="${MESH_REFINE}" 'BEGIN {
        c = n * r; i = int(c + 0.5)
        if (i < 1 || c - i > 1e-9 || i - c > 1e-9) exit 1
        print i
    }'
}
if ! CELLS_X="$(scale_cells "${BASE_CELLS_X}")" || ! CELLS_Y="$(scale_cells "${BASE_CELLS_Y}")"; then
    echo "Error: MESH_REFINE=${MESH_REFINE} does not give whole cell numbers" \
         "(${BASE_CELLS_X} x ${BASE_CELLS_Y} cells)" >&2
    exit 1
fi

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
WORK_BMD="$(mktemp "${TMPDIR:-/tmp}/blockMeshDict.XXXXXX")"
WORK_CD="$(mktemp "${TMPDIR:-/tmp}/controlDict.XXXXXX")"
trap 'for w in "${WORK}" "${WORK_SFD}" "${WORK_BMD}" "${WORK_CD}"; do rm -f "${w}" "${w}.new"; done' EXIT
cp "${TP}" "${WORK}"
cp "${SFD}" "${WORK_SFD}"
cp "${BMD}" "${WORK_BMD}"
cp "${CD}" "${WORK_CD}"

TARGET="${TP}"
set_entry "${WORK}" activeStep                            "${STEP}"
set_entry "${WORK}" viscosityVariant                      "$([[ ${STEP} -ge 4 ]] && echo "${VARIANT}" || echo none)"
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

TARGET="${BMD}"
set_entry "${WORK_BMD}" lengthX                           "${LENGTH_X}"
set_entry "${WORK_BMD}" cellsX                            "${CELLS_X}"
set_entry "${WORK_BMD}" cellsY                            "${CELLS_Y}"
set_entry "${WORK_BMD}" gradingY                          "${GRADING_Y}"

TARGET="${CD}"
set_entry "${WORK_CD}" endTime                            "${END_TIME}"
set_entry "${WORK_CD}" writeInterval                      "${WRITE_INTERVAL}"

install_if_changed "${WORK}" "${TP}"
install_if_changed "${WORK_SFD}" "${SFD}"
install_if_changed "${WORK_BMD}" "${BMD}"
install_if_changed "${WORK_CD}" "${CD}"

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
    "temperatureDependent" "${T_DEPENDENT}" \
    "domain length [m]"    "${LENGTH_X}" \
    "cells (x, y)"         "${CELLS_X} x ${CELLS_Y} (MESH_REFINE ${MESH_REFINE})" \
    "y grading"            "${GRADING_Y}" \
    "endTime [s]"          "${END_TIME}" \
    "writeInterval [s]"    "${WRITE_INTERVAL}"
