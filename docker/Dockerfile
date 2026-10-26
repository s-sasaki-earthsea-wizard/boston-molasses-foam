# -----------------------------------------------------------------------------
# OpenFOAM v2312 (OpenCFD) development/runtime image for interThermalFoam
#
# Base: official "-dev" image = OpenFOAM runtime + compilers, wmake, headers
# and solver sources ($FOAM_SOLVERS/multiphase/VoF is used by the solver).
# /openfoam/profile.rc loads the OpenFOAM environment; it is sourced by
# entrypoint.sh and by interactive shells (/etc/bash.bashrc).
# -----------------------------------------------------------------------------
FROM opencfd/openfoam-dev:2312

ENV PROJECT_ROOT=/work

# Project overrides, loaded for interactive and non-interactive shells:
#  - user-compiled binaries/libraries go to the mounted project (build/) so
#    that they survive `docker compose run --rm`; the default
#    $HOME/OpenFOAM/... is neither writable nor persistent here
#  - start in the project root if a shell lands in $HOME
RUN cat >> /openfoam/profile.rc <<'EOF'

# --- interThermalFoam project overrides --------------------------------------
if [ -n "${PROJECT_ROOT:-}" ] && [ -n "${WM_PROJECT_DIR:-}" ]
then
    export FOAM_USER_APPBIN="${PROJECT_ROOT}/build/bin"
    export FOAM_USER_LIBBIN="${PROJECT_ROOT}/build/lib"
    export PATH="${FOAM_USER_APPBIN}:${PATH}"
    export LD_LIBRARY_PATH="${FOAM_USER_LIBBIN}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

    if [ "${PWD}" = "${HOME}" ] && [ -d "${PROJECT_ROOT}" ]
    then
        cd "${PROJECT_ROOT}" || true
    fi
fi
EOF

COPY --chmod=0755 entrypoint.sh /usr/local/bin/foam-entrypoint

WORKDIR ${PROJECT_ROOT}

ENTRYPOINT ["/usr/local/bin/foam-entrypoint"]
