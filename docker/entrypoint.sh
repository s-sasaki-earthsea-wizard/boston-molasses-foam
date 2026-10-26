#!/bin/bash
# -----------------------------------------------------------------------------
# Container entrypoint (replaces the upstream /openfoam/run)
#
#   1. map an arbitrary --user UID/GID to the user "openfoam" (nss_wrapper)
#   2. load the OpenFOAM environment and the project overrides
#      (/openfoam/profile.rc, extended in the Dockerfile)
#   3. exec the command, or an interactive bash when none is given
#
# Why not /openfoam/run: it sources the environment from /bin/sh (dash), where
# `. file args` does not replace the positional parameters. The container
# command then reaches OpenFOAM's etc/config.sh/setup as FOAM_SETTINGS and
# fragments of it are eval'ed in $HOME before the real command runs.
# -----------------------------------------------------------------------------

uid="$(id -u)"
gid="$(id -g)"

nss_lib="/usr/lib/$(uname -m)-linux-gnu/libnss_wrapper.so"

if [ "${uid}" != 0 ] && [ -e "${nss_lib}" ]
then
    export USER=openfoam
    export HOME="/tmp/home-${USER}"
    mkdir -p "${HOME}"

    export NSS_WRAPPER_PASSWD=/tmp/passwd.nss
    export NSS_WRAPPER_GROUP=/tmp/group.nss
    echo "${USER}:x:${uid}:${gid}:${USER}:${HOME}:/bin/bash" > "${NSS_WRAPPER_PASSWD}"
    echo "${USER}:x:${gid}:" > "${NSS_WRAPPER_GROUP}"
    export LD_PRELOAD="${nss_lib}${LD_PRELOAD:+:${LD_PRELOAD}}"
fi

# Keep the command away from OpenFOAM's etc/bashrc argument parsing
args=("$@")
set --

if [ -f /openfoam/profile.rc ]
then
    . /openfoam/profile.rc
fi

if [ "${#args[@]}" -eq 0 ]
then
    exec /bin/bash -i
fi

exec "${args[@]}"
