# -----------------------------------------------------------------------------
# docker.mk - Docker image and container helpers
#
# COMPOSE_OVERRIDE is a local, uncommitted compose file (e.g. a mount for the
# disk that symlinked case-* directories point to). It is added to every
# compose call when it exists: with an explicit -f, docker compose does not
# load docker-compose.override.yml by itself.
# -----------------------------------------------------------------------------

COMPOSE_FILE     ?= docker/docker-compose.yml
COMPOSE_OVERRIDE ?= docker/docker-compose.override.yml
SERVICE          ?= openfoam
HOST_UID         ?= $(shell id -u)
HOST_GID         ?= $(shell id -g)

COMPOSE_FILES    := $(COMPOSE_FILE) $(wildcard $(COMPOSE_OVERRIDE))

# Files written into the mounted project are owned by the invoking user
COMPOSE          := HOST_UID=$(HOST_UID) HOST_GID=$(HOST_GID) docker compose $(addprefix -f ,$(COMPOSE_FILES))

# Non-interactive container command (no TTY, removed after exit)
FOAM_RUN         := $(COMPOSE) run --rm -T $(SERVICE)

# $(call foam,<command>)
#   Run a shell command from the project root inside the container with the
#   OpenFOAM environment loaded. pipefail keeps `solver | tee log` failures.
#   Note: <command> must not contain commas (make argument separator).
foam = $(FOAM_RUN) bash -o pipefail -c 'cd "$$PROJECT_ROOT" && $(1)'

.PHONY: docker-build docker-shell

docker-build: ## Build the OpenFOAM v2312 Docker image
	$(COMPOSE) build

docker-shell: ## Open an interactive bash shell in the container
	$(COMPOSE) run --rm $(SERVICE)
