# =============================================================================
#  Boston Molasses Flood (1919) - 2D thermo-viscous VOF dam-break
#  OpenFOAM v2312 + Docker
#
#  Every OpenFOAM command is wrapped in `docker compose run --rm openfoam ...`
#  (see make/docker.mk); animations are rendered on the host with uv and
#  ffmpeg (see make/post.mk). Run `make help` for the list of targets.
# =============================================================================

SHELL         := /bin/bash
.SHELLFLAGS   := -eu -o pipefail -c
.DEFAULT_GOAL := help

# Paths relative to the project root (mounted at /work in the container)
CASE_DIR      := case
SOLVER_NAME   := interThermalFoam
SOLVER_DIR    := solver/$(SOLVER_NAME)
BUILD_DIR     := build

include make/docker.mk
include make/solver.mk
include make/case.mk
include make/post.mk

.PHONY: help clean clean-all

help: ## Show this help
	@printf 'Usage: make <target>\n\nTargets:\n'
	@awk 'BEGIN { FS = ":.*## " } /^[a-zA-Z0-9_-]+:.*## / { printf "  %-14s %s\n", $$1, $$2 }' $(MAKEFILE_LIST)

clean: ## Remove result time directories, solver logs and postProcessing
	bash scripts/clean.sh $(CASE_DIR)

clean-all: ## clean + mesh, initialised 0/, solver binary and wmake objects
	bash scripts/clean.sh --all $(CASE_DIR)
