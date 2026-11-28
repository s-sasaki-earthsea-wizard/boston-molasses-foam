# -----------------------------------------------------------------------------
# cooling.mk - the 1D cooling column of #9 (host: uv; no OpenFOAM)
#
#   cooling-runs  run every scenario of cooling/scenario.py (both days, six
#                 layer thicknesses, the paths and one-at-a-time changes) and
#                 write the summary to COOLING_SUMMARY and the figures
#                 docs/media/cooling-*.svg; the revision written into both is
#                 the last commit of cooling/, marked -dirty if cooling/ has
#                 changes that are not committed
#
# Example: make cooling-runs COOLING_JOBS=4
# -----------------------------------------------------------------------------

COOLING_JOBS     ?= $(shell nproc 2>/dev/null || echo 1)
COOLING_SUMMARY  := docs/cooling/summary.json
COOLING_FIGURES  := docs/media
COOLING_REVISION  = $(shell git log -1 --format=%h -- cooling)$(if $(shell git status --porcelain -- cooling),-dirty)

.PHONY: cooling-runs

cooling-runs: ## Run the cooling-column scenarios (#9): summary to docs/cooling, figures to docs/media
	uv run python -m cooling.scenario --json $(COOLING_SUMMARY) --figures $(COOLING_FIGURES) \
		--revision $(COOLING_REVISION) --jobs $(COOLING_JOBS)
