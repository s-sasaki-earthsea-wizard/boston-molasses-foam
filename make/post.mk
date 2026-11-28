# -----------------------------------------------------------------------------
# post.mk - animation runs and their rendering
#
#   run-anim    copy the inputs of CASE_DIR into ANIM_DIR, set a dense output
#               (ANIM_END, ANIM_WRITE) and run step ANIM_STEP there
#   anim        render frames, mp4 and GIF of ANIM_DIR into ANIM_OUT
#               (python -m post.animate on the host via uv; needs ffmpeg)
#   readme-gif  anim + copy the GIF to docs/media for the README
#   metrics     front, volume and bounds of METRICS_DIRS (python -m post.metrics)
#   readme-fronts  plot the fronts of the README reference runs to docs/media
#               (fails before any output unless every run finished)
#
# Example: make run-anim anim ANIM_STEP=4
#          make anim ANIM_PANELS=flow,depth,froude ANIM_ARGS='--t-max 10'
#          make metrics METRICS_DIRS="case-step4-low case-step4-high"
# -----------------------------------------------------------------------------

ANIM_STEP    ?= 2
ANIM_DIR     ?= case-anim-step$(ANIM_STEP)
ANIM_END     ?= 15
ANIM_WRITE   ?= 0.04
ANIM_PANELS  ?= flow,temperature,depth
ANIM_ARGS    ?=
ANIM_OUT     := $(ANIM_DIR)/postProcessing/anim
README_GIF   := docs/media/dambreak-step$(ANIM_STEP).gif

METRICS_DIRS      ?= $(CASE_DIR)
METRICS_ARGS      ?=
README_FRONTS     := docs/media/fronts-step4.svg
README_FRONT_DIRS ?= case-step1 case-step4-low case-step4-central case-step4-high

# Cell centres (0/C) let post/ recover the grid from the mesh
# $(call centres_cmd,<case dir>)
centres_cmd = cd $(1) && postProcess -func writeCellCentres -time 0 2>&1 | tee log.writeCellCentres
# $(call need_centres,<case dir>): write 0/C unless it exists
need_centres = $(if $(wildcard $(1)/0/C.gz $(1)/0/C),,$(call foam,$(call centres_cmd,$(1)));)

.PHONY: run-anim anim readme-gif metrics readme-fronts

run-anim: ## Run step ANIM_STEP in ANIM_DIR with dense output (default: Step 2, 15 s, every 0.04 s)
	bash scripts/prepare_anim_case.sh $(CASE_DIR) $(ANIM_DIR) $(ANIM_END) $(ANIM_WRITE)
	$(MAKE) --no-print-directory run-step$(ANIM_STEP) CASE_DIR=$(ANIM_DIR)

anim: ## Render frames, mp4 and GIF of ANIM_DIR into ANIM_DIR/postProcessing/anim
	@test -f $(ANIM_DIR)/0/.stamp || { echo "Error: $(ANIM_DIR) has no results; run 'make run-anim' first" >&2; exit 1; }
	$(call need_centres,$(ANIM_DIR))
	uv run python -m post.animate $(ANIM_DIR) --panels $(ANIM_PANELS) --out $(ANIM_OUT) $(ANIM_ARGS)

readme-gif: anim ## Copy the GIF of ANIM_DIR to docs/media for the README
	mkdir -p docs/media
	cp $(ANIM_OUT)/animation.gif $(README_GIF)

metrics: ## Front, volume and bounds of METRICS_DIRS (default CASE_DIR; host: uv)
	@for d in $(METRICS_DIRS); do test -f $$d/0/.stamp || { echo "Error: $$d has no results" >&2; exit 1; }; done
	$(foreach d,$(METRICS_DIRS),$(call need_centres,$(d)))
	uv run python -m post.metrics $(METRICS_DIRS) $(METRICS_ARGS)

readme-fronts: ## Plot the fronts of README_FRONT_DIRS to docs/media/fronts-step4.svg (finished runs only)
	$(MAKE) --no-print-directory metrics METRICS_DIRS="$(README_FRONT_DIRS)" METRICS_ARGS="--require-finished --plot $(README_FRONTS)"
