# -----------------------------------------------------------------------------
# post.mk - animation runs and their rendering
#
#   run-anim    copy the inputs of CASE_DIR into ANIM_DIR, set a dense output
#               (ANIM_END, ANIM_WRITE) and run step ANIM_STEP there
#   anim        render frames, mp4 and GIF of ANIM_DIR into ANIM_OUT
#               (python -m post.animate on the host via uv; needs ffmpeg)
#   readme-gif  anim + copy the GIF to docs/media for the README
#
# Example: make run-anim anim ANIM_STEP=4
#          make anim ANIM_PANELS=flow,depth,froude ANIM_ARGS='--t-max 10'
# -----------------------------------------------------------------------------

ANIM_STEP    ?= 2
ANIM_DIR     ?= case-anim-step$(ANIM_STEP)
ANIM_END     ?= 15
ANIM_WRITE   ?= 0.04
ANIM_PANELS  ?= flow,temperature,depth
ANIM_ARGS    ?=
ANIM_OUT     := $(ANIM_DIR)/postProcessing/anim
README_GIF   := docs/media/dambreak-step$(ANIM_STEP).gif

# Cell centres (0/C) let post.animate recover the grid from the mesh
CENTRES_CMD  := cd $(ANIM_DIR) && postProcess -func writeCellCentres -time 0 2>&1 | tee log.writeCellCentres

.PHONY: run-anim anim readme-gif

run-anim: ## Run step ANIM_STEP in ANIM_DIR with dense output (default: Step 2, 15 s, every 0.04 s)
	bash scripts/prepare_anim_case.sh $(CASE_DIR) $(ANIM_DIR) $(ANIM_END) $(ANIM_WRITE)
	$(MAKE) --no-print-directory run-step$(ANIM_STEP) CASE_DIR=$(ANIM_DIR)

anim: ## Render frames, mp4 and GIF of ANIM_DIR into ANIM_DIR/postProcessing/anim
	@test -f $(ANIM_DIR)/0/.stamp || { echo "Error: $(ANIM_DIR) has no results; run 'make run-anim' first" >&2; exit 1; }
	@if [ ! -e $(ANIM_DIR)/0/C.gz ] && [ ! -e $(ANIM_DIR)/0/C ]; then $(call foam,$(CENTRES_CMD)); fi
	uv run python -m post.animate $(ANIM_DIR) --panels $(ANIM_PANELS) --out $(ANIM_OUT) $(ANIM_ARGS)

readme-gif: anim ## Copy the GIF of ANIM_DIR to docs/media for the README
	mkdir -p docs/media
	cp $(ANIM_OUT)/animation.gif $(README_GIF)
