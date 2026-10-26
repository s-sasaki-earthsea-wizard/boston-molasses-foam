# -----------------------------------------------------------------------------
# solver.mk - build the custom solver interThermalFoam with wmake
#
# The binary is written to $(BUILD_DIR)/bin (FOAM_USER_APPBIN is redirected
# there by docker/Dockerfile), wmake objects to $(SOLVER_DIR)/Make/<options>.
# -----------------------------------------------------------------------------

SOLVER_BIN     := $(BUILD_DIR)/bin/$(SOLVER_NAME)
SOLVER_SOURCES := $(wildcard $(SOLVER_DIR)/*.C $(SOLVER_DIR)/*.H) \
                  $(SOLVER_DIR)/Make/files $(SOLVER_DIR)/Make/options

.PHONY: build-solver

build-solver: ## Compile interThermalFoam with wmake (in the container)
	$(call foam,wmake $(SOLVER_DIR))

# Prerequisite of the run targets: rebuild only when sources changed
$(SOLVER_BIN): $(SOLVER_SOURCES)
	$(call foam,wmake $(SOLVER_DIR))
