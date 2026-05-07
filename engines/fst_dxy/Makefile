CC      = gcc
CFLAGS  = -O3 -march=native -fopenmp -Wall
LDFLAGS = -lz -lm -fopenmp

# -----------------------------------------------------------------------------
# unified_ancestry/engines/fst_dxy/Makefile
#
# This Makefile builds the candidate-agnostic Fst/dXY/theta engine that lives
# in unified_ancestry/ and is reusable by any downstream module.
#
# NOTE: rare_sfs_pairwise and export_q_residual_dosage were previously built
# here but now live in MODULE_5B_inversion_followup/engines/ alongside their
# R consumers (plot_ld_panels.R, 07_plot_rare_sfs_heatmap.R). See
# MODULE_5B_inversion_followup/engines/Makefile.
# -----------------------------------------------------------------------------

all: region_popstats

region_popstats: region_popstats.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f region_popstats

test: all
	@echo "Compile test passed: 1 binary (region_popstats)"

.PHONY: all clean test
