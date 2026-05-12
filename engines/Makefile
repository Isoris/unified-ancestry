# =============================================================================
# engines/Makefile — Build all C/C++ engines for unified_ancestry
#
# Targets:
#   make             # builds all 5 binaries
#   make instant_q
#   make region_popstats
#   make rare_sfs_pairwise
#   make export_q_residual_dosage
#   make hobs_windower
#   make clean       # removes all binaries
#   make install     # installs binaries to $(PREFIX)/bin (default: ~/.local)
# =============================================================================

CXX      = g++
CC       = gcc
CXXFLAGS = -O3 -march=native -fopenmp -std=c++17 -Wall
CFLAGS   = -O3 -march=native -fopenmp -Wall
LDFLAGS  = -lz -lm -fopenmp

PREFIX  ?= $(HOME)/.local

CPP_BINS = instant_q
C_BINS   = region_popstats rare_sfs_pairwise export_q_residual_dosage hobs_windower \
           codon_stats pi_NS fdM region_test ultrabootstrap fst_outlier_scan xpehh \
           homolog_index_build homolog_index_query homolog_atlas_server
ALL_BINS = $(CPP_BINS) $(C_BINS)

all: $(ALL_BINS)

instant_q: instant_q.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

region_popstats: region_popstats.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

rare_sfs_pairwise: rare_sfs_pairwise.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

export_q_residual_dosage: export_q_residual_dosage.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

hobs_windower: hobs_windower.c
	$(CC) -O3 -o $@ $< -lz -lm

codon_stats: codon_stats.c
	$(CC) $(CFLAGS) -o $@ $< -lm -fopenmp

pi_NS: pi_NS.c
	$(CC) $(CFLAGS) -o $@ $< -lm -fopenmp

fdM: fdM.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

region_test: region_test.c
	$(CC) $(CFLAGS) -o $@ $< -lm -fopenmp

ultrabootstrap: ultrabootstrap.c
	$(CC) $(CFLAGS) -o $@ $< -lm -fopenmp

fst_outlier_scan: fst_outlier_scan.c
	$(CC) -O3 -march=native -Wall -o $@ $< -lm

xpehh: xpehh.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

homolog_index_build: homolog_index_build.c homolog_index.h
	$(CC) -O3 -march=native -Wall -o $@ $<

homolog_index_query: homolog_index_query.c homolog_index.h
	$(CC) -O3 -march=native -Wall -o $@ $<

homolog_atlas_server: homolog_atlas_server.c homolog_index.h
	$(CC) -O3 -march=native -Wall -o $@ $< -lpthread

install: $(ALL_BINS)
	mkdir -p $(PREFIX)/bin
	cp $(ALL_BINS) $(PREFIX)/bin/

clean:
	rm -f $(ALL_BINS)

test: all
	@echo "Compile test passed: $(words $(ALL_BINS)) binaries built"
	@for b in $(ALL_BINS); do echo "  - $$b"; done

.PHONY: all install clean test
