# UNIFIED_ANCESTRY_C_BINARY_ADDITIONS_SPEC — additions to instant_q / region_popstats / hobs_windower

**Status**: SPEC ONLY. Not yet implemented. Awaiting audit.
**Source**: derived from the requirements of `REGIME_ANNOTATION_SPEC.md`,
`FISH_ANCESTRY_SCROLLER_SPEC.md`, and `COPY_ORIGIN_PAINTING_SPEC.md`.
**Position in pipeline**: upstream — adds capabilities the four
annotation/visualisation specs need. Without these additions, the
specs cannot run on real data (only on hand-computed mock TSVs).

## Scope and constraints

This spec covers **additions to the three existing C binaries** in the
unified ancestry module (MODULE_2B):

- `instant_q.cpp` — Engine B C++ fixed-F EM ancestry caller.
- `region_popstats.c` — region-level Hudson Fst, dXY, MI, Hp, θπ,
  Tajima's D, with `--downsample`, `--type 0/1/2`, `--range`,
  `--windows`.
- `hobs_windower.c` — windowed observed heterozygosity.

It does NOT propose new binaries. It does NOT change the core
algorithms of any existing binary. It adds columns, flags, and one
projection mode to support the spec requirements.

The cross-cutting principle: **the C binaries stay annotation-
agnostic.** They consume BEAGLE GLs and group labels; they emit
numeric columns. SnpEff/SIFT/VESM annotation joins, regime
interpretation labels, and visual encoding all happen downstream in
R/Python.

## What needs to be added — by binary

### `region_popstats.c` — adds 5 new outputs + 1 input mode

**Driver**: `REGIME_ANNOTATION_SPEC.md` Layer 4 (underdominance) and
Layer 3 (POD load aggregation) both consume BEAGLE-GL-posterior-based
per-arrangement statistics. The user-facing spec mandates the
naming convention (`HWE_FIS` not bare `FIS`, `arrangement_FST_like`
not bare `FST` between-arrangement) to resolve the I = Inversion
vs I = Inbreeding collision in inversion-context manuscripts.

#### Addition 1 — per-regime HWE statistics

Treat each candidate regime as a single biallelic locus where the
"alleles" are the two arrangements (A = STD, B = INV). Each fish's
karyotype state (HOM_A / HET / HOM_B) is provided by the caller via
a `--karyotype` file. Under BEAGLE GL posteriors, the per-state
counts are expected counts, not hard-call counts.

**New columns** in the per-regime output row:

| column                       | type   | meaning                                          |
|------------------------------|--------|--------------------------------------------------|
| `HWE_Hobs`                   | float  | observed HET frequency = `E[HET_count] / N`      |
| `HWE_Hexp`                   | float  | `2 · p_A · p_B` where p_A = `(2·E[HOM_A] + E[HET]) / (2N)` |
| `HWE_FIS`                    | float  | `1 − HWE_Hobs / HWE_Hexp`                        |
| `HWE_pvalue`                 | float  | chi-square HWE test p-value (or exact for small N) |
| `HWE_HET_excess_deficit`     | str    | `HET_excess` / `HWE_like` / `HET_deficit`        |

**Sign convention** (mandatory):

```
HWE_FIS < 0   →  HET_excess     (Hobs > Hexp)
HWE_FIS = 0   →  HWE_like
HWE_FIS > 0   →  HET_deficit    (Hobs < Hexp)
```

The categorical thresholds for HWE_HET_excess_deficit:

| HWE_FIS range  | label          |
|----------------|----------------|
| ≤ −0.10        | `HET_excess`   |
| −0.10 to +0.10 | `HWE_like`     |
| ≥ +0.10        | `HET_deficit`  |

(Thresholds are defaults; audit chat to confirm against cohort
empirical distribution.)

#### Addition 2 — between-arrangement differentiation, disambiguated name

The existing binary already computes Hudson Fst between arbitrary
group labels (output column pattern: `Fst_<group1>_<group2>`). What
is missing is the **naming convention enforcement** for the
inversion-context use case.

**New flag**: `--output_naming disambiguated`

When this flag is set, the binary emits an additional column:

| column                  | type  | meaning                                           |
|-------------------------|-------|---------------------------------------------------|
| `arrangement_FST_like`  | float | Hudson Fst between HOM_A and HOM_B samples within the regime |

Computed by the same Hudson estimator as existing Fst output; only
the column name differs. The duplicate output (existing `Fst_b1_b3`
or whatever group naming the caller uses + the new
`arrangement_FST_like`) is intentional — it lets downstream code
read whichever name it expects, and lets the manuscript phrase
"FAT_arrangement" or "arrangement_FST_like" without re-running the
binary.

When `--output_naming disambiguated` is set, the binary also adds
a `naming_convention_version` metadata field in the header:
`naming_convention_version=disambiguated_v1`.

#### Addition 3 — expected alt-dosage sum per group

**Driver**: Layer 3a/3b POD load aggregation in
`REGIME_ANNOTATION_SPEC.md`. The C binary must NOT know about
SnpEff / SIFT / VESM annotations. It DOES emit the per-variant,
per-group expected alt-allele dosage that R can then multiply by
a per-variant deleterious-score vector.

**New flag**: `--emit_per_variant_group_dosage`

When set, emits a long-format TSV alongside the standard window
output, with columns:

| column                | type  | meaning                                     |
|-----------------------|-------|---------------------------------------------|
| `chrom`               | str   |                                             |
| `pos_bp`              | int   | variant position                            |
| `ref`                 | str   | reference allele                            |
| `alt`                 | str   | alternative allele                          |
| `group`               | str   | group label as passed via `--karyotype`     |
| `n_samples_in_group`  | int   | size of this group                          |
| `E_alt_dosage_sum`    | float | sum across samples in group of `E[alt_dosage | GL]` |
| `E_alt_freq`          | float | `E_alt_dosage_sum / (2 · n_samples_in_group)` |
| `mean_GL_peak`        | float | mean across samples of `max(P_AA, P_AB, P_BB)` |

R-side, the SnpEff-weighted burden becomes:

```r
burden_per_group <- per_variant_group_dosage %>%
  inner_join(snpeff_scores, by = c("chrom", "pos_bp", "ref", "alt")) %>%
  group_by(group) %>%
  summarise(weighted_burden = sum(E_alt_dosage_sum * deleterious_score))
```

This keeps the C binary annotation-agnostic and lets R own the
annotation-version dependency (SnpEff databases update; the C
binary doesn't need rebuilds when they do).

#### Addition 4 — polarisation column for per-variant POD filter

**Driver**: Layer 3a polarisation filter
`|freq_HOM_STD − freq_HOM_INV| ≥ 0.6`. When
`--emit_per_variant_group_dosage` is set with at least two groups
named exactly `HOM_STD` and `HOM_INV` (or `HOM_A` and `HOM_B` —
configurable via `--polarisation_groups HOM_A,HOM_B`), the binary
also emits:

| column            | type  | meaning                                            |
|-------------------|-------|----------------------------------------------------|
| `polarisation`    | float | `|E_alt_freq(group1) − E_alt_freq(group2)|`        |
| `het_intermediate`| bool  | true iff HET group's E_alt_freq is between the two homozygous groups' E_alt_freq |
| `gl_peak_pass`    | bool  | true iff ≥ 80% of samples in each group have GL peak ≥ 0.85 |

These are pre-computed filter columns so R doesn't have to walk the
data twice.

#### Addition 5 — BEAGLE GL input mode (verify, document)

**This may already exist.** The spec assumes the binary consumes
BEAGLE-phased GLs throughout. The audit chat should:

- Verify `region_popstats.c` accepts BEAGLE GL input (the `.beagle.gz`
  format with `marker`, `allele1`, `allele2`, then triplets of
  `Ind0_AA Ind0_AB Ind0_BB` per sample).
- If it currently only takes hard genotypes or `.geno` format,
  add `--input_format beagle_gl` mode.
- Document the mandatory input format in the binary's `--help`.

Without GL input mode, the entire annotation layer collapses to
hard-call-based estimates, which the user-facing specs explicitly
forbid for break-point regions.

### `instant_q.cpp` — adds `--reference_F` projection mode

**Driver**: `FISH_ANCESTRY_SCROLLER_SPEC.md` label-switching
problem. Per-RF NGSadmix-style runs produce Q matrices whose K
columns are permuted relative to the global run. The scroller spec
documents a downstream alignment pipeline (F-based correlation,
Q-based fallback in flanks, regime-aware neighbour smoothing) that
fixes this post-hoc.

**The better fix is upstream**: anchor the local EM's column
ordering to the global F at initialisation. This is the same
principle as ADMIXTURE's `-P` mode but used as a *constraint on
initial conditions*, not as a fully fixed reference.

#### Addition — `--reference_F <global.fopt>` flag

When set, instant_q:

1. Reads the global F matrix (SNP × K from a previous global run).
2. Restricts to SNPs that intersect the current local run's input.
3. Initialises the local F at those SNPs to the global values
   (with small perturbation for EM stability), and the local F at
   non-overlapping SNPs to a random or uniform initial state.
4. Runs the EM as usual.

The output Q and F columns will be **ordered to match the global K
ordering** with high probability. Local F values at overlapping
SNPs will drift from global values during EM (because biology), but
the column identity is anchored.

This dramatically reduces label-switching at source. The downstream
F-based alignment pipeline (Hungarian matching on K×K correlation
matrix) still runs as a sanity check and emits an alignment-
confidence score, but most RFs will arrive already aligned.

**New output column** in the per-RF metadata:

| column                     | type   | meaning                                       |
|----------------------------|--------|-----------------------------------------------|
| `reference_F_used`         | bool   | whether `--reference_F` was active            |
| `n_anchored_snps`          | int    | SNPs intersecting global F                    |
| `mean_anchored_correlation`| float  | mean correlation of local F to global F at anchored SNPs after EM |

The downstream alignment pipeline can use
`mean_anchored_correlation` as a fast-path: when ≥ 0.95, skip the
K! permutation enumeration entirely.

**Backward compatibility**: when `--reference_F` is absent,
instant_q behaves exactly as before. No existing call sites break.

### `hobs_windower.c` — adds `--groups` partitioning

**Driver**: `FISH_ANCESTRY_SCROLLER_SPEC.md` per-brick
heterozygosity and `REGIME_ANNOTATION_SPEC.md` per-arrangement-group
exposed-burden context. The annotation specs need Hobs computed
separately per arrangement group (HOM_A / HET / HOM_B) at the same
windows, not just at the cohort level.

#### Addition — `--groups <karyotype.tsv>` flag

When set, the binary takes a sample-to-group mapping (same format
as `region_popstats.c`'s `--karyotype`) and emits Hobs separately
per group at each window.

Output columns extend from the current `Hobs` to:

| column                | type   | meaning                                          |
|-----------------------|--------|--------------------------------------------------|
| `Hobs_<group_label>`  | float  | observed HET in the named group                  |
| `n_<group_label>`     | int    | sample count in the named group at this window   |

For example with `--groups karyotype.tsv` where the file defines
HOM_A / HET / HOM_B, the output gains columns:
`Hobs_HOM_A`, `Hobs_HET`, `Hobs_HOM_B`, `n_HOM_A`, `n_HET`,
`n_HOM_B`.

Total cohort Hobs column (`Hobs`) remains for backward compatibility.

## Cross-cutting requirement — BEAGLE GL input format consistency

All three binaries must accept BEAGLE GL input under one consistent
flag name. Recommendation: `--input_format beagle_gl` across all
three. Documented mandatory format:

```
marker    allele1   allele2   Ind0_AA   Ind0_AB   Ind0_BB   Ind1_AA   ...
chr1:1234  A         T         0.95      0.04      0.01      0.02      ...
```

If any binary currently requires a different format, this is a
blocking compatibility issue — the audit chat must verify all three
binaries can read the same upstream BEAGLE output.

## Cross-cutting requirement — naming convention enforcement

The user-facing specs (REGIME_ANNOTATION_SPEC §"Naming convention",
`HANDOFF.md` §11 Critical wording rules) mandate:

- `HWE_FIS` (not bare `FIS`) for the within-cohort heterozygosity
  excess/deficit statistic.
- `arrangement_FST_like` or `FAT_arrangement` (not bare `FST`) for
  between-arrangement differentiation.

These column names are **enforced at the C-binary output layer** —
the strings written into the TSV header. Once they live in the C
binary's output, downstream code can't accidentally rename them
back. The audit chat should:

- grep all three binaries' source for any output of `"FIS"`, `"Fst"`,
  `"F_IS"`, `"F-ST"` as standalone column names (not as part of
  `Fst_<group1>_<group2>`).
- Verify the disambiguated names appear when `--output_naming
  disambiguated` is set.

## What this spec does NOT cover

- New C binaries. The copy-origin painting module (per
  `COPY_ORIGIN_PAINTING_SPEC.md`) needs a BAM-walker for paralogue-
  informative variant counting — that's a different kind of binary
  (mosdepth-style) and out of scope here.
- Brick construction. R/Python orchestration, not a C binary's job.
- The label-switching alignment pipeline itself (`--reference_F`
  reduces but does not eliminate the need for it). The post-hoc
  alignment in the scroller spec still runs.
- Centromere proximity, regime annotation labels, POD scoring, or
  any of Layer 1 / 2 / 2b / 3 / 4 interpretation logic. All
  downstream.
- Visualisation. The atlas pages consume precomputed TSVs.
- SnpEff / SIFT / VESM annotation joins. R-side, downstream of
  Addition 3 (`--emit_per_variant_group_dosage`).

## Audit questions for the next chat

1. **Existence of BEAGLE GL input mode in each binary.** The audit
   chat should compile each binary and run a 10-sample test to
   confirm `.beagle.gz` is consumed correctly by all three. If
   any binary lacks this, it's the highest-priority addition.
2. **Was F_IS removed from the C migration?** UserMemories notes
   that R-side F_IS / Ashman_D dosage functions were removed when
   the C binaries were finalized. The audit chat should verify
   whether F_IS was re-implemented in C (under the old or any name)
   or simply dropped. If dropped, Addition 1 (HWE_FIS) is genuinely
   new C code, not just a rename.
3. **Existing Fst output column naming.** From earlier project
   conversations, the existing C binary emits Hudson Fst with
   column names like `Fst_<group1>_<group2>` (e.g. `Fst_b1_b3`).
   Confirm the exact naming pattern in the current source. The
   `--output_naming disambiguated` flag should add
   `arrangement_FST_like` alongside the existing column, not
   replace it (backward compatibility).
4. **Thresholds for HWE_HET_excess_deficit categorical.** ±0.10 is
   a starting value. Compute the empirical HWE_FIS distribution
   across the 226-sample cohort's neutral regions to decide if
   ±0.10 is appropriate or if a tighter / looser threshold fits.
5. **Performance of `--emit_per_variant_group_dosage`.** A regime
   may contain 10–100k variants × 226 samples × 3 groups → up to
   ~70M output rows for a single LG. Verify the binary can stream
   this output without memory blowup, and that gzip compression is
   used on disk.
6. **`--reference_F` EM initialization.** The recommended approach
   is to seed local F at overlapping SNPs from the global F with a
   small Gaussian perturbation (σ ≈ 0.02) to avoid EM collapse to
   the exact reference. Audit chat should:
   - Verify the perturbation magnitude doesn't bias estimates.
   - Test on a synthetic dataset where the truth Q is known per
     RF and confirm column ordering is recovered.
7. **`--groups` interaction with windowing.** `hobs_windower.c`
   currently emits per-window cohort Hobs. With `--groups`, do
   sample counts vary across windows due to GL coverage gaps?
   Confirm `n_<group>` is per-window, not a constant column.
8. **Per-regime vs per-window output for HWE statistics.** Addition
   1 outputs HWE statistics **at the regime level** (one row per
   regime), not per window. Confirm `region_popstats.c` has a
   per-region output mode in addition to its per-window mode. If
   not, the addition needs a new `--mode regime` flag and a
   `--regions <BED>` input.
9. **Naming-convention enforcement scope.** The mandate is for
   spec-owned columns. Existing published statistics (Hudson Fst,
   Tajima's D, MI) keep their canonical names. The audit chat
   should produce a list of which output columns are renamed under
   `--output_naming disambiguated` and which are unchanged, so the
   intent is unambiguous.
10. **Output schema versioning.** When new columns are added, R-side
    consumers may break if they expect a fixed column count. The
    binary should emit a `schema_version` header line at the top of
    every output file. Audit chat should standardize this across
    all three binaries.
11. **Compilation and test integration on LANTA.** The unified
    ancestry module is built on the LANTA HPC. Audit chat should:
    - Confirm the build script (Makefile or CMake) is updated for
      each addition.
    - Add unit tests with known-truth fixtures (small VCF + known
      HWE_FIS / Fst / Hobs values) for each new output.
    - Add a SLURM launcher template for running the new modes.

## Status checklist for implementation (DO NOT BEGIN YET)

When the audit chat has approved this spec, implementation order:

- [ ] Inventory check: compile all three binaries on LANTA, confirm
      current `--help` output, identify which BEAGLE GL input mode
      each already supports.
- [ ] Addition 5 (BEAGLE GL input mode) first — it's the prerequisite
      for everything else if any binary currently lacks it.
- [ ] `region_popstats.c` Addition 1 (HWE_FIS + family). Unit test
      with hand-computed cohort.
- [ ] `region_popstats.c` Addition 2 (`arrangement_FST_like` under
      `--output_naming disambiguated`). Verify backward compatibility
      with existing pipelines.
- [ ] `region_popstats.c` Addition 3 (`--emit_per_variant_group_dosage`).
      Test memory footprint on a 1k-variant × 226-sample × 3-group
      synthetic.
- [ ] `region_popstats.c` Addition 4 (polarisation / het_intermediate
      / gl_peak_pass pre-filter columns).
- [ ] `hobs_windower.c` `--groups` addition. Unit test.
- [ ] `instant_q.cpp` `--reference_F` projection mode. The most
      algorithmically interesting addition; test on synthetic data
      with known column-permutation truth.
- [ ] Cross-binary schema_version harmonization.
- [ ] Wire all new modes into the unified-ancestry orchestrator
      (Snakefile / SLURM launchers).
- [ ] Update the R-side dispatcher (load_bridge.R / dispatcher) to
      consume the new columns and verify naming convention.
- [ ] Update the wiki and the binary `--help` text for each addition.

## Open questions for Quentin

These are questions the audit chat cannot answer alone; Quentin's
input is needed:

1. **Should `arrangement_FST_like` replace or accompany the existing
   `Fst_<group1>_<group2>` column?** The spec recommends accompany
   (for backward compatibility), but if Quentin prefers a cleaner
   one-name-one-column policy, the binary can do that with a flag
   like `--output_naming disambiguated_only`.

2. **`HOM_STD` / `HOM_INV` vs `HOM_A` / `HOM_B` group labels.** The
   user-facing specs use both interchangeably. The C binaries
   accept whatever labels the karyotype file provides, but
   `--polarisation_groups` needs a canonical default. Recommend
   `HOM_A,HOM_B` because it's polarity-neutral (no implicit "STD
   is the reference"), but Quentin may prefer the more biological
   `HOM_STD,HOM_INV`.

3. **Output file naming**. The new per-variant per-group dosage
   output (Addition 3) is a separate long-format TSV. Naming
   suggestion: `<input_prefix>.per_variant_group_dosage.tsv.gz`.
   Quentin to confirm or override.

4. **`--reference_F` perturbation magnitude.** σ ≈ 0.02 is a
   guess. Quentin's experience with the instant_q EM convergence
   behaviour on this cohort should inform the default.

5. **Does Quentin actually have BEAGLE GLs ready** for all 226
   samples across all 28 LGs? If not, the prerequisite work is to
   complete that pipeline first — every Addition in this spec
   assumes GL input is available.
