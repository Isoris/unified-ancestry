# MODULE: Hobs/HWE Secondary Inversion Confirmation

## One-line summary
Subset-aware, multiscale Hobs/HWE sliding-window framework for secondary inversion confirmation in a structured hatchery dataset.

## Citation credits
- **Hobs sliding-window logic**: Claire Mérot ([angsd_pipeline](https://github.com/clairemerot/angsd_pipeline))
- **ANGSD HWE**: Korneliussen TS, Albrechtsen A, Nielsen R (2014) BMC Bioinformatics
- **ANGSD thetaStat**: single-sample heterozygosity and local theta
- **Patched ANGSD (-doHetFreq)**: [Isoris/angsd_fixed_HWE](https://github.com/Isoris/angsd_fixed_HWE)
- **ngsF-HMM ROH/FROH**: Vieira FG et al. (2016) Bioinformatics

## What this module is NOT
- NOT the primary inversion detection method
- NOT the main heterozygosity estimate for the fish
- NOT a replacement for theta/ROH/FROH
- NOT valid when run on the full pooled hatchery panel without subsetting

## What this module IS
A secondary, group-aware confirmation layer that asks: within ancestry-aware subsets, is there localized heterozygote excess/deficit or HWE distortion at candidate inversion intervals?

## Pipeline

```
Step 1: Build subset BAM lists (from K=8 ancestry clusters ± pruning)
Step 2: ANGSD -doHWE -doHetFreq 1 per subset × chromosome
Step 3: hobs_windower (C binary) → site-level Hobs/F + 7-scale windows
Step 4: Candidate interval overlay + distortion classification
```

## Window scales
5kb, 10kb, 50kb, 100kb, 250kb, 500kb, 1Mb

## Interpretation patterns
1. **High mean F + high median F + many outlier sites** → broad localized genotype distortion, stronger inversion-like support
2. **High mean F but normal median F** → signal driven by a few extreme loci; inspect carefully
3. **Low Hobs across mean and median** → regional heterozygote deficit
4. **Signal disappears after pruning related individuals** → likely family-structure / relatedness confounding
5. **Signal appears only in full pooled panel but not in ancestry-aware subsets** → likely structure / Wahlund effect, not strong inversion support

## Methods text (for manuscript)
"HWE/Hobs-based scans were treated as a secondary, subset-based diagnostic rather than a primary genome-wide inversion detection signal, because strong hatchery structure, relatedness, and potential Wahlund effects can generate apparent heterozygote deficits or excesses unrelated to true local structural polymorphism. Accordingly, HWE/Hobs summaries were computed within ancestry-aware and, where applicable, pruned subsets, and interpreted only in conjunction with independent evidence such as local diversity, ROH, and other structural signals."

## Limitations text (for manuscript)
"HWE/Hobs-derived regional scans are sensitive to sample composition and do not by themselves distinguish inversion polymorphism from other causes of Hardy-Weinberg distortion, including family structure, ancestry mixture, mapping bias, and local paralogy; therefore these summaries were used only as a secondary confirmation layer."
