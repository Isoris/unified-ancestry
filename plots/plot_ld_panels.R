#!/usr/bin/env Rscript

# =============================================================================
# plot_ld_panels.R — Inversion LD Analysis Panel System
#
# Produces two figure types from your mockups:
#
# FIGURE TYPE 1 (Image 1): "What to Show" — 6+3 core/optional panels
#   Panel 1: LD in All Samples                (raw r² heatmap)
#   Panel 2: LD without First-Degree Relatives (pruned r²)
#   Panel 3: LD after Ancestry/Q Correction    (residual r²)
#   Panel 4: HOM_INV Only                      (subset r²)
#   Panel 5: Inversion-Class Contrast          (split: upper=all, lower=all−HOM_INV)
#   Panel 6: Local PCA of HOM_INV             (scatter, colored by family/Q)
#   Optional A: HOM_STD Only
#   Optional B: Raw All − No First-Degree      (difference heatmap)
#   Optional C: Raw All − Q Corrected          (difference heatmap)
#
# FIGURE TYPE 2 (Image 2): "Latent Group Decomposition" — 7 panels
#   Panel A: Pooled LD Heatmap (all samples)
#   Panel B: Per-group LD heatmaps (1 per latent group, tiled)
#   Panel C: Composite RGB overlay (3-group visualization)
#   Panel D: Dominant Group per Cell
#   Panel E: Dominance Strength (Purity)
#   Panel F: Contribution Entropy
#   Panel G: Block Summary Table
#
# INPUTS:
#   --dosage_file <chr.dosage.tsv.gz>   Per-chr dosage matrix (markers × samples)
#   --sites_file <chr.sites.tsv.gz>     Site positions
#   --groups <file>                     sample_id<TAB>group_id (inversion genotypes OR ancestry)
#   --pruned_list <file>                Samples to EXCLUDE for Panel 2 (first-degree relatives)
#   --q_residual_dosage <file>          Q-corrected dosage (for Panel 3). If absent, skip.
#   --candidate_id <id>                 Candidate label
#   --chrom <chr>                       Chromosome
#   --start <bp> --end <bp>             Region
#   --flank <bp>                        Flank size (default: 200000)
#   --thin <N>                          Use every Nth SNP (default: 5 for speed)
#   --outdir <dir>
#   --type 1|2|both                     Which figure type (default: both)
#
# Requires: data.table, ggplot2. Optional: patchwork, viridis, RColorBrewer.
#
# The LD metric is r² = cor(dosage_i, dosage_j)² for SNP pairs i,j.
# For N SNPs in the region, the heatmap is N×N (thinned by --thin).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))
has_viridis   <- suppressWarnings(require(viridis, quietly = TRUE))

# ── CLI ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
dosage_file <- NULL; sites_file <- NULL; groups_file <- NULL
pruned_list_file <- NULL; q_residual_file <- NULL
candidate_id <- "cand"; chrom <- ""; start_bp <- 0; end_bp <- 0
flank_bp <- 200000; thin <- 5; outdir <- "plots/ld_panels"
fig_type <- "both"

i <- 1L
while (i <= length(args)) {
  switch(args[i],
    "--dosage_file"  = { dosage_file <- args[i+1]; i <- i+2L },
    "--sites_file"   = { sites_file <- args[i+1]; i <- i+2L },
    "--groups"       = { groups_file <- args[i+1]; i <- i+2L },
    "--pruned_list"  = { pruned_list_file <- args[i+1]; i <- i+2L },
    "--q_residual_dosage" = { q_residual_file <- args[i+1]; i <- i+2L },
    "--candidate_id" = { candidate_id <- args[i+1]; i <- i+2L },
    "--chrom"        = { chrom <- args[i+1]; i <- i+2L },
    "--start"        = { start_bp <- as.numeric(args[i+1]); i <- i+2L },
    "--end"          = { end_bp <- as.numeric(args[i+1]); i <- i+2L },
    "--flank"        = { flank_bp <- as.numeric(args[i+1]); i <- i+2L },
    "--thin"         = { thin <- as.integer(args[i+1]); i <- i+2L },
    "--outdir"       = { outdir <- args[i+1]; i <- i+2L },
    "--type"         = { fig_type <- args[i+1]; i <- i+2L },
    { i <- i+1L }
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

# ── Color scales ─────────────────────────────────────────────────────────────

# Warm LD scale (yellow-orange-red-black, like Image 2 Panel A)
ld_colors_warm <- colorRampPalette(c(
  "#000004", "#1B0C41", "#4A0C6B", "#781C6D", "#A52C60",
  "#CF4446", "#ED6925", "#FB9B06", "#F7D13D", "#FCFFA4"
))(256)

# Blue-white-red for difference maps (Image 1, Panels 5/B/C)
diff_colors <- colorRampPalette(c(
  "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
))(256)

# Group palette (10 groups, vivid, like Image 2)
GROUP_COLORS <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#FFD700", "#984EA3",
  "#00CED1", "#FF7F00", "#228B22", "#D2691E", "#ADFF2F"
)

# ── Load data ────────────────────────────────────────────────────────────────

message("[LD] Loading dosage...")
dos_dt <- fread(dosage_file)
sites <- fread(sites_file)

# Region filter
reg_start <- max(0, start_bp - flank_bp)
reg_end <- end_bp + flank_bp

if ("pos" %in% names(sites)) {
  keep <- which(sites$pos >= reg_start & sites$pos <= reg_end)
} else if ("position" %in% names(sites)) {
  keep <- which(sites$position >= reg_start & sites$position <= reg_end)
  setnames(sites, "position", "pos")
} else {
  stop("sites file needs 'pos' or 'position' column")
}

# Thin
if (thin > 1) keep <- keep[seq(1, length(keep), by = thin)]
message("[LD] Region: ", reg_start, "-", reg_end, " → ", length(keep), " SNPs (thin=", thin, ")")

sites_sub <- sites[keep]
sample_cols <- setdiff(names(dos_dt), c("marker", "chrom", "pos", "position"))

# Dosage matrix: SNPs × samples
dos_mat <- as.matrix(dos_dt[keep, ..sample_cols])
rownames(dos_mat) <- sites_sub$pos
storage.mode(dos_mat) <- "double"
pos_vec <- sites_sub$pos
n_snps <- nrow(dos_mat)
n_samples <- ncol(dos_mat)
sample_names <- sample_cols

message("[LD] ", n_snps, " SNPs × ", n_samples, " samples")

# ── Load groups ──────────────────────────────────────────────────────────────

grp <- NULL
if (!is.null(groups_file) && file.exists(groups_file)) {
  grp_raw <- fread(groups_file, header = FALSE)
  if (ncol(grp_raw) >= 2) {
    setnames(grp_raw, c("sample_id", "group_id"))
    # Strip BAM extensions
    grp_raw[, sample_id := sub(".*/(.*?)(\\.sorted)?\\.markdup\\.bam$", "\\1",
                                sub("\\.bam$", "", sample_id))]
    grp <- setNames(grp_raw$group_id, grp_raw$sample_id)
  }
}

# ── Load pruned list ──
pruned_samples <- NULL
if (!is.null(pruned_list_file) && file.exists(pruned_list_file)) {
  pruned_samples <- readLines(pruned_list_file)
  pruned_samples <- sub(".*/(.*?)(\\.sorted)?\\.markdup\\.bam$", "\\1",
                         sub("\\.bam$", "", pruned_samples))
  pruned_samples <- pruned_samples[nzchar(pruned_samples)]
}

# ── Load Q-residual dosage ──
dos_residual <- NULL
if (!is.null(q_residual_file) && file.exists(q_residual_file)) {
  res_dt <- fread(q_residual_file)
  res_cols <- intersect(sample_cols, names(res_dt))
  dos_residual <- as.matrix(res_dt[keep, ..res_cols])
  storage.mode(dos_residual) <- "double"
  message("[LD] Q-residual dosage loaded: ", nrow(dos_residual), " × ", ncol(dos_residual))
}

# =============================================================================
# CORE FUNCTION: compute r² matrix
# =============================================================================

compute_r2_matrix <- function(dos, impute_median = TRUE) {
  # dos: SNPs × samples
  if (impute_median) {
    for (mi in seq_len(nrow(dos))) {
      na_idx <- is.na(dos[mi, ])
      if (any(na_idx)) dos[mi, na_idx] <- median(dos[mi, !na_idx], na.rm = TRUE)
    }
  }
  r <- cor(t(dos), use = "pairwise.complete.obs")
  r[!is.finite(r)] <- 0
  r * r  # r²
}

# =============================================================================
# HEATMAP RENDERING (base R — fast, publication quality)
# =============================================================================

# Triangle heatmap: plot only lower triangle
plot_ld_heatmap <- function(r2, pos, title, color_scale = ld_colors_warm,
                             zlim = c(0, 1), inv_start = NULL, inv_end = NULL,
                             out_file = NULL, width = 7, height = 6.5) {
  n <- nrow(r2)
  pos_mb <- pos / 1e6

  if (!is.null(out_file)) pdf(out_file, width = width, height = height)

  par(mar = c(4, 4, 3, 1.5), family = "sans")
  image(pos_mb, pos_mb, r2,
        col = color_scale, zlim = zlim,
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = title, cex.main = 1.1, cex.lab = 0.9, cex.axis = 0.8,
        useRaster = TRUE)

  # Candidate interval lines
  if (!is.null(inv_start) && !is.null(inv_end)) {
    abline(v = c(inv_start, inv_end) / 1e6, col = "#FFFFFF80", lwd = 1.5, lty = 2)
    abline(h = c(inv_start, inv_end) / 1e6, col = "#FFFFFF80", lwd = 1.5, lty = 2)
  }

  if (!is.null(out_file)) dev.off()
}

# Difference heatmap (blue-white-red)
plot_diff_heatmap <- function(diff_mat, pos, title, upper_label = "", lower_label = "",
                               out_file = NULL, width = 7, height = 6.5) {
  n <- nrow(diff_mat)
  pos_mb <- pos / 1e6
  zmax <- quantile(abs(diff_mat[is.finite(diff_mat)]), 0.98, na.rm = TRUE)

  if (!is.null(out_file)) pdf(out_file, width = width, height = height)

  par(mar = c(4, 4, 3, 1.5))
  image(pos_mb, pos_mb, diff_mat,
        col = diff_colors, zlim = c(-zmax, zmax),
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = title, cex.main = 1.0, useRaster = TRUE)

  # Labels for split panels
  if (nzchar(upper_label)) {
    text(min(pos_mb) + diff(range(pos_mb)) * 0.8,
         max(pos_mb) - diff(range(pos_mb)) * 0.1,
         upper_label, col = "grey40", cex = 0.8, font = 2)
  }
  if (nzchar(lower_label)) {
    text(min(pos_mb) + diff(range(pos_mb)) * 0.2,
         min(pos_mb) + diff(range(pos_mb)) * 0.1,
         lower_label, col = "grey40", cex = 0.8, font = 2)
  }

  if (!is.null(out_file)) dev.off()
}

# Split heatmap: upper triangle = mat_A, lower triangle = mat_B
plot_split_heatmap <- function(mat_upper, mat_lower, pos, title,
                                upper_label = "All samples", lower_label = "Contrast",
                                col_upper = ld_colors_warm, col_lower = diff_colors,
                                out_file = NULL, width = 7, height = 6.5) {
  n <- nrow(mat_upper)
  pos_mb <- pos / 1e6

  # Combine: upper tri from mat_upper, lower tri from mat_lower
  combined <- matrix(0, n, n)
  combined[upper.tri(combined, diag = TRUE)] <- mat_upper[upper.tri(mat_upper, diag = TRUE)]
  combined[lower.tri(combined)] <- mat_lower[lower.tri(mat_lower)]

  if (!is.null(out_file)) pdf(out_file, width = width, height = height)

  # Plot upper and lower separately for different color scales
  par(mar = c(4, 4, 3, 1.5))

  # Background: lower triangle
  lower_only <- matrix(NA, n, n)
  lower_only[lower.tri(lower_only)] <- mat_lower[lower.tri(mat_lower)]
  zmax_l <- quantile(abs(mat_lower[is.finite(mat_lower)]), 0.98, na.rm = TRUE)
  image(pos_mb, pos_mb, lower_only,
        col = col_lower, zlim = c(-zmax_l, zmax_l),
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = title, cex.main = 1.0, useRaster = TRUE)

  # Overlay: upper triangle
  upper_only <- matrix(NA, n, n)
  upper_only[upper.tri(upper_only, diag = TRUE)] <- mat_upper[upper.tri(mat_upper, diag = TRUE)]
  par(new = TRUE)
  image(pos_mb, pos_mb, upper_only,
        col = col_upper, zlim = c(0, 1),
        xlab = "", ylab = "", main = "", axes = FALSE, useRaster = TRUE)

  # Diagonal line
  abline(0, 1, col = "white", lwd = 1)

  # Labels
  text(min(pos_mb) + diff(range(pos_mb)) * 0.75,
       max(pos_mb) - diff(range(pos_mb)) * 0.05,
       upper_label, col = "white", cex = 0.8, font = 2)
  text(min(pos_mb) + diff(range(pos_mb)) * 0.2,
       min(pos_mb) + diff(range(pos_mb)) * 0.05,
       lower_label, col = "grey30", cex = 0.8, font = 2)

  if (!is.null(out_file)) dev.off()
}

# =============================================================================
# FIGURE TYPE 1: "What to Show" — Q-correction panels
# =============================================================================

build_figure_type1 <- function() {
  message("[LD] Building Figure Type 1: Q-correction panels...")

  # Panel 1: All samples
  r2_all <- compute_r2_matrix(dos_mat)
  plot_ld_heatmap(r2_all, pos_vec,
    title = paste0("1. LD in All Samples (", candidate_id, ")"),
    inv_start = start_bp, inv_end = end_bp,
    out_file = file.path(outdir, paste0("panel1_ld_all_", candidate_id, ".pdf")))
  message("[LD]   Panel 1 done")

  # Panel 2: No first-degree relatives
  if (!is.null(pruned_samples)) {
    keep_cols <- which(!sample_names %in% pruned_samples)
    if (length(keep_cols) > 10) {
      r2_pruned <- compute_r2_matrix(dos_mat[, keep_cols, drop = FALSE])
      plot_ld_heatmap(r2_pruned, pos_vec,
        title = paste0("2. LD without First-Degree Relatives (n=", length(keep_cols), ")"),
        inv_start = start_bp, inv_end = end_bp,
        out_file = file.path(outdir, paste0("panel2_ld_no_firstdeg_", candidate_id, ".pdf")))

      # Optional B: Difference (all − pruned)
      diff_rel <- r2_all - r2_pruned
      plot_diff_heatmap(diff_rel, pos_vec,
        title = "B. Raw All − No First-Degree (family contribution)",
        out_file = file.path(outdir, paste0("panelB_diff_relatedness_", candidate_id, ".pdf")))
    }
  }
  message("[LD]   Panel 2 done")

  # Panel 3: Q-corrected
  if (!is.null(dos_residual)) {
    r2_qcorr <- compute_r2_matrix(dos_residual)
    plot_ld_heatmap(r2_qcorr, pos_vec,
      title = "3. LD after Ancestry/Q Correction",
      inv_start = start_bp, inv_end = end_bp,
      out_file = file.path(outdir, paste0("panel3_ld_q_corrected_", candidate_id, ".pdf")))

    # Optional C: Difference (all − Q-corrected)
    diff_q <- r2_all - r2_qcorr
    plot_diff_heatmap(diff_q, pos_vec,
      title = "C. Raw All − Q Corrected (ancestry contribution)",
      out_file = file.path(outdir, paste0("panelC_diff_q_correction_", candidate_id, ".pdf")))
  }
  message("[LD]   Panel 3 done")

  # Panels 4, 5, A: Require inversion genotype groups
  if (!is.null(grp)) {
    grp_matched <- grp[sample_names]
    grp_matched[is.na(grp_matched)] <- "UNKNOWN"

    # Identify inversion classes
    inv_labels <- c("HOM_INV", "HOMO_2", "INV", "FULL_B")
    std_labels <- c("HOM_STD", "HOMO_1", "STD", "FULL_A")
    het_labels <- c("HET", "HALF")

    idx_inv <- which(grp_matched %in% inv_labels)
    idx_std <- which(grp_matched %in% std_labels)
    idx_het <- which(grp_matched %in% het_labels)

    # Panel 4: HOM_INV only
    if (length(idx_inv) >= 5) {
      r2_inv <- compute_r2_matrix(dos_mat[, idx_inv, drop = FALSE])
      plot_ld_heatmap(r2_inv, pos_vec,
        title = paste0("4. HOM_INV Only (n=", length(idx_inv), ")"),
        inv_start = start_bp, inv_end = end_bp,
        out_file = file.path(outdir, paste0("panel4_ld_hom_inv_", candidate_id, ".pdf")))
    }
    message("[LD]   Panel 4 done")

    # Optional A: HOM_STD only
    if (length(idx_std) >= 5) {
      r2_std <- compute_r2_matrix(dos_mat[, idx_std, drop = FALSE])
      plot_ld_heatmap(r2_std, pos_vec,
        title = paste0("A. HOM_STD Only (n=", length(idx_std), ")"),
        inv_start = start_bp, inv_end = end_bp,
        out_file = file.path(outdir, paste0("panelA_ld_hom_std_", candidate_id, ".pdf")))
    }

    # Panel 5: Split — upper = all, lower = all − HOM_INV
    if (length(idx_inv) >= 5) {
      non_inv <- setdiff(seq_len(n_samples), idx_inv)
      if (length(non_inv) >= 5) {
        r2_no_inv <- compute_r2_matrix(dos_mat[, non_inv, drop = FALSE])
        diff_removal <- r2_all - r2_no_inv  # what HOM_INV contributed

        plot_split_heatmap(r2_all, diff_removal, pos_vec,
          title = paste0("5. Inversion-Class Contrast (", candidate_id, ")"),
          upper_label = "All samples",
          lower_label = paste0("All − (no HOM_INV)"),
          out_file = file.path(outdir, paste0("panel5_contrast_removal_", candidate_id, ".pdf")))
      }
    }
    message("[LD]   Panel 5 done")

    # Panel 6: Local PCA of HOM_INV colored by family/Q
    if (length(idx_inv) >= 10) {
      inv_dos <- dos_mat[, idx_inv, drop = FALSE]
      inv_dos_clean <- inv_dos
      for (mi in seq_len(nrow(inv_dos_clean))) {
        na_idx <- is.na(inv_dos_clean[mi, ])
        if (any(na_idx)) inv_dos_clean[mi, na_idx] <- median(inv_dos_clean[mi, !na_idx], na.rm = TRUE)
      }

      # PCA on SNP dosage
      pca <- prcomp(t(inv_dos_clean), center = TRUE, scale. = FALSE)
      pc_dt <- data.table(
        sample = sample_names[idx_inv],
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2]
      )

      # Color by ancestry cluster if available
      if (!is.null(grp)) {
        pc_dt[, ancestry := grp[sample]]
        pc_dt[is.na(ancestry), ancestry := "unknown"]
      } else {
        pc_dt[, ancestry := "all"]
      }

      anc_pal <- setNames(GROUP_COLORS[seq_len(length(unique(pc_dt$ancestry)))],
                           unique(pc_dt$ancestry))

      p6 <- ggplot(pc_dt, aes(x = PC1, y = PC2, color = ancestry)) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(values = anc_pal) +
        theme_bw(base_size = 10) +
        labs(title = paste0("6. Local PCA of HOM_INV (", candidate_id, ")"),
             subtitle = "Colored by family/ancestry cluster",
             color = "Group")

      ggsave(file.path(outdir, paste0("panel6_pca_hom_inv_", candidate_id, ".pdf")),
             p6, width = 7, height = 6)
    }
    message("[LD]   Panel 6 done")
  }

  message("[LD] Figure Type 1 complete")
}

# =============================================================================
# FIGURE TYPE 2: "Latent Group Decomposition"
# =============================================================================

build_figure_type2 <- function() {
  if (is.null(grp)) {
    message("[LD] Figure Type 2 requires --groups file"); return(invisible())
  }

  message("[LD] Building Figure Type 2: Latent group decomposition...")

  grp_matched <- grp[sample_names]
  grp_matched[is.na(grp_matched)] <- "UNKNOWN"
  groups_present <- sort(unique(grp_matched[grp_matched != "UNKNOWN"]))
  n_groups <- length(groups_present)

  message("[LD]   ", n_groups, " groups: ", paste(groups_present, collapse = ", "))

  # Panel A: Pooled LD
  r2_all <- compute_r2_matrix(dos_mat)
  plot_ld_heatmap(r2_all, pos_vec,
    title = "A. Pooled LD Heatmap (All Samples)",
    inv_start = start_bp, inv_end = end_bp,
    out_file = file.path(outdir, paste0("decomp_A_pooled_", candidate_id, ".pdf")))
  message("[LD]   Panel A done")

  # Panel B: Per-group LD heatmaps
  r2_per_group <- list()
  for (gi in seq_along(groups_present)) {
    g <- groups_present[gi]
    idx <- which(grp_matched == g)
    if (length(idx) < 3) {
      r2_per_group[[g]] <- matrix(0, n_snps, n_snps)
      next
    }
    r2_g <- compute_r2_matrix(dos_mat[, idx, drop = FALSE])
    r2_per_group[[g]] <- r2_g

    plot_ld_heatmap(r2_g, pos_vec,
      title = paste0("Group ", g, " (n=", length(idx), ")"),
      inv_start = start_bp, inv_end = end_bp,
      out_file = file.path(outdir, paste0("decomp_B_group_", g, "_", candidate_id, ".pdf")),
      width = 5, height = 4.5)
  }
  message("[LD]   Panel B done (", n_groups, " groups)")

  # Panel D: Dominant Group per Cell
  # For each SNP pair, which group has the highest r²?
  dom_group <- matrix(0L, n_snps, n_snps)
  dom_strength <- matrix(0, n_snps, n_snps)
  entropy_mat <- matrix(0, n_snps, n_snps)

  for (si in seq_len(n_snps)) {
    for (sj in si:n_snps) {
      vals <- vapply(groups_present, function(g) r2_per_group[[g]][si, sj], numeric(1))
      total <- sum(vals, na.rm = TRUE)
      if (total > 1e-10) {
        dom_group[si, sj] <- dom_group[sj, si] <- which.max(vals)
        dom_strength[si, sj] <- dom_strength[sj, si] <- max(vals) / total
        # Shannon entropy of proportions
        props <- vals / total
        props <- props[props > 0]
        entropy_mat[si, sj] <- entropy_mat[sj, si] <- -sum(props * log2(props))
      }
    }
  }

  # Plot D: Dominant group (categorical heatmap)
  pos_mb <- pos_vec / 1e6
  dom_group_f <- dom_group
  dom_group_f[dom_group_f == 0] <- NA

  pdf(file.path(outdir, paste0("decomp_D_dominant_group_", candidate_id, ".pdf")),
      width = 7, height = 6.5)
  par(mar = c(4, 4, 3, 5))
  pal_dom <- GROUP_COLORS[seq_len(n_groups)]
  image(pos_mb, pos_mb, dom_group_f,
        col = pal_dom, zlim = c(0.5, n_groups + 0.5),
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = "D. Dominant Group per Cell", useRaster = TRUE)
  if (!is.null(start_bp)) {
    abline(v = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
    abline(h = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
  }
  # Legend
  legend("topright", legend = groups_present, fill = pal_dom,
         cex = 0.6, bg = "white", title = "Group")
  dev.off()
  message("[LD]   Panel D done")

  # Plot E: Dominance Strength (Purity)
  purity_colors <- colorRampPalette(c("#000004", "#3B0F70", "#8C2981",
                                       "#DE4968", "#FE9F6D", "#FCFDBF"))(256)
  pdf(file.path(outdir, paste0("decomp_E_purity_", candidate_id, ".pdf")),
      width = 7, height = 6.5)
  par(mar = c(4, 4, 3, 1.5))
  image(pos_mb, pos_mb, dom_strength,
        col = purity_colors, zlim = c(0, 1),
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = "E. Dominance Strength (Purity)", useRaster = TRUE)
  if (!is.null(start_bp)) {
    abline(v = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
    abline(h = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
  }
  dev.off()
  message("[LD]   Panel E done")

  # Plot F: Contribution Entropy
  entropy_colors <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6",
                                        "#BDD7E7", "#EFF3FF",
                                        "#FEE5D9", "#FCAE91", "#FB6A4A",
                                        "#CB181D"))(256)
  max_entropy <- log2(n_groups)
  pdf(file.path(outdir, paste0("decomp_F_entropy_", candidate_id, ".pdf")),
      width = 7, height = 6.5)
  par(mar = c(4, 4, 3, 1.5))
  image(pos_mb, pos_mb, entropy_mat,
        col = entropy_colors, zlim = c(0, max_entropy),
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = "F. Contribution Entropy", useRaster = TRUE)
  if (!is.null(start_bp)) {
    abline(v = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
    abline(h = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
  }
  dev.off()
  message("[LD]   Panel F done")

  # Panel G: Block Summary Table
  # Bin SNPs into blocks by position and summarize dominant group
  block_size_bp <- (end_bp - start_bp) / 8  # ~8 blocks across interval
  block_starts <- seq(start_bp, end_bp, by = block_size_bp)

  block_summary <- data.table(
    block_start_Mb = round(block_starts / 1e6, 2),
    block_end_Mb = round(pmin(block_starts + block_size_bp, end_bp) / 1e6, 2)
  )

  block_summary[, dominant_group := vapply(seq_len(.N), function(bi) {
    bs <- block_starts[bi]; be <- min(bs + block_size_bp, end_bp)
    snp_idx <- which(pos_vec >= bs & pos_vec <= be)
    if (length(snp_idx) < 2) return("N/A")

    # Mean r² per group in this block
    grp_means <- vapply(seq_along(groups_present), function(gi) {
      sub <- r2_per_group[[groups_present[gi]]][snp_idx, snp_idx]
      mean(sub[upper.tri(sub)], na.rm = TRUE)
    }, numeric(1))
    groups_present[which.max(grp_means)]
  }, character(1))]

  block_summary[, purity := vapply(seq_len(.N), function(bi) {
    bs <- block_starts[bi]; be <- min(bs + block_size_bp, end_bp)
    snp_idx <- which(pos_vec >= bs & pos_vec <= be)
    if (length(snp_idx) < 2) return(0)
    sub_pur <- dom_strength[snp_idx, snp_idx]
    round(mean(sub_pur[upper.tri(sub_pur)], na.rm = TRUE), 2)
  }, numeric(1))]

  fwrite(block_summary,
         file.path(outdir, paste0("decomp_G_block_summary_", candidate_id, ".tsv")),
         sep = "\t")
  message("[LD]   Panel G done")

  # Panel C: RGB Composite (3-group overlay)
  # Assign first 3 groups to R, G, B channels
  if (n_groups >= 3) {
    message("[LD]   Panel C: RGB composite...")

    r_mat <- r2_per_group[[groups_present[1]]]
    g_mat <- r2_per_group[[groups_present[2]]]
    b_mat <- r2_per_group[[groups_present[3]]]

    # Normalize each to [0,1]
    norm01 <- function(m) {
      mx <- quantile(m[is.finite(m)], 0.98, na.rm = TRUE)
      if (mx <= 0) mx <- 1
      pmin(m / mx, 1)
    }
    r_n <- norm01(r_mat); g_n <- norm01(g_mat); b_n <- norm01(b_mat)

    # Write as PNG using R's png device with rgb() coloring
    png_path <- file.path(outdir, paste0("decomp_C_rgb_composite_", candidate_id, ".png"))
    png(png_path, width = 800, height = 800, res = 150)
    par(mar = c(4, 4, 3, 1))
    plot(NA, xlim = range(pos_mb), ylim = range(pos_mb),
         xlab = "Position (Mb)", ylab = "Position (Mb)",
         main = paste0("C. RGB Composite: R=", groups_present[1],
                        " G=", groups_present[2], " B=", groups_present[3]),
         cex.main = 0.9, asp = 1)
    # Rasterize: build color matrix
    cols <- matrix(rgb(r_n, g_n, b_n), nrow = n_snps)
    rasterImage(as.raster(cols),
                min(pos_mb), min(pos_mb), max(pos_mb), max(pos_mb),
                interpolate = FALSE)
    if (!is.null(start_bp)) {
      abline(v = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
      abline(h = c(start_bp, end_bp) / 1e6, col = "white", lwd = 1, lty = 2)
    }
    dev.off()
    message("[LD]   Panel C done")
  }

  message("[LD] Figure Type 2 complete")
}

# =============================================================================
# DISPATCH
# =============================================================================

if (fig_type %in% c("1", "both")) build_figure_type1()
if (fig_type %in% c("2", "both")) build_figure_type2()

message("[LD] All done. Output: ", outdir)
