#!/usr/bin/env Rscript

# =============================================================================
# 07_plot_rare_sfs_heatmap.R — Pairwise Rare Allele Sharing Heatmaps
#
# Bergström et al. (2020) Science Fig. 2 style:
# One heatmap per MAC bin (doubleton, tripleton, quadrupleton, ...),
# showing pairwise sharing counts. Rows/columns = samples, ordered
# by group assignment with black/white group sidebar.
#
# Input:
#   --sfs_dir <dir>        Directory with .bin2.tsv, .bin3.tsv, etc.
#                          (output of rare_sfs_pairwise)
#   --prefix <prefix>      File prefix (default: rare_sfs)
#   --groups <file>        sample_id<tab>group_id (same as C binary input)
#   --bins 2,3,4,5         Which bins to plot
#   --cluster              Hierarchical clustering within groups
#   --log_scale            Log10 transform counts (like Bergström)
#   --outdir <dir>         Output directory
#
# Uses base R heatmap or ComplexHeatmap if available.
#
# Citation: Bergström A et al. (2020) Science 367:eaay5012
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

has_CH <- suppressWarnings(require(ComplexHeatmap, quietly = TRUE))
has_circlize <- suppressWarnings(require(circlize, quietly = TRUE))

# ── CLI ──
args <- commandArgs(trailingOnly = TRUE)
sfs_dir <- "."
prefix <- "rare_sfs"
groups_file <- NULL
bins_str <- "2,3,4,5"
do_cluster <- TRUE
log_scale <- TRUE
outdir <- "plots/rare_sfs"

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--sfs_dir"   = { sfs_dir <- args[i+1]; i <- i+2 },
    "--prefix"    = { prefix <- args[i+1]; i <- i+2 },
    "--groups"    = { groups_file <- args[i+1]; i <- i+2 },
    "--bins"      = { bins_str <- args[i+1]; i <- i+2 },
    "--cluster"   = { do_cluster <- TRUE; i <- i+1 },
    "--no_cluster"= { do_cluster <- FALSE; i <- i+1 },
    "--log_scale" = { log_scale <- TRUE; i <- i+1 },
    "--linear"    = { log_scale <- FALSE; i <- i+1 },
    "--outdir"    = { outdir <- args[i+1]; i <- i+2 },
    { i <- i+1 }
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
bins <- as.integer(strsplit(bins_str, ",")[[1]])

# ── Load groups ──
grp_dt <- NULL
if (!is.null(groups_file) && file.exists(groups_file)) {
  grp_dt <- fread(groups_file, header = FALSE)
  if (ncol(grp_dt) >= 2) {
    setnames(grp_dt, c("sample_id", "group_id"))
    # Strip BAM extensions
    grp_dt[, sample_id := sub(".*/(.*?)(\\.sorted)?\\.markdup\\.bam$", "\\1",
                               sub("\\.bam$", "", sample_id))]
  } else {
    grp_dt <- NULL
  }
}

# ── Group colors (alternating dark/light for Bergström sidebar look) ──
make_group_palette <- function(groups) {
  ug <- unique(groups)
  n <- length(ug)
  # Alternating grayscale + color accents
  base_colors <- c(
    "#2166AC", "#D6604D", "#4DAF4A", "#FF7F00", "#984EA3",
    "#E7298A", "#66C2A5", "#E6AB02", "#A6761D", "#666666",
    "#1B9E77", "#7570B3", "#E41A1C", "#377EB8", "#4DAF4A",
    "#F781BF", "#A65628", "#FFFF33", "#999999", "#B3DE69"
  )
  setNames(base_colors[seq_len(n)], ug)
}

# ── Bin labels ──
bin_labels <- c("2" = "Doubleton (MAC=2)", "3" = "Tripleton (MAC=3)",
                "4" = "Quadrupleton (MAC=4)", "5" = "Quintupleton (MAC=5)",
                "6" = "MAC=6", "7" = "MAC=7", "8" = "MAC=8", "9" = "MAC=9",
                "10" = "MAC=10")

# ── Process each bin ──
for (b in bins) {
  f_sample <- file.path(sfs_dir, paste0(prefix, ".bin", b, ".tsv"))
  f_group <- file.path(sfs_dir, paste0(prefix, ".group_summary.bin", b, ".tsv"))

  if (!file.exists(f_sample)) {
    message("[rare_sfs_plot] Missing: ", f_sample, " — skip bin ", b)
    next
  }

  # ── Load sample-level matrix ──
  dt <- fread(f_sample)
  sample_names <- dt$sample
  mat <- as.matrix(dt[, -1])
  rownames(mat) <- colnames(mat) <- sample_names

  if (log_scale) {
    mat <- log10(mat + 1)
  }

  label <- bin_labels[as.character(b)]
  if (is.na(label)) label <- paste0("MAC=", b)

  # ── Determine sample order ──
  sample_order <- sample_names
  group_vec <- NULL

  if (!is.null(grp_dt)) {
    # Match groups to matrix rows
    grp_matched <- grp_dt[match(sample_names, grp_dt$sample_id)]
    group_vec <- grp_matched$group_id
    group_vec[is.na(group_vec)] <- "ungrouped"
    names(group_vec) <- sample_names

    # Order: by group, then optionally cluster within group
    if (do_cluster) {
      new_order <- character(0)
      for (g in unique(group_vec[!is.na(group_vec)])) {
        idx <- which(group_vec == g)
        if (length(idx) > 2) {
          sub_mat <- mat[idx, idx]
          d <- as.dist(max(sub_mat) - sub_mat + 1)
          hc <- hclust(d, method = "average")
          new_order <- c(new_order, sample_names[idx[hc$order]])
        } else {
          new_order <- c(new_order, sample_names[idx])
        }
      }
      sample_order <- new_order
    } else {
      sample_order <- sample_names[order(group_vec)]
    }
  } else if (do_cluster) {
    d <- as.dist(max(mat) - mat + 1)
    hc <- hclust(d, method = "average")
    sample_order <- sample_names[hc$order]
  }

  mat_ordered <- mat[sample_order, sample_order]

  # ── Plot with ComplexHeatmap if available ──
  if (has_CH && has_circlize) {
    message("[rare_sfs_plot] ComplexHeatmap: bin ", b)

    # Color function
    if (log_scale) {
      col_fun <- circlize::colorRamp2(
        c(0, quantile(mat_ordered[mat_ordered > 0], 0.5, na.rm = TRUE),
          max(mat_ordered, na.rm = TRUE)),
        c("#000033", "#2166AC", "#FDDBC7")
      )
    } else {
      col_fun <- circlize::colorRamp2(
        c(0, median(mat_ordered, na.rm = TRUE), max(mat_ordered, na.rm = TRUE)),
        c("#000033", "#2166AC", "#FDDBC7")
      )
    }

    # Group sidebar annotation
    ha <- NULL
    if (!is.null(group_vec)) {
      grp_pal <- make_group_palette(group_vec)
      ordered_groups <- group_vec[sample_order]
      ha <- ComplexHeatmap::HeatmapAnnotation(
        Group = ordered_groups,
        col = list(Group = grp_pal),
        show_legend = TRUE,
        annotation_name_side = "left",
        simple_anno_size = unit(4, "mm")
      )
    }

    ht <- ComplexHeatmap::Heatmap(
      mat_ordered,
      name = if (log_scale) "log10(count+1)" else "count",
      col = col_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = nrow(mat_ordered) <= 50,
      show_column_names = nrow(mat_ordered) <= 50,
      row_names_gp = grid::gpar(fontsize = 5),
      column_names_gp = grid::gpar(fontsize = 5),
      top_annotation = ha,
      left_annotation = if (!is.null(group_vec)) {
        ComplexHeatmap::rowAnnotation(
          Group = group_vec[sample_order],
          col = list(Group = grp_pal),
          show_legend = FALSE,
          simple_anno_size = unit(4, "mm")
        )
      } else NULL,
      column_title = paste0(label, " — Pairwise Rare Allele Sharing"),
      column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
      heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8))
    )

    pdf_path <- file.path(outdir, paste0("rare_sfs_heatmap_bin", b, ".pdf"))
    sz <- max(8, nrow(mat_ordered) * 0.05 + 3)
    pdf(pdf_path, width = sz, height = sz)
    ComplexHeatmap::draw(ht)
    dev.off()
    message("[rare_sfs_plot] Wrote: ", pdf_path)

  } else {
    # ── Fallback: base R heatmap ──
    message("[rare_sfs_plot] base heatmap: bin ", b)

    pdf_path <- file.path(outdir, paste0("rare_sfs_heatmap_bin", b, ".pdf"))
    sz <- max(8, nrow(mat_ordered) * 0.04 + 3)
    pdf(pdf_path, width = sz, height = sz)

    heatmap(mat_ordered,
            Rowv = NA, Colv = NA,
            col = colorRampPalette(c("#000033", "#2166AC", "#FDDBC7", "#B2182B"))(100),
            scale = "none",
            margins = c(6, 6),
            cexRow = if (nrow(mat_ordered) <= 50) 0.5 else 0,
            cexCol = if (nrow(mat_ordered) <= 50) 0.5 else 0,
            main = paste0(label, "\nPairwise Rare Allele Sharing"),
            labRow = if (nrow(mat_ordered) > 50) "" else rownames(mat_ordered),
            labCol = if (nrow(mat_ordered) > 50) "" else colnames(mat_ordered))

    dev.off()
    message("[rare_sfs_plot] Wrote: ", pdf_path)
  }

  # ── Group-level summary heatmap ──
  if (file.exists(f_group)) {
    gdt <- fread(f_group)
    gnames <- gdt$group
    gmat <- as.matrix(gdt[, -1])
    rownames(gmat) <- colnames(gmat) <- gnames

    if (log_scale) gmat <- log10(gmat + 1)

    pdf_path <- file.path(outdir, paste0("rare_sfs_group_heatmap_bin", b, ".pdf"))
    ng <- nrow(gmat)
    pdf(pdf_path, width = max(6, ng * 0.6 + 2), height = max(5, ng * 0.5 + 2))

    heatmap(gmat,
            Rowv = NA, Colv = NA,
            col = colorRampPalette(c("#000033", "#2166AC", "#FDDBC7", "#B2182B"))(100),
            scale = "none",
            margins = c(8, 8),
            cexRow = 0.8, cexCol = 0.8,
            main = paste0(label, "\nGroup-Level Mean Pairwise Sharing"))

    dev.off()
    message("[rare_sfs_plot] Wrote: ", pdf_path, " (group summary)")
  }
}

message("[rare_sfs_plot] Done. Output: ", outdir)
