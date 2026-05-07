#!/usr/bin/env Rscript

# =============================================================================
# plot_local_Q_diagnostics.R — Diagnostic figure suite for Engine B output
#
# Generates the multi-level diagnostic figure set:
#   Level 1: Per-window tracks along the chromosome (delta12, H, ENA)
#   Level 2: Distributions by inversion segment (boxplots, Wilcoxon)
#   Level 3: Sample × window heatmaps (ordered by class)
#   Level 4: Cell-level maps (per-SNP Q probabilities, zoomed)
#
# Usage:
#   Rscript plot_local_Q_diagnostics.R \
#     --summary <chr.local_Q_summary.tsv.gz> \
#     --samples <chr.local_Q_samples.tsv.gz> \
#     --classes <snake2_band_assignments.tsv.gz> \
#     --candidate <candidate_region.tsv> \
#     --sample_list <samples.txt> \
#     [--cov <pcangsd.cov>] \
#     [--chr C_gar_LG01] \
#     [--outdir plots/] \
#     [--format png|pdf] \
#     [--K 8]
#
# Requires: data.table, ggplot2, patchwork, viridis, scales, ggpubr
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
})

# If ggpubr available, use for significance brackets
HAS_GGPUBR <- requireNamespace("ggpubr", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# =============================================================================
# PALETTE
# =============================================================================
Q_PALETTE <- c(
  Q1 = "#4E79A7", Q2 = "#F28E2B", Q3 = "#E15759", Q4 = "#76B7B2",
  Q5 = "#59A14F", Q6 = "#EDC948", Q7 = "#B07AA1", Q8 = "#FF9DA7",
  Q9 = "#9C755F", Q10 = "#BAB0AC", Q11 = "#86BCB6", Q12 = "#D37295"
)

CLASS_PALETTE <- c(
  HOM_STD     = "#4E79A7",
  HET         = "#F28E2B",
  HOM_INV     = "#E15759",
  Recombinant = "#59A14F"
)

SEGMENT_PALETTE <- c(
  "Left flank"    = "#E8E8E8",
  "Inv. left half" = "#D4E5F7",
  "Inversion core" = "#FFCCCC",
  "Inv. right half" = "#D4E5F7",
  "Right flank"    = "#E8E8E8"
)

# =============================================================================
# PARSE ARGS
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
summary_file <- NULL; samples_file <- NULL; classes_file <- NULL
candidate_file <- NULL; sample_list_file <- NULL; cov_file <- NULL
chr_name <- NULL; outdir <- "plots"; fmt <- "png"; K <- 8

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a == "--summary")     { summary_file     <- args[i+1]; i <- i+2 }
  else if (a == "--samples") { samples_file    <- args[i+1]; i <- i+2 }
  else if (a == "--classes") { classes_file    <- args[i+1]; i <- i+2 }
  else if (a == "--candidate") { candidate_file <- args[i+1]; i <- i+2 }
  else if (a == "--sample_list") { sample_list_file <- args[i+1]; i <- i+2 }
  else if (a == "--cov")     { cov_file        <- args[i+1]; i <- i+2 }
  else if (a == "--chr")     { chr_name        <- args[i+1]; i <- i+2 }
  else if (a == "--outdir")  { outdir          <- args[i+1]; i <- i+2 }
  else if (a == "--format")  { fmt             <- args[i+1]; i <- i+2 }
  else if (a == "--K")       { K               <- as.integer(args[i+1]); i <- i+2 }
  else { i <- i+1 }
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# LOAD DATA
# =============================================================================

message("[plot] Loading data...")

# Summary (per-window aggregated)
if (!is.null(summary_file) && file.exists(summary_file)) {
  dt_sum <- fread(summary_file)
  if (!is.null(chr_name)) dt_sum <- dt_sum[chrom == chr_name]
  dt_sum[, mid_bp := (start_bp + end_bp) / 2]
  dt_sum[, mid_Mb := mid_bp / 1e6]
  message("[plot] Summary: ", nrow(dt_sum), " windows")
} else {
  dt_sum <- NULL
  message("[plot] No summary file")
}

# Per-sample (window × sample)
if (!is.null(samples_file) && file.exists(samples_file)) {
  dt_samp <- fread(samples_file)
  if (!is.null(chr_name)) dt_samp <- dt_samp[chrom == chr_name]
  dt_samp[, mid_bp := (start_bp + end_bp) / 2]
  dt_samp[, mid_Mb := mid_bp / 1e6]
  message("[plot] Per-sample: ", nrow(dt_samp), " rows, ",
          uniqueN(dt_samp$sample_id %||% dt_samp$sample_idx), " samples")
} else {
  dt_samp <- NULL
  message("[plot] No per-sample file")
}

# Sample list (for ID mapping)
sample_ids <- NULL
if (!is.null(sample_list_file) && file.exists(sample_list_file)) {
  raw <- readLines(sample_list_file)
  raw <- raw[nzchar(trimws(raw))]
  # Strip BAM path if present
  sample_ids <- basename(raw)
  sample_ids <- sub("\\.sorted\\.markdup\\.bam$", "", sample_ids)
  sample_ids <- sub("\\.markdup\\.bam$", "", sample_ids)
  sample_ids <- sub("\\.sorted\\.bam$", "", sample_ids)
  sample_ids <- sub("\\.bam$", "", sample_ids)
  sample_ids <- sub("\\.cram$", "", sample_ids)
  message("[plot] Sample list: ", length(sample_ids), " IDs")

  # Map sample indices to IDs in dt_samp
  if (!is.null(dt_samp) && "sample_idx" %in% names(dt_samp)) {
    idx_col <- dt_samp$sample_idx
    if (is.numeric(idx_col) && max(idx_col, na.rm = TRUE) <= length(sample_ids)) {
      dt_samp[, sample_id := sample_ids[sample_idx]]
    }
  }
}

# PCAngsd covariance (for sample ordering)
cov_order <- NULL
if (!is.null(cov_file) && file.exists(cov_file)) {
  cov_mat <- as.matrix(fread(cov_file, header = FALSE))
  n <- nrow(cov_mat)
  # Power iteration for PC1
  v <- rep(1/sqrt(n), n)
  for (iter in 1:100) {
    w <- cov_mat %*% v
    v <- as.numeric(w / sqrt(sum(w^2)))
  }
  cov_order <- order(v)
  message("[plot] Cov order: ", n, " samples, PC1 range [",
          round(v[cov_order[1]], 3), ", ", round(v[cov_order[n]], 3), "]")
}

# Band assignments / classes
dt_class <- NULL
if (!is.null(classes_file) && file.exists(classes_file)) {
  dt_class <- fread(classes_file)
  message("[plot] Classes: ", nrow(dt_class), " assignments")
}

# Candidate region
candidate <- NULL
if (!is.null(candidate_file) && file.exists(candidate_file)) {
  cand_raw <- fread(candidate_file)
  if (!is.null(chr_name)) cand_raw <- cand_raw[chrom == chr_name]
  if (nrow(cand_raw) > 0) {
    candidate <- list(
      chrom = cand_raw$chrom[1],
      start = cand_raw$start_bp[1] %||% cand_raw$start[1],
      end   = cand_raw$end_bp[1]   %||% cand_raw$end[1]
    )
    message("[plot] Candidate: ", candidate$chrom, ":",
            candidate$start, "-", candidate$end)
  }
}

# =============================================================================
# THEME
# =============================================================================
theme_diag <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2,
                                 hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(color = "grey40", size = base_size - 1,
                                     hjust = 0.5, margin = margin(b = 6)),
      axis.title = element_text(size = base_size - 1),
      axis.text = element_text(size = base_size - 2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position = "bottom",
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

# =============================================================================
# HELPER: assign inversion segments
# =============================================================================
assign_segment <- function(pos_Mb, inv_start_Mb, inv_end_Mb) {
  inv_mid <- (inv_start_Mb + inv_end_Mb) / 2
  flank_margin <- (inv_end_Mb - inv_start_Mb) * 0.5

  fcase(
    pos_Mb < inv_start_Mb - flank_margin * 0.1, "Left flank",
    pos_Mb < inv_mid, "Inv. left half",
    pos_Mb < inv_end_Mb, "Inv. right half",
    pos_Mb < inv_end_Mb + flank_margin * 0.1, "Inv. right half",
    default = "Right flank"
  )
}

save_plot <- function(gg, name, w = 12, h = 8) {
  path <- file.path(outdir, paste0(name, ".", fmt))
  ggsave(path, gg, width = w, height = h, dpi = 300, bg = "white")
  message("[plot] Saved: ", path)
}

# =============================================================================
# LEVEL 1: PER-WINDOW TRACKS ALONG CHROMOSOME
# =============================================================================

if (!is.null(dt_sum) && nrow(dt_sum) > 0) {
  message("[plot] Level 1: Per-window tracks...")

  # Inversion core shading
  inv_rect <- if (!is.null(candidate)) {
    annotate("rect",
             xmin = candidate$start / 1e6, xmax = candidate$end / 1e6,
             ymin = -Inf, ymax = Inf,
             fill = "#FFCCCC", alpha = 0.3)
  } else {
    NULL
  }

  # Smoothed ribbon (loess)
  add_smooth_ribbon <- function(gg, yvar, color = "navy") {
    gg + geom_ribbon(
      stat = "smooth", method = "loess", span = 0.05,
      aes(ymin = after_stat(ymin), ymax = after_stat(ymax)),
      fill = color, alpha = 0.15
    ) + geom_smooth(
      method = "loess", span = 0.05,
      color = color, linewidth = 0.7, se = FALSE
    )
  }

  # Delta12 track
  p_d12 <- ggplot(dt_sum, aes(x = mid_Mb, y = mean_delta12)) +
    inv_rect +
    geom_point(size = 0.3, alpha = 0.3, color = "navy") +
    theme_diag() +
    labs(x = NULL, y = expression(Delta[12]~"(dominance)")) +
    scale_y_continuous(limits = c(0, 1)) +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = 0.9,
             label = "clear →", hjust = 0, size = 2.5, color = "grey40") +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = 0.15,
             label = "ambig. →", hjust = 0, size = 2.5, color = "grey40")
  p_d12 <- add_smooth_ribbon(p_d12, "mean_delta12", "navy")

  # Entropy track
  p_H <- ggplot(dt_sum, aes(x = mid_Mb, y = mean_entropy)) +
    inv_rect +
    geom_point(size = 0.3, alpha = 0.3, color = "#6A0DAD") +
    theme_diag() +
    labs(x = NULL, y = "H (entropy, nats)") +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = max(dt_sum$mean_entropy, na.rm=T) * 0.95,
             label = "mixed →", hjust = 0, size = 2.5, color = "grey40") +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = min(dt_sum$mean_entropy, na.rm=T) * 1.1,
             label = "focused →", hjust = 0, size = 2.5, color = "grey40")
  p_H <- add_smooth_ribbon(p_H, "mean_entropy", "#6A0DAD")

  # ENA track
  p_ena <- ggplot(dt_sum, aes(x = mid_Mb, y = mean_ena)) +
    inv_rect +
    geom_point(size = 0.3, alpha = 0.3, color = "#006400") +
    theme_diag() +
    labs(x = "Genomic position (Mb)", y = "ENA (exp(H))") +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = max(dt_sum$mean_ena, na.rm=T) * 0.95,
             label = "more groups →", hjust = 0, size = 2.5, color = "grey40") +
    annotate("text", x = max(dt_sum$mid_Mb) * 1.02, y = 1.05,
             label = "~1 group →", hjust = 0, size = 2.5, color = "grey40")
  p_ena <- add_smooth_ribbon(p_ena, "mean_ena", "#006400")

  # Optional: max_Q track (mean across samples)
  p_maxq <- NULL
  if (!is.null(dt_samp) && "max_q" %in% names(dt_samp)) {
    maxq_sum <- dt_samp[, .(mean_maxq = mean(max_q, na.rm = TRUE)),
                         by = .(window_id, mid_Mb)]
    p_maxq <- ggplot(maxq_sum, aes(x = mid_Mb, y = mean_maxq)) +
      inv_rect +
      geom_point(size = 0.3, alpha = 0.3, color = "#8B4513") +
      theme_diag() +
      labs(x = NULL, y = "Mean max_Q")
    p_maxq <- add_smooth_ribbon(p_maxq, "mean_maxq", "#8B4513")
  }

  # Optional: assigned_pop bar strip
  p_pop <- NULL
  if (!is.null(dt_samp) && "assigned_pop" %in% names(dt_samp)) {
    # Mode assignment per window
    pop_mode <- dt_samp[, .(pop = as.character(
      names(sort(table(assigned_pop), decreasing = TRUE))[1]
    )), by = .(window_id, mid_Mb)]

    p_pop <- ggplot(pop_mode, aes(x = mid_Mb, y = 0.5, fill = pop)) +
      inv_rect +
      geom_tile(width = diff(range(pop_mode$mid_Mb)) / nrow(pop_mode),
                height = 1) +
      scale_fill_manual(values = setNames(
        Q_PALETTE[1:K], as.character(1:K)
      ), name = "Assigned\npop") +
      theme_diag() +
      theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
            panel.grid = element_blank()) +
      labs(x = NULL) +
      scale_y_continuous(limits = c(0, 1))
  }

  # Assemble Level 1
  level1_plots <- list(p_d12, p_H, p_ena)
  if (!is.null(p_maxq)) level1_plots <- c(level1_plots, list(p_maxq))
  if (!is.null(p_pop))  level1_plots <- c(list(p_pop), level1_plots)

  p_level1 <- wrap_plots(level1_plots, ncol = 1) +
    plot_annotation(
      title = paste0("LEVEL 1: PER-WINDOW TRACKS — ", chr_name %||% ""),
      subtitle = "Local ancestry/structure summaries across genomic position",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", size = 10, hjust = 0.5)
      )
    )

  save_plot(p_level1, paste0("level1_tracks_", chr_name %||% "all"),
            w = 14, h = 2.5 * length(level1_plots) + 1)
}

# =============================================================================
# LEVEL 2: DISTRIBUTIONS BY INVERSION SEGMENT (boxplots + significance)
# =============================================================================

if (!is.null(dt_samp) && !is.null(candidate) && nrow(dt_samp) > 0) {
  message("[plot] Level 2: Segment distributions...")

  inv_s <- candidate$start / 1e6
  inv_e <- candidate$end / 1e6
  inv_mid <- (inv_s + inv_e) / 2

  dt_samp[, segment := fcase(
    mid_Mb < inv_s,                    "Left flank",
    mid_Mb >= inv_s & mid_Mb < inv_mid, "Inv. left\nhalf",
    mid_Mb >= inv_mid & mid_Mb <= inv_e, "Core",
    mid_Mb > inv_e & mid_Mb <= inv_e + (inv_e - inv_s) * 0.5, "Inv. right\nhalf",
    default = "Right flank"
  )]

  dt_samp[, segment := factor(segment, levels = c(
    "Left flank", "Inv. left\nhalf", "Core", "Inv. right\nhalf", "Right flank"
  ))]

  # Merge class assignments if available
  if (!is.null(dt_class) && "sample_id" %in% names(dt_class)) {
    # Find a class column
    class_col <- intersect(names(dt_class),
                            c("class", "band", "genotype", "assignment"))[1]
    if (!is.null(class_col)) {
      if (!"sample_id" %in% names(dt_samp)) {
        dt_samp[, sample_id := paste0("ind", sample_idx)]
      }
      dt_samp <- merge(dt_samp,
                        dt_class[, c("sample_id", class_col), with = FALSE],
                        by = "sample_id", all.x = TRUE)
      setnames(dt_samp, class_col, "inv_class", skip_absent = TRUE)
    }
  }

  make_segment_boxplot <- function(data, yvar, ylabel, ylim = NULL) {
    p <- ggplot(data, aes(x = segment, y = get(yvar))) +
      geom_boxplot(aes(fill = segment), outlier.size = 0.5, outlier.alpha = 0.3,
                   width = 0.7, linewidth = 0.3) +
      geom_jitter(width = 0.15, size = 0.2, alpha = 0.15) +
      scale_fill_manual(values = c(
        "Left flank" = "#D4D4D4", "Inv. left\nhalf" = "#A8C8E8",
        "Core" = "#FFAAAA", "Inv. right\nhalf" = "#A8C8E8",
        "Right flank" = "#D4D4D4"
      ), guide = "none") +
      theme_diag(base_size = 9) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7)) +
      labs(x = NULL, y = ylabel)

    if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)

    # Add Wilcoxon significance (if ggpubr available)
    if (HAS_GGPUBR) {
      comparisons <- list(
        c("Left flank", "Core"),
        c("Inv. left\nhalf", "Core"),
        c("Core", "Right flank")
      )
      p <- p + ggpubr::stat_compare_means(
        comparisons = comparisons, method = "wilcox.test",
        label = "p.signif", size = 2.5, tip.length = 0.01,
        step.increase = 0.05
      )
    }
    p
  }

  # Aggregate to window level for boxplots
  win_agg <- dt_samp[, .(
    mean_delta12 = mean(delta12, na.rm = TRUE),
    mean_entropy = mean(entropy, na.rm = TRUE),
    mean_ena = mean(ena, na.rm = TRUE)
  ), by = .(window_id, segment)]

  p_box_d12 <- make_segment_boxplot(win_agg, "mean_delta12",
                                     expression(Delta[12]~"(dominance)"), c(0, 1))
  p_box_H   <- make_segment_boxplot(win_agg, "mean_entropy",
                                     "H (entropy, nats)")
  p_box_ena <- make_segment_boxplot(win_agg, "mean_ena",
                                     "ENA (exp(H))")

  p_level2 <- (p_box_d12 | p_box_H | p_box_ena) +
    plot_annotation(
      title = "LEVEL 2: DISTRIBUTIONS BY INVERSION SEGMENT",
      subtitle = "How metrics differ between spatial regions",
      caption = "Boxes = IQR, line = median, points = windows",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", size = 9, hjust = 0.5)
      )
    )

  save_plot(p_level2, paste0("level2_segments_", chr_name %||% "all"),
            w = 12, h = 5)
}

# =============================================================================
# LEVEL 3: SAMPLE × WINDOW HEATMAPS
# =============================================================================

if (!is.null(dt_samp) && nrow(dt_samp) > 0) {
  message("[plot] Level 3: Heatmaps...")

  # Get unique windows and samples
  windows <- sort(unique(dt_samp$mid_Mb))
  sample_col <- if ("sample_id" %in% names(dt_samp)) "sample_id" else "sample_idx"
  samples <- unique(dt_samp[[sample_col]])

  # Determine sample order (cov-based PC1 or by mean delta12)
  if (!is.null(cov_order) && length(cov_order) == length(samples)) {
    sample_order <- samples[cov_order]
  } else if ("inv_class" %in% names(dt_samp)) {
    # Order by class, then by mean delta12 within class
    ord <- dt_samp[, .(md12 = mean(delta12, na.rm = TRUE),
                        cls = inv_class[1]),
                    by = sample_col]
    ord <- ord[order(cls, md12)]
    sample_order <- ord[[sample_col]]
  } else {
    ord <- dt_samp[, .(md12 = mean(delta12, na.rm = TRUE)), by = sample_col]
    ord <- ord[order(md12)]
    sample_order <- ord[[sample_col]]
  }

  dt_samp[, sample_f := factor(get(sample_col), levels = sample_order)]

  make_heatmap <- function(data, fill_var, fill_label, palette, limits = NULL) {
    p <- ggplot(data, aes(x = mid_Mb, y = sample_f, fill = get(fill_var))) +
      geom_raster() +
      theme_diag(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            legend.key.width = unit(0.8, "cm"),
            legend.key.height = unit(0.25, "cm")) +
      labs(x = "Genomic position (Mb)", y = "Samples (ordered)", fill = fill_label)

    if (is.character(palette) && length(palette) == 1) {
      p <- p + scale_fill_viridis_c(option = palette, limits = limits, na.value = "grey80")
    } else {
      p <- p + scale_fill_gradientn(colours = palette, limits = limits, na.value = "grey80")
    }

    # Add candidate shading
    if (!is.null(candidate)) {
      p <- p + annotate("rect",
                         xmin = candidate$start/1e6, xmax = candidate$end/1e6,
                         ymin = -Inf, ymax = Inf,
                         fill = NA, color = "red", linewidth = 0.5, linetype = "dashed")
    }
    p
  }

  # Delta12 heatmap
  p_heat_d12 <- make_heatmap(dt_samp, "delta12",
                              expression(Delta[12]),
                              c("white", "#08306B"), c(0, 1))

  # Entropy heatmap
  p_heat_H <- make_heatmap(dt_samp, "entropy", "H (nats)",
                            c("white", "#4A0072"),
                            c(0, max(dt_samp$entropy, na.rm = TRUE)))

  # ENA heatmap
  p_heat_ena <- make_heatmap(dt_samp, "ena", "ENA",
                              c("white", "#006400"),
                              c(1, max(dt_samp$ena, na.rm = TRUE)))

  p_level3 <- (p_heat_d12 | p_heat_H | p_heat_ena) +
    plot_annotation(
      title = "LEVEL 3: SAMPLE × WINDOW HEATMAPS",
      subtitle = "Detailed patterns across individuals and windows",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(color = "grey40", size = 9, hjust = 0.5)
      )
    )

  save_plot(p_level3, paste0("level3_heatmaps_", chr_name %||% "all"),
            w = 15, h = 7)
}

# =============================================================================
# LEVEL 4: CELL-LEVEL MAPS (dominant group argmax + metric heatmaps)
# =============================================================================

if (!is.null(dt_samp) && "assigned_pop" %in% names(dt_samp) && nrow(dt_samp) > 0) {
  message("[plot] Level 4: Cell-level maps...")

  # Zoom to candidate region or use full chromosome
  if (!is.null(candidate)) {
    margin_bp <- (candidate$end - candidate$start) * 0.3
    zoom_dt <- dt_samp[mid_bp >= (candidate$start - margin_bp) &
                         mid_bp <= (candidate$end + margin_bp)]
  } else {
    zoom_dt <- dt_samp
  }

  if (nrow(zoom_dt) > 0) {
    sample_col <- if ("sample_id" %in% names(zoom_dt)) "sample_id" else "sample_idx"

    # Reorder samples
    if (exists("sample_order")) {
      zoom_dt[, sample_f := factor(get(sample_col), levels = sample_order)]
    } else {
      zoom_dt[, sample_f := factor(get(sample_col))]
    }

    # Dominant group (argmax)
    zoom_dt[, pop_f := factor(assigned_pop, levels = 1:K)]

    p_cell_pop <- ggplot(zoom_dt, aes(x = mid_Mb, y = sample_f, fill = pop_f)) +
      geom_raster() +
      scale_fill_manual(
        values = setNames(Q_PALETTE[1:K], as.character(1:K)),
        name = "Group", drop = FALSE
      ) +
      theme_diag(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank()) +
      labs(x = "Genomic position (Mb)", y = "Samples",
           title = "Dominant group\n(argmax)")

    # Delta12, H, ENA heatmaps for zoom
    p_cell_d12 <- ggplot(zoom_dt, aes(x = mid_Mb, y = sample_f, fill = delta12)) +
      geom_raster() +
      scale_fill_gradientn(colours = c("white", "#08306B"),
                           limits = c(0, 1), name = expression(Delta[12])) +
      theme_diag(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank()) +
      labs(x = "Genomic position (Mb)", y = NULL, title = expression(Delta[12]~"(dominance)"))

    p_cell_H <- ggplot(zoom_dt, aes(x = mid_Mb, y = sample_f, fill = entropy)) +
      geom_raster() +
      scale_fill_gradientn(colours = c("white", "#4A0072"),
                           name = "H") +
      theme_diag(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank()) +
      labs(x = "Genomic position (Mb)", y = NULL, title = "H (entropy, nats)")

    p_cell_ena <- ggplot(zoom_dt, aes(x = mid_Mb, y = sample_f, fill = ena)) +
      geom_raster() +
      scale_fill_gradientn(colours = c("white", "#006400"),
                           name = "ENA") +
      theme_diag(base_size = 8) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            panel.grid = element_blank()) +
      labs(x = "Genomic position (Mb)", y = NULL, title = "ENA (exp(H))")

    p_level4 <- (p_cell_pop | p_cell_d12 | p_cell_H | p_cell_ena) +
      plot_annotation(
        title = "CELL-LEVEL MAPS",
        subtitle = if (!is.null(candidate))
          paste0("Zoomed region: ", round(candidate$start/1e6, 2), " – ",
                  round(candidate$end/1e6, 2), " Mb")
        else "Full chromosome",
        theme = theme(
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(color = "grey40", size = 9, hjust = 0.5)
        )
      )

    save_plot(p_level4, paste0("level4_cellmaps_", chr_name %||% "all"),
              w = 16, h = 7)
  }
}

# =============================================================================
# COMBINED FIGURE (all 4 levels, single composite)
# =============================================================================

if (!is.null(dt_sum) && !is.null(dt_samp) && nrow(dt_samp) > 0) {
  message("[plot] Building combined figure...")

  # Simplified versions for composite
  # Just use the individual plots already created and combine
  # This creates a large composite matching the schema in Image 2

  # We already saved individual panels; also save a note about assembly
  message("[plot] Individual panels saved. For composite figure matching the")
  message("       uploaded schema, use Inkscape to assemble:")
  message("         level1_tracks + level2_segments + level3_heatmaps + level4_cellmaps")
  message("       Or run plot_composite_figure.R (to be created).")
}

# =============================================================================
# HOW TO READ annotation block
# =============================================================================

if (!is.null(dt_samp) && nrow(dt_samp) > 0) {
  message("[plot] Generating interpretation guide...")

  guide_text <- paste(
    "HOW TO READ:",
    "",
    "• High Δ₁₂, low H, ENA ≈ 1",
    "  = one group dominates (clear assignment)",
    "",
    "• Low Δ₁₂, high H, high ENA",
    "  = mixed / ambiguous ancestry",
    "",
    "• Compare segments to see where",
    "  structure is clean vs. mixed",
    "",
    "• Heatmaps reveal patchwork",
    "  and recombinants",
    "",
    "• Cell-level maps show per-sample",
    "  ancestry switching along the chromosome",
    sep = "\n"
  )

  writeLines(guide_text, file.path(outdir, "interpretation_guide.txt"))
}

message("[plot] All done. Outputs in: ", outdir)
