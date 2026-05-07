#!/usr/bin/env Rscript

# =============================================================================
# 05_plot_hobs_hwe.R — Part 12: Hobs/HWE Plotting Suite
#
# Genome-wide tracks per subset, heatmaps, outlier burden,
# subset comparisons, candidate overlay stacked figures.
#
# Reads the output of hobs_windower (Step 3) and candidate_overlay (Step 4).
#
# Modes (--mode):
#   genome_tracks     Genome-wide Hobs/F line tracks per subset, one panel per chr
#   heatmap           Subset × window heatmap (mean Hobs or mean F)
#   outlier_burden    Stacked bar of outlier fractions per chr per subset
#   subset_compare    Paired comparison of two subsets across all chr
#   candidate_stack   Per-candidate stacked figure (Hobs + F + outlier + Fst)
#   all               Run all modes
#
# Usage:
#   Rscript 05_plot_hobs_hwe.R --config 00_hobs_hwe_config.sh --mode all
#   Rscript 05_plot_hobs_hwe.R --config ... --mode candidate_stack --scale 50kb
#   Rscript 05_plot_hobs_hwe.R --config ... --mode genome_tracks --subset Q1_full
#
# Part 12 of the Hobs/HWE Secondary Confirmation Module spec.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))

# ── CLI ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
config_file <- NULL
mode <- "all"
scale_label <- "50kb"
subset_filter <- NULL
candidate_file <- NULL
fst_dir <- NULL
outdir <- NULL

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--config"    = { config_file <- args[i+1]; i <- i+2 },
    "--mode"      = { mode <- args[i+1]; i <- i+2 },
    "--scale"     = { scale_label <- args[i+1]; i <- i+2 },
    "--subset"    = { subset_filter <- args[i+1]; i <- i+2 },
    "--candidates"= { candidate_file <- args[i+1]; i <- i+2 },
    "--fst_dir"   = { fst_dir <- args[i+1]; i <- i+2 },
    "--outdir"    = { outdir <- args[i+1]; i <- i+2 },
    { i <- i+1 }
  )
}

# ── Load config ──────────────────────────────────────────────────────────────
# Parse shell config to R variables
parse_shell_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#|^\\s*$", lines)]
  lines <- lines[grepl("^export\\s+|^[A-Z_]+=", lines)]
  env <- list()
  for (l in lines) {
    l <- sub("^export\\s+", "", l)
    eq <- regexpr("=", l, fixed = TRUE)
    if (eq < 2) next
    key <- substr(l, 1, eq - 1)
    val <- substr(l, eq + 1, nchar(l))
    val <- gsub("^[\"']|[\"']$", "", val)
    val <- gsub("\\$\\{[^}]+:-([^}]+)\\}", "\\1", val)
    env[[key]] <- val
  }
  env
}

if (!is.null(config_file) && file.exists(config_file)) {
  cfg <- parse_shell_config(config_file)
} else {
  # Try default location
  default <- file.path(dirname(sys.frame(1)$ofile %||% "."), "00_hobs_hwe_config.sh")
  if (file.exists(default)) cfg <- parse_shell_config(default)
  else stop("Config not found. Use --config <path>")
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a

HOBS_OUTDIR   <- cfg$HOBS_OUTDIR %||% "hobs_hwe_confirmation"
WINDOWS_DIR   <- file.path(HOBS_OUTDIR, "03_hobs_windows")
OVERLAY_DIR   <- file.path(HOBS_OUTDIR, "04_candidate_overlay")
MANIFEST_PATH <- file.path(HOBS_OUTDIR, "subsets", "subset_manifest.tsv")
REF_FAI       <- cfg$REF_FAI %||% cfg$REF %||% ""
outdir        <- outdir %||% file.path(HOBS_OUTDIR, "05_plots")
fst_dir       <- fst_dir %||% cfg$STATS_CACHE_DIR %||% ""

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Load manifest ────────────────────────────────────────────────────────────
manifest <- NULL
if (file.exists(MANIFEST_PATH)) {
  manifest <- fread(MANIFEST_PATH)
  if (!is.null(subset_filter)) manifest <- manifest[subset_id == subset_filter]
}

# ── Load chromosome info ─────────────────────────────────────────────────────
chroms <- NULL
if (nzchar(REF_FAI) && file.exists(REF_FAI)) {
  fai <- fread(REF_FAI, header = FALSE, select = 1:2)
  setnames(fai, c("chrom", "length"))
  chroms <- fai$chrom
}

# ── Ancestry palette (from config or default) ───────────────────────────────
PALETTE <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#86BCB6", "#D37295"
)

# ── Helper: load windows for one subset × one chr × one scale ────────────────
load_windows <- function(subset_id, chr, scale) {
  path <- file.path(WINDOWS_DIR, subset_id, paste0(chr, ".win", scale, ".tsv"))
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  dt[, subset_id := subset_id]
  dt
}

# ── Helper: load all windows for one subset × one scale ──────────────────────
load_all_windows <- function(subset_id, scale, chroms_vec) {
  dts <- lapply(chroms_vec, function(chr) load_windows(subset_id, chr, scale))
  rbindlist(dts[!vapply(dts, is.null, logical(1))], fill = TRUE)
}

# ── Helper: load all subsets × all chr × one scale ───────────────────────────
load_all <- function(subset_ids, scale, chroms_vec) {
  dts <- lapply(subset_ids, function(sid) load_all_windows(sid, scale, chroms_vec))
  rbindlist(dts[!vapply(dts, is.null, logical(1))], fill = TRUE)
}

fmt_mb <- function(bp) sprintf("%.1f", bp / 1e6)

# ── Candidate shading rect layer ─────────────────────────────────────────────
candidate_rect <- function(cand_row) {
  annotate("rect",
    xmin = as.numeric(cand_row$start_bp) / 1e6,
    xmax = as.numeric(cand_row$end_bp) / 1e6,
    ymin = -Inf, ymax = Inf, alpha = 0.15, fill = "grey50"
  )
}

# =============================================================================
# MODE 1: GENOME-WIDE TRACKS
# =============================================================================

plot_genome_tracks <- function(subset_ids, scale, chroms_vec) {
  message("[plot] Genome-wide tracks: scale=", scale)

  for (sid in subset_ids) {
    dt <- load_all_windows(sid, scale, chroms_vec)
    if (is.null(dt) || nrow(dt) == 0) { message("  skip ", sid); next }

    dt[, chrom := factor(chrom, levels = chroms_vec)]

    # ── Hobs track ──
    p_hobs <- ggplot(dt, aes(x = window_center / 1e6, y = mean_Hobs)) +
      geom_line(linewidth = 0.25, color = "#4E79A7", alpha = 0.8) +
      facet_wrap(~chrom, scales = "free_x", ncol = 4) +
      theme_bw(base_size = 9) +
      theme(strip.text = element_text(size = 7),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 6)) +
      labs(title = paste0("Hobs genome track — ", sid, " (", scale, ")"),
           x = "Position (Mb)", y = "Mean Hobs")

    ggsave(file.path(outdir, paste0("genome_Hobs_", sid, "_", scale, ".pdf")),
           p_hobs, width = 14, height = ceiling(length(chroms_vec) / 4) * 2.5 + 1)

    # ── F track ──
    p_f <- ggplot(dt, aes(x = window_center / 1e6, y = mean_F)) +
      geom_line(linewidth = 0.25, color = "#E15759", alpha = 0.8) +
      facet_wrap(~chrom, scales = "free_x", ncol = 4) +
      theme_bw(base_size = 9) +
      theme(strip.text = element_text(size = 7),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 6)) +
      labs(title = paste0("F genome track — ", sid, " (", scale, ")"),
           x = "Position (Mb)", y = "Mean F")

    ggsave(file.path(outdir, paste0("genome_F_", sid, "_", scale, ".pdf")),
           p_f, width = 14, height = ceiling(length(chroms_vec) / 4) * 2.5 + 1)

    message("  ", sid, " done (", nrow(dt), " windows)")
  }
}

# =============================================================================
# MODE 2: HEATMAP (subset × window)
# =============================================================================

plot_heatmap <- function(subset_ids, scale, chroms_vec) {
  message("[plot] Heatmap: scale=", scale)

  for (chr in chroms_vec) {
    dts <- lapply(subset_ids, function(sid) load_windows(sid, chr, scale))
    dts <- dts[!vapply(dts, is.null, logical(1))]
    if (length(dts) == 0) next
    dt <- rbindlist(dts, fill = TRUE)

    # Bin positions to match across subsets
    dt[, pos_bin := round(window_center / 1e5) * 1e5]

    # Cast to wide: rows = subsets, cols = position bins
    hobs_wide <- dcast(dt, subset_id ~ pos_bin, value.var = "mean_Hobs",
                        fun.aggregate = mean)
    f_wide <- dcast(dt, subset_id ~ pos_bin, value.var = "mean_F",
                     fun.aggregate = mean)

    # ggplot heatmap via geom_tile
    p_heat <- ggplot(dt, aes(x = window_center / 1e6, y = subset_id, fill = mean_Hobs)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma", name = "Hobs") +
      theme_minimal(base_size = 9) +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 6)) +
      labs(title = paste0("Hobs heatmap — ", chr, " (", scale, ")"),
           x = "Position (Mb)", y = "Subset")

    ht <- max(3, length(subset_ids) * 0.5 + 2)
    ggsave(file.path(outdir, paste0("heatmap_Hobs_", chr, "_", scale, ".pdf")),
           p_heat, width = 12, height = ht)

    # F heatmap
    p_heat_f <- ggplot(dt, aes(x = window_center / 1e6, y = subset_id, fill = mean_F)) +
      geom_tile() +
      scale_fill_viridis_c(option = "inferno", name = "F") +
      theme_minimal(base_size = 9) +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 6)) +
      labs(title = paste0("F heatmap — ", chr, " (", scale, ")"),
           x = "Position (Mb)", y = "Subset")

    ggsave(file.path(outdir, paste0("heatmap_F_", chr, "_", scale, ".pdf")),
           p_heat_f, width = 12, height = ht)
  }
}

# =============================================================================
# MODE 3: OUTLIER BURDEN
# =============================================================================

plot_outlier_burden <- function(subset_ids, scale, chroms_vec) {
  message("[plot] Outlier burden: scale=", scale)

  dt_all <- load_all(subset_ids, scale, chroms_vec)
  if (nrow(dt_all) == 0) return(invisible())

  # Per subset × chromosome: mean fraction of outlier windows
  burden <- dt_all[, .(
    frac_lo_Hobs = mean(frac_low_Hobs_outlier, na.rm = TRUE),
    frac_hi_Hobs = mean(frac_high_Hobs_outlier, na.rm = TRUE),
    frac_lo_F    = mean(frac_low_F_outlier, na.rm = TRUE),
    frac_hi_F    = mean(frac_high_F_outlier, na.rm = TRUE),
    n_windows    = .N
  ), by = .(subset_id, chrom)]

  burden[, chrom := factor(chrom, levels = chroms_vec)]

  # Melt for stacked bar
  burden_m <- melt(burden,
    id.vars = c("subset_id", "chrom", "n_windows"),
    measure.vars = c("frac_lo_Hobs", "frac_hi_Hobs", "frac_lo_F", "frac_hi_F"),
    variable.name = "outlier_type", value.name = "fraction"
  )

  burden_m[, outlier_type := factor(outlier_type, levels = c(
    "frac_lo_Hobs", "frac_hi_Hobs", "frac_lo_F", "frac_hi_F"
  ), labels = c("Low Hobs", "High Hobs", "Low F", "High F"))]

  p <- ggplot(burden_m, aes(x = chrom, y = fraction, fill = outlier_type)) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~subset_id, ncol = 2) +
    scale_fill_manual(values = c(
      "Low Hobs" = "#2166AC", "High Hobs" = "#B2182B",
      "Low F" = "#4393C3", "High F" = "#D6604D"
    )) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          strip.text = element_text(size = 8)) +
    labs(title = paste0("Outlier burden per chromosome (", scale, ")"),
         x = "Chromosome", y = "Mean fraction outlier sites", fill = "Type")

  ht <- max(6, ceiling(length(subset_ids) / 2) * 3 + 2)
  ggsave(file.path(outdir, paste0("outlier_burden_", scale, ".pdf")),
         p, width = 14, height = ht)
  message("  done: ", nrow(burden), " subset × chr combos")
}

# =============================================================================
# MODE 4: SUBSET COMPARISON
# =============================================================================

plot_subset_compare <- function(subset_ids, scale, chroms_vec) {
  message("[plot] Subset comparison: scale=", scale)

  if (length(subset_ids) < 2) {
    message("  Need ≥2 subsets for comparison — skip"); return(invisible())
  }

  dt_all <- load_all(subset_ids, scale, chroms_vec)
  if (nrow(dt_all) == 0) return(invisible())

  dt_all[, chrom := factor(chrom, levels = chroms_vec)]

  # Per-subset genome-wide summary
  summary_dt <- dt_all[, .(
    median_Hobs = median(mean_Hobs, na.rm = TRUE),
    mad_Hobs    = mad(mean_Hobs, na.rm = TRUE),
    median_F    = median(mean_F, na.rm = TRUE),
    mad_F       = mad(mean_F, na.rm = TRUE),
    n_windows   = .N
  ), by = subset_id]

  fwrite(summary_dt, file.path(outdir, paste0("subset_summary_", scale, ".tsv")),
         sep = "\t")

  # Per-chromosome paired overlay: pick first two subsets as primary contrast
  s1 <- subset_ids[1]; s2 <- subset_ids[2]
  d1 <- dt_all[subset_id == s1]
  d2 <- dt_all[subset_id == s2]

  # Merge by approximate position
  d1[, pos_bin := round(window_center / 5e4) * 5e4]
  d2[, pos_bin := round(window_center / 5e4) * 5e4]
  merged <- merge(d1[, .(chrom, pos_bin, Hobs_1 = mean_Hobs, F_1 = mean_F)],
                  d2[, .(chrom, pos_bin, Hobs_2 = mean_Hobs, F_2 = mean_F)],
                  by = c("chrom", "pos_bin"))
  merged[, delta_Hobs := Hobs_1 - Hobs_2]
  merged[, delta_F := F_1 - F_2]

  p_delta <- ggplot(merged, aes(x = pos_bin / 1e6, y = delta_Hobs)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 0.3, color = "#4E79A7", alpha = 0.7) +
    facet_wrap(~chrom, scales = "free_x", ncol = 4) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7),
          axis.text.x = element_text(size = 6)) +
    labs(title = paste0("ΔHobs: ", s1, " − ", s2, " (", scale, ")"),
         subtitle = "Positive = higher Hobs in first subset",
         x = "Position (Mb)", y = "ΔHobs")

  ggsave(file.path(outdir, paste0("subset_compare_deltaHobs_", s1, "_vs_", s2,
                                    "_", scale, ".pdf")),
         p_delta, width = 14, height = ceiling(length(chroms_vec) / 4) * 2.5 + 1)

  message("  done")
}

# =============================================================================
# MODE 5: CANDIDATE OVERLAY STACKED FIGURE
# =============================================================================

plot_candidate_stack <- function(subset_ids, scale, chroms_vec) {
  message("[plot] Candidate stacked figures: scale=", scale)

  if (!has_patchwork) {
    message("  patchwork not available — skip stacked figures")
    return(invisible())
  }

  # Load candidate overlay if available
  overlay_path <- file.path(OVERLAY_DIR, paste0("candidate_hobs_overlay_", scale, ".tsv"))
  if (!file.exists(overlay_path)) {
    message("  No overlay file: ", overlay_path, " — skip")
    return(invisible())
  }
  overlay <- fread(overlay_path)

  # Load candidate table
  cand_file <- candidate_file
  if (is.null(cand_file) || !file.exists(cand_file %||% "")) {
    # Try to infer from overlay
    cands <- unique(overlay[, .(chrom, start, end, candidate_id)])
  } else {
    cands <- fread(cand_file)
    if ("start_bp" %in% names(cands)) setnames(cands, "start_bp", "start")
    if ("end_bp" %in% names(cands)) setnames(cands, "end_bp", "end")
  }

  cand_dir <- file.path(outdir, "candidate_stacks")
  dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)

  for (ci in seq_len(nrow(cands))) {
    cr <- cands[ci]
    cid <- cr$candidate_id
    c_chr <- cr$chrom
    c_start <- as.numeric(cr$start)
    c_end <- as.numeric(cr$end)
    flank_bp <- max(200000, (c_end - c_start) * 0.3)
    reg_start <- max(0, c_start - flank_bp)
    reg_end <- c_end + flank_bp

    # Collect windows for all subsets in this region
    win_list <- list()
    for (sid in subset_ids) {
      w <- load_windows(sid, c_chr, scale)
      if (is.null(w)) next
      w <- w[window_center >= reg_start & window_center <= reg_end]
      if (nrow(w) == 0) next
      win_list[[length(win_list) + 1]] <- w
    }
    if (length(win_list) == 0) next
    win_dt <- rbindlist(win_list, fill = TRUE)

    n_subsets <- uniqueN(win_dt$subset_id)
    pal <- setNames(PALETTE[seq_len(n_subsets)], unique(win_dt$subset_id))

    # ── Panel A: Hobs tracks ──
    pA <- ggplot(win_dt, aes(x = window_center / 1e6, y = mean_Hobs,
                              color = subset_id)) +
      annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
               ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      scale_color_manual(values = pal) +
      theme_bw(base_size = 9) +
      theme(legend.position = "none") +
      labs(title = paste0("Candidate ", cid, " — ", c_chr, ":",
                          fmt_mb(c_start), "–", fmt_mb(c_end), " Mb"),
           subtitle = paste0("Scale: ", scale, " | Subsets: ", n_subsets),
           x = NULL, y = "Mean Hobs")

    # ── Panel B: F tracks ──
    pB <- ggplot(win_dt, aes(x = window_center / 1e6, y = mean_F,
                              color = subset_id)) +
      annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
               ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      scale_color_manual(values = pal) +
      theme_bw(base_size = 9) +
      theme(legend.position = "none") +
      labs(x = NULL, y = "Mean F")

    # ── Panel C: Outlier fraction ──
    pC <- ggplot(win_dt, aes(x = window_center / 1e6, y = frac_high_F_outlier,
                              color = subset_id)) +
      annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
               ymin = -Inf, ymax = Inf, alpha = 0.12, fill = "grey50") +
      geom_line(linewidth = 0.5, alpha = 0.8) +
      geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      scale_color_manual(values = pal) +
      theme_bw(base_size = 9) +
      theme(legend.position = "none") +
      labs(x = NULL, y = "Frac high-F outlier")

    # ── Panel D: Overlay summary bar ──
    ov_sub <- overlay[candidate_id == cid]
    if (nrow(ov_sub) > 0) {
      pD <- ggplot(ov_sub, aes(x = subset_id, y = as.numeric(mean_F_interval),
                                fill = pattern_class)) +
        geom_col(width = 0.6) +
        coord_flip() +
        scale_fill_manual(values = c(
          "broad_localized_distortion" = "#D6604D",
          "few_extreme_loci"           = "#F4A582",
          "regional_het_deficit"       = "#92C5DE",
          "weak_signal"                = "#E0E0E0",
          "no_signal"                  = "#F5F5F5",
          "insufficient_data"          = "#BDBDBD"
        ), drop = FALSE) +
        theme_bw(base_size = 9) +
        labs(x = NULL, y = "Mean F in interval", fill = "Pattern")
    } else {
      pD <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5,
        label = "No overlay data", size = 3, color = "grey60")
    }

    # ── Assemble ──
    combined <- pA / pB / pC / pD + plot_layout(heights = c(3, 3, 2, 2))

    ggsave(file.path(cand_dir,
           paste0("candidate_", cid, "_", c_chr, "_hobs_stack_", scale, ".pdf")),
           combined, width = 10, height = 12)

    message("  candidate ", cid, " done")
  }
}

# =============================================================================
# DISPATCH
# =============================================================================

subset_ids <- if (!is.null(manifest)) manifest$subset_id else character(0)
if (length(subset_ids) == 0) {
  # Try to auto-discover from directory listing
  if (dir.exists(WINDOWS_DIR)) {
    subset_ids <- list.dirs(WINDOWS_DIR, recursive = FALSE, full.names = FALSE)
  }
}

if (is.null(chroms) || length(chroms) == 0) {
  # Auto-discover from first subset's files
  if (length(subset_ids) > 0) {
    tsv_files <- list.files(file.path(WINDOWS_DIR, subset_ids[1]),
                            pattern = paste0("\\.win", scale_label, "\\.tsv$"))
    chroms <- sub(paste0("\\.win", scale_label, "\\.tsv$"), "", tsv_files)
    chroms <- sort(chroms)
  }
}

message("[plot] Subsets: ", paste(subset_ids, collapse = ", "))
message("[plot] Chromosomes: ", length(chroms))
message("[plot] Scale: ", scale_label)
message("[plot] Mode: ", mode)

run_mode <- function(m) {
  switch(m,
    "genome_tracks"  = plot_genome_tracks(subset_ids, scale_label, chroms),
    "heatmap"        = plot_heatmap(subset_ids, scale_label, chroms),
    "outlier_burden" = plot_outlier_burden(subset_ids, scale_label, chroms),
    "subset_compare" = plot_subset_compare(subset_ids, scale_label, chroms),
    "candidate_stack"= plot_candidate_stack(subset_ids, scale_label, chroms),
    message("[plot] Unknown mode: ", m)
  )
}

if (mode == "all") {
  for (m in c("genome_tracks", "heatmap", "outlier_burden",
              "subset_compare", "candidate_stack")) {
    tryCatch(run_mode(m), error = function(e) message("[WARN] ", m, ": ", e$message))
  }
} else {
  run_mode(mode)
}

message("[plot] Complete. Output: ", outdir)
