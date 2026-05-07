#!/usr/bin/env Rscript

# =============================================================================
# 06_plot_fst_dxy_tracks.R — Windowed Fst/dXY/MI/theta track plots
#
# Reads region_popstats.c output (per-chr .popstats.tsv files).
# Produces Mérot-style multi-track figures:
#   - Pairwise Fst curves per group contrast (colored lines)
#   - dXY / dA / MI tracks
#   - Theta_pi per group
#   - Pooled Hp / Tajima_D
#   - Candidate interval shading + optional sub-bloc ribbons
#
# Input file format (from region_popstats):
#   --popstats_dir <dir>      Directory of <chr>.popstats.tsv files
#   --candidates <file>       Candidate table (candidate_id, chrom, start_bp, end_bp)
#   --blocs <file>            Optional: sub-bloc table (candidate_id, bloc_id,
#                              start_bp, end_bp, bloc_type [flank/outer/inner/core/recomb])
#   --groups_order <file>     Optional: group display order + colors
#
# Modes:
#   genome_fst        All chromosomes, Fst panel per pairwise comparison
#   candidate_stack   Per-candidate stacked multi-metric figure
#   single_chr        One chromosome, all metrics
#
# Usage:
#   Rscript 06_plot_fst_dxy_tracks.R \
#     --popstats_dir region_stats_cache/ \
#     --candidates candidates.tsv \
#     --mode candidate_stack \
#     --outdir plots/fst_dxy/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))

# ── CLI ──
args <- commandArgs(trailingOnly = TRUE)
popstats_dir <- "region_stats_cache"
candidates_file <- NULL
blocs_file <- NULL
groups_order_file <- NULL
mode <- "candidate_stack"
outdir <- "plots/fst_dxy"
scale <- NULL  # auto from data

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "--popstats_dir" = { popstats_dir <- args[i+1]; i <- i+2 },
    "--candidates"   = { candidates_file <- args[i+1]; i <- i+2 },
    "--blocs"        = { blocs_file <- args[i+1]; i <- i+2 },
    "--groups_order" = { groups_order_file <- args[i+1]; i <- i+2 },
    "--mode"         = { mode <- args[i+1]; i <- i+2 },
    "--outdir"       = { outdir <- args[i+1]; i <- i+2 },
    { i <- i+1 }
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── Palette ──
PALETTE <- c(
  "#2166AC", "#D6604D", "#4DAF4A", "#FF7F00", "#984EA3",
  "#E7298A", "#66C2A5", "#E6AB02", "#A6761D", "#666666",
  "#1B9E77", "#7570B3"
)

fmt_mb <- function(bp) sprintf("%.2f", bp / 1e6)

# ── Load popstats for a chromosome ──
load_popstats <- function(chr) {
  f <- file.path(popstats_dir, paste0(chr, ".popstats.tsv"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, skip = 1)  # skip metadata comment
  dt[, mid_bp := (start + end) / 2]
  dt
}

# ── Detect pairwise columns ──
detect_pairwise <- function(dt, prefix = "Fst") {
  grep(paste0("^", prefix, "_"), names(dt), value = TRUE)
}

# ── Sub-bloc ribbon layer ──
bloc_ribbon <- function(blocs_dt, cand_id) {
  if (is.null(blocs_dt)) return(NULL)
  b <- blocs_dt[candidate_id == cand_id]
  if (nrow(b) == 0) return(NULL)

  bloc_colors <- c(
    flank = "#F0F0F0", outer = "#D9E8F7", inner = "#B3CDE3",
    core = "#6BAED6", recomb = "#FDD0A2", breakpoint = "#FDAE6B"
  )

  list(
    geom_rect(data = b,
      aes(xmin = start_bp / 1e6, xmax = end_bp / 1e6,
          ymin = -Inf, ymax = Inf, fill = bloc_type),
      alpha = 0.2, inherit.aes = FALSE),
    scale_fill_manual(values = bloc_colors, name = "Bloc", guide = "none")
  )
}

# ── Candidate interval rect ──
cand_rect <- function(c_start, c_end) {
  annotate("rect", xmin = c_start / 1e6, xmax = c_end / 1e6,
           ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey30")
}

# =============================================================================
# MODE: candidate_stack — per-candidate multi-metric stacked figure
# =============================================================================

plot_candidate_stack <- function() {
  if (is.null(candidates_file) || !file.exists(candidates_file)) {
    message("[fst_plots] No candidates file"); return(invisible())
  }
  if (!has_patchwork) {
    message("[fst_plots] patchwork required for stacked figures"); return(invisible())
  }

  cands <- fread(candidates_file)
  blocs <- if (!is.null(blocs_file) && file.exists(blocs_file)) fread(blocs_file) else NULL

  for (ci in seq_len(nrow(cands))) {
    cr <- cands[ci]
    cid <- cr$candidate_id
    chr <- cr$chrom
    c_start <- as.numeric(cr$start_bp)
    c_end <- as.numeric(cr$end_bp)
    inv_len <- c_end - c_start
    flank <- max(200000, inv_len * 0.3)
    reg_start <- max(0, c_start - flank)
    reg_end <- c_end + flank

    dt <- load_popstats(chr)
    if (is.null(dt) || nrow(dt) == 0) next
    dt <- dt[mid_bp >= reg_start & mid_bp <= reg_end]
    if (nrow(dt) < 3) next

    # ── Panel A: Fst tracks ──
    fst_cols <- detect_pairwise(dt, "Fst")
    if (length(fst_cols) > 0) {
      fst_m <- melt(dt, id.vars = c("mid_bp"), measure.vars = fst_cols,
                     variable.name = "contrast", value.name = "Fst")
      fst_m[, contrast := sub("^Fst_", "", contrast)]

      pal_fst <- setNames(PALETTE[seq_along(fst_cols)], sub("^Fst_", "", fst_cols))

      pA <- ggplot(fst_m, aes(x = mid_bp / 1e6, y = Fst, color = contrast)) +
        cand_rect(c_start, c_end) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        scale_color_manual(values = pal_fst) +
        theme_bw(base_size = 9) +
        theme(legend.position = "top", legend.key.size = unit(0.3, "cm"),
              legend.text = element_text(size = 7)) +
        labs(title = paste0("Candidate ", cid, " — ", chr, ":",
                            fmt_mb(c_start), "–", fmt_mb(c_end), " Mb"),
             x = NULL, y = expression(F[ST]), color = "Contrast")
    } else {
      pA <- ggplot() + theme_void() + annotate("text", x = .5, y = .5,
        label = "No Fst columns", color = "grey60", size = 3)
    }

    # ── Panel B: dXY tracks ──
    dxy_cols <- detect_pairwise(dt, "dXY")
    if (length(dxy_cols) > 0) {
      dxy_m <- melt(dt, id.vars = "mid_bp", measure.vars = dxy_cols,
                     variable.name = "contrast", value.name = "dXY")
      dxy_m[, contrast := sub("^dXY_", "", contrast)]
      pal_dxy <- setNames(PALETTE[seq_along(dxy_cols)], sub("^dXY_", "", dxy_cols))

      pB <- ggplot(dxy_m, aes(x = mid_bp / 1e6, y = dXY, color = contrast)) +
        cand_rect(c_start, c_end) +
        geom_line(linewidth = 0.6, alpha = 0.85) +
        scale_color_manual(values = pal_dxy) +
        theme_bw(base_size = 9) +
        theme(legend.position = "none") +
        labs(x = NULL, y = expression(d[XY]), color = NULL)
    } else {
      pB <- ggplot() + theme_void()
    }

    # ── Panel C: theta_pi per group ──
    tpi_cols <- grep("^theta_pi_", names(dt), value = TRUE)
    tpi_cols <- setdiff(tpi_cols, "theta_pi_all")  # keep only per-group
    if (length(tpi_cols) > 0) {
      tpi_m <- melt(dt, id.vars = "mid_bp", measure.vars = tpi_cols,
                     variable.name = "group", value.name = "theta_pi")
      tpi_m[, group := sub("^theta_pi_", "", group)]
      pal_tpi <- setNames(PALETTE[seq_along(tpi_cols)], sub("^theta_pi_", "", tpi_cols))

      pC <- ggplot(tpi_m, aes(x = mid_bp / 1e6, y = theta_pi, color = group)) +
        cand_rect(c_start, c_end) +
        geom_line(linewidth = 0.5, alpha = 0.8) +
        scale_color_manual(values = pal_tpi) +
        theme_bw(base_size = 9) +
        theme(legend.position = "none") +
        labs(x = NULL, y = expression(theta[pi]), color = NULL)
    } else {
      pC <- ggplot() + theme_void()
    }

    # ── Panel D: Tajima_D + Hp ──
    pD <- ggplot(dt, aes(x = mid_bp / 1e6)) +
      cand_rect(c_start, c_end) +
      geom_line(aes(y = Tajima_D), color = "#E15759", linewidth = 0.5, alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      theme_bw(base_size = 9) +
      labs(x = NULL, y = "Tajima's D")

    # ── Panel E: MI ──
    mi_cols <- detect_pairwise(dt, "MInorm")
    if (length(mi_cols) > 0) {
      mi_m <- melt(dt, id.vars = "mid_bp", measure.vars = mi_cols,
                    variable.name = "contrast", value.name = "MInorm")
      mi_m[, contrast := sub("^MInorm_", "", contrast)]
      pal_mi <- setNames(PALETTE[seq_along(mi_cols)], sub("^MInorm_", "", mi_cols))

      pE <- ggplot(mi_m, aes(x = mid_bp / 1e6, y = MInorm, color = contrast)) +
        cand_rect(c_start, c_end) +
        geom_line(linewidth = 0.5, alpha = 0.8) +
        scale_color_manual(values = pal_mi) +
        theme_bw(base_size = 9) +
        theme(legend.position = "none") +
        labs(x = "Position (Mb)", y = "MI (normalized)", color = NULL)
    } else {
      pE <- ggplot(dt, aes(x = mid_bp / 1e6, y = Hp)) +
        cand_rect(c_start, c_end) +
        geom_line(color = "#4E79A7", linewidth = 0.5) +
        theme_bw(base_size = 9) +
        labs(x = "Position (Mb)", y = "Hp")
    }

    # ── Assemble ──
    combined <- pA / pB / pC / pD / pE +
      plot_layout(heights = c(3, 2, 2, 2, 2))

    fname <- paste0("fst_stack_candidate_", cid, "_", chr, ".pdf")
    ggsave(file.path(outdir, fname), combined, width = 10, height = 14)
    message("[fst_plots] ", fname)
  }
}

# =============================================================================
# MODE: genome_fst — all chromosomes, one Fst panel per pairwise contrast
# =============================================================================

plot_genome_fst <- function() {
  files <- list.files(popstats_dir, pattern = "\\.popstats\\.tsv$", full.names = TRUE)
  if (length(files) == 0) { message("[fst_plots] No popstats files"); return(invisible()) }

  all_dt <- rbindlist(lapply(files, function(f) {
    dt <- fread(f, skip = 1)
    dt[, mid_bp := (start + end) / 2]
    dt
  }), fill = TRUE)

  fst_cols <- detect_pairwise(all_dt, "Fst")
  if (length(fst_cols) == 0) { message("[fst_plots] No Fst columns"); return(invisible()) }

  chroms <- sort(unique(all_dt$chrom))
  all_dt[, chrom := factor(chrom, levels = chroms)]

  for (fc in fst_cols) {
    contrast_name <- sub("^Fst_", "", fc)
    p <- ggplot(all_dt, aes(x = mid_bp / 1e6, y = .data[[fc]])) +
      geom_line(linewidth = 0.2, color = "#2166AC", alpha = 0.7) +
      facet_wrap(~chrom, scales = "free_x", ncol = 4) +
      theme_bw(base_size = 9) +
      theme(strip.text = element_text(size = 7),
            axis.text.x = element_text(size = 6)) +
      labs(title = paste0("Hudson Fst: ", contrast_name),
           x = "Position (Mb)", y = expression(F[ST]))

    ht <- ceiling(length(chroms) / 4) * 2.5 + 1
    ggsave(file.path(outdir, paste0("genome_Fst_", contrast_name, ".pdf")),
           p, width = 14, height = ht)
    message("[fst_plots] genome Fst: ", contrast_name)
  }
}

# =============================================================================
# DISPATCH
# =============================================================================

message("[fst_plots] Mode: ", mode)
switch(mode,
  "candidate_stack" = plot_candidate_stack(),
  "genome_fst"      = plot_genome_fst(),
  message("[fst_plots] Unknown mode: ", mode)
)
message("[fst_plots] Done. Output: ", outdir)
