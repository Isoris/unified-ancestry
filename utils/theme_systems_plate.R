#!/usr/bin/env Rscript

# =============================================================================
# theme_systems_plate.R — White Modular Systems Figure Theme
#
# Publication-quality scientific figure system for population genomics.
# Implements a complete visual language: theme, palettes, helpers for
# scatter/heatmap/track/boxplot/table/annotation, patchwork composers.
#
# Style: white printable background, modular panel architecture, strong
# hierarchy, information-dense but organized, publication-safe, restrained
# color accents. Feels like a scientific systems plate / modular analysis
# dashboard.
#
# Usage:
#   source("theme_systems_plate.R")
#   p <- ggplot(df, aes(x, y)) + geom_point() + theme_plate()
#   fig <- compose_systems_plate(list(pA, pB, pC), title = "Figure 1")
#
# Requires: ggplot2. Optional: patchwork, grid, gridExtra, gridtext.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

has_patchwork <- suppressWarnings(require(patchwork, quietly = TRUE))
has_gridExtra <- suppressWarnings(require(gridExtra, quietly = TRUE))
has_gridtext  <- suppressWarnings(require(gridtext, quietly = TRUE))

# =============================================================================
# §1  FOUNDATIONAL CONSTANTS
# =============================================================================

# ── Background tones ──
PLATE_BG        <- "#FFFFFF"       # main figure background: pure white
PANEL_BG        <- "#FAFAFA"       # panel interior: barely warm
PANEL_BORDER    <- "#D4D4D4"       # panel frame: light warm gray
STRIP_BG        <- "#EDEDED"       # facet strip background
STRIP_TEXT_COL  <- "#2D2D2D"       # facet strip text
MODULE_FRAME    <- "#C8C8C8"       # heavier frame for module grouping
GRID_COL        <- "#E8E8E8"       # major gridlines
GRID_MINOR_COL  <- "#F2F2F2"       # minor gridlines (rarely used)

# ── Text hierarchy ──
TEXT_TITLE       <- "#1A1A1A"      # near-black, not pure black
TEXT_SUBTITLE    <- "#5C5C5C"      # medium gray
TEXT_BODY        <- "#3A3A3A"      # body / axis titles
TEXT_LABEL       <- "#4A4A4A"      # axis tick labels
TEXT_CAPTION     <- "#8C8C8C"      # captions, notes, small annotations
TEXT_MUTED       <- "#ABABAB"      # watermarks, very secondary

# ── Font sizes (points) ──
SIZE_MAIN_TITLE  <- 14
SIZE_SUBTITLE    <- 10
SIZE_PANEL_LABEL <- 12             # A, B, C panel letters
SIZE_STRIP       <- 9              # facet strips
SIZE_AXIS_TITLE  <- 9.5
SIZE_AXIS_TEXT   <- 8
SIZE_LEGEND_TITLE <- 8.5
SIZE_LEGEND_TEXT <- 7.5
SIZE_CAPTION     <- 7
SIZE_ANNOTATION  <- 7

# ── Spacing (pt) ──
MARGIN_PLOT      <- margin(8, 8, 8, 8)
MARGIN_PANEL     <- margin(4, 4, 4, 4)

# ── Line weights ──
LINE_AXIS        <- 0.3
LINE_TICK        <- 0.25
LINE_PANEL_BORDER <- 0.4
LINE_DATA_DEFAULT <- 0.7
LINE_DATA_THIN   <- 0.4
LINE_RIBBON_ALPHA <- 0.15

# =============================================================================
# §2  CORE THEME
# =============================================================================

#' Systems Plate Theme for ggplot2
#'
#' @param base_size Base font size (default 9.5)
#' @param grid "x", "y", "xy", "none" — which major gridlines to show
#' @param panel_border Logical: draw panel border rectangle
#' @param strip_style "default", "header", "clean"
theme_plate <- function(base_size = 9.5, grid = "y", panel_border = FALSE,
                         strip_style = "default") {
  t <- theme_minimal(base_size = base_size) %+replace%
    theme(
      # ── Plot-level ──
      plot.background = element_rect(fill = PLATE_BG, color = NA),
      plot.title = element_text(
        size = SIZE_MAIN_TITLE, face = "bold", color = TEXT_TITLE,
        hjust = 0, margin = margin(0, 0, 4, 0)
      ),
      plot.subtitle = element_text(
        size = SIZE_SUBTITLE, color = TEXT_SUBTITLE,
        hjust = 0, margin = margin(0, 0, 8, 0)
      ),
      plot.caption = element_text(
        size = SIZE_CAPTION, color = TEXT_CAPTION,
        hjust = 0, margin = margin(6, 0, 0, 0)
      ),
      plot.margin = MARGIN_PLOT,
      plot.title.position = "plot",
      plot.caption.position = "plot",

      # ── Panel ──
      panel.background = element_rect(fill = PANEL_BG, color = NA),
      panel.border = if (panel_border) {
        element_rect(fill = NA, color = PANEL_BORDER, linewidth = LINE_PANEL_BORDER)
      } else {
        element_blank()
      },
      panel.spacing = unit(10, "pt"),

      # ── Grid ──
      panel.grid.major.x = if (grepl("x", grid)) {
        element_line(color = GRID_COL, linewidth = 0.2)
      } else element_blank(),
      panel.grid.major.y = if (grepl("y", grid)) {
        element_line(color = GRID_COL, linewidth = 0.2)
      } else element_blank(),
      panel.grid.minor = element_blank(),

      # ── Axes ──
      axis.title.x = element_text(
        size = SIZE_AXIS_TITLE, color = TEXT_BODY,
        margin = margin(6, 0, 0, 0)
      ),
      axis.title.y = element_text(
        size = SIZE_AXIS_TITLE, color = TEXT_BODY,
        margin = margin(0, 6, 0, 0)
      ),
      axis.text = element_text(size = SIZE_AXIS_TEXT, color = TEXT_LABEL),
      axis.ticks = element_line(color = PANEL_BORDER, linewidth = LINE_TICK),
      axis.ticks.length = unit(2.5, "pt"),
      axis.line = element_line(color = "#BFBFBF", linewidth = LINE_AXIS),

      # ── Strip (facets) ──
      strip.background = if (strip_style == "header") {
        element_rect(fill = STRIP_BG, color = PANEL_BORDER, linewidth = 0.3)
      } else if (strip_style == "clean") {
        element_rect(fill = "transparent", color = NA)
      } else {
        element_rect(fill = STRIP_BG, color = NA)
      },
      strip.text = element_text(
        size = SIZE_STRIP, face = "bold", color = STRIP_TEXT_COL,
        margin = margin(3, 3, 3, 3)
      ),

      # ── Legend ──
      legend.background = element_rect(fill = PLATE_BG, color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.key.size = unit(12, "pt"),
      legend.title = element_text(
        size = SIZE_LEGEND_TITLE, face = "bold", color = TEXT_BODY
      ),
      legend.text = element_text(size = SIZE_LEGEND_TEXT, color = TEXT_LABEL),
      legend.spacing.y = unit(2, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    )
  t
}

#' Compact variant — tighter for multi-panel patchwork
theme_plate_compact <- function(...) {
  theme_plate(base_size = 8.5, ...) %+replace%
    theme(
      plot.margin = margin(4, 4, 4, 4),
      plot.title = element_text(size = 10, face = "bold", color = TEXT_TITLE, hjust = 0),
      plot.subtitle = element_text(size = 8, color = TEXT_SUBTITLE, hjust = 0),
      panel.spacing = unit(6, "pt"),
      legend.key.size = unit(10, "pt")
    )
}

# =============================================================================
# §3  COLOR PALETTES
# =============================================================================

# ── §3A Metric color families (consistent across all figures) ──

pal_delta12    <- c(lo = "#D1E5F0", mid = "#4393C3", hi = "#08519C")    # blue family
pal_entropy    <- c(lo = "#E8E0F0", mid = "#9E7CB5", hi = "#54278F")    # purple family
pal_ena        <- c(lo = "#D5EFCF", mid = "#74C476", hi = "#006D2C")    # green family
pal_fst        <- c(lo = "#FEE5D9", mid = "#FB6A4A", hi = "#A50F15")    # red family
pal_theta      <- c(lo = "#DEEBF7", mid = "#6BAED6", hi = "#08306B")    # steel blue
pal_hobs       <- c(lo = "#FFF7BC", mid = "#FEC44F", hi = "#D95F0E")    # amber family
pal_tajima     <- c(neg = "#2166AC", zero = "#F7F7F7", pos = "#B2182B") # diverging

# ── §3B Categorical class colors ──

pal_inversion_class <- c(
  HOM_STD     = "#4393C3",    # cool blue
  HET         = "#7CAE7A",    # muted sage green
  HOM_INV     = "#D6604D",    # warm brick red
  AMBIGUOUS   = "#C0C0C0",    # neutral gray
  RECOMBINANT = "#F4A460"     # sandy
)

pal_ancestry_k8 <- c(
  Q1 = "#3B6FA0", Q2 = "#CF6839", Q3 = "#C44E52", Q4 = "#6BA08E",
  Q5 = "#5A8F4A", Q6 = "#C9A83E", Q7 = "#8B6DAD", Q8 = "#E8919C"
)

# ── §3C Six accent variants (same structure, different accent hue) ──

plate_accents <- list(
  editorial_blue = list(
    primary = "#2C5F8A", secondary = "#7BA7CC", tertiary = "#BDD4E7",
    accent = "#D6604D", neutral = "#7A7A7A", light = "#E8F0F8"
  ),
  deep_indigo = list(
    primary = "#3F3B6C", secondary = "#8E87C5", tertiary = "#CBC6E8",
    accent = "#D4744E", neutral = "#6E6E6E", light = "#EEEDF6"
  ),
  slate_teal = list(
    primary = "#2E6B62", secondary = "#6FAF9D", tertiary = "#B8DED2",
    accent = "#C75B3F", neutral = "#787878", light = "#E5F4F0"
  ),
  steel_lavender = list(
    primary = "#5C5480", secondary = "#9A92C4", tertiary = "#C9C4E0",
    accent = "#B8704A", neutral = "#757575", light = "#EDEBF4"
  ),
  warm_technical = list(
    primary = "#8C5E3C", secondary = "#C49573", tertiary = "#E0C9B5",
    accent = "#3D7A8A", neutral = "#7D7D7D", light = "#F5EDE5"
  ),
  neutral_scientific = list(
    primary = "#505050", secondary = "#8A8A8A", tertiary = "#C5C5C5",
    accent = "#3C7EA8", neutral = "#6E6E6E", light = "#F0F0F0"
  )
)

# Active accent (default: editorial_blue)
.plate_accent <- plate_accents$editorial_blue

#' Set the active accent palette
#' @param name One of: editorial_blue, deep_indigo, slate_teal,
#'   steel_lavender, warm_technical, neutral_scientific
set_plate_accent <- function(name) {
  if (name %in% names(plate_accents)) {
    .plate_accent <<- plate_accents[[name]]
    message("[plate] Accent set: ", name)
  } else {
    stop("Unknown accent: ", name, ". Options: ",
         paste(names(plate_accents), collapse = ", "))
  }
}

# ── §3D Continuous color scales (for heatmaps) ──

# LD heatmap: white → warm amber → deep brown (clean, not garish)
scale_fill_ld <- function(...) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FFF5EB", "#FDD49E", "#FDBB84",
               "#FC8D59", "#E34A33", "#B30000", "#4A0000"),
    ...
  )
}

# Difference heatmap: blue ← white → red
scale_fill_diff <- function(...) {
  scale_fill_gradientn(
    colors = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
               "#FDDBC7", "#EF8A62", "#B2182B"),
    ...
  )
}

# Similarity: white → slate blue
scale_fill_sim <- function(...) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#ECF0F6", "#C6DBEF", "#9ECAE1",
               "#6BAED6", "#3182BD", "#08519C"),
    ...
  )
}

# Hobs / frequency: white → warm gold
scale_fill_hobs <- function(...) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FFF7BC", "#FEE391", "#FEC44F",
               "#FE9929", "#D95F0E", "#993404"),
    ...
  )
}

# Entropy: white → indigo
scale_fill_entropy <- function(...) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#EFEDF5", "#DADAEB", "#BCBDDC",
               "#9E9AC8", "#756BB1", "#54278F"),
    ...
  )
}

# Purity / dominance: white → dark teal
scale_fill_purity <- function(...) {
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#E5F5F9", "#CCECE6", "#99D8C9",
               "#66C2A4", "#2CA25F", "#006D2C"),
    ...
  )
}

# =============================================================================
# §4  PANEL LABEL HELPER
# =============================================================================

#' Add a bold panel letter (A, B, C...) to a ggplot
#' @param p A ggplot object
#' @param label The panel letter (e.g. "A")
#' @param x Relative x position (0-1, default 0.02)
#' @param y Relative y position (0-1, default 0.98)
add_panel_label <- function(p, label, x = -0.02, y = 1.04) {
  p + labs(tag = label) +
    theme(plot.tag = element_text(
      size = SIZE_PANEL_LABEL, face = "bold", color = TEXT_TITLE,
      hjust = 0, vjust = 1
    ),
    plot.tag.position = c(x, y))
}

# =============================================================================
# §5  HEATMAP HELPERS
# =============================================================================

#' Systems-plate heatmap via ggplot2 geom_tile
#'
#' @param mat Named matrix (rows × cols)
#' @param title Heatmap title
#' @param fill_scale A scale_fill_* function
#' @param show_labels Show row/col labels if n < max_labels
#' @param max_labels Threshold for showing labels
#' @param annotation_row Named vector for row sidebar (optional)
#' @param annotation_col Named vector for column sidebar (optional)
plate_heatmap <- function(mat, title = "", fill_scale = scale_fill_sim(),
                           show_labels = TRUE, max_labels = 60,
                           annotation_row = NULL, annotation_col = NULL) {
  n_row <- nrow(mat); n_col <- ncol(mat)
  rn <- rownames(mat) %||% paste0("r", seq_len(n_row))
  cn <- colnames(mat) %||% paste0("c", seq_len(n_col))

  dt <- data.table(
    row = rep(rn, times = n_col),
    col = rep(cn, each = n_row),
    value = as.vector(mat)
  )
  dt[, row := factor(row, levels = rev(rn))]
  dt[, col := factor(col, levels = cn)]

  p <- ggplot(dt, aes(x = col, y = row, fill = value)) +
    geom_tile(color = NA, linewidth = 0) +
    fill_scale +
    theme_plate(grid = "none", panel_border = FALSE) +
    theme(
      axis.text.x = if (n_col <= max_labels && show_labels) {
        element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)
      } else element_blank(),
      axis.text.y = if (n_row <= max_labels && show_labels) {
        element_text(size = 6)
      } else element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      legend.key.height = unit(20, "pt"),
      legend.key.width = unit(8, "pt")
    ) +
    labs(title = title) +
    coord_fixed(ratio = 1)

  p
}

# =============================================================================
# §6  TRACK PLOT HELPERS (chromosome profiles)
# =============================================================================

#' Genomic track plot — line with optional ribbon and region shading
#'
#' @param dt data.table with pos_mb, value columns
#' @param ylab Y-axis label
#' @param color Line color (default: primary accent)
#' @param ribbon_col Ribbon column name for confidence band (optional)
#' @param region_start Region start in Mb (optional, for shading)
#' @param region_end Region end in Mb
plate_track <- function(dt, ylab = "", color = NULL, value_col = "value",
                         ribbon_lo = NULL, ribbon_hi = NULL,
                         region_start = NULL, region_end = NULL) {
  color <- color %||% .plate_accent$primary

  p <- ggplot(dt, aes(x = pos_mb, y = .data[[value_col]])) +
    theme_plate_compact(grid = "y", panel_border = TRUE)

  # Region shading
  if (!is.null(region_start) && !is.null(region_end)) {
    p <- p + annotate("rect",
      xmin = region_start, xmax = region_end,
      ymin = -Inf, ymax = Inf,
      fill = .plate_accent$light, color = NA
    )
  }

  # Ribbon
  if (!is.null(ribbon_lo) && !is.null(ribbon_hi)) {
    p <- p + geom_ribbon(
      aes(ymin = .data[[ribbon_lo]], ymax = .data[[ribbon_hi]]),
      fill = color, alpha = LINE_RIBBON_ALPHA
    )
  }

  # Line
  p <- p + geom_line(linewidth = LINE_DATA_DEFAULT, color = color, alpha = 0.85) +
    labs(x = "Position (Mb)", y = ylab)

  p
}

#' Multi-line track (e.g. Fst per pairwise comparison)
plate_track_multi <- function(dt, ylab = "", color_col = "contrast",
                               value_col = "value", palette = NULL,
                               region_start = NULL, region_end = NULL) {
  p <- ggplot(dt, aes(x = pos_mb, y = .data[[value_col]],
                        color = .data[[color_col]])) +
    theme_plate_compact(grid = "y", panel_border = TRUE)

  if (!is.null(region_start) && !is.null(region_end)) {
    p <- p + annotate("rect",
      xmin = region_start, xmax = region_end,
      ymin = -Inf, ymax = Inf,
      fill = .plate_accent$light, color = NA
    )
  }

  p <- p +
    geom_line(linewidth = LINE_DATA_DEFAULT, alpha = 0.85) +
    labs(x = "Position (Mb)", y = ylab, color = NULL)

  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }

  p + theme(legend.position = "top",
            legend.key.size = unit(10, "pt"),
            legend.text = element_text(size = 7))
}

# =============================================================================
# §7  BOXPLOT / VIOLIN HELPERS
# =============================================================================

#' Systems-plate boxplot
plate_boxplot <- function(dt, x_col, y_col, fill_col = NULL,
                           palette = NULL, ylab = "", xlab = "") {
  fill_col <- fill_col %||% x_col

  p <- ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]],
                        fill = .data[[fill_col]])) +
    geom_boxplot(
      width = 0.6, outlier.size = 0.6, outlier.alpha = 0.4,
      color = "#808080", linewidth = 0.35
    ) +
    theme_plate(grid = "y") +
    labs(x = xlab, y = ylab, fill = NULL)

  if (!is.null(palette)) {
    p <- p + scale_fill_manual(values = palette)
  } else {
    p <- p + scale_fill_manual(values = unname(pal_ancestry_k8))
  }

  p + theme(legend.position = "none")
}

#' Compact violin with inner boxplot
plate_violin <- function(dt, x_col, y_col, fill_col = NULL,
                          palette = NULL, ylab = "", xlab = "") {
  fill_col <- fill_col %||% x_col

  p <- ggplot(dt, aes(x = .data[[x_col]], y = .data[[y_col]],
                        fill = .data[[fill_col]])) +
    geom_violin(
      scale = "width", width = 0.75, color = "#A0A0A0", linewidth = 0.25,
      alpha = 0.7, trim = TRUE
    ) +
    geom_boxplot(
      width = 0.12, outlier.shape = NA, color = "#505050", linewidth = 0.3,
      fill = "white", alpha = 0.8
    ) +
    theme_plate(grid = "y") +
    labs(x = xlab, y = ylab, fill = NULL)

  if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)
  p + theme(legend.position = "none")
}

# =============================================================================
# §8  SCATTER HELPERS
# =============================================================================

#' Systems-plate scatter plot (e.g. PCA)
plate_scatter <- function(dt, x_col, y_col, color_col = NULL,
                           palette = NULL, xlab = "", ylab = "",
                           point_size = 1.8, point_alpha = 0.75) {
  aes_base <- aes(x = .data[[x_col]], y = .data[[y_col]])
  if (!is.null(color_col)) {
    aes_base <- aes(x = .data[[x_col]], y = .data[[y_col]],
                     color = .data[[color_col]])
  }

  p <- ggplot(dt, aes_base) +
    geom_point(size = point_size, alpha = point_alpha) +
    theme_plate(grid = "xy", panel_border = TRUE) +
    labs(x = xlab, y = ylab, color = NULL)

  if (!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }

  p
}

# =============================================================================
# §9  TABLE / SUMMARY BOX HELPERS
# =============================================================================

#' Create a styled table grob for embedding in patchwork
#'
#' @param dt data.table or data.frame
#' @param title Table title (shown as header row)
#' @param header_fill Fill color for header row
#' @param body_fill Fill color for body rows (alternating)
#' @return A tableGrob suitable for patchwork::wrap_elements()
plate_table <- function(dt, title = NULL, header_fill = "#E8E8E8",
                         body_fill = c("#FFFFFF", "#F6F6F6")) {
  if (!has_gridExtra) {
    message("[plate] gridExtra required for plate_table")
    return(grid::textGrob("(table requires gridExtra)", gp = gpar(fontsize = 8)))
  }

  n_rows <- nrow(dt); n_cols <- ncol(dt)

  # Theme for gridExtra tableGrob
  tt <- gridExtra::ttheme_minimal(
    base_size = 7.5,
    core = list(
      bg_params = list(
        fill = rep(body_fill, length.out = n_rows),
        col = "#D4D4D4", lwd = 0.3
      ),
      fg_params = list(fontsize = 7, col = TEXT_BODY)
    ),
    colhead = list(
      bg_params = list(fill = header_fill, col = "#D4D4D4", lwd = 0.3),
      fg_params = list(fontsize = 7.5, fontface = "bold", col = TEXT_TITLE)
    ),
    rowhead = list(
      fg_params = list(fontsize = 7, fontface = "bold", col = TEXT_BODY)
    )
  )

  g <- gridExtra::tableGrob(dt, rows = NULL, theme = tt)

  # Add title above if provided
  if (!is.null(title)) {
    title_grob <- grid::textGrob(
      title,
      gp = gpar(fontsize = 9, fontface = "bold", col = TEXT_TITLE),
      just = "left", x = unit(0.02, "npc")
    )
    g <- gtable::gtable_add_rows(g, heights = unit(1.2, "lines"), pos = 0)
    g <- gtable::gtable_add_grob(g, title_grob, t = 1, l = 1, r = ncol(g))
  }

  g
}

#' Key-metric summary box (e.g. "Fst = 0.42, dXY = 0.03, n = 85")
#'
#' @param metrics Named list of key metrics
#' @param title Box title
#' @return A grob
plate_metric_box <- function(metrics, title = "Summary") {
  lines <- paste0("**", names(metrics), "**: ", unlist(metrics))
  text <- paste(c(title, lines), collapse = "\n")

  if (has_gridtext) {
    # Rich text box
    grid::grobTree(
      grid::rectGrob(
        gp = gpar(fill = "#F5F5F5", col = PANEL_BORDER, lwd = 0.5),
        width = unit(0.95, "npc"), height = unit(0.9, "npc")
      ),
      gridtext::richtext_grob(
        text = paste0(
          "<b style='font-size:9pt;color:", TEXT_TITLE, "'>", title, "</b><br>",
          paste0("<span style='font-size:7.5pt;color:", TEXT_BODY, "'>",
                 paste0(names(metrics), ": ", unlist(metrics)), "</span>",
                 collapse = "<br>")
        ),
        x = unit(0.1, "npc"), y = unit(0.85, "npc"),
        hjust = 0, vjust = 1, padding = unit(c(4, 4, 4, 4), "pt")
      )
    )
  } else {
    # Fallback plain text
    grid::grobTree(
      grid::rectGrob(
        gp = gpar(fill = "#F5F5F5", col = PANEL_BORDER, lwd = 0.5)
      ),
      grid::textGrob(
        paste(c(title, paste0(names(metrics), ": ", unlist(metrics))),
              collapse = "\n"),
        x = 0.08, y = 0.85, hjust = 0, vjust = 1,
        gp = gpar(fontsize = 7.5, col = TEXT_BODY, lineheight = 1.3)
      )
    )
  }
}

#' Interpretation / key-finding callout box
plate_callout <- function(text, type = "finding") {
  bg <- switch(type,
    finding = "#EBF4FA",    # pale blue
    warning = "#FFF3E0",    # pale amber
    note    = "#F5F5F5",    # pale gray
    "#F5F5F5"
  )
  border <- switch(type,
    finding = "#3B6FA0",
    warning = "#E8903E",
    note    = "#A0A0A0",
    "#A0A0A0"
  )

  grid::grobTree(
    grid::rectGrob(
      gp = gpar(fill = bg, col = border, lwd = 0.6),
      width = unit(0.96, "npc"), height = unit(0.92, "npc")
    ),
    grid::textGrob(
      text, x = 0.06, y = 0.5, hjust = 0, vjust = 0.5,
      gp = gpar(fontsize = 7.5, col = TEXT_BODY, lineheight = 1.25)
    )
  )
}

# =============================================================================
# §10  ANNOTATION BAR HELPERS
# =============================================================================

#' Horizontal annotation strip (e.g. group identity bar above a heatmap)
plate_annotation_bar <- function(categories, palette, height = 0.8) {
  n <- length(categories)
  dt <- data.table(x = seq_len(n), cat = factor(categories, levels = unique(categories)))

  ggplot(dt, aes(x = x, y = 1, fill = cat)) +
    geom_tile(height = height, color = NA) +
    scale_fill_manual(values = palette) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(0, 0, 0, 0)) +
    coord_cartesian(expand = FALSE)
}

# =============================================================================
# §11  PATCHWORK COMPOSITION HELPERS
# =============================================================================

#' Compose a systems-plate figure from panels
#'
#' @param panels List of ggplot objects or grobs
#' @param labels Character vector of panel labels (A, B, C, ...)
#' @param title Figure title
#' @param subtitle Figure subtitle
#' @param layout Patchwork layout string or design
#' @param widths Relative column widths
#' @param heights Relative row heights
compose_systems_plate <- function(panels, labels = LETTERS[seq_along(panels)],
                                   title = NULL, subtitle = NULL,
                                   layout = NULL, widths = NULL, heights = NULL) {
  if (!has_patchwork) {
    message("[plate] patchwork required for compose_systems_plate")
    return(panels[[1]])
  }

  # Add panel labels
  for (i in seq_along(panels)) {
    if (inherits(panels[[i]], "gg")) {
      panels[[i]] <- add_panel_label(panels[[i]], labels[i])
    } else {
      # Wrap grob
      panels[[i]] <- patchwork::wrap_elements(panels[[i]])
    }
  }

  # Assemble
  combined <- panels[[1]]
  for (i in 2:length(panels)) {
    combined <- combined + panels[[i]]
  }

  if (!is.null(layout)) {
    combined <- combined + plot_layout(design = layout)
  } else {
    combined <- combined + plot_layout(
      ncol = min(length(panels), 3),
      widths = widths, heights = heights
    )
  }

  # Title annotation
  combined <- combined +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(
          size = SIZE_MAIN_TITLE + 2, face = "bold", color = TEXT_TITLE,
          hjust = 0, margin = margin(0, 0, 2, 0)
        ),
        plot.subtitle = element_text(
          size = SIZE_SUBTITLE + 1, color = TEXT_SUBTITLE,
          hjust = 0, margin = margin(0, 0, 10, 0)
        ),
        plot.background = element_rect(fill = PLATE_BG, color = NA)
      )
    )

  combined
}

#' Save a systems-plate figure at publication resolution
save_plate <- function(fig, filename, width = 10, height = 8, dpi = 300) {
  ggsave(filename, fig, width = width, height = height, dpi = dpi,
         bg = PLATE_BG, device = cairo_pdf)
  message("[plate] Saved: ", filename, " (", width, "×", height, " in)")
}

#' Save as both PDF and PNG
save_plate_dual <- function(fig, prefix, width = 10, height = 8, dpi = 300) {
  ggsave(paste0(prefix, ".pdf"), fig, width = width, height = height,
         bg = PLATE_BG, device = cairo_pdf)
  ggsave(paste0(prefix, ".png"), fig, width = width, height = height,
         dpi = dpi, bg = PLATE_BG)
  message("[plate] Saved: ", prefix, ".pdf + .png")
}

# =============================================================================
# §12  CONVENIENCE: IMAGE()-BASED HEATMAP FOR BASE R
# =============================================================================

#' Systems-plate heatmap using base R image() — for SNP×SNP matrices
#'
#' Matches the theme system but uses base R for speed with large matrices.
#' Produces publication-quality output compatible with the plate theme.
plate_image_heatmap <- function(mat, pos_mb = NULL, title = "",
                                 colors = NULL, zlim = NULL,
                                 region_start_mb = NULL, region_end_mb = NULL) {
  if (is.null(colors)) {
    colors <- colorRampPalette(c(
      "#FFFFFF", "#FFF5EB", "#FDD49E", "#FDBB84",
      "#FC8D59", "#E34A33", "#B30000", "#4A0000"
    ))(256)
  }
  if (is.null(zlim)) zlim <- c(0, quantile(mat[is.finite(mat)], 0.98, na.rm = TRUE))
  if (is.null(pos_mb)) pos_mb <- seq_len(nrow(mat))

  par(
    mar = c(3.5, 3.5, 2.5, 1),
    family = "sans",
    bg = PLATE_BG,
    fg = TEXT_BODY,
    col.axis = TEXT_LABEL,
    col.lab = TEXT_BODY,
    col.main = TEXT_TITLE,
    cex.main = 0.95,
    cex.axis = 0.7,
    cex.lab = 0.8
  )

  image(pos_mb, pos_mb, mat, col = colors, zlim = zlim,
        xlab = "Position (Mb)", ylab = "Position (Mb)",
        main = title, useRaster = TRUE, axes = FALSE)

  # Clean axes
  axis(1, col = "#BFBFBF", col.ticks = "#BFBFBF", lwd = LINE_AXIS, lwd.ticks = LINE_TICK)
  axis(2, col = "#BFBFBF", col.ticks = "#BFBFBF", lwd = LINE_AXIS, lwd.ticks = LINE_TICK)
  box(col = PANEL_BORDER, lwd = LINE_PANEL_BORDER)

  # Region shading lines
  if (!is.null(region_start_mb) && !is.null(region_end_mb)) {
    abline(v = c(region_start_mb, region_end_mb),
           col = paste0(.plate_accent$primary, "60"), lwd = 1, lty = 2)
    abline(h = c(region_start_mb, region_end_mb),
           col = paste0(.plate_accent$primary, "60"), lwd = 1, lty = 2)
  }
}

# =============================================================================
# §13  READY MESSAGE
# =============================================================================

message("[plate] White Modular Systems Figure Theme loaded")
message("[plate]   theme_plate(), theme_plate_compact()")
message("[plate]   plate_heatmap(), plate_track(), plate_scatter()")
message("[plate]   plate_boxplot(), plate_violin()")
message("[plate]   plate_table(), plate_metric_box(), plate_callout()")
message("[plate]   compose_systems_plate(), save_plate()")
message("[plate]   set_plate_accent('editorial_blue'|'deep_indigo'|...)")
message("[plate]   Active accent: editorial_blue")
