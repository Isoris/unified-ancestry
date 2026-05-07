#!/usr/bin/env Rscript
# =============================================================================
# instant_q.R — R wrapper for the Instant Q Engine (v12.1 REWIRED)
#
# v12.1 changes:
#   1. Sources sample_map.R for robust Ind↔CGA handling
#   2. get_Q() always returns sample_id with CGA names
#   3. Precompute cache files include sample_id column with CGA names
#   4. configure_instant_q() accepts sample_map_file param
#   5. get_region_stats() included here (lightweight version)
#      For the full dispatcher, use region_stats_dispatcher.R
#
# Usage:
#   source("instant_q.R")
#   configure_instant_q(config_file = "00_ancestry_config.sh")
#   q <- get_Q("C_gar_LG01", start = 35e6, end = 41e6)
#   summary <- get_Q_summary("C_gar_LG01")
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
EPS_Q <- 1e-10

# =============================================================================
# CONFIGURATION
# =============================================================================

.iq_env <- new.env(parent = emptyenv())

configure_instant_q <- function(
    config_file    = NULL,
    binary         = NULL,
    beagle_dir     = NULL,
    fopt_file      = NULL,
    qinit_file     = NULL,
    cache_dir      = NULL,
    sample_list    = NULL,
    sample_map_file = NULL,
    window_size    = 100L,
    window_step    = 20L,
    em_iter        = 20L,
    K              = NULL,
    ncores         = 4L
) {
  # Parse shell config if provided
  if (!is.null(config_file) && file.exists(config_file)) {
    cfg <- parse_shell_config(config_file)
    beagle_dir      <- beagle_dir      %||% cfg$BEAGLE_DIR %||% cfg$STEP1_BEAGLE_DIR
    fopt_file       <- fopt_file       %||% cfg$BEST_FOPT
    qinit_file      <- qinit_file      %||% cfg$BEST_QOPT
    cache_dir       <- cache_dir       %||% cfg$LOCAL_Q_DIR
    sample_list     <- sample_list     %||% cfg$SAMPLE_LIST %||% cfg$STEP1_SAMPLE_LIST
    binary          <- binary          %||% cfg$INSTANT_Q_BIN
    sample_map_file <- sample_map_file %||% cfg$SAMPLE_MAP_R
    # Canonical K — what the inversion precomp RDS flattens in. Resolved
    # from config's CANONICAL_K (preferred) or DEFAULT_K (fallback).
    .canon_from_cfg <- cfg$CANONICAL_K %||% cfg$DEFAULT_K
  } else {
    .canon_from_cfg <- NULL
  }

  .iq_env$binary      <- binary %||% find_instant_q_binary()
  .iq_env$beagle_dir  <- beagle_dir
  .iq_env$fopt_file   <- fopt_file
  .iq_env$qinit_file  <- qinit_file
  .iq_env$cache_dir   <- cache_dir %||% "local_Q"
  .iq_env$sample_list <- sample_list
  .iq_env$window_size <- as.integer(window_size)
  .iq_env$window_step <- as.integer(window_step)
  .iq_env$em_iter     <- as.integer(em_iter)
  .iq_env$ncores      <- as.integer(ncores)
  # Canonical K: resolved from config's CANONICAL_K (preferred) or the K arg
  # passed to configure_instant_q, or a hard-coded 8 default.
  .iq_env$canonical_K <- as.integer(K %||% .canon_from_cfg %||% 8L)
  # 2026-04-17 chat-16: sample-set identity is now a registered group_id
  # from sample_registry (typically "all_226"), not a content hash. The
  # group name is read from the SAMPLE_GROUP env var set by the launcher
  # or defaults to "all_226". See registries/DATABASE_DESIGN.md for the
  # rationale — this replaces the chat-15 .compute_sample_set_tag() which
  # produced R-vs-bash hash mismatches, silently no-op'ing every downstream
  # merge.
  .iq_env$sample_group <- Sys.getenv("SAMPLE_GROUP", "all_226")
  if (!nzchar(.iq_env$sample_group)) .iq_env$sample_group <- "all_226"
  # Back-compat: keep sample_set_tag field name so existing callers that
  # pass it through still work. New code should use sample_group.
  .iq_env$sample_set_tag <- .iq_env$sample_group

  # ── Load sample_map.R for robust name resolution ──
  .iq_env$smap <- NULL
  if (!is.null(sample_map_file) && file.exists(sample_map_file)) {
    tryCatch({
      source(sample_map_file, local = TRUE)
      samples_ind <- Sys.getenv("SAMPLES_IND", "")
      .iq_env$smap <- load_sample_map(
        if (nzchar(samples_ind)) samples_ind else NULL
      )
    }, error = function(e) {
      message("[instant_q.R] sample_map.R load failed: ", conditionMessage(e))
    })
  } else if (exists("smap", envir = .GlobalEnv)) {
    # Already loaded by load_bridge.R
    .iq_env$smap <- get("smap", envir = .GlobalEnv)
  }

  # ── Load sample names (prefer smap, fallback to raw readLines) ──
  if (!is.null(.iq_env$smap)) {
    .iq_env$samples <- .iq_env$smap$real_names  # Always CGA
    message("[instant_q.R] Sample names from sample_map: ",
            .iq_env$samples[1], " .. ", tail(.iq_env$samples, 1))
  } else if (!is.null(sample_list) && file.exists(sample_list)) {
    .iq_env$samples <- readLines(sample_list)
    .iq_env$samples <- .iq_env$samples[nzchar(.iq_env$samples)]
  }

  # Auto-detect K from qinit
  if (!is.null(qinit_file) && file.exists(qinit_file)) {
    first_line <- readLines(qinit_file, n = 1)
    .iq_env$K <- length(strsplit(trimws(first_line), "\\s+")[[1]])
  }
  if (!is.null(K)) .iq_env$K <- as.integer(K)

  dir.create(.iq_env$cache_dir, recursive = TRUE, showWarnings = FALSE)

  message("[instant_q.R] Configured: K=", .iq_env$K %||% "auto",
          ", cache=", .iq_env$cache_dir,
          ", binary=", .iq_env$binary %||% "NOT FOUND",
          ", smap=", if (!is.null(.iq_env$smap)) "YES" else "NO")
  invisible(TRUE)
}

parse_shell_config <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines)]
  lines <- lines[grepl("=", lines)]
  cfg <- list()
  for (ln in lines) {
    ln <- sub("^\\s*export\\s+", "", ln)
    eq <- regexpr("=", ln)
    if (eq > 0) {
      key <- trimws(substr(ln, 1, eq - 1))
      val <- trimws(substr(ln, eq + 1, nchar(ln)))
      val <- gsub('^"|"$', '', val)
      val <- gsub("^'|'$", '', val)
      cfg[[key]] <- val
    }
  }
  cfg
}

find_instant_q_binary <- function() {
  candidates <- c(
    Sys.getenv("INSTANT_Q_BIN", ""),
    file.path(dirname(sys.frame(1)$ofile %||% "."), "src", "instant_q"),
    file.path(dirname(sys.frame(1)$ofile %||% "."), "instant_q"),
    "instant_q",
    file.path(Sys.getenv("HOME"), ".local", "bin", "instant_q")
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(p)) return(normalizePath(p))
  }
  w <- Sys.which("instant_q")
  if (nzchar(w)) return(w)
  NULL
}

# =============================================================================
# SAMPLE NAME RESOLUTION (v12.1)
# =============================================================================

#' Ensure sample_ids are CGA names. Uses smap if available, else returns as-is.
ensure_cga_names <- function(names_vec) {
  if (is.null(names_vec) || length(names_vec) == 0) return(names_vec)
  sm <- .iq_env$smap
  if (!is.null(sm) && sm$detect(names_vec[1]) == "ind") {
    return(sm$to_real_vec(names_vec))
  }
  names_vec
}

#' Get the canonical CGA sample names (226)
get_sample_ids_iq <- function() {
  if (!is.null(.iq_env$smap)) return(.iq_env$smap$real_names)
  .iq_env$samples %||% character(0)
}

# =============================================================================
# BEAGLE FILE RESOLUTION
# =============================================================================

resolve_beagle <- function(chr) {
  bd <- .iq_env$beagle_dir
  if (is.null(bd)) stop("[instant_q.R] beagle_dir not configured")
  patterns <- c(
    paste0("*", chr, "*.beagle.gz"),
    paste0("*", chr, "*.beagle"),
    paste0("catfish.", chr, ".rf.txt.thin_500.beagle.gz"),
    paste0("catfish.", chr, "*.beagle.gz")
  )
  for (pat in patterns) {
    files <- Sys.glob(file.path(bd, pat))
    if (length(files) > 0) return(files[1])
  }
  stop("[instant_q.R] No BEAGLE found for ", chr, " in ", bd)
}

# =============================================================================
# CACHE MANAGEMENT
# =============================================================================

# ── Sample-set fingerprint helpers ───────────────────────────────────────────
# Cache files are tagged with a short sample-set identifier so that running
# the same chrom/K against a different sample subset doesn't silently
# overwrite the prior result. Tag format: "N<count>_<6-char-sha>", e.g.
# "N226_a3f2c1". Computed once at configure time from the sorted sample list.
#
# If sample_list is missing or empty, the tag falls back to "unknown". Legacy
# files (no tag in the filename) are read with an "N226_legacy" assumption —
# the instant_q engine today always processes whichever samples are columns
# in the BEAGLE file, and production has been 226-sample throughout, so this
# is the safe assumption for pre-2026-04-17 caches.
# 2026-04-17 chat-16: DEPRECATED. The chat-15 implementation of this function
# hashed the sample list with digest::digest(ids, algo="sha1") which serializes
# the R object before hashing — producing a hash that did NOT match the bash
# launcher's `sort | sha1sum | cut -c1-6` over the same list. Result: every
# cache read via this tag silently returned NA (no file found), every
# merge_local_Q_into_invlikeness() silently no-op'd, every reg$compute$ancestry_*
# method silently returned NULL.
#
# Replacement: sample-set identity is now a registered group_id from the
# sample_registry (read from env SAMPLE_GROUP, default "all_226"). The group
# name is a single source of truth — no hash to disagree on. See
# registries/DATABASE_DESIGN.md for the design.
#
# This function is kept as a shim that returns the current sample_group so
# old callers that reference it don't break mid-migration.
.compute_sample_set_tag <- function(sample_list_path) {
  sg <- Sys.getenv("SAMPLE_GROUP", "all_226")
  if (!nzchar(sg)) sg <- "all_226"
  sg
}

# Public accessor in case callers want to pass the tag into downstream file
# naming (e.g. stats_cache files).
get_sample_set_tag <- function() {
  .iq_env$sample_set_tag %||% "unknown"
}

# Chat-16: preferred accessor — returns the registered group_id (e.g. "all_226")
# that this cache was computed over. Used by results_registry put_* methods
# when they need to stamp who_1.group_id on a new row.
get_active_sample_group <- function() {
  .iq_env$sample_group %||% Sys.getenv("SAMPLE_GROUP", "all_226")
}
# 2026-04-17 rewrite: the cache is sharded by K under <LOCAL_Q_DIR>/K<NN>/ so
# a full K=2..20 sweep can coexist. The canonical K (flattened into the
# inversion precomp RDS) defaults to .iq_env$canonical_K which configure_instant_q
# reads from CANONICAL_K or DEFAULT_K in the ancestry config.
#
# Backward compatibility: if a caller passes K=NULL we still resolve to
# canonical. On read, we also try the legacy flat layout <LOCAL_Q_DIR>/<chr>.*
# if the K-sharded path is missing — this lets pre-2026-04-17 caches keep
# working without a rerun. New writes always go to the sharded path.
.resolve_K <- function(K) {
  if (is.null(K)) {
    k <- .iq_env$canonical_K %||% 8L
    return(as.integer(k))
  }
  as.integer(K)
}
.k_subdir <- function(K) {
  K <- .resolve_K(K)
  file.path(.iq_env$cache_dir, sprintf("K%02d", K))
}
cache_path_summary <- function(chr, K = NULL, sample_set = NULL) {
  set_tag <- sample_set %||% get_sample_set_tag()
  file.path(.k_subdir(K),
            paste0(chr, ".", set_tag, ".local_Q_summary.tsv.gz"))
}
cache_path_samples <- function(chr, K = NULL, sample_set = NULL) {
  set_tag <- sample_set %||% get_sample_set_tag()
  file.path(.k_subdir(K),
            paste0(chr, ".", set_tag, ".local_Q_samples.tsv.gz"))
}
cache_path_meta <- function(chr, K = NULL, sample_set = NULL) {
  set_tag <- sample_set %||% get_sample_set_tag()
  file.path(.k_subdir(K),
            paste0(chr, ".", set_tag, ".local_Q_meta.tsv"))
}
# Legacy flat-layout paths — canonical K + untagged sample set (pre-chat-15).
# Only used for read fallback; new writes never land here.
.legacy_path_summary <- function(chr) {
  file.path(.iq_env$cache_dir, paste0(chr, ".local_Q_summary.tsv.gz"))
}
.legacy_path_samples <- function(chr) {
  file.path(.iq_env$cache_dir, paste0(chr, ".local_Q_samples.tsv.gz"))
}
# Pre-sample-set-tag path inside a K shard (post chat-15 K-aware rewrite,
# pre chat-15 sample-set-tag addition). Also tolerated on read.
.k_shard_untagged_summary <- function(chr, K = NULL) {
  file.path(.k_subdir(K), paste0(chr, ".local_Q_summary.tsv.gz"))
}
.k_shard_untagged_samples <- function(chr, K = NULL) {
  file.path(.k_subdir(K), paste0(chr, ".local_Q_samples.tsv.gz"))
}
# Resolve to whichever layout has the file. Preference order:
#   1. K-sharded + sample-set tag (current layout)
#   2. K-sharded + untagged    (transient intermediate)
#   3. Legacy flat (canonical K only)
# Each step also tries the .tsv variant (no .gz).
resolve_cache_summary <- function(chr, K = NULL, sample_set = NULL) {
  for (p in c(cache_path_summary(chr, K, sample_set),
              sub("\\.gz$", "", cache_path_summary(chr, K, sample_set)),
              .k_shard_untagged_summary(chr, K),
              sub("\\.gz$", "", .k_shard_untagged_summary(chr, K)))) {
    if (file.exists(p)) return(p)
  }
  if (is.null(K) || .resolve_K(K) == (.iq_env$canonical_K %||% 8L)) {
    for (p in c(.legacy_path_summary(chr),
                sub("\\.gz$", "", .legacy_path_summary(chr)))) {
      if (file.exists(p)) return(p)
    }
  }
  NA_character_
}
resolve_cache_samples <- function(chr, K = NULL, sample_set = NULL) {
  for (p in c(cache_path_samples(chr, K, sample_set),
              sub("\\.gz$", "", cache_path_samples(chr, K, sample_set)),
              .k_shard_untagged_samples(chr, K),
              sub("\\.gz$", "", .k_shard_untagged_samples(chr, K)))) {
    if (file.exists(p)) return(p)
  }
  if (is.null(K) || .resolve_K(K) == (.iq_env$canonical_K %||% 8L)) {
    for (p in c(.legacy_path_samples(chr),
                sub("\\.gz$", "", .legacy_path_samples(chr)))) {
      if (file.exists(p)) return(p)
    }
  }
  NA_character_
}
cache_exists <- function(chr, K = NULL, sample_set = NULL) {
  !is.na(resolve_cache_summary(chr, K, sample_set))
}
update_manifest <- function(chr, K = NULL, sample_set = NULL) {
  K <- .resolve_K(K)
  set_tag <- sample_set %||% get_sample_set_tag()
  manifest_file <- file.path(.iq_env$cache_dir, "manifest.tsv")
  row <- data.table(
    chrom = chr,
    K     = K,
    sample_set = set_tag,
    summary_file = cache_path_summary(chr, K, set_tag),
    samples_file = cache_path_samples(chr, K, set_tag),
    timestamp = as.character(Sys.time())
  )
  if (file.exists(manifest_file)) {
    old <- fread(manifest_file)
    if (!"K" %in% names(old)) old[, K := .iq_env$canonical_K %||% 8L]
    if (!"sample_set" %in% names(old)) old[, sample_set := "N226_legacy"]
    old <- old[!(chrom == chr & K == row$K & sample_set == row$sample_set)]
    fwrite(rbind(old, row, fill = TRUE), manifest_file, sep = "\t")
  } else {
    fwrite(row, manifest_file, sep = "\t")
  }
}

# =============================================================================
# get_Q() — Single region query
# =============================================================================

get_Q <- function(chr, start, end, K = NULL, force = FALSE) {
  bin <- .iq_env$binary
  if (is.null(bin) || !file.exists(bin)) {
    message("[instant_q.R] C++ binary not found, falling back to R implementation")
    return(get_Q_R(chr, start, end))
  }

  K_eff <- .resolve_K(K)
  # Resolve per-K qinit/fopt the same way get_Q_precompute does
  canon <- .iq_env$canonical_K %||% 8L
  qinit_use <- .iq_env$qinit_file
  fopt_use  <- .iq_env$fopt_file
  if (K_eff != canon) {
    pad_canon <- sprintf("K%02d", canon)
    pad_K     <- sprintf("K%02d", K_eff)
    if (!is.null(qinit_use)) qinit_use <- sub(pad_canon, pad_K, qinit_use, fixed = TRUE)
    if (!is.null(fopt_use))  fopt_use  <- sub(pad_canon, pad_K, fopt_use,  fixed = TRUE)
  }

  beagle <- resolve_beagle(chr)
  cmd <- sprintf(
    "%s --beagle %s --fopt %s --qinit %s --chr %s --start %d --end %d --em_iter %d",
    shQuote(bin), shQuote(beagle), shQuote(fopt_use),
    shQuote(qinit_use), shQuote(chr),
    as.integer(start), as.integer(end), .iq_env$em_iter
  )

  message("[instant_q.R] ", cmd)
  result <- fread(cmd = cmd, header = TRUE, sep = "\t")

  # v12.1: Always assign CGA sample names
  samples <- get_sample_ids_iq()
  if (nrow(result) > 0 && length(samples) > 0) {
    n <- min(nrow(result), length(samples))
    if ("sample_id" %in% names(result)) {
      result$sample_id[1:n] <- samples[1:n]
    } else {
      result[, sample_id := c(samples[1:n],
                               rep(NA_character_, nrow(result) - n))]
    }
  }

  result
}

# =============================================================================
# get_Q_summary()
# =============================================================================

get_Q_summary <- function(chr, K = NULL, force = FALSE) {
  cf <- resolve_cache_summary(chr, K)
  if (!is.na(cf) && !force) {
    dt <- fread(cf)
    message("[instant_q.R] Loaded cached summary for ", chr,
            " K=", .resolve_K(K), ": ", nrow(dt), " windows (", cf, ")")
    return(dt)
  }
  message("[instant_q.R] No cache for ", chr, " K=", .resolve_K(K),
          " — triggering precompute")
  get_Q_precompute(chr, K = K, force = force)
  cf2 <- resolve_cache_summary(chr, K)
  if (!is.na(cf2)) return(fread(cf2))
  warning("[instant_q.R] Precompute failed for ", chr, " K=", .resolve_K(K))
  data.table()
}

# =============================================================================
# get_Q_precompute()
# =============================================================================

get_Q_precompute <- function(chr = NULL, K = NULL, sample_output = FALSE,
                               force = FALSE,
                               qinit_file = NULL, fopt_file = NULL) {
  bin <- .iq_env$binary
  if (is.null(bin) || !file.exists(bin)) {
    message("[instant_q.R] C++ binary not found, falling back to R precompute")
    return(get_Q_precompute_R(chr, force))
  }

  K_eff <- .resolve_K(K)
  if (!is.null(chr) && cache_exists(chr, K_eff) && !force) {
    message("[instant_q.R] Cache exists for ", chr, " K=", K_eff,
            " (use force=TRUE to recompute)")
    return(invisible(NULL))
  }

  # Resolve per-K qinit/fopt. If the caller passed explicit paths, honour them.
  # Otherwise derive from the configured canonical paths by swapping K<NN> in
  # the filename (matches the BEST_QOPT / BEST_FOPT naming in 00_ancestry_config.sh).
  .swap_K_in_filename <- function(path, K_eff) {
    if (is.null(path) || !nzchar(path)) return(path)
    canon <- .iq_env$canonical_K %||% 8L
    if (K_eff == canon) return(path)
    pad_canon <- sprintf("K%02d", canon)
    pad_K     <- sprintf("K%02d", K_eff)
    sub(pad_canon, pad_K, path, fixed = TRUE)
  }
  qinit_use <- qinit_file %||% .swap_K_in_filename(.iq_env$qinit_file, K_eff)
  fopt_use  <- fopt_file  %||% .swap_K_in_filename(.iq_env$fopt_file,  K_eff)
  if (!is.null(qinit_use) && !file.exists(qinit_use)) {
    warning("[instant_q.R] qinit for K=", K_eff, " not found: ", qinit_use,
            " — precompute will fail")
  }

  beagle <- resolve_beagle(chr %||% stop("chr required"))
  tmpdir <- tempfile(pattern = "iq_")
  dir.create(tmpdir, showWarnings = FALSE)

  cmd <- sprintf(
    "%s --beagle %s --fopt %s --qinit %s --precompute --outdir %s --chr %s --window_size %d --window_step %d --em_iter %d --ncores %d%s",
    shQuote(bin), shQuote(beagle), shQuote(fopt_use),
    shQuote(qinit_use), shQuote(tmpdir), shQuote(chr),
    .iq_env$window_size, .iq_env$window_step, .iq_env$em_iter,
    .iq_env$ncores, if (sample_output) " --sample_output" else ""
  )

  message("[instant_q.R] ", cmd)
  ret <- system(cmd)
  if (ret != 0) {
    warning("[instant_q.R] C++ engine returned ", ret)
    return(invisible(NULL))
  }

  # K-sharded output dir
  out_shard <- .k_subdir(K_eff)
  dir.create(out_shard, recursive = TRUE, showWarnings = FALSE)

  summary_raw <- file.path(tmpdir, paste0(chr, ".local_Q_summary.tsv"))
  if (file.exists(summary_raw)) {
    dt <- fread(summary_raw)
    fwrite(dt, cache_path_summary(chr, K_eff), sep = "\t")
    message("[instant_q.R] Cached summary K=", K_eff, ": ", nrow(dt), " windows")
  }

  if (sample_output) {
    sample_raw <- file.path(tmpdir, paste0(chr, ".local_Q_samples.tsv"))
    if (file.exists(sample_raw)) {
      dt <- fread(sample_raw)
      samples <- get_sample_ids_iq()
      if ("sample_idx" %in% names(dt) && length(samples) > 0) {
        dt[, sample_id := samples[sample_idx + 1L]]
      } else if (!"sample_id" %in% names(dt) && length(samples) > 0) {
        n_samp <- length(samples)
        if (nrow(dt) %% n_samp == 0) {
          dt[, sample_id := rep(samples, nrow(dt) %/% n_samp)]
        }
      }
      fwrite(dt, cache_path_samples(chr, K_eff), sep = "\t")
    }
  }

  meta_raw <- file.path(tmpdir, paste0(chr, ".local_Q_meta.tsv"))
  if (file.exists(meta_raw)) {
    file.copy(meta_raw, cache_path_meta(chr, K_eff), overwrite = TRUE)
  }

  update_manifest(chr, K_eff)
  unlink(tmpdir, recursive = TRUE)
  invisible(NULL)
}

# =============================================================================
# R FALLBACK
# =============================================================================

get_Q_R <- function(chr, start, end) {
  beagle <- resolve_beagle(chr)
  message("[instant_q.R] R fallback: loading BEAGLE for region ", chr, ":", start, "-", end)

  bgl <- fread(beagle, nThread = 4L)
  markers <- bgl[[1]]
  positions <- as.integer(sub(".*_", "", markers))

  keep <- which(positions >= start & positions <= end)
  if (length(keep) < 10) {
    warning("[instant_q.R] Only ", length(keep), " sites in region")
    return(data.table())
  }

  gl_mat <- as.matrix(bgl[keep, 4:ncol(bgl)])
  n_sites <- nrow(gl_mat)
  n_ind <- ncol(gl_mat) %/% 3

  q_init <- as.matrix(fread(.iq_env$qinit_file, header = FALSE))
  K <- ncol(q_init)
  n_use <- min(n_ind, nrow(q_init))

  f_mat <- NULL
  if (!is.null(.iq_env$fopt_file) && file.exists(.iq_env$fopt_file)) {
    f_all <- as.matrix(fread(.iq_env$fopt_file, header = FALSE))
    if (max(keep) <= nrow(f_all)) {
      f_mat <- f_all[keep, , drop = FALSE]
    }
  }
  if (is.null(f_mat)) f_mat <- matrix(0.25, n_sites, K)

  Q <- compute_window_Q_R(gl_mat[, seq_len(3 * n_use), drop = FALSE],
                            f_mat, q_init[seq_len(n_use), , drop = FALSE],
                            n_iter = .iq_env$em_iter %||% 20L)

  # v12.1: use CGA names
  samples <- get_sample_ids_iq()
  if (length(samples) == 0) samples <- paste0("ind", seq_len(n_use))
  samples <- samples[seq_len(n_use)]

  rows <- lapply(seq_len(n_use), function(i) {
    qv <- Q[i, ]
    qs <- sort(qv, decreasing = TRUE)
    H <- -sum(pmax(qv, EPS_Q) * log(pmax(qv, EPS_Q)))
    row <- as.list(setNames(round(qv, 4), paste0("Q", seq_len(K))))
    row$sample_id <- samples[i]  # CGA name
    row$max_q <- round(qs[1], 4)
    row$delta12 <- round(qs[1] - qs[2], 4)
    row$entropy <- round(H, 4)
    row$ena <- round(exp(H), 4)
    row$assigned_pop <- which.max(qv)
    row
  })
  rbindlist(rows, fill = TRUE)
}

compute_window_Q_R <- function(gl_win, f_win, q_prior, n_iter = 20L) {
  n_sites <- nrow(gl_win)
  n_ind <- nrow(q_prior)
  K <- ncol(q_prior)
  Q <- q_prior
  prev_ll <- -Inf

  for (iter in seq_len(n_iter)) {
    prodA <- matrix(0, n_ind, K)
    prodB <- matrix(0, n_ind, K)
    ll <- 0

    for (j in seq_len(n_sites)) {
      fj <- pmin(pmax(f_win[j, ], EPS_Q), 1 - EPS_Q)
      for (i in seq_len(n_ind)) {
        p <- sum(Q[i, ] * fj)
        p <- min(max(p, EPS_Q), 1 - EPS_Q)
        pp0 <- (1 - p)^2     * gl_win[j, 3*(i-1) + 1]
        pp1 <- 2 * p * (1-p) * gl_win[j, 3*(i-1) + 2]
        pp2 <- p^2           * gl_win[j, 3*(i-1) + 3]
        denom <- pp0 + pp1 + pp2
        if (denom < EPS_Q) next
        ll <- ll + log(denom)
        exp_g <- (pp1 + 2 * pp2) / denom
        prodA[i, ] <- prodA[i, ] + exp_g / (1 - p) * fj
        prodB[i, ] <- prodB[i, ] + (2 - exp_g) / p * (1 - fj)
      }
    }
    for (i in seq_len(n_ind)) {
      Q[i, ] <- Q[i, ] * (prodA[i, ] + prodB[i, ])
      Q[i, ] <- pmax(Q[i, ], EPS_Q)
      Q[i, ] <- Q[i, ] / sum(Q[i, ])
    }
    if (abs(ll - prev_ll) < 1e-4 * abs(prev_ll) && iter > 3) break
    prev_ll <- ll
  }
  Q
}

get_Q_precompute_R <- function(chr, force = FALSE) {
  if (cache_exists(chr) && !force) return(invisible(NULL))

  beagle <- resolve_beagle(chr)
  message("[instant_q.R] R precompute for ", chr, " (slow — compile C++ for speed)")

  bgl <- fread(beagle, nThread = 4L)
  markers <- bgl[[1]]
  positions <- as.integer(sub(".*_", "", markers))
  gl_mat <- as.matrix(bgl[, 4:ncol(bgl)])
  n_sites <- nrow(gl_mat)
  n_ind <- ncol(gl_mat) %/% 3

  q_init <- as.matrix(fread(.iq_env$qinit_file, header = FALSE))
  K <- ncol(q_init)
  n_use <- min(n_ind, nrow(q_init))

  f_all <- NULL
  if (!is.null(.iq_env$fopt_file) && file.exists(.iq_env$fopt_file)) {
    f_all <- as.matrix(fread(.iq_env$fopt_file, header = FALSE))
  }

  ws <- .iq_env$window_size; wt <- .iq_env$window_step
  n_win <- max(1, (n_sites - ws) %/% wt + 1)
  message("[instant_q.R] ", n_win, " windows, ", n_sites, " sites")

  summaries <- vector("list", n_win)
  t0 <- proc.time()

  for (wi in seq_len(n_win)) {
    s_idx <- (wi - 1) * wt + 1
    e_idx <- min(s_idx + ws - 1, n_sites)
    if (e_idx - s_idx + 1 < 10) next

    gl_sub <- gl_mat[s_idx:e_idx, seq_len(3 * n_use), drop = FALSE]
    f_sub <- if (!is.null(f_all) && e_idx <= nrow(f_all)) {
      f_all[s_idx:e_idx, , drop = FALSE]
    } else {
      matrix(0.25, e_idx - s_idx + 1, K)
    }

    Q <- compute_window_Q_R(gl_sub, f_sub,
                             q_init[seq_len(n_use), , drop = FALSE],
                             n_iter = .iq_env$em_iter)

    d12 <- vapply(seq_len(n_use), function(i) {
      qs <- sort(Q[i, ], decreasing = TRUE)
      qs[1] - qs[2]
    }, numeric(1))

    H_vec <- vapply(seq_len(n_use), function(i) {
      qv <- pmax(Q[i, ], EPS_Q)
      -sum(qv * log(qv))
    }, numeric(1))

    summaries[[wi]] <- data.table(
      window_id = paste0(chr, "_w", wi),
      chrom = chr,
      start_bp = positions[s_idx],
      end_bp = positions[e_idx],
      n_sites = e_idx - s_idx + 1L,
      mean_delta12 = round(mean(d12), 4),
      mean_entropy = round(mean(H_vec), 4),
      mean_ena = round(mean(exp(H_vec)), 4),
      sd_delta12 = round(sd(d12), 4)
    )

    if (wi %% 200 == 0) {
      elapsed <- round((proc.time() - t0)[3], 1)
      message("[instant_q.R]   window ", wi, "/", n_win,
              " (", round(wi / max(elapsed, 0.1), 1), " win/s)")
    }
  }

  result <- rbindlist(summaries[!vapply(summaries, is.null, logical(1))])
  dir.create(.iq_env$cache_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(result, cache_path_summary(chr), sep = "\t")
  update_manifest(chr)

  elapsed <- round((proc.time() - t0)[3], 1)
  message("[instant_q.R] ", chr, " done: ", nrow(result), " windows in ", elapsed, "s")
  invisible(NULL)
}

# =============================================================================
# LIGHTWEIGHT get_region_stats() (bundled for standalone use)
# =============================================================================
# For the full dispatcher with C binary routing, use region_stats_dispatcher.R

get_region_stats <- function(chr, start, end,
                              what = "Q",
                              groups = NULL) {
  result <- list()

  if ("Q" %in% what) {
    result$Q <- get_Q(chr, start, end)
  }

  need_dosage <- any(c("Fst", "Hobs", "dXY", "theta_pi", "HWE",
                        "dosage_matrix") %in% what)
  dos <- NULL

  if (need_dosage) {
    beagle <- resolve_beagle(chr)
    bgl <- fread(beagle, nThread = 4L)
    markers <- bgl[[1]]
    positions <- as.integer(sub(".*_", "", markers))
    keep <- which(positions >= start & positions <= end)

    if (length(keep) > 0) {
      gl_raw <- as.matrix(bgl[keep, 4:ncol(bgl)])
      n_sites <- nrow(gl_raw)
      n_ind <- ncol(gl_raw) %/% 3

      dos <- matrix(0, n_sites, n_ind)
      for (i in seq_len(n_ind)) {
        dos[, i] <- gl_raw[, 3*(i-1) + 2] + 2 * gl_raw[, 3*(i-1) + 3]
      }

      # v12.1: CGA column names
      colnames(dos) <- get_sample_ids_iq()[seq_len(n_ind)]
    }
  }

  if ("dosage_matrix" %in% what) {
    result$dosage_matrix <- dos
  }

  if ("Hobs" %in% what && !is.null(dos)) {
    het_mat <- (dos >= 0.3 & dos <= 1.7)
    per_sample_het <- colMeans(het_mat, na.rm = TRUE)
    result$Hobs <- data.table(
      sample_id = colnames(dos),
      Hobs = round(per_sample_het, 4)
    )
    result$Hobs_mean <- round(mean(per_sample_het, na.rm = TRUE), 4)
    result$Hobs_sd <- round(sd(per_sample_het, na.rm = TRUE), 4)
    result$Hobs_CV <- round(result$Hobs_sd / max(result$Hobs_mean, EPS_Q), 4)
  }

  if ("theta_pi" %in% what && !is.null(dos)) {
    n_ind <- ncol(dos)
    p_hat <- rowMeans(dos, na.rm = TRUE) / 2
    per_site_pi <- 2 * p_hat * (1 - p_hat) * n_ind / (n_ind - 1)
    result$theta_pi <- round(mean(per_site_pi, na.rm = TRUE), 6)
  }

  if ("Fst" %in% what && !is.null(dos) && !is.null(groups) && length(groups) >= 2) {
    all_samples <- colnames(dos)
    grp_names <- names(groups)
    fst_results <- list()
    for (a in seq_len(length(grp_names) - 1)) {
      for (b in (a + 1):length(grp_names)) {
        ga <- grp_names[a]; gb <- grp_names[b]
        idx_a <- which(all_samples %in% groups[[ga]])
        idx_b <- which(all_samples %in% groups[[gb]])
        if (length(idx_a) < 2 || length(idx_b) < 2) next
        p1 <- rowMeans(dos[, idx_a, drop = FALSE], na.rm = TRUE) / 2
        p2 <- rowMeans(dos[, idx_b, drop = FALSE], na.rm = TRUE) / 2
        n1 <- length(idx_a); n2 <- length(idx_b)
        num <- (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
        den <- p1*(1-p2) + p2*(1-p1)
        valid <- den > EPS_Q
        fst_val <- if (sum(valid) > 0) sum(num[valid]) / sum(den[valid]) else NA_real_
        fst_results[[paste0("Fst_", ga, "_", gb)]] <- round(max(0, fst_val), 4)
      }
    }
    result$Fst <- fst_results
  }

  if ("dXY" %in% what && !is.null(dos) && !is.null(groups) && length(groups) >= 2) {
    all_samples <- colnames(dos)
    grp_names <- names(groups)
    dxy_results <- list()
    for (a in seq_len(length(grp_names) - 1)) {
      for (b in (a + 1):length(grp_names)) {
        ga <- grp_names[a]; gb <- grp_names[b]
        idx_a <- which(all_samples %in% groups[[ga]])
        idx_b <- which(all_samples %in% groups[[gb]])
        if (length(idx_a) < 2 || length(idx_b) < 2) next
        p1 <- rowMeans(dos[, idx_a, drop = FALSE], na.rm = TRUE) / 2
        p2 <- rowMeans(dos[, idx_b, drop = FALSE], na.rm = TRUE) / 2
        dxy_per_site <- p1 * (1 - p2) + p2 * (1 - p1)
        dxy_results[[paste0("dXY_", ga, "_", gb)]] <- round(mean(dxy_per_site, na.rm = TRUE), 6)
      }
    }
    result$dXY <- dxy_results
  }

  result
}

# =============================================================================
# MERGE INTO INV_LIKENESS
# =============================================================================

merge_local_Q_into_invlikeness <- function(window_dt, K = NULL,
                                              cache_dir = NULL) {
  # 2026-04-17: K-aware merge. Reads the K-sharded cache (new layout) or the
  # flat layout (legacy, canonical-K only). Output columns are suffixed with
  # _K<NN> so multiple K levels can coexist in the same precomp RDS if a
  # later caller wants to merge more than one.
  #
  # New columns added:
  #   localQ_delta12_K<NN>  — mean |Q_max1 - Q_max2| per window
  #   localQ_entropy_K<NN>  — mean Shannon entropy of Q vector per window
  #   localQ_ena_K<NN>      — mean effective number of ancestries per window
  cache_dir <- cache_dir %||% .iq_env$cache_dir %||% "local_Q"
  K_eff <- .resolve_K(K)
  pad_K <- sprintf("K%02d", K_eff)

  manifest_file <- file.path(cache_dir, "manifest.tsv")
  if (!file.exists(manifest_file)) {
    # Not fatal: might still find summary files if the user seeded the cache
    # manually. Warn then keep going.
    message("[instant_q.R] No manifest at ", manifest_file,
            " — will scan for per-chrom files")
  } else {
    manifest <- fread(manifest_file)
    message("[instant_q.R] Merging K=", K_eff, " from ", nrow(manifest),
            " manifest rows (cache=", cache_dir, ")")
  }

  col_delta12 <- paste0("localQ_delta12_", pad_K)
  col_entropy <- paste0("localQ_entropy_", pad_K)
  col_ena     <- paste0("localQ_ena_",     pad_K)

  # Save .iq_env$cache_dir for resolve_cache_summary() calls; temporarily
  # point at the requested cache_dir.
  saved_cache_dir <- .iq_env$cache_dir
  .iq_env$cache_dir <- cache_dir
  on.exit(.iq_env$cache_dir <- saved_cache_dir, add = TRUE)

  n_merged <- 0L
  for (chr in unique(window_dt$chrom)) {
    sf <- resolve_cache_summary(chr, K_eff)
    if (is.na(sf)) {
      message("[instant_q.R]   ", chr, " K=", K_eff, ": no cache found, skipped")
      next
    }
    q_summ <- fread(sf)
    if (nrow(q_summ) == 0) next
    if ("window_id" %in% names(q_summ) && "global_window_id" %in% names(window_dt)) {
      q_summ[, global_window_id := window_id]
      q_sub <- q_summ[, list(global_window_id,
                               V_delta12 = mean_delta12,
                               V_entropy = mean_entropy,
                               V_ena     = mean_ena)]
      setnames(q_sub, c("V_delta12", "V_entropy", "V_ena"),
               c(col_delta12, col_entropy, col_ena))
      window_dt <- merge(window_dt, q_sub, by = "global_window_id",
                            all.x = TRUE, suffixes = c("", ".q"))
      n_merged <- n_merged + 1L
    }
  }
  n_filled <- if (col_delta12 %in% names(window_dt)) {
    sum(is.finite(window_dt[[col_delta12]]), na.rm = TRUE)
  } else 0L
  message("[instant_q.R] Merge K=", K_eff, ": ", n_merged, " chrom(s) merged, ",
          n_filled, "/", nrow(window_dt), " windows filled")
  window_dt
}

list_cached_Q <- function(cache_dir = NULL, K = NULL) {
  cache_dir <- cache_dir %||% .iq_env$cache_dir %||% "local_Q"
  mf <- file.path(cache_dir, "manifest.tsv")
  if (file.exists(mf)) {
    m <- fread(mf)
    if (!"K" %in% names(m)) m[, K := .iq_env$canonical_K %||% 8L]
    if (!is.null(K)) m <- m[K == .resolve_K(K)]
    return(m)
  }
  # Fallback: scan. Look in both sharded (K<NN>/) and legacy (flat) layouts.
  out <- list()
  K_dirs <- list.files(cache_dir, pattern = "^K[0-9]{2}$", full.names = TRUE)
  for (kd in K_dirs) {
    k_int <- as.integer(sub("^K", "", basename(kd)))
    files <- list.files(kd, pattern = "\\.local_Q_summary\\.", full.names = FALSE)
    if (length(files)) out[[basename(kd)]] <- data.table(
      chrom = sub("\\.local_Q_summary\\..*$", "", files), K = k_int
    )
  }
  flat <- list.files(cache_dir, pattern = "\\.local_Q_summary\\.", full.names = FALSE)
  if (length(flat)) out[["legacy"]] <- data.table(
    chrom = sub("\\.local_Q_summary\\..*$", "", flat),
    K = .iq_env$canonical_K %||% 8L
  )
  if (!length(out)) return(data.table(chrom = character(0), K = integer(0)))
  r <- rbindlist(out, fill = TRUE)
  if (!is.null(K)) r <- r[K == .resolve_K(K)]
  r
}
