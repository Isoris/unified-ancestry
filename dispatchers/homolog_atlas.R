# =============================================================================
# homolog_atlas.R — R wrapper for the homolog index + codon_stats binaries.
#
# Use case: the unified-ancestry atlas page wants "click a gene → see all
# chains / paralogs / orthologs / dN/dS instantly". This file wraps the C
# binaries so the atlas code only sees tibbles.
#
# Exposed functions:
#
#   atlas_lookup(index, gene_id, source = "all", min_identity = NULL,
#                min_bitscore = NULL, limit = NULL)
#     → tibble of hits from homolog_index_query (one row per hit).
#
#   atlas_orthologs(index, gene_id, ...)
#     → cross-organism hits (miniprot source only).
#
#   atlas_paralogs(index, gene_id, ...)
#     → within-organism hits (diamond source only, self-hit removed).
#
#   codon_stats_run(pairs = NULL, fasta = NULL, method = "ng86",
#                   kappa = NULL, ncores = 1)
#     → tibble of per-pair dN/dS / 4D divergence (one row per CDS pair).
#
#   atlas_build(query, out_dir, targets, diamond_files, diamond_subjects,
#               mode = "chain", threads = 8, ...)
#     → invokes dispatchers/build_homolog_atlas.sh; returns index path.
#
#   atlas_health_check(index)
#     → quick sanity report (path, size_mb, binary_responds).
#
#   atlas_engines_dir(override = NULL)
#     → resolves the engines/ directory (env var, option, sibling fallback).
#
# Optional deps (graceful fallback): jsonlite, tibble.
# Required: base R only.
# =============================================================================

# ── Path resolution ─────────────────────────────────────────────────────────

#' Locate the unified-ancestry engines/ directory.
#'
#' Resolution order:
#'   1. \code{override} argument (if non-NULL, returned as-is)
#'   2. \code{getOption("unified.ancestry.engines_dir")}
#'   3. \code{Sys.getenv("UNIFIED_ANCESTRY_ENGINES_DIR")}
#'   4. Sibling \code{../engines} relative to the sourced script
#'   5. \code{engines}, \code{../engines}, \code{../../engines} from getwd()
#'
#' @param override optional explicit path
#' @return absolute path to engines directory
#' @export
atlas_engines_dir <- function(override = NULL) {
  if (!is.null(override) && nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  opt <- getOption("unified.ancestry.engines_dir", default = NULL)
  if (!is.null(opt) && dir.exists(opt)) return(normalizePath(opt))
  env <- Sys.getenv("UNIFIED_ANCESTRY_ENGINES_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))
  # Script-sibling fallback (only works if file was source()'d)
  this <- tryCatch(
    {
      frames <- sys.frames()
      ofile <- NULL
      for (i in rev(seq_along(frames))) {
        f <- frames[[i]]$ofile
        if (!is.null(f) && nzchar(f)) { ofile <- f; break }
      }
      ofile
    },
    error = function(e) NULL
  )
  if (!is.null(this)) {
    cand <- normalizePath(file.path(dirname(this), "..", "engines"), mustWork = FALSE)
    if (dir.exists(cand)) return(cand)
  }
  for (cand in c("engines", "../engines", "../../engines")) {
    if (dir.exists(cand)) return(normalizePath(cand))
  }
  stop("Cannot locate unified-ancestry engines/. Set UNIFIED_ANCESTRY_ENGINES_DIR ",
       "or options(unified.ancestry.engines_dir = '...').")
}

# ── Internal: invoke a C engine ─────────────────────────────────────────────

.run_engine <- function(name, args, engines_dir = NULL) {
  bin <- file.path(atlas_engines_dir(engines_dir), name)
  if (!file.exists(bin)) {
    stop(sprintf("Engine '%s' not built at %s. Run: make -C %s",
                 name, bin, dirname(bin)))
  }
  res <- suppressWarnings(system2(bin, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status")
  if (!is.null(status) && status != 0) {
    msg <- if (length(res) > 0) paste(res, collapse = "\n") else "(no output)"
    stop(sprintf("%s exited with status %d:\n%s", name, status, msg))
  }
  res
}

.as_tibble_safe <- function(df) {
  if (requireNamespace("tibble", quietly = TRUE)) tibble::as_tibble(df) else df
}

# ── Lookup ──────────────────────────────────────────────────────────────────

.empty_hits <- function() {
  .as_tibble_safe(data.frame(
    query = character(0), target = character(0),
    source_label = character(0), source_type = character(0),
    identity = numeric(0), bitscore = numeric(0), evalue = numeric(0),
    qstart = integer(0), qend = integer(0),
    tstart = integer(0), tend = integer(0),
    alnlen = integer(0), n_segments = integer(0), strand = character(0),
    stringsAsFactors = FALSE
  ))
}

#' Query the homolog index for all hits of a gene.
#'
#' @param index path to .holindx
#' @param gene_id gene id as it appears in the query proteome
#' @param source one of "all" / "diamond" / "miniprot"
#' @param min_identity NULL or fractional 0..1
#' @param min_bitscore NULL or numeric
#' @param limit NULL or integer N (returns top-N by bitscore desc)
#' @param engines_dir override engines/ location
#' @return tibble with columns documented in
#'   engines/schemas/homolog_index.query.output.schema.json
#' @export
atlas_lookup <- function(index, gene_id,
                          source = c("all", "diamond", "miniprot"),
                          min_identity = NULL,
                          min_bitscore = NULL,
                          limit = NULL,
                          engines_dir = NULL) {
  source <- match.arg(source)
  if (!file.exists(index)) stop("Index not found: ", index)
  args <- c("--index", index, "--gene", gene_id, "--format", "json")
  if (source != "all")        args <- c(args, "--source", source)
  if (!is.null(min_identity)) args <- c(args, "--min_identity", format(min_identity, scientific = FALSE))
  if (!is.null(min_bitscore)) args <- c(args, "--min_bitscore", format(min_bitscore, scientific = FALSE))
  if (!is.null(limit))        args <- c(args, "--limit", as.integer(limit))

  raw <- .run_engine("homolog_index_query", args, engines_dir = engines_dir)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' required for atlas_lookup; install.packages('jsonlite')")
  }
  payload <- jsonlite::fromJSON(paste(raw, collapse = "\n"), simplifyVector = TRUE)
  if (!isTRUE(payload$found)) return(.empty_hits())
  hits <- payload$hits
  if (!is.data.frame(hits) || nrow(hits) == 0) return(.empty_hits())
  hits$query <- payload$query
  hits <- hits[, c("query", "target", "source_label", "source_type",
                   "identity", "bitscore", "evalue",
                   "qstart", "qend", "tstart", "tend",
                   "alnlen", "n_segments", "strand")]
  .as_tibble_safe(hits)
}

#' Cross-organism hits only (miniprot source).
#' @inheritParams atlas_lookup
#' @export
atlas_orthologs <- function(index, gene_id, ...) {
  atlas_lookup(index, gene_id, source = "miniprot", ...)
}

#' Within-organism hits only (diamond source), self-hit removed.
#' @inheritParams atlas_lookup
#' @export
atlas_paralogs <- function(index, gene_id, ...) {
  hits <- atlas_lookup(index, gene_id, source = "diamond", ...)
  hits[hits$target != hits$query, , drop = FALSE]
}

# ── codon_stats wrapper ─────────────────────────────────────────────────────

.empty_codon_stats <- function() {
  num <- numeric(0); int <- integer(0)
  .as_tibble_safe(data.frame(
    pair_id = character(0),
    n_aligned_bp = int, n_codons_total = int, n_codons_used = int,
    n_codons_gap = int, n_codons_stop = int,
    n_4d_sites = int, n_4d_diffs = int, n_4d_ts = int, n_4d_tv = int,
    S = num, N = num, sd = num, nd = num,
    sd_ts = num, sd_tv = num, nd_ts = num, nd_tv = num,
    pS = num, pN = num, dS = num, dN = num,
    omega = num, p4d = num, d4d = num,
    stringsAsFactors = FALSE
  ))
}

#' Run codon_stats on a TSV of aligned CDS pairs.
#'
#' @param pairs path to TSV (pair_id<TAB>seqA<TAB>seqB), or NULL
#' @param fasta path to multi-FASTA with consecutive A/B records, or NULL
#' @param method "ng86" (Nei-Gojobori, JC) or "yn00" (Yang-Nielsen κ-weighted, K80)
#' @param kappa optional numeric for YN00; if NULL, estimated from data
#' @param ncores OpenMP threads across pairs
#' @param engines_dir override engines/ location
#' @return tibble with one row per pair (columns in
#'   engines/schemas/codon_stats.output.schema.json)
#' @export
codon_stats_run <- function(pairs = NULL, fasta = NULL,
                             method = c("ng86", "yn00"),
                             kappa = NULL, ncores = 1,
                             engines_dir = NULL) {
  method <- match.arg(method)
  if (is.null(pairs) && is.null(fasta))
    stop("Need pairs (TSV) or fasta (multi-FASTA).")
  if (!is.null(pairs) && !file.exists(pairs)) stop("pairs file not found: ", pairs)
  if (!is.null(fasta) && !file.exists(fasta)) stop("fasta file not found: ", fasta)

  args <- character(0)
  if (!is.null(pairs)) args <- c(args, "--pairs", pairs)
  if (!is.null(fasta)) args <- c(args, "--fasta", fasta)
  args <- c(args, "--method", method, "--ncores", as.character(ncores))
  if (!is.null(kappa)) args <- c(args, "--kappa", format(kappa, scientific = FALSE))

  raw <- .run_engine("codon_stats", args, engines_dir = engines_dir)
  comment_idx <- grep("^#", raw)
  rows <- raw[setdiff(seq_along(raw), comment_idx)]
  if (length(rows) <= 1) return(.empty_codon_stats())

  tmp <- tempfile(fileext = ".tsv")
  writeLines(rows, tmp)
  df <- utils::read.table(tmp, header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("NA", "Inf", "-Inf"),
                           check.names = FALSE)
  unlink(tmp)
  .as_tibble_safe(df)
}

# ── Atlas build orchestration ───────────────────────────────────────────────

#' Build a homolog atlas index via the shell dispatcher.
#'
#' Wraps dispatchers/build_homolog_atlas.sh.
#'
#' @param query path to query proteome (FASTA)
#' @param out_dir output directory (created if missing)
#' @param targets list of list(path=, label=) entries (miniprot)
#' @param diamond_files list of list(path=, label=) entries (pre-computed DIAMOND TSV)
#' @param diamond_subjects list of list(path=, label=) entries (run DIAMOND blastp on the fly)
#' @param mode "strict" / "chain" (default) / "loose"
#' @param threads integer
#' @param keep_intermediate keep miniprot/diamond intermediate files
#' @param config optional path to JSON config (overrides matching args)
#' @param dry_run only echo the planned commands
#' @return invisible(path to homolog_atlas.holindx)
#' @export
atlas_build <- function(query, out_dir,
                         targets = list(), diamond_files = list(),
                         diamond_subjects = list(),
                         mode = c("chain", "strict", "loose"),
                         threads = 8,
                         keep_intermediate = FALSE,
                         config = NULL,
                         dry_run = FALSE,
                         engines_dir = NULL) {
  mode <- match.arg(mode)
  script <- normalizePath(file.path(atlas_engines_dir(engines_dir), "..", "dispatchers",
                                     "build_homolog_atlas.sh"), mustWork = FALSE)
  if (!file.exists(script))
    stop("build_homolog_atlas.sh not found at ", script)

  args <- c("--query", query, "--out_dir", out_dir,
            "--mode", mode, "--threads", as.character(threads))
  for (t in targets)          args <- c(args, "--target",             paste(t$path, t$label, sep = "="))
  for (d in diamond_files)    args <- c(args, "--diamond",            paste(d$path, d$label, sep = "="))
  for (s in diamond_subjects) args <- c(args, "--run_diamond_blastp", paste(s$path, s$label, sep = "="))
  if (keep_intermediate) args <- c(args, "--keep_intermediate")
  if (!is.null(config))  args <- c(args, "--config", config)
  if (dry_run)           args <- c(args, "--dry_run")

  status <- system2(script, args = args)
  if (status != 0) stop("atlas_build failed with status ", status)
  index_path <- file.path(out_dir, "homolog_atlas.holindx")
  if (dry_run) return(invisible(index_path))
  if (!file.exists(index_path)) stop("Atlas index not produced at ", index_path)
  invisible(normalizePath(index_path))
}

# ── Health check ────────────────────────────────────────────────────────────

#' Quick sanity report on a built index.
#'
#' Verifies the binary is built, the file exists, and a no-such-gene query
#' returns the expected JSON. Useful as an atlas-page bootstrap check.
#'
#' @param index path to .holindx
#' @return list of diagnostics
#' @export
atlas_health_check <- function(index, engines_dir = NULL) {
  ed <- atlas_engines_dir(engines_dir)
  bin <- file.path(ed, "homolog_index_query")
  if (!file.exists(bin)) stop("homolog_index_query not built at ", bin)
  if (!file.exists(index)) stop("Index not found: ", index)
  size_mb <- file.size(index) / 1024^2
  raw <- suppressWarnings(system2(bin, c("--index", index,
                                           "--gene", "__atlas_health_check__",
                                           "--format", "json"),
                                    stdout = TRUE, stderr = TRUE))
  list(
    index_path        = normalizePath(index),
    engines_dir       = ed,
    binary_path       = bin,
    size_mb           = round(size_mb, 3),
    binary_responds   = length(raw) > 0,
    last_response     = if (length(raw) > 0) raw[1] else NA_character_
  )
}
