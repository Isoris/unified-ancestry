#!/usr/bin/env Rscript
# =============================================================================
# region_stats_dispatcher.R — Unified Stats Dispatcher (v12.1 REWIRED)
#
# v12.1 changes:
#   1. resolve_groups() can pull from sample_registry via reg$get_group()
#   2. All dosage matrices have CGA colnames (via smap/instant_q)
#   3. configure_dispatcher() loads registry if available
#   4. group names like "inv_LG05_3e6_8e6_HET" resolve from registry
#
# Central interface for ALL population genetics statistics from any region.
# Routes to the correct backend:
#   - Q (Group 1): Engine B C++ via instant_q.R
#   - Dosage-derived popstats (Groups 2-3): C binary region_popstats
#   - Hobs/HWE (Group 4): ANGSD -doHWE via hobs_windower
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
EPS <- 1e-10

# =============================================================================
# CONFIGURATION
# =============================================================================

.disp_env <- new.env(parent = emptyenv())

configure_dispatcher <- function(config_file = NULL, ...) {
  # Load instant_q config first
  iq_path <- file.path(dirname(sys.frame(1)$ofile %||% "."),
                        "wrappers", "instant_q.R")
  if (file.exists(iq_path)) {
    source(iq_path, local = FALSE)
    configure_instant_q(config_file = config_file, ...)
  }

  if (!is.null(config_file) && file.exists(config_file)) {
    cfg <- parse_shell_config(config_file)
    .disp_env$beagle_dir    <- cfg$BEAGLE_DIR %||% cfg$STEP1_BEAGLE_DIR
    .disp_env$sample_list   <- cfg$SAMPLE_LIST %||% cfg$STEP1_SAMPLE_LIST
    .disp_env$pruned_list   <- cfg$STEP1_PRUNED_LIST
    .disp_env$cache_dir     <- cfg$STATS_CACHE_DIR %||% cfg$LOCAL_Q_DIR %||% "region_stats_cache"
    .disp_env$snake1_dir    <- cfg$SNAKE1_DIR
    .disp_env$decomp_dir    <- cfg$DECOMPOSITION_DIR
    .disp_env$invdir        <- cfg$INVDIR
    .disp_env$best_qopt     <- cfg$BEST_QOPT
    .disp_env$best_fopt     <- cfg$BEST_FOPT
    .disp_env$K             <- as.integer(cfg$DEFAULT_K %||% "8")
    .disp_env$config_file   <- config_file

    # C binary paths
    .disp_env$popstats_bin  <- cfg$POPSTATS_BIN %||%
      file.path(dirname(config_file), "engines", "region_popstats")
    .disp_env$hobs_dir      <- cfg$HOBS_OUTDIR %||%
      file.path(cfg$BASE %||% ".", "hobs_hwe_confirmation")

    # v12.1: registry
    .disp_env$registry_dir  <- cfg$REGISTRY_DIR %||%
      file.path(cfg$BASE %||% ".", "sample_registry")
  }

  # v12.1: Try to load registry
  # Chat-18: prefer the full 4-table loader (registries/api/R/registry_loader.R)
  # over the old flat sample-only one. Either works via top-level aliases
  # (reg$has_group / reg$get_group / reg$list_groups), but the full one also
  # exposes reg$samples / reg$intervals / reg$evidence / reg$results for any
  # downstream caller that needs them. If a script sourced load_bridge.R
  # first, reg is already in the global env and branch 1 wins — we only
  # reach branch 2 in standalone invocations.
  .disp_env$reg <- NULL
  if (exists("reg", envir = .GlobalEnv) && !is.null(get("reg", envir = .GlobalEnv))) {
    .disp_env$reg <- get("reg", envir = .GlobalEnv)
    message("[dispatcher] Using registry from global env")
  } else {
    # Full-reg only (2026-04-24 migration). Flat sample_registry.R fallback
    # removed; v10 full loader is mandatory. Search order: explicit env
    # override first, then canonical v10 paths.
    reg_r <- Sys.getenv("SAMPLE_REGISTRY_R", "")
    loader_flavor <- "none"

    if (nzchar(reg_r) && file.exists(reg_r)) {
      loader_flavor <- "override"
    } else {
      base <- Sys.getenv("BASE", "")
      full_cands <- c(
        file.path(.disp_env$registry_dir %||% ".", "..", "registries",
                  "api", "R", "registry_loader.R"),
        file.path(base, "registries", "api", "R", "registry_loader.R"),
        file.path(base, "inversion-popgen-toolkit", "registries", "api", "R",
                  "registry_loader.R"),
        "registries/api/R/registry_loader.R"
      )
      for (p in full_cands) {
        if (file.exists(p)) { reg_r <- p; loader_flavor <- "full"; break }
      }
    }

    if (nzchar(reg_r) && file.exists(reg_r)) {
      tryCatch({
        source(reg_r, local = TRUE)
        # Full loader auto-resolves registries_root via $REGISTRIES / $BASE
        .disp_env$reg <- load_registry()
        n_groups <- tryCatch(nrow(.disp_env$reg$list_groups()),
                             error = function(e) NA_integer_)
        message("[dispatcher] Registry loaded (", loader_flavor,
                "): ", n_groups, " groups")
      }, error = function(e) {
        message("[dispatcher] Registry load failed: ", conditionMessage(e))
      })
    }
  }

  dir.create(.disp_env$cache_dir %||% "region_stats_cache",
             recursive = TRUE, showWarnings = FALSE)
  message("[dispatcher] Configured (v12.1: registry-aware)")
  invisible(TRUE)
}

# =============================================================================
# GROUP AUTO-RESOLUTION (v12.1 — registry-first)
# =============================================================================

resolve_groups <- function(groups_from, chr = NULL) {
  if (is.list(groups_from)) return(groups_from)

  # ── v12.1: Try registry first ──
  # If groups_from looks like a registry group ID (or comma-separated list)
  reg <- .disp_env$reg
  if (!is.null(reg) && is.character(groups_from) && length(groups_from) == 1) {

    # Check if it's a comma-separated list of group IDs
    candidate_ids <- trimws(strsplit(groups_from, ",")[[1]])

    # Attempt registry lookup
    registry_groups <- list()
    for (gid in candidate_ids) {
      if (reg$has_group(gid)) {
        members <- reg$get_group(gid)
        if (length(members) > 0) {
          registry_groups[[gid]] <- members
        }
      }
    }

    if (length(registry_groups) >= 1) {
      message("[dispatcher] Resolved ", length(registry_groups),
              " groups from registry: ", paste(names(registry_groups), collapse = ", "))
      return(registry_groups)
    }

    # Try pattern match: "inv_*" → find all groups matching pattern
    if (grepl("^inv_|^ghsl_|^ancestry_", groups_from)) {
      all_groups <- reg$list_groups()
      if (nrow(all_groups) > 0) {
        matching <- all_groups[grepl(groups_from, group_id, fixed = FALSE)]
        if (nrow(matching) > 0) {
          grps <- list()
          for (i in seq_len(nrow(matching))) {
            gid <- matching$group_id[i]
            members <- reg$get_group(gid)
            if (length(members) > 0) grps[[gid]] <- members
          }
          if (length(grps) >= 1) {
            message("[dispatcher] Pattern-matched ", length(grps), " groups from registry")
            return(grps)
          }
        }
      }
    }
  }

  # ── Fallback: named resolution strategies ──
  if (groups_from == "precomp_bands") {
    band_file <- file.path(.disp_env$snake1_dir %||% "",
                            "precomp", paste0(chr, ".precomp.rds"))
    if (!is.null(chr) && file.exists(band_file)) {
      pc <- readRDS(band_file)
      if ("band_assignments" %in% names(pc)) {
        ba <- pc$band_assignments
        groups <- split(names(ba), ba)
        names(groups) <- paste0("band", seq_along(groups))
        return(groups)
      }
    }
    message("[dispatcher] Cannot resolve precomp_bands for ", chr)
    return(NULL)

  } else if (groups_from == "C01i_decomposition") {
    decomp_dir <- .disp_env$decomp_dir %||% ""
    if (nzchar(decomp_dir)) {
      decomp_files <- list.files(decomp_dir, pattern = "decomposition.*\\.tsv",
                                  full.names = TRUE)
      if (!is.null(chr)) decomp_files <- grep(chr, decomp_files, value = TRUE)
      if (length(decomp_files) > 0) {
        dt <- fread(decomp_files[1])
        class_col <- intersect(names(dt), c("genotype", "class", "inv_class",
                                              "assignment"))[1]
        if (!is.null(class_col) && "sample_id" %in% names(dt)) {
          return(split(dt$sample_id, dt[[class_col]]))
        }
      }
    }

    # v12.1: try registry fallback for decomposition groups
    if (!is.null(reg)) {
      all_groups <- reg$list_groups()
      if (nrow(all_groups) > 0) {
        inv_groups <- all_groups[grepl("^inv_.*_(HOM_REF|HET|HOM_INV)$", group_id)]
        if (!is.null(chr)) inv_groups <- inv_groups[grepl(chr, inv_id) | grepl(chr, group_id)]
        if (nrow(inv_groups) > 0) {
          grps <- list()
          for (i in seq_len(nrow(inv_groups))) {
            gid <- inv_groups$group_id[i]
            members <- reg$get_group(gid)
            if (length(members) > 0) grps[[gid]] <- members
          }
          if (length(grps) >= 1) return(grps)
        }
      }
    }

    message("[dispatcher] Cannot resolve C01i_decomposition")
    return(NULL)

  } else if (groups_from == "ancestry_clusters") {
    # v12.1: try registry first
    if (!is.null(reg)) {
      K <- .disp_env$K %||% 8
      grps <- list()
      for (qi in seq_len(K)) {
        gid <- paste0("ancestry_K", K, "_Q", qi)
        if (reg$has_group(gid)) {
          members <- reg$get_group(gid)
          if (length(members) > 0) grps[[paste0("Q", qi)]] <- members
        }
      }
      if (length(grps) >= 2) {
        message("[dispatcher] Ancestry clusters from registry: ", length(grps), " groups")
        return(grps)
      }
    }

    # Fallback to Q matrix
    qopt <- .disp_env$best_qopt
    slist <- .disp_env$sample_list
    if (!is.null(qopt) && file.exists(qopt) && !is.null(slist) && file.exists(slist)) {
      q <- as.matrix(fread(qopt, header = FALSE))
      samples <- readLines(slist)
      samples <- samples[nzchar(samples)]
      assigned <- apply(q, 1, which.max)
      n_use <- min(length(samples), nrow(q))
      groups <- split(samples[1:n_use], assigned[1:n_use])
      names(groups) <- paste0("Q", names(groups))
      return(groups)
    }
    return(NULL)
  }

  # Try as file path
  if (file.exists(groups_from)) {
    dt <- fread(groups_from)
    grp_col <- intersect(names(dt), c("group", "class", "genotype", "band"))[1]
    id_col <- intersect(names(dt), c("sample_id", "sample", "id"))[1]
    if (!is.null(grp_col) && !is.null(id_col)) {
      return(split(dt[[id_col]], dt[[grp_col]]))
    }
  }

  message("[dispatcher] Unknown groups_from: ", groups_from)
  NULL
}

# =============================================================================
# C BINARY DELEGATION: region_popstats
# =============================================================================

route_to_c_popstats <- function(chr, start = NULL, end = NULL,
                                 groups = NULL,
                                 fixed_win = "50000:10000",
                                 downsample = 1, type = 2) {
  bin <- .disp_env$popstats_bin
  if (is.null(bin) || !file.exists(bin)) {
    message("[dispatcher] region_popstats binary not found: ", bin %||% "NULL")
    return(NULL)
  }

  bgl <- list.files(.disp_env$beagle_dir, pattern = paste0(".*", chr, ".*beagle"),
                     full.names = TRUE)[1]
  if (is.na(bgl) || !file.exists(bgl)) {
    message("[dispatcher] No BEAGLE for ", chr)
    return(NULL)
  }

  tmpdir <- tempdir()
  out_file <- file.path(tmpdir, paste0(chr, ".popstats_tmp.tsv"))

  groups_arg <- ""
  tmp_files <- character(0)
  if (!is.null(groups) && length(groups) > 0) {
    specs <- character(0)
    for (gname in names(groups)) {
      gfile <- file.path(tmpdir, paste0("grp_", gname, ".txt"))
      writeLines(groups[[gname]], gfile)
      tmp_files <- c(tmp_files, gfile)
      specs <- c(specs, paste0(gname, ":", gfile))
    }
    groups_arg <- paste0("--groups ", paste(specs, collapse = ","))
  }

  range_arg <- ""
  if (!is.null(start) && !is.null(end)) {
    range_arg <- paste0("--range ", as.integer(start), ":", as.integer(end))
  }

  cmd <- paste(
    shQuote(bin),
    "--beagle", shQuote(bgl),
    "--sample_list", shQuote(.disp_env$sample_list),
    "--chr", shQuote(chr),
    "--fixed_win", fixed_win,
    "--downsample", downsample,
    "--type", type,
    "--out", shQuote(out_file),
    groups_arg,
    range_arg
  )

  ret <- system(cmd, intern = FALSE)
  on.exit(file.remove(c(out_file, tmp_files)[file.exists(c(out_file, tmp_files))]),
          add = TRUE)

  if (ret != 0 || !file.exists(out_file)) {
    message("[dispatcher] region_popstats failed for ", chr)
    return(NULL)
  }

  dt <- fread(out_file, skip = 1)
  dt
}

# =============================================================================
# HOBS/HWE DATA LOADER
# =============================================================================

load_hobs_windows <- function(chr, subset_id, scale = "50kb") {
  hobs_dir <- .disp_env$hobs_dir
  path <- file.path(hobs_dir, "03_hobs_windows", subset_id,
                     paste0(chr, ".win", scale, ".tsv"))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

# =============================================================================
# MAIN DISPATCHER: get_region_stats()
# =============================================================================

get_region_stats <- function(chr, start, end,
                              what = "Q",
                              groups = NULL,
                              groups_from = NULL,
                              popstats_opts = list()) {
  result <- list()

  need_groups <- any(c("Fst", "dXY", "dA", "MI", "family_likeness") %in% what)
  if (need_groups && is.null(groups) && !is.null(groups_from)) {
    groups <- resolve_groups(groups_from, chr)
  }

  # ── Q (Engine B) ──
  if ("Q" %in% what) {
    if (exists("get_Q", mode = "function")) {
      result$Q <- get_Q(chr, start, end)
    }
  }

  # ── C-binary popstats ──
  c_stats <- c("Fst", "dXY", "dA", "MI", "theta_pi", "theta_W", "Tajima_D", "Hp")
  need_c <- any(c_stats %in% what)
  if (need_c) {
    fw <- popstats_opts$fixed_win %||% "50000:10000"
    ds <- popstats_opts$downsample %||% 1
    tp <- popstats_opts$type %||% 2

    c_dt <- route_to_c_popstats(chr, start, end, groups = groups,
                                 fixed_win = fw, downsample = ds, type = tp)
    if (!is.null(c_dt) && nrow(c_dt) > 0) {
      result$popstats <- c_dt

      if ("theta_pi" %in% what && "theta_pi_all" %in% names(c_dt))
        result$theta_pi <- round(mean(c_dt$theta_pi_all, na.rm = TRUE), 6)
      if ("theta_W" %in% what && "theta_W_all" %in% names(c_dt))
        result$theta_W <- round(mean(c_dt$theta_W_all, na.rm = TRUE), 6)
      if ("Tajima_D" %in% what && "Tajima_D" %in% names(c_dt))
        result$Tajima_D <- round(mean(c_dt$Tajima_D, na.rm = TRUE), 4)

      fst_cols <- grep("^Fst_", names(c_dt), value = TRUE)
      dxy_cols <- grep("^dXY_", names(c_dt), value = TRUE)
      if ("Fst" %in% what && length(fst_cols) > 0)
        result$Fst <- lapply(setNames(fst_cols, fst_cols), function(col)
          round(mean(c_dt[[col]], na.rm = TRUE), 6))
      if ("dXY" %in% what && length(dxy_cols) > 0)
        result$dXY <- lapply(setNames(dxy_cols, dxy_cols), function(col)
          round(mean(c_dt[[col]], na.rm = TRUE), 6))
    }
  }

  # ── Hobs/HWE ──
  if ("Hobs" %in% what || "HWE" %in% what) {
    hobs_results <- list()
    subsets_dir <- file.path(.disp_env$hobs_dir %||% "", "03_hobs_windows")
    if (dir.exists(subsets_dir)) {
      subset_dirs <- list.dirs(subsets_dir, recursive = FALSE, full.names = FALSE)
      for (sid in subset_dirs) {
        hw <- load_hobs_windows(chr, sid, scale = "50kb")
        if (!is.null(hw) && nrow(hw) > 0) {
          in_region <- hw[window_center >= start & window_center <= end]
          if (nrow(in_region) > 0) {
            hobs_results[[sid]] <- list(
              mean_Hobs = round(mean(in_region$mean_Hobs, na.rm = TRUE), 6),
              n_windows = nrow(in_region)
            )
          }
        }
      }
    }
    if (length(hobs_results) > 0) result$Hobs <- hobs_results
  }

  result
}

# =============================================================================
# BATCH MODE
# =============================================================================

batch_candidate_stats <- function(candidates_file, what = "ALL",
                                   groups = NULL, groups_from = NULL,
                                   outdir = NULL, popstats_opts = list()) {
  outdir <- outdir %||% file.path(.disp_env$cache_dir, "candidate_stats")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  cands <- fread(candidates_file)
  chr_col <- intersect(names(cands), c("chrom", "chr"))[1]
  start_col <- intersect(names(cands), c("start_bp", "start"))[1]
  end_col <- intersect(names(cands), c("end_bp", "end"))[1]
  id_col <- intersect(names(cands), c("candidate_id", "interval_id"))[1]

  if (is.na(chr_col) || is.na(start_col) || is.na(end_col))
    stop("[dispatcher] Candidate file needs chrom, start, end columns")

  if (identical(what, "ALL"))
    what <- c("Q", "Fst", "dXY", "theta_pi", "theta_W", "Tajima_D", "Hobs")

  message("[dispatcher] Batch: ", nrow(cands), " candidates")

  all_rows <- list()
  for (ci in seq_len(nrow(cands))) {
    chr <- cands[[chr_col]][ci]
    s <- cands[[start_col]][ci]
    e <- cands[[end_col]][ci]
    cid <- if (!is.na(id_col)) cands[[id_col]][ci] else paste0(chr, "_", s, "_", e)

    grps <- groups
    if (is.null(grps) && !is.null(groups_from))
      grps <- resolve_groups(groups_from, chr)

    stats <- get_region_stats(chr, s, e, what = what, groups = grps,
                               popstats_opts = popstats_opts)

    flat <- list(candidate_id = cid, chrom = chr, start_bp = s, end_bp = e)
    for (nm in names(stats)) {
      val <- stats[[nm]]
      if (is.data.table(val)) next
      if (is.list(val) && !is.data.table(val)) {
        for (sub_nm in names(val)) {
          sv <- val[[sub_nm]]
          if (is.list(sv)) {
            for (ssn in names(sv)) flat[[paste0(sub_nm, "_", ssn)]] <- sv[[ssn]]
          } else if (length(sv) == 1) flat[[sub_nm]] <- sv
        }
      } else if (length(val) == 1) flat[[nm]] <- val
    }
    all_rows[[ci]] <- flat
    if (ci %% 10 == 0) message("[dispatcher] ", ci, "/", nrow(cands))
  }

  dt <- rbindlist(all_rows, fill = TRUE)
  out_file <- file.path(outdir, "candidate_stats.tsv")
  fwrite(dt, out_file, sep = "\t")
  message("[dispatcher] Saved: ", out_file)
  dt
}
