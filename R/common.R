## ============================================================================
## common.R -- shared definitions for all R steps of the pipeline
## ============================================================================
## Sourced by every script in R/. Provides paths, the panel definition, marker
## roles, the sample sheet, the positivity threshold estimator, the statistical
## helpers and the autofluorescence rules.
##
## Configuration (environment variables):
##   SFP_PANEL    panel definition: "ntnk" (default), "mye" or the path to an
##                R file that defines PANEL (see config/panel_ntnk.R)
##   SFP_SAMPLES  sample sheet (CSV, see docs/SAMPLE_SHEET.md)
##   SFP_DESIGN   study design for the statistics (R file, see
##                config/design_template.R)
## ============================================================================

suppressPackageStartupMessages({ library(data.table) })

`%||%` <- function(a, b) if (is.null(a) || !length(a) || (length(a) == 1 && is.na(a))) b else a

## Locate the repository root from the calling script ("--file=" argument of Rscript).
.script_dir <- function() {
  ## Rscript encodes spaces in the script path as "~+~"
  f <- gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
  if (is.na(f)) "R" else dirname(normalizePath(f))
}
REPO <- normalizePath(file.path(.script_dir(), ".."), mustWork = FALSE)
if (!dir.exists(file.path(REPO, "config"))) REPO <- normalizePath(file.path(.script_dir(), "..", ".."), mustWork = FALSE)

dir_ok <- function(p) { if (!dir.exists(p)) dir.create(p, recursive = TRUE); p }

## ---- Panel -------------------------------------------------------------------
## A panel is a table with the columns channel, marker and role:
##   "lineage"  lineage markers; they span the clustering space
##   "state"    activation and function markers; kept out of the clustering and
##              evaluated per population (differential state)
##   "qc"       quality control only (viability dye)
## DISCARDED lists markers without signal above noise (see 05_thresholds.R);
## they are neither clustered nor tested.
.panel_arg <- Sys.getenv("SFP_PANEL", "ntnk")
.panel_file <- if (.panel_arg %in% c("ntnk", "mye")) file.path(REPO, "config", paste0("panel_", .panel_arg, ".R")) else .panel_arg
stopifnot(file.exists(.panel_file))
source(.panel_file, local = FALSE)          # defines PANEL, PANEL_ID and DISCARDED
stopifnot(all(c("channel", "marker", "role") %in% names(PANEL)))
IS_MYE <- identical(PANEL_ID, "mye")

LINEAGE_MARKERS <- setdiff(PANEL$marker[PANEL$role == "lineage"], DISCARDED)
STATE_MARKERS   <- setdiff(PANEL$marker[PANEL$role == "state"], DISCARDED)
ALL_MARKERS     <- PANEL$marker[PANEL$role %in% c("lineage", "state")]   # including discarded markers
SCATTER         <- c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W")

## Marker names from a flowFrame: the $PnS description without the
## ": dye - Area" suffix, mapped to the canonical names of the panel.
marker_names <- function(ff) {
  p <- flowCore::pData(flowCore::parameters(ff))
  m <- trimws(sub(" *:.*$", "", ifelse(is.na(p$desc), p$name, p$desc)))
  keep <- is.na(p$desc) | p$name %in% c(SCATTER, "TIME", "Time")
  m[keep] <- p$name[keep]
  hit <- match(p$name, PANEL$channel)
  m[!is.na(hit)] <- PANEL$marker[hit[!is.na(hit)]]
  make.unique(m)
}

## ---- Excluded clusters -----------------------------------------------------------
## Cluster numbers removed as debris after the preliminary clustering
## (see 09_fragment_cluster.R). Set per run with SFP_EXCLUDE (comma-separated).
EXCLUDED_CLUSTERS <- trimws(strsplit(Sys.getenv("SFP_EXCLUDE", ""), ",")[[1]])
EXCLUDED_CLUSTERS <- EXCLUDED_CLUSTERS[nzchar(EXCLUDED_CLUSTERS)]

## ---- Sample sheet ------------------------------------------------------------
## One row per FCS file. Required columns: file (base name of the FCS file),
## sample_id (text; leading zeros are kept) and run (measurement run, the batch
## variable). Any further columns (grouping variables, covariates) are named in
## the design file and are never interpreted by this module.
read_samples <- function(path = Sys.getenv("SFP_SAMPLES")) {
  if (!nzchar(path) || !file.exists(path)) stop("sample sheet not found; set SFP_SAMPLES")
  s <- fread(path, colClasses = "character", na.strings = c("", "NA"))
  miss <- setdiff(c("file", "sample_id", "run"), names(s))
  if (length(miss)) stop("sample sheet lacks column(s): ", paste(miss, collapse = ", "))
  if (anyDuplicated(s$file)) stop("duplicated file names in the sample sheet")
  s
}

## Attach sample-sheet columns to a table with a column `file`. Unknown files
## stop the run instead of silently producing missing values.
attach_samples <- function(dt, cols, samples = read_samples()) {
  miss <- setdiff(unique(dt$file), samples$file)
  if (length(miss)) stop(length(miss), " file(s) missing from the sample sheet, e.g. ", miss[1])
  for (cc in cols) set(dt, j = cc, value = samples[[cc]][match(dt$file, samples$file)])
  invisible(dt)
}

## ---- Study design -----------------------------------------------------------------
## Loaded only by the statistics scripts. See config/design_template.R.
load_design <- function(path = Sys.getenv("SFP_DESIGN")) {
  if (!nzchar(path) || !file.exists(path)) stop("design file not found; set SFP_DESIGN")
  e <- new.env()
  sys.source(file.path(REPO, "config", "design_helpers.R"), envir = e)
  sys.source(path, envir = e)
  d <- list(comparisons = e$COMPARISONS,
            protected_covariates = e$PROTECTED_COVARIATES %||% character(0),
            correlate_with = e$CORRELATE_WITH %||% character(0),
            one_sided = e$ONE_SIDED %||% list(markers = character(0), comparisons = character(0), alternative = "less"),
            extra_adjustment = e$EXTRA_ADJUSTMENT %||% character(0),
            factorial = e$FACTORIAL)
  vars <- unique(c(unlist(lapply(d$comparisons, function(v) c(v$var, names(v$str)))),
                   d$factorial$outer$var, d$factorial$inner$var))
  d$group_columns <- vars
  if (length(d$extra_adjustment) && is.null(names(d$extra_adjustment))) stop("EXTRA_ADJUSTMENT must be a named vector")
  d
}

## ---- arcsinh -------------------------------------------------------------------------
arcsinh_trans <- function(x, cofactor = 6000) asinh(x / cofactor)

## ---- Positivity threshold -----------------------------------------------------------
## The threshold is set by the distance from the negative mode, not by
## bimodality. A two-component mixture does not find rare positive populations;
## it splits the dominant negative peak instead. For broad continua (for
## example CD127) mixture models from different implementations also disagree,
## so the threshold would depend on the software rather than on the data.
##
## Procedure: the mode of a kernel density estimate (restricted to the 0.5th to
## 99.5th percentile), a robust spread from the left half only (the positive
## cells lie to the right of the mode), threshold = mode + k * spread (k = 3).
## Quality measure "reach": distance of the 99.5th percentile from the mode, in
## units of the spread. A marker without signal reaches a few units, one with a
## clear positive population reaches double digits.
## All values are used (no subsampling): an unseeded subsample made the
## threshold differ between repeated runs.
threshold_negative_mode <- function(x, n = Inf, k = 3) {
  x <- x[is.finite(x)]
  if (!length(x)) return(list(threshold = NA_real_, neg_mode = NA_real_,
                              neg_sd = NA_real_, reach = NA_real_, components = NA_integer_))
  xs <- if (is.finite(n) && length(x) > n) sample(x, n) else x
  gr <- range(quantile(xs, c(.005, .995), na.rm = TRUE))
  d  <- tryCatch(density(xs, from = gr[1], to = gr[2], n = 512, adjust = 1),
                 error = function(e) NULL)
  mode <- if (is.null(d)) median(xs) else d$x[which.max(d$y)]
  left <- x[x <= mode]
  ns <- if (length(left) < 100) NA_real_ else 1.4826 * median(abs(left - mode))
  if (!is.finite(ns) || ns <= 0) ns <- mad(x, na.rm = TRUE)
  list(threshold = mode + k * ns, neg_mode = mode, neg_sd = ns,
       reach = (unname(quantile(x, .995, na.rm = TRUE)) - mode) / max(ns, 1e-9),
       components = NA_integer_)
}

rate_marker <- function(reach, frac_pos) {
  ifelse(is.na(reach), "not assessable",
  ifelse(reach < 3 & frac_pos < 0.01,
         "no negative population or no signal -- do not use as +/-",
  ifelse(reach < 3,  "no signal above noise",
  ifelse(reach < 6,  "weak signal",
  ifelse(frac_pos < 0.002, "very rare -- qualitative only", "clear signal")))))
}

## ---- Statistics ------------------------------------------------------------------------
## Minimum number of samples per group for a test, minimum number of cells of
## the reference population (lymphocytes or leukocytes) for a sample to be
## included.
MIN_N <- 5
MIN_L <- 100

## Two-sided Mann-Whitney U test with R's default: exact p values when both
## groups have fewer than 50 values. From R 4.6.0 on the exact distribution is
## also used in the presence of ties (shift algorithm of Streitberg & Roehmel,
## EDV in Medizin und Biologie 1987;18:12-19; see the R 4.6.0 help of wilcox.test);
## older R versions fall back to the normal approximation with continuity
## correction when ties are present. The pipeline therefore requires R >= 4.6.0.
mw <- function(a, b, alternative = "two.sided") {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) < MIN_N || length(b) < MIN_N)
    return(list(p = NA_real_, n1 = length(a), n2 = length(b), tested = FALSE))
  t <- suppressWarnings(wilcox.test(a, b, alternative = alternative))
  list(p = unname(t$p.value), n1 = length(a), n2 = length(b), tested = TRUE)
}
if (getRversion() < "4.6.0") warning("R < 4.6.0: Mann-Whitney p values with ties are not exact; results will differ")

## Cliff's delta = P(b > a) - P(b < a) with a percentile bootstrap interval.
cliffs_delta <- function(a, b) mean(outer(b, a, ">")) - mean(outer(b, a, "<"))
cliffs_delta_ci <- function(a, b, B = 2000L, seed = 1L) {
  set.seed(seed)
  d <- replicate(B, cliffs_delta(sample(a, replace = TRUE), sample(b, replace = TRUE)))
  c(delta = cliffs_delta(a, b), lower = unname(quantile(d, 0.025)), upper = unname(quantile(d, 0.975)))
}

## Rank residuals after one or more covariates (adjusted comparison: the
## Mann-Whitney test is applied to the residuals of the ranks regressed on the
## covariates). na.exclude keeps the full length with NA where a covariate is
## missing; without it, residuals() is shorter and assignment would recycle
## values onto the wrong samples. Residuals are rounded to 9 decimals: when a
## run contributes a single sample to the tested subset, its residual is exactly
## zero in theory but carries floating-point noise below 1e-9 in practice;
## rounding restores the tie so that R and Python rank it identically.
resid_rank <- function(value, ...) {
  f <- list(...)
  tryCatch({
    d <- data.frame(r = rank(value, na.last = "keep"), f)
    round(as.numeric(residuals(lm(r ~ ., data = d, na.action = na.exclude))), 9)
  }, error = function(e) rep(NA_real_, length(value)))
}

## ---- Autofluorescence load ------------------------------------------------------------
## NK/T panel: a small fraction of cells carries a far-red autofluorescence that
## the unmixing distributes onto the dyes of the 665-710 nm window and that
## creates false positivity there: FOXP3 (BB700), CD25 (SBUV665), CD223 (RR688),
## CD49a (PerCP-eF710), CXCR5 (RB705). Evidence from the raw detector data of
## the same stained cells: the excess at 670-710 nm appears under UV, violet,
## blue and red excitation at the same time, which no single antibody dye does
## (see R/unmixing/af_signature.R). Affected cells stay in their lineage and in
## every count; only these channels are masked in the expression statistics.
AF_MARKERS <- c("FOXP3", "CD25", "CD223", "CD49a", "CXCR5")
## min_extra sets the strictness of the rule; 2 is the default, 1 and 3 are
## used for the sensitivity analysis.
mask_af <- function(Zm, af) {
  sp <- intersect(AF_MARKERS, colnames(Zm))
  if (any(af) && length(sp)) Zm[af, sp] <- NA_real_
  Zm
}
af_load <- function(Z, thr, min_extra = 2L) {
  pp <- function(m) Z[[m]] > thr[[m]]
  n <- as.integer(pp("CD223")) + as.integer(pp("CD49a")) + as.integer(pp("CXCR5"))
  pp("FOXP3") & pp("CD25") & n >= as.integer(min_extra)
}

## MYE/LYM panel: the unmixing of this panel has no autofluorescence channel,
## and autofluorescence projects onto CD45, CD1c, CD206, androgen receptor,
## CD14, CD64, CD33 and BDCA-2 (see python/unmixing/af_projection.py). A cell is
## flagged when at least min_extra (default 3) of the seven sinks are positive
## at the same time, a pattern that no real lineage shows (cDC2 = CD1c alone,
## macrophage = CD206/CD64/CD14 without AR/BDCA-2/CD1c).
if (IS_MYE) {
  AF_SINKS <- c("CD1c", "MRC1 (CD206)", "Androgen R", "FCGR1 (CD64)", "BDCA-2", "CD33", "CD14")
  AF_MARKERS <- c("CD45", AF_SINKS)
  af_load <- function(Z, thr, min_extra = 3L) {
    pp <- function(m) !is.na(Z[[m]]) & Z[[m]] > thr[[m]]
    n <- Reduce(`+`, lapply(intersect(AF_SINKS, names(Z)), function(m) as.integer(pp(m))))
    n >= as.integer(min_extra)
  }
}

message("common.R loaded | panel ", PANEL_ID, " | ", length(LINEAGE_MARKERS), " lineage markers, ",
        length(STATE_MARKERS), " state markers")
