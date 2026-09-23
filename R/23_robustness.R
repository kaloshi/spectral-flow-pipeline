## ============================================================================
## 23_robustness.R -- every expression contrast staged by cell cut-off and by
## sample statistic
## ============================================================================
##   Rscript R/23_robustness.R <folder> [<output folder>]
## Input:  annotated_cells.csv.gz, clusters.csv.gz, transformed_arcsinh.csv.gz,
##         thresholds.csv, 20_tests.xlsx (sheet Samples)
## Output: 23_robustness.xlsx
##
## With few cells, a sample value mostly measures its own sampling noise: at
## 20 cells the uncertainty of a sample median can approach the spread between
## samples. No sample-level expression finding is therefore reported without
## staging. Every contrast of the expression layer is computed
##   * at three minimum cell numbers (default 30, with 20 and 50 as neighbours
##     so that the sensitivity stays visible), and
##   * with two sample statistics (median = default; mean = more stable with
##     few cells).
##
## Why no geometric mean: the geometric mean is the usual MFI summary in flow
## cytometry because fluorescence is roughly log-normal. It is already contained
## here: for x >> cofactor, arcsinh(x / c) ~ ln(2x / c), so the arithmetic mean
## of arcsinh values is the geometric mean in the bright range, with defined
## behaviour at zero and for the negative values that exist after unmixing (a
## geometric mean of raw values is undefined there and would force a selection
## of positive cells). Because Mann-Whitney works on ranks, the choice changes
## little.
##
## The z-score is unchanged: it is a monotone standardisation across samples
## and does not change the rank test. Only the first step varies, i.e. how the
## cells of a sample are summarised into one number.
##
## Output: one row per contrast with six p values and a verdict:
##   "robust"                 p < 0.05 at all three cut-offs and both statistics
##   "depends on cut-off"     only at some cut-offs
##   "depends on statistic"   only with the median or only with the mean
##   "not significant"        nowhere
## In addition n per group and Cliff's delta, to show where power is lost
## rather than an effect.
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(writexl) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
args <- commandArgs(TRUE); OUT <- args[1]
DEST <- if (length(args) >= 2) args[2] else OUT
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)
CUTS <- as.integer(strsplit(Sys.getenv("SFP_CUTS", "20,30,50"), ",")[[1]])
MIN_N <- as.integer(Sys.getenv("SFP_MIN_N", "5"))

C  <- fread(file.path(OUT, "annotated_cells.csv.gz"), colClasses = c(sample_id = "character"))
K  <- fread(file.path(OUT, "clusters.csv.gz"), select = "row_all")
Z  <- fread(file.path(OUT, "transformed_arcsinh.csv.gz"))[K$row_all]
ST <- fread(file.path(OUT, "thresholds.csv")); thr <- setNames(ST$threshold, ST$marker)
P  <- as.data.table(readxl::read_excel(file.path(OUT, "20_tests.xlsx"), sheet = "Samples", col_types = "text"))
P  <- P[included == "TRUE"]
MARKER <- setdiff(intersect(PANEL$marker, names(Z)), c(DISCARDED, "CD45", PANEL$marker[PANEL$role == "qc"]))
stopifnot(nrow(C) == nrow(Z))

## Autofluorescence mask as in R/20_tests.R
af <- af_load(Z, thr)
Zm <- as.matrix(Z[, ..MARKER])
AFM <- intersect(AF_MARKERS, MARKER)
if (any(af) && length(AFM)) Zm[af, AFM] <- NA_real_
cat("Cells:", format(nrow(Z), big.mark = " "), "| markers:", length(MARKER),
    "| AF-masked:", sum(af), "cells\n")

## Populations as in R/20_tests.R: lineage and sublineage, without the T/NK
## container in the MYE/LYM panel. The list must match R/20_tests.R exactly.
OUTSIDE <- if (IS_MYE) c("T/NK", "NK/T", "T/NK (other panel)", "doublets", "artefact", "unassigned") else
                       c("artefact", "unassigned")
E <- rbind(data.table(row = seq_len(nrow(C)), file = C$file, level = "lineage",
                      parent = NA_character_, unit = C$lineage),
           data.table(row = seq_len(nrow(C)), file = C$file, level = "sublineage",
                      parent = C$lineage, unit = C$sublineage))
E <- E[!unit %in% OUTSIDE & !(parent %in% OUTSIDE & !is.na(parent))]
## Populations without subdivision appear as lineage and as a sublineage of the
## same name; the sublineage is dropped when it is its own parent.
E <- E[!(level == "sublineage" & unit == parent)]
E <- unique(E, by = c("row", "level", "unit"))

## Two groups per comparison (design file)
groups_of <- function(v, p) {
  s <- if (is.null(v$str)) p else p[get(names(v$str)) == unname(v$str)]
  list(s[get(v$var) == v$a, file], s[get(v$var) == v$b, file])
}
COMP <- setNames(DESIGN$comparisons, vapply(DESIGN$comparisons, `[[`, "", "name"))

delta <- function(a, b) mean(outer(b, a, ">")) - mean(outer(b, a, "<"))
res <- list(); t0 <- Sys.time()
combos <- unique(E[, .(level, parent, unit)])
cat("Populations:", nrow(combos), "| markers:", length(MARKER), "| contrasts:", length(COMP), "\n")
for (i in seq_len(nrow(combos))) {
  k <- combos[i]
  idx <- E[level == k$level & unit == k$unit &
             (is.na(k$parent) | parent == k$parent), .(row, file)]
  if (!nrow(idx)) next
  for (m in MARKER) {
    v <- Zm[idx$row, m]
    D <- data.table(file = idx$file, value = v)[is.finite(value)]
    if (!nrow(D)) next
    stat <- D[, .(n = .N, median = median(value), mean = mean(value)), by = file]
    for (vn in names(COMP)) {
      gr <- groups_of(COMP[[vn]], P); ga <- gr[[1]]; gb <- gr[[2]]
      row <- list(level = k$level, parent = k$parent, unit = k$unit, marker = m, comparison = vn)
      for (cz in CUTS) {
        s <- stat[n >= cz]
        a <- s[file %in% ga]; b <- s[file %in% gb]
        for (nm in c("median", "mean")) {
          sp <- paste0("p_", nm, "_", cz); sd_ <- paste0("d_", nm, "_", cz)
          if (nrow(a) < MIN_N || nrow(b) < MIN_N) { row[[sp]] <- NA_real_; row[[sd_]] <- NA_real_; next }
          row[[sp]] <- suppressWarnings(wilcox.test(a[[nm]], b[[nm]])$p.value)
          row[[sd_]] <- delta(a[[nm]], b[[nm]])
        }
        row[[paste0("n_", cz)]] <- paste0(nrow(a), "/", nrow(b))
      }
      res[[length(res) + 1L]] <- as.data.table(row)
    }
  }
  if (i %% 5 == 0) cat(".", sep = "")
}
cat("\n")
R <- rbindlist(res, fill = TRUE)
pcols <- grep("^p_", names(R), value = TRUE)
R[, n_tested := rowSums(!is.na(.SD)), .SDcols = pcols]
R[, n_significant := rowSums(.SD < 0.05, na.rm = TRUE), .SDcols = pcols]
pm <- grep("^p_median_", names(R), value = TRUE); pt <- grep("^p_mean_", names(R), value = TRUE)
R[, sig_median := rowSums(.SD < 0.05, na.rm = TRUE), .SDcols = pm]
R[, sig_mean := rowSums(.SD < 0.05, na.rm = TRUE), .SDcols = pt]
R[, verdict := fifelse(n_tested == 0, "not testable",
               fifelse(n_significant == n_tested & n_tested >= 4, "robust",
               fifelse(n_significant == 0, "not significant",
               fifelse((sig_median == 0) != (sig_mean == 0), "depends on statistic",
                       "depends on cut-off"))))]
setorderv(R, paste0("p_median_", CUTS[1]), na.last = TRUE)
cat("\n=== Verdicts across all contrasts ===\n")
print(R[, .N, by = verdict][order(-N)], row.names = FALSE)

write_xlsx(list(All = as.data.frame(R), Robust = as.data.frame(R[verdict == "robust"]),
                Depends_on_cutoff = as.data.frame(R[verdict == "depends on cut-off"][order(get(paste0("p_median_", CUTS[1])))]),
                Notes = data.frame(
                  item = c("z-score", "median versus mean", "cut-off", "robust", "depends on cut-off"),
                  meaning = c("unchanged; the standardisation across samples is monotone and does not change the rank test",
                    "first step: how the cells of a sample become one number. Median is the default, the mean is more stable with few cells",
                    "minimum number of cells per sample and population; 30 is the default",
                    "p < 0.05 at all cut-offs and with both statistics",
                    "only at some cut-offs -- usually because single samples with very few cells drop out"))),
           file.path(DEST, "23_robustness.xlsx"))
cat(sprintf("\nwritten: %s (%d contrasts, %.1f min)\n", file.path(DEST, "23_robustness.xlsx"),
            nrow(R), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
