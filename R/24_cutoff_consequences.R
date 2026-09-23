## ============================================================================
## 24_cutoff_consequences.R -- what a higher cell cut-off costs
## ============================================================================
##   Rscript R/24_cutoff_consequences.R <folder> [<output folder>]
## Input:  annotated_cells.csv.gz, 20_tests.xlsx (sheet Samples)
## Output: 24_cutoff_consequences.xlsx
##
## The cut-off is the minimum number of cells per sample and population for the
## expression layer (MIN_CELLS, default 30), not the inclusion of a sample by
## its reference cell count (MIN_L), which is unaffected.
##
## Per population and cut-off:
##   * how many samples remain, per cell of the two-by-two design (if set)
##   * how many of the comparisons remain testable (>= MIN_N per group)
##   * whether the population drops out completely
## and across all populations the sum of testable contrasts.
##
## The cut-off is not neutral: samples with few cells are systematically
## different samples (low yield). The median reference cell count of the
## retained samples is therefore reported as well.
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(writexl) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
args <- commandArgs(TRUE); OUT <- args[1]
DEST <- if (length(args) >= 2) args[2] else OUT
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)
CUTS <- as.integer(strsplit(Sys.getenv("SFP_CUTS", "20,30,50,100,200"), ",")[[1]])
MIN_N <- as.integer(Sys.getenv("SFP_MIN_N", "5"))
DEFAULT_CUT <- as.integer(Sys.getenv("SFP_MIN_CELLS", "30"))

C <- fread(file.path(OUT, "annotated_cells.csv.gz"), colClasses = c(sample_id = "character"))
P <- as.data.table(readxl::read_excel(file.path(OUT, "20_tests.xlsx"), sheet = "Samples", col_types = "text"))
P <- P[included == "TRUE"][, n_reference := as.numeric(n_reference)]
OUTSIDE <- c("unassigned", "doublets", "artefact")
if (IS_MYE) OUTSIDE <- c(OUTSIDE, grep("^T/NK|^NK/T", unique(c(C$lineage, C$sublineage)), value = TRUE))

E <- rbind(C[, .(file, level = "lineage", unit = lineage)],
           C[, .(file, level = "sublineage", unit = sublineage)])
E <- E[!unit %in% OUTSIDE]
N <- E[, .(n = .N), by = .(file, level, unit)]
N <- merge(N, P[, c("file", DESIGN$group_columns, "n_reference"), with = FALSE], by = "file")

## Number of comparisons with at least MIN_N samples in both groups
testable <- function(d) {
  sum(vapply(DESIGN$comparisons, function(v) {
    s <- if (is.null(v$str)) d else d[get(names(v$str)) == unname(v$str)]
    sum(s[[v$var]] == v$a) >= MIN_N && sum(s[[v$var]] == v$b) >= MIN_N }, logical(1)))
}
FAC <- DESIGN$factorial
cell_counts <- function(d) {
  if (is.null(FAC)) return(list())
  out <- list()
  for (o in FAC$outer$levels) for (i in FAC$inner$levels)
    out[[paste0("n_", o, ":", i)]] <- sum(d[[FAC$outer$var]] == o & d[[FAC$inner$var]] == i)
  out
}
res <- list()
for (cz in CUTS) {
  S <- N[n >= cz]
  for (k in unique(N[, .(level, unit)])[, paste(level, unit, sep = "|")]) {
    tl <- strsplit(k, "|", fixed = TRUE)[[1]]
    d <- S[level == tl[1] & unit == tl[2]]
    all_ <- N[level == tl[1] & unit == tl[2]]
    res[[length(res) + 1L]] <- as.data.table(c(list(
      cutoff = cz, level = tl[1], unit = tl[2],
      samples = nrow(d), samples_at_default = nrow(all_[n >= DEFAULT_CUT]), samples_total = nrow(all_)),
      cell_counts(d),
      list(contrasts = if (nrow(d)) testable(d) else 0L,
           reference_median = if (nrow(d)) median(d$n_reference) else NA_real_)))
  }
}
R <- rbindlist(res)
R[, loss_vs_default := samples_at_default - samples]

U <- R[, .(populations = uniqueN(paste(level, unit)),
           populations_with_contrast = uniqueN(paste(level, unit)[contrasts > 0]),
           testable_contrasts = sum(contrasts),
           samples_sum = sum(samples),
           reference_median = round(median(reference_median, na.rm = TRUE))), by = cutoff]
cat("\n=== Consequences of the cut-off across all populations ===\n"); print(U, row.names = FALSE)
cat("\n'testable contrasts' = population x comparison with at least ", MIN_N,
    " samples per group. Every marker is tested on each of them.\n", sep = "")

W <- dcast(R, level + unit ~ cutoff, value.var = c("samples", "contrasts"))
setorderv(W, c("level", paste0("samples_", DEFAULT_CUT)), c(1, -1))

write_xlsx(list(Summary = as.data.frame(U), Per_population = as.data.frame(W),
                Complete = as.data.frame(R)),
           file.path(DEST, "24_cutoff_consequences.xlsx"))
cat("\nwritten:", file.path(DEST, "24_cutoff_consequences.xlsx"), "\n")
