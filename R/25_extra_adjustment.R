## ============================================================================
## 25_extra_adjustment.R -- do additional sample covariates explain the hits?
## ============================================================================
## Usage:  Rscript R/25_extra_adjustment.R <folder> <output.xlsx>
## Input:  <folder>/20_tests.xlsx; design file with EXTRA_ADJUSTMENT
##
## For every nominal hit (all abundance tests with p < 0.05; expression tests of
## the z-median layer with p < 0.05) the same rank test is repeated on rank
## residuals: after each additional covariate alone, after all of them, and
## after all of them together with run and log cell yield. Method and rounding
## as in R/20_tests.R (resid_rank: na.exclude, 9 decimals).
## For every hit in addition: is the readout itself associated with a
## covariate (Kruskal-Wallis for factors, Spearman for numeric covariates), and
## are the two groups of the comparison unbalanced in it (Fisher's exact test
## with simulated p for factors, Mann-Whitney for numeric covariates)?
## ============================================================================
suppressPackageStartupMessages({ library(readxl); library(data.table); library(writexl) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
COV <- DESIGN$extra_adjustment
if (!length(COV)) stop("EXTRA_ADJUSTMENT is empty in the design file")
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; DEST <- args[2]
TAB <- file.path(OUT, "20_tests.xlsx")
SEED <- as.integer(Sys.getenv("SFP_SEED", "1234"))

## ---- Data --------------------------------------------------------------------------
P <- as.data.table(read_excel(TAB, sheet = "Samples", col_types = "text"))
H <- as.data.table(read_excel(TAB, sheet = "Counts", guess_max = 100000))
X <- as.data.table(read_excel(TAB, sheet = "Expression_per_sample", guess_max = 100000))
attach_samples(P, names(COV))
for (cc in names(COV)) {
  if (COV[[cc]] == "numeric") set(P, j = cc, value = suppressWarnings(as.numeric(P[[cc]])))
  else set(P, j = cc, value = factor(P[[cc]]))
}
if (any(vapply(names(COV), function(cc) anyNA(P[[cc]][P$included == "TRUE"]) && COV[[cc]] == "factor", logical(1))))
  stop("missing values in a factor covariate")
P[, `:=`(n_reference = as.numeric(n_reference))]

## ---- Values as in R/20_tests.R ---------------------------------------------------------
Hl <- melt(H, id.vars = c("file", "level", "parent", "unit"),
           measure.vars = c("pct_reference", "pct_CD3", "pct_parent", "clr"),
           variable.name = "metric", value.name = "value")[is.finite(value)]
Hl[, marker := NA_character_]
Xl <- rbind(X[, .(file, level, parent, unit, marker, metric = "z_median", value = z)],
            X[, .(file, level, parent, unit, marker, metric = "pct_positive", value = pct_positive)])[is.finite(value)]
W <- rbind(Hl[, .(file, level, parent, unit, marker, metric = as.character(metric), value)], Xl)
W <- merge(W, P[included == "TRUE", c("file", "run", DESIGN$group_columns, names(COV), "n_reference"), with = FALSE], by = "file")

## ---- Hits --------------------------------------------------------------------------
TA <- as.data.table(read_excel(TAB, sheet = "Tests_abundance", guess_max = 100000))
TE <- as.data.table(read_excel(TAB, sheet = "Tests_expression", guess_max = 100000))
cols <- c("type", "level", "parent", "unit", "marker", "metric", "comparison", "n_a", "n_b", "cliffs_delta", "p_two_sided")
SA <- TA[p_two_sided < 0.05][, `:=`(type = "abundance", marker = NA_character_)]
SE <- TE[metric == "z_median" & p_two_sided < 0.05][, type := "expression"]
HITS <- rbind(SA[, ..cols], SE[, ..cols])
COMP <- setNames(DESIGN$comparisons, vapply(DESIGN$comparisons, `[[`, "", "name"))
subset_of <- function(d, name) {
  v <- COMP[[name]]; if (is.null(v)) stop("unknown comparison: ", name)
  list(d = if (is.null(v$str)) d else d[get(names(v$str)) == unname(v$str)], var = v$var, a = v$a, b = v$b)
}
pmw <- function(x, g, a, b) mw(x[g == a], x[g == b])$p
covs <- function(d, which) setNames(lapply(which, function(cc) if (COV[[cc]] == "factor") droplevels(d[[cc]]) else d[[cc]]), which)

set.seed(SEED)
rows <- rbindlist(lapply(seq_len(nrow(HITS)), function(i) {
  z <- HITS[i]
  d <- W[level == z$level & unit == z$unit & metric == z$metric &
         (is.na(z$parent) | parent == z$parent) &
         (is.na(z$marker) | (!is.na(marker) & marker == z$marker))]
  t <- subset_of(d, z$comparison); d <- t$d
  if (!nrow(d)) return(NULL)
  g <- d[[t$var]]
  ok <- g %in% c(t$a, t$b)
  d <- d[ok]; g <- g[ok]
  if (length(unique(g)) < 2) return(NULL)
  out <- list(n_total = nrow(d), p_raw_recomputed = pmw(d$value, g, t$a, t$b))
  for (cc in names(COV)) out[[paste0("p_", cc, "_adj")]] <- pmw(do.call(resid_rank, c(list(d$value), covs(d, cc))), g, t$a, t$b)
  out[["p_all_extra_adj"]] <- pmw(do.call(resid_rank, c(list(d$value), covs(d, names(COV)))), g, t$a, t$b)
  out[["p_full_adj"]] <- pmw(do.call(resid_rank, c(list(d$value, run = factor(d$run), yield = log(d$n_reference)),
                                                    covs(d, names(COV)))), g, t$a, t$b)
  for (cc in names(COV)) {
    x <- d[[cc]]
    out[[paste0("p_value_vs_", cc)]] <- if (COV[[cc]] == "factor")
      tryCatch(kruskal.test(d$value, droplevels(x))$p.value, error = function(e) NA_real_) else
      tryCatch(suppressWarnings(cor.test(d$value, x, method = "spearman"))$p.value, error = function(e) NA_real_)
    out[[paste0("p_", cc, "_unbalanced")]] <- if (COV[[cc]] == "factor")
      tryCatch(fisher.test(table(droplevels(x), g), simulate.p.value = TRUE, B = 20000)$p.value, error = function(e) NA_real_) else
      tryCatch(mw(x[g == t$a], x[g == t$b])$p, error = function(e) NA_real_)
  }
  cbind(z, as.data.table(out))
}), fill = TRUE)

## ---- Distribution of the covariates across the groups ----------------------------------
PE <- P[included == "TRUE"]
MAIN <- Filter(function(v) is.null(v$str), DESIGN$comparisons)
set.seed(SEED)
distribution <- rbindlist(lapply(names(COV), function(cc) rbindlist(lapply(MAIN, function(v) {
  g <- PE[[v$var]]
  p <- if (COV[[cc]] == "factor") fisher.test(table(PE[[cc]], g), simulate.p.value = TRUE, B = 20000)$p.value else
         mw(PE[[cc]][g == v$a], PE[[cc]][g == v$b])$p
  data.table(covariate = cc, grouping = v$var, p = p) }))))
cross <- lapply(names(COV)[COV == "factor"], function(cc)
  dcast(PE[, .N, by = c(cc, DESIGN$group_columns)], as.formula(paste(cc, "~", paste(DESIGN$group_columns, collapse = " + "))),
        value.var = "N", fill = 0))

## ---- Output ------------------------------------------------------------------------
holds <- function(x) sum(is.finite(x) & x < 0.05)
summ <- c(list(hits = nrow(rows),
               raw_reproduced = sum(abs(round(rows$p_raw_recomputed, 6) - round(rows$p_two_sided, 6)) < 1e-9, na.rm = TRUE)),
          setNames(lapply(names(COV), function(cc) holds(rows[[paste0("p_", cc, "_adj")]])), paste0("holds_after_", names(COV))),
          list(holds_after_all_extra = holds(rows$p_all_extra_adj), holds_after_all = holds(rows$p_full_adj)),
          setNames(lapply(names(COV), function(cc) holds(rows[[paste0("p_value_vs_", cc)]])), paste0("value_depends_on_", names(COV))),
          setNames(lapply(names(COV), function(cc) holds(rows[[paste0("p_", cc, "_unbalanced")]])), paste0("unbalanced_", names(COV))))
summ <- as.data.table(summ)
sheets <- list(Summary = summ, Hits = rows, Distribution = distribution)
for (j in seq_along(cross)) sheets[[paste0("Cross_table_", j)]] <- cross[[j]]
write_xlsx(sheets, DEST)
print(summ)
cat("\nwritten:", DEST, "\n")
