## ============================================================================
## 26_cell_count_bias.R -- do the nominal hits depend on cell numbers?
## ============================================================================
## Usage:  Rscript R/26_cell_count_bias.R <folder> <output.xlsx>
## Input:  <folder>/20_tests.xlsx; design file
##
## Two different questions that must not be mixed:
##   abundance:   the cell number of the population IS the numerator of the
##                value (fraction = N / denominator). An association of value
##                and N is enforced by construction; an adjustment for it would
##                be circular and is only reported for information (column
##                circular = yes). Informative is the DENOMINATOR (reference
##                cells, CD3 cells, parent population), i.e. the yield.
##   expression:  the number of cells in the median (n_marker) does not enter
##                the value. An association of value and cell number would be a
##                genuine sign of a cell-number bias (few cells -> unstable
##                median). Here the adjustment is meaningful.
## Method as in R/20_tests.R: exact Mann-Whitney (mw()), rank residuals
## (resid_rank), Spearman.
## ============================================================================
suppressPackageStartupMessages({ library(readxl); library(data.table); library(writexl) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; DEST <- args[2]
TAB <- file.path(OUT, "20_tests.xlsx")

P <- as.data.table(read_excel(TAB, sheet = "Samples", col_types = "text"))
H <- as.data.table(read_excel(TAB, sheet = "Counts", guess_max = 100000))
X <- as.data.table(read_excel(TAB, sheet = "Expression_per_sample", guess_max = 100000))
PE <- P[included == "TRUE"][, n_reference := as.numeric(n_reference)]
H <- H[file %in% PE$file]
X <- X[file %in% PE$file]
COMP <- setNames(DESIGN$comparisons, vapply(DESIGN$comparisons, `[[`, "", "name"))
FAC <- DESIGN$factorial
if (!is.null(FAC)) PE[, cell := paste0(get(FAC$outer$var), ":", get(FAC$inner$var))]
CELLS <- if (is.null(FAC)) character(0) else
  as.vector(t(outer(FAC$outer$levels, FAC$inner$levels, paste, sep = ":")))

pmwv <- function(a, b) mw(a, b)$p
delta <- function(a, b) { a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (!length(a) || !length(b)) NA_real_ else cliffs_delta(a, b) }
groups_ab <- function(d, v, col) {
  s <- if (is.null(v$str)) d else d[get(names(v$str)) == unname(v$str)]
  list(s[get(v$var) == v$a][[col]], s[get(v$var) == v$b][[col]])
}

## ---------------------------------------------------------------------------
## A) Cell numbers themselves: do the groups differ?
## ---------------------------------------------------------------------------
group_row <- function(dt, what, level = NA_character_, parent = NA_character_) {
  ## dt: file + count, exactly once per sample
  d <- merge(dt, PE[, c("file", DESIGN$group_columns, if (length(CELLS)) "cell"), with = FALSE], by = "file")
  out <- list(level = level, parent = parent, unit = what)
  for (cl in CELLS) { x <- d[cell == cl, count]; out[[paste0("n_", cl)]] <- length(x); out[[paste0("median_", cl)]] <- median(x) }
  if (length(CELLS)) out$p_four_cells <- tryCatch(kruskal.test(d$count, factor(d$cell))$p.value, error = function(e) NA_real_)
  for (vn in names(COMP)) { g <- groups_ab(d, COMP[[vn]], "count"); out[[paste0("p_", vn)]] <- pmwv(g[[1]], g[[2]]) }
  for (vn in names(COMP)) { g <- groups_ab(d, COMP[[vn]], "count"); out[[paste0("delta_", vn)]] <- delta(g[[1]], g[[2]]) }
  c(out, list(median_all = median(d$count), min_cells = min(d$count),
              frac_below_30 = mean(d$count < 30), frac_zero = mean(d$count == 0)))
}
denominators <- rbindlist(list(
  group_row(PE[, .(file, count = n_reference)], "DENOMINATOR: yield (n_reference)"),
  group_row(unique(H[, .(file, count = n_CD3)])[is.finite(count)], "DENOMINATOR: n_CD3")))
per_unit <- rbindlist(lapply(split(H, by = c("level", "parent", "unit"), drop = TRUE), function(h)
  group_row(h[, .(file, count = N)], h$unit[1], h$level[1], h$parent[1])))
cell_counts <- rbind(denominators, per_unit)
if (length(CELLS)) setorder(cell_counts, p_four_cells)

## ---------------------------------------------------------------------------
## B) per hit: cell-number bias
## ---------------------------------------------------------------------------
Hl <- melt(H, id.vars = c("file", "level", "parent", "unit", "N", "n_reference", "n_CD3", "n_parent"),
           measure.vars = c("pct_reference", "pct_CD3", "pct_parent", "clr"),
           variable.name = "metric", value.name = "value")[is.finite(value)]
Hl[, metric := as.character(metric)]
Hl[, `:=`(marker = NA_character_, cells = as.numeric(N),
          denominator = fcase(metric == "pct_reference", as.numeric(n_reference),
                              metric == "pct_CD3",       as.numeric(n_CD3),
                              metric == "pct_parent",    as.numeric(n_parent),
                              metric == "clr",           as.numeric(n_reference)),
          denominator_kind = fcase(metric == "pct_reference", "reference cells of the sample",
                                   metric == "pct_CD3",       "CD3 cells of the sample",
                                   metric == "pct_parent",    "parent population",
                                   metric == "clr",           "reference cells of the sample (reference)"),
          cell_count_kind = "N of the population (= numerator)", circular = "yes")]
Xl <- rbind(
  X[, .(file, level, parent, unit, marker, metric = "z_median", value = z,
        cells = as.numeric(n_marker), denominator = as.numeric(n_reference))],
  X[, .(file, level, parent, unit, marker, metric = "pct_positive", value = pct_positive,
        cells = as.numeric(n_marker), denominator = as.numeric(n_reference))])[is.finite(value)]
Xl[, `:=`(denominator_kind = "reference cells of the sample", cell_count_kind = "cells in the median (n_marker)", circular = "no")]
cols <- c("file", "level", "parent", "unit", "marker", "metric", "value", "cells", "denominator",
          "denominator_kind", "cell_count_kind", "circular")
W <- rbind(Hl[, ..cols], Xl[, ..cols])
W <- merge(W, PE[, c("file", "run", DESIGN$group_columns, "n_reference"), with = FALSE], by = "file")

TA <- as.data.table(read_excel(TAB, sheet = "Tests_abundance", guess_max = 100000))
TE <- as.data.table(read_excel(TAB, sheet = "Tests_expression", guess_max = 100000))
hc <- c("type", "level", "parent", "unit", "marker", "metric", "comparison", "n_a", "n_b", "cliffs_delta", "p_two_sided")
SA <- TA[p_two_sided < 0.05][, `:=`(type = "abundance", marker = NA_character_)]
SE <- TE[metric == "z_median" & p_two_sided < 0.05][, type := "expression"]
B <- rbind(SA[, ..hc], SE[, ..hc])

pmw <- function(x, g, a, b) pmwv(x[g == a], x[g == b])
sp_rho <- function(x, y) {
  r <- tryCatch(suppressWarnings(cor.test(x, y, method = "spearman")), error = function(e) NULL)
  if (is.null(r)) c(NA_real_, NA_real_) else c(unname(r$estimate), r$p.value)
}
fetch <- function(z) {
  d <- W[level == z$level & unit == z$unit & metric == z$metric &
         (is.na(z$parent) | parent == z$parent) &
         (is.na(z$marker) | (!is.na(marker) & marker == z$marker))]
  v <- COMP[[z$comparison]]; if (is.null(v)) stop("unknown comparison: ", z$comparison)
  if (!is.null(v$str)) d <- d[get(names(v$str)) == unname(v$str)]
  if (!nrow(d)) return(NULL)
  g <- d[[v$var]]; ok <- g %in% c(v$a, v$b)
  list(d = d[ok], g = g[ok], a = v$a, b = v$b)
}

rows <- rbindlist(lapply(seq_len(nrow(B)), function(i) {
  z <- B[i]; h <- fetch(z); if (is.null(h)) return(NULL)
  d <- h$d; g <- h$g
  if (length(unique(g)) < 2) return(NULL)
  za <- d$cells[g == h$a]; zb <- d$cells[g == h$b]
  na_ <- d$denominator[g == h$a]; nb_ <- d$denominator[g == h$b]
  r_z <- sp_rho(d$value, d$cells); r_n <- sp_rho(d$value, d$denominator)
  e_z <- resid_rank(d$value, zz = log10(d$cells + 1))
  e_n <- resid_rank(d$value, nn = log10(d$denominator + 1))
  e_a <- resid_rank(d$value, yield = log(d$n_reference))
  data.table(z, n_total = nrow(d),
    p_raw_recomputed = pmw(d$value, g, h$a, h$b),
    cell_count_kind = d$cell_count_kind[1], circular = d$circular[1],
    median_cells_a = median(za), median_cells_b = median(zb), min_cells = min(d$cells),
    p_cells_between_groups = pmwv(za, zb), delta_cells = delta(za, zb),
    frac_zero_cells = mean(d$cells == 0), frac_below_30_cells = mean(d$cells < 30),
    rho_value_cells = r_z[1], p_value_cells = r_z[2],
    p_cells_adj = pmw(e_z, g, h$a, h$b),
    denominator_kind = d$denominator_kind[1],
    median_denominator_a = median(na_), median_denominator_b = median(nb_),
    p_denominator_between_groups = pmwv(na_, nb_), delta_denominator = delta(na_, nb_),
    rho_value_denominator = r_n[1], p_value_denominator = r_n[2],
    p_denominator_adj = pmw(e_n, g, h$a, h$b),
    p_yield_adj = pmw(e_a, g, h$a, h$b))
}), fill = TRUE)

## ---------------------------------------------------------------------------
## C) Sensitivity: does the hit hold when only cell-rich samples count?
##    abundance  -> threshold on the YIELD (reference cells), not on N (that
##                  would select on the value itself)
##    expression -> threshold on the cells in the median
## ---------------------------------------------------------------------------
THRESHOLDS <- list(abundance = c(500, 1000), expression = c(50, 100))
sens <- rbindlist(lapply(seq_len(nrow(B)), function(i) {
  z <- B[i]; h <- fetch(z); if (is.null(h)) return(NULL)
  rbindlist(lapply(THRESHOLDS[[z$type]], function(s) {
    v <- if (z$type == "abundance") h$d$n_reference else h$d$cells
    k <- v >= s; d <- h$d[k]; g <- h$g[k]
    data.table(z, threshold_on = if (z$type == "abundance") "yield n_reference" else "cells in the median",
      threshold = s, n_a_remaining = sum(g == h$a), n_b_remaining = sum(g == h$b),
      p_remaining = pmwv(d$value[g == h$a], d$value[g == h$b]),
      delta_remaining = delta(d$value[g == h$a], d$value[g == h$b]))
  }))
}), fill = TRUE)

## ---------------------------------------------------------------------------
holds <- function(x) sum(is.finite(x) & x < 0.05)
summ <- data.table(hits = nrow(rows),
  raw_reproduced = sum(abs(round(rows$p_raw_recomputed, 6) - round(rows$p_two_sided, 6)) < 1e-9, na.rm = TRUE),
  cells_unequal_between_groups = holds(rows$p_cells_between_groups),
  denominator_unequal_between_groups = holds(rows$p_denominator_between_groups),
  value_depends_on_cells = holds(rows$p_value_cells),
  value_depends_on_denominator = holds(rows$p_value_denominator),
  expr_value_depends_on_cells = holds(rows[type == "expression", p_value_cells]),
  expr_holds_after_cells_adj = holds(rows[type == "expression", p_cells_adj]),
  expr_hits = nrow(rows[type == "expression"]),
  holds_after_denominator_adj = holds(rows$p_denominator_adj),
  holds_after_yield_adj = holds(rows$p_yield_adj),
  cell_counts_unequal_four_cells = if (length(CELLS)) sum(cell_counts$p_four_cells < 0.05, na.rm = TRUE) else NA_integer_,
  cell_count_rows = nrow(cell_counts))

write_xlsx(list(Summary = summ, Cell_counts = cell_counts, Hits = rows, Sensitivity = sens), DEST)
print(t(summ))
cat("\nwritten:", DEST, "\n")
