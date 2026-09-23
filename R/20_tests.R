## 20_tests.R -- group comparisons, correlations and yield tests on the
## annotated populations: abundance (three denominators and CLR) and expression
## of all non-defining markers (median arcsinh per sample and population,
## z-scored across samples; in addition % positive with the frozen thresholds).
##
## Usage: Rscript R/20_tests.R <folder>
## Input:  <folder>/annotated_cells.csv.gz, transformed_arcsinh.csv.gz, thresholds.csv
## Output: 20_tests.xlsx, 20_abundance_long.csv.gz, 20_expression_long.csv.gz
##
## Rules: a sample enters when it has at least MIN_L cells of the reference
## population (lymphocytes in the NK/T panel, leukocytes in the MYE/LYM panel);
## a test needs at least MIN_N samples per group; two-sided Mann-Whitney p
## values without correction (multiplicity is handled downstream, see
## R/23_robustness.R and METHODS.md); one-sided p values only for a
## pre-specified directional hypothesis (design file).
## Secondary, reported separately: p after adjustment for run and for run plus
## log cell yield (Mann-Whitney on rank residuals), Spearman correlation with
## numeric sample-sheet columns and with the cell yield (raw and run-partial).
suppressPackageStartupMessages({library(data.table); library(writexl)})
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
args <- commandArgs(TRUE); OUT <- args[1]
## Minimum number of cells per sample and population for the expression layer.
## With 20 cells the sample median carries an own uncertainty of roughly the
## size of the between-sample spread; 30 is the default. Samples with few cells
## are systematically different samples (low yield), so the cut-off is not
## neutral; R/24_cutoff_sensitivity.R stages 20/30/50.
MIN_CELLS <- as.integer(Sys.getenv("SFP_MIN_CELLS", "30"))
## Sensitivity option: SFP_FSC_FRACTION=high restricts everything to the large
## FSC fraction.
FSC_MODE <- Sys.getenv("SFP_FSC_FRACTION", "all")
SUF <- paste0(if (FSC_MODE == "high") "_fsc_high" else "",
              if (Sys.getenv("SFP_AF_RULE", "2") != "2") paste0("_af", Sys.getenv("SFP_AF_RULE", "2")) else "")

C <- fread(file.path(OUT, "annotated_cells.csv.gz"), colClasses = c(sample_id = "character"))
Z <- fread(file.path(OUT, "transformed_arcsinh.csv.gz"))[C$row_all]
if (FSC_MODE == "high") { keep <- C$fsc_fraction == "high"; C <- C[keep]; Z <- Z[keep]; message("Sensitivity: large FSC fraction only, ", nrow(C), " cells") }
ST <- fread(file.path(OUT, "thresholds.csv")); thr <- setNames(ST$threshold, ST$marker)
MARKER <- setdiff(PANEL$marker, c(DISCARDED, "Viability", "CD45", PANEL$marker[PANEL$role == "qc"]))   # qc markers are never tested
stopifnot(all(MARKER %in% names(Z)))

## ---- Sample sheet: grouping variables and numeric covariates ----------------------
SAMPLES <- read_samples()
COR_COLS <- DESIGN$correlate_with

## ---- Remove artefacts, cut-off on the reference population ---------------------
n_art <- sum(C$lineage == "artefact")
keep <- C$lineage != "artefact"
C <- C[keep]; Z <- Z[keep]; C[, row := .I]   ## filter Z in parallel, otherwise rows shift
stopifnot(nrow(Z) == nrow(C))
ref <- C[, .(n_reference = .N), by = file]
## sample_id always from the sample sheet (text, leading zeros kept); run from
## the cell table, checked against the sample sheet.
P <- unique(C[, .(file, run)])
stopifnot(!anyDuplicated(P$file))
attach_samples(P, c("sample_id", DESIGN$group_columns, COR_COLS), SAMPLES)
if (any(P$run != SAMPLES$run[match(P$file, SAMPLES$file)])) stop("run in the cell table differs from the sample sheet")
setcolorder(P, c("file", "sample_id", "run"))
for (cc in COR_COLS) set(P, j = cc, value = suppressWarnings(as.numeric(P[[cc]])))
P <- merge(P, ref, by = "file")
P[, included := n_reference >= MIN_L]
message(sprintf("Artefact cells removed: %d | samples: %d, with >= %d reference cells: %d",
                n_art, nrow(P), MIN_L, sum(P$included)))

## ---- Frequencies (zero-filled) --------------------------------------------------------
T_LIN <- c("CD4", "CD8", "gd T", "CD8dim/CD8- TCRgd- T cells")
PARENTS <- c("CD4", "CD8", "ILC")
## MYE/LYM panel: T/NK is one lineage without subsets (it is evaluated in the
## NK/T panel); the subsets lie under the B-cell and myeloid lineages.
if (IS_MYE) {
  T_LIN <- c("NK/T", "T/NK")
  PARENTS <- setdiff(unique(as.character(C$lineage)), c(T_LIN, "doublets", "unassigned"))
  ## Plasma cells are a lineage with exactly one subset; as a parent they would
  ## only produce constants (100 %, CLR 0), so they are tested at lineage level.
  PARENTS <- setdiff(PARENTS, "Plasma cells")
  ## T/NK, doublets and unassigned cells stay in the cell table and in the
  ## leukocyte denominator, but are not evaluated as populations.
  LINEAGES_OUT <- c("NK/T", "T/NK", "T/NK (other panel)", "doublets", "unassigned")
  message("Populations: MYE/LYM panel (subsets under the B-cell and myeloid lineages)")
}
n_cd3 <- C[lineage %in% T_LIN, .(n_CD3 = .N), by = file]
H1 <- C[, .N, by = .(file, unit = lineage)]
H1 <- merge(CJ(file = ref$file, unit = sort(unique(C$lineage))), H1, by = c("file", "unit"), all.x = TRUE)
H1[is.na(N), N := 0L]; H1[, `:=`(level = "lineage", parent = NA_character_)]
sub_grid <- unique(C[lineage %in% PARENTS, .(parent = lineage, unit = sublineage)])
grid2 <- CJ(file = ref$file, k = seq_len(nrow(sub_grid)))
grid2 <- cbind(grid2[, .(file)], sub_grid[grid2$k])
H2 <- C[lineage %in% PARENTS, .N, by = .(file, parent = lineage, unit = sublineage)]
H2 <- merge(grid2, H2, by = c("file", "parent", "unit"), all.x = TRUE)
H2[is.na(N), N := 0L]; H2[, level := "sublineage"]
H <- rbind(H1, H2, fill = TRUE)
H <- merge(H, ref, by = "file")
H <- merge(H, n_cd3, by = "file", all.x = TRUE); H[is.na(n_CD3), n_CD3 := 0L]
H <- merge(H, H1[, .(file, parent = unit, n_parent = N)], by = c("file", "parent"), all.x = TRUE)
H[, is_T := (level == "lineage" & unit %in% T_LIN) | (level == "sublineage" & parent %in% c("CD4", "CD8"))]
H[, pct_reference := 100 * N / n_reference]
H[, pct_CD3       := fifelse(is_T, 100 * N / pmax(n_CD3, 1), NA_real_)]
H[, pct_parent    := fifelse(level == "sublineage", 100 * N / pmax(n_parent, 1), NA_real_)]
## CLR within the composition (lineages: all lineages of the sample;
## sublineages: within the parent), pseudocount 0.5 against zeros. Rounded to 12
## decimals: identical compositions must give identical CLR values (ties);
## otherwise the summation order produces 1e-16 differences that rank tests
## read as non-ties.
H[, clr := { x <- log(N + 0.5); round(x - mean(x), 12) }, by = .(file, level, parent)]
Hl <- melt(H, id.vars = c("file", "level", "parent", "unit", "N"),
           measure.vars = c("pct_reference", "pct_CD3", "pct_parent", "clr"),
           variable.name = "metric", value.name = "value")[is.finite(value)]
Hl <- merge(Hl, P, by = "file")[included == TRUE]

## ---- Expression per sample x population ------------------------------------------
UNITS <- rbind(C[, .(row, file, level = "lineage", parent = NA_character_, unit = lineage)],
               C[lineage %in% PARENTS, .(row, file, level = "sublineage", parent = lineage, unit = sublineage)])
UNITS <- UNITS[unit != "unassigned"]
if (exists("LINEAGES_OUT")) {
  before <- uniqueN(UNITS$unit)
  UNITS <- UNITS[!(unit %in% LINEAGES_OUT)]   # on both levels
  message("Not evaluated as populations: ", paste(LINEAGES_OUT, collapse = ", "),
          " (populations ", before, " -> ", uniqueN(UNITS$unit), ")")
}
n_e <- UNITS[, .(n = .N), by = .(file, level, parent, unit)][n >= MIN_CELLS]
UNITS <- merge(UNITS, n_e, by = c("file", "level", "parent", "unit"))
Zm <- as.matrix(Z[, ..MARKER]); SW <- thr[MARKER]
## Autofluorescence load (see common.R): for flagged cells only the channels
## affected by autofluorescence are set to NA. The cell stays in the abundance
## and in the expression of all other markers. The minimum cell number applies
## again per marker.
## SFP_AF_RULE (sensitivity analysis, R/unmixing/af_sensitivity.R):
##   "2"      NK/T default: mask when FOXP3+ CD25+ and >= 2 of CD223/CD49a/CXCR5
##   "1"/"3"  weaker or stricter rule for the additional markers
##   "off"    no treatment
##   "cells"  flagged cells removed from the expression (not from the abundance)
## The default per panel must match the flag set during annotation: NK/T "2",
## MYE/LYM "3" (>= 3 of 7 autofluorescence sinks).
AF_DEFAULT <- if (IS_MYE) "3" else "2"
AF_RULE <- Sys.getenv("SFP_AF_RULE", AF_DEFAULT)
AFM <- intersect(AF_MARKERS, MARKER)
n_af <- 0L
if ("af_load" %in% names(C) && AF_RULE != "off") {
  af_v <- if (AF_RULE == "cells") C$af_load else af_load(Z, thr, as.integer(AF_RULE))
  n_af <- sum(af_v)
  if (AF_RULE == "cells") { Zm[af_v, ] <- NA_real_ } else { Zm[af_v, AFM] <- NA_real_ }
  message(sprintf("Autofluorescence rule '%s': %d cells (%.2f %%), masked %s",
                  AF_RULE, n_af, 100 * n_af / nrow(C),
                  if (AF_RULE == "cells") "all markers" else paste(AFM, collapse = ", ")))
}
Xmed <- UNITS[, {
  m  <- Zm[row, , drop = FALSE]
  nn <- colSums(!is.na(m))
  md <- suppressWarnings(apply(m, 2, median, na.rm = TRUE))
  pz <- 100 * colMeans(sweep(m, 2, SW, ">"), na.rm = TRUE)
  md[nn < MIN_CELLS] <- NA_real_; pz[nn < MIN_CELLS] <- NA_real_
  c(list(n = .N), as.list(md), setNames(as.list(pz), paste0("pos_", MARKER)),
    setNames(as.list(as.numeric(nn)), paste0("nmk_", MARKER)))
}, by = .(file, level, parent, unit)]
X1 <- melt(Xmed, id.vars = c("file", "level", "parent", "unit", "n"), measure.vars = MARKER,
           variable.name = "marker", value.name = "median_arcsinh")
X2 <- melt(Xmed, id.vars = c("file", "level", "parent", "unit", "n"), measure.vars = paste0("pos_", MARKER),
           variable.name = "marker", value.name = "pct_positive")
X2[, marker := sub("^pos_", "", marker)]
X3 <- melt(Xmed, id.vars = c("file", "level", "parent", "unit", "n"), measure.vars = paste0("nmk_", MARKER),
           variable.name = "marker", value.name = "n_marker")
X3[, marker := sub("^nmk_", "", marker)]
X1[, marker := as.character(marker)]
X <- merge(X1, X2, by = c("file", "level", "parent", "unit", "n", "marker"))
X <- merge(X, X3, by = c("file", "level", "parent", "unit", "n", "marker"))
X <- merge(X, P, by = "file")[included == TRUE]
## na.rm is needed: the autofluorescence mask can leave single samples NA;
## without na.rm, sd() returns NA and the z-score of the whole
## population-marker group would disappear (the test would drop silently).
X[, z := { s <- sd(median_arcsinh, na.rm = TRUE)
           if (is.finite(s) && s > 0) (median_arcsinh - mean(median_arcsinh, na.rm = TRUE)) / s else NA_real_ },
  by = .(level, parent, unit, marker)]

## Defining markers per population (flagged, not removed): a group difference
## in a marker that defines the population reflects a shift of the gate or
## cluster boundary, not a biological finding.
T_DEF <- c("CD3", "CD4", "CD8a", "TCRgd", "Lineage")
ILC_DEF <- c("CD3", "Lineage", "CD94", "c-kit", "CD127", "NKp44", "CRTH2", "CD56", "CD4")
DEF <- list(
  "lineage|B cells" = c("Lineage", "CXCR5", "HLA-DR", "CD3"), "lineage|Myeloid" = c("Lineage", "HLA-DR", "CD3"),
  "lineage|gd T" = T_DEF, "lineage|CD4" = T_DEF, "lineage|CD8" = T_DEF, "lineage|CD8dim/CD8- TCRgd- T cells" = T_DEF,
  "lineage|NK" = c("CD3", "CD56", "CD94", "TCRgd", "CD4"), "lineage|ILC" = ILC_DEF,
  "sublineage|Treg" = c(T_DEF, "FOXP3", "CD25"), "sublineage|CD4 CD45RA+RO-" = c(T_DEF, "CD45RA", "CD45RO"),
  "sublineage|CD4 Trm CD103+" = c(T_DEF, "CD103", "CD69"),
  "sublineage|CD4 memory (mostly CD69+)" = c(T_DEF, "CD45RA", "CD45RO", "CD103", "CD69", "FOXP3", "CD25"),
  "sublineage|CD8 CD45RA+RO-" = c(T_DEF, "CD45RA", "CD45RO"), "sublineage|CD8 Trm" = c(T_DEF, "CD103", "CD69"),
  "sublineage|CD8 CD69+ CD103-" = c(T_DEF, "CD103", "CD69"),
  "sublineage|CD8 CD69-" = c(T_DEF, "CD45RA", "CD45RO", "CD103", "CD69"),
  "sublineage|ILC3 NKp44+" = ILC_DEF, "sublineage|ILC3 NKp44-" = ILC_DEF, "sublineage|ILC1-like (c-kit- CRTH2-)" = ILC_DEF)
if (IS_MYE) {
  B_DEF   <- c("CD19", "CD27", "CD38", "IgD", "IgM", "IgA", "HLA-DR")
  ## Maturity is defined by CD206 or VSIG4, pDC by CD123 and BDCA-2; CD163 and
  ## CD11c are not used by the rules and define no population. The boundary
  ## between B cells and plasma cells runs through CD38, HLA-DR and CD27.
  MNP_DEF <- c("HLA-DR", "CD14", "MRC1 (CD206)", "VSIG4", "FCGR1 (CD64)", "CD68")
  DC_DEF  <- c("HLA-DR", "CD1c", "CD141", "IL-3Ra", "BDCA-2", "CD14")
  ## CD16 co-defines the non-classical monocyte, CD64 the cDC1 (by exclusion).
  ## "myeloid, not further classifiable" is defined by the absence of all
  ## macrophage and DC markers.
  MNP_DEF <- c(MNP_DEF, "CD16")
  DC_DEF  <- c(DC_DEF, "FCGR1 (CD64)")
  REST_DEF <- unique(c(MNP_DEF, DC_DEF))
  MYE_DEF <- list(
    "lineage|B cells" = c("CD19", "HLA-DR", "CD38", "CD27"), "lineage|Plasma cells" = B_DEF, "lineage|Macrophages" = MNP_DEF,
    "lineage|Dendritic cells" = DC_DEF, "lineage|Monocytes" = c("CD14", "CD16", "HLA-DR", "FCGR1 (CD64)"),
    "lineage|other myeloid" = REST_DEF)
  for (sl in unique(c(C$sublineage))) {
    if (grepl("^B |^Plasma cell|^Plasmablast", sl)) MYE_DEF[[paste0("sublineage|", sl)]] <- B_DEF
    else if (grepl("^Macrophage|^Monocyte", sl)) MYE_DEF[[paste0("sublineage|", sl)]] <- MNP_DEF
    else if (grepl("^cDC|^pDC", sl)) MYE_DEF[[paste0("sublineage|", sl)]] <- DC_DEF
    else if (grepl("HLA-DR-negative", sl)) MYE_DEF[[paste0("sublineage|", sl)]] <- c("HLA-DR")
    else if (grepl("not further classifiable", sl)) MYE_DEF[[paste0("sublineage|", sl)]] <- REST_DEF
  }
  DEF <- MYE_DEF
}
DEFt <- rbindlist(lapply(names(DEF), function(k) {
  s <- strsplit(k, "|", fixed = TRUE)[[1]]
  data.table(level = s[1], unit = s[2], marker = DEF[[k]], defining = TRUE)
}))
X <- merge(X, DEFt, by = c("level", "unit", "marker"), all.x = TRUE)
X[is.na(defining), defining := FALSE]

## ---- Baseline per marker and sample --------------------------------------------------
## Sample-wise staining intensity can leak through the fixed unmixing matrix
## into neighbouring channels. A single brightness index per sample does not
## capture this, because the coupling is channel-specific and changes sign.
## The baseline of the same marker is used instead: per sample, the median of
## this channel on a population that cannot carry it (an internal FMO-like
## reference from the data): B cells in the NK/T panel, T/NK cells in the
## MYE/LYM panel (SFP_BASELINE_POP). A marker gets a baseline only when the
## population does not carry it (pooled median below the frozen threshold) and
## at least MIN_CELLS cells per sample are present.
## The baseline does not enter the reported p values; it is used only for the
## additional column p_baseline_adj.
BASELINE_POP <- Sys.getenv("SFP_BASELINE_POP", if (IS_MYE) "T/NK (other panel)" else "B cells")
BASELINE <- NULL
i_b <- which(C$lineage == BASELINE_POP | (if ("sublineage" %in% names(C)) C$sublineage == BASELINE_POP else FALSE))
if (length(i_b) >= 500) {
  cand <- intersect(MARKER, names(thr))
  fit <- cand[vapply(cand, function(m) {
    med <- median(Z[[m]][i_b], na.rm = TRUE); is.finite(med) && med < thr[[m]] }, logical(1))]
  if (length(fit)) {
    B <- rbindlist(lapply(fit, function(m) {
      d <- data.table(file = C$file[i_b], w = Z[[m]][i_b])
      d <- d[, .(n = .N, baseline = median(w, na.rm = TRUE)), by = file][n >= MIN_CELLS]
      d[, .(file, marker = m, baseline)] }))
    BASELINE <- B
    message(sprintf("Baseline per marker: population '%s', %d of %d markers usable, %d sample values",
                    BASELINE_POP, length(fit), length(cand), nrow(B)))
  }
}
if (is.null(BASELINE)) message("Baseline per marker: population '", BASELINE_POP,
                               "' not found or too small -- p_baseline_adj stays NA")

KEEP_P <- c("sample_id", "run", DESIGN$group_columns, COR_COLS, "n_reference")
Xl <- rbind(X[, c(list(file = file, level = level, parent = parent, unit = unit, marker = marker, defining = defining,
                       n = n, metric = "z_median", value = z), .SD), .SDcols = KEEP_P],
            X[, c(list(file = file, level = level, parent = parent, unit = unit, marker = marker, defining = defining,
                       n = n, metric = "pct_positive", value = pct_positive), .SD), .SDcols = KEEP_P])[is.finite(value)]
if (!is.null(BASELINE)) Xl <- merge(Xl, BASELINE, by = c("file", "marker"), all.x = TRUE) else Xl[, baseline := NA_real_]

## ---- Tests -------------------------------------------------------------------------
COMPARISONS <- DESIGN$comparisons
## Cliff's delta, > 0: group b higher
cliff <- function(a, b) mean(sign(outer(b, a, "-")))
test_one <- function(d, v, keys) {
  has_baseline <- "baseline" %in% names(d)
  if (!has_baseline) d <- copy(d)[, baseline := NA_real_]
  if (!is.null(v$str)) d <- d[get(names(v$str)) == unname(v$str)]
  d <- d[get(v$var) %in% c(v$a, v$b)]
  d[, {
    g <- get(v$var); a <- value[g == v$a]; b <- value[g == v$b]
    r <- mw(a, b); r1 <- mw(a, b, alternative = DESIGN$one_sided$alternative)
    e1 <- resid_rank(value, run = factor(run))
    e2 <- resid_rank(value, run = factor(run), yield = log(n_reference))
    e3 <- if (has_baseline && sum(is.finite(baseline)) >= 10)
            resid_rank(value, run = factor(run), yield = log(n_reference), baseline = baseline)
          else rep(NA_real_, length(value))
    list(n_a = r$n1, n_b = r$n2,
         median_a = suppressWarnings(median(a)), median_b = suppressWarnings(median(b)),
         difference = suppressWarnings(median(b) - median(a)),
         cliffs_delta = if (r$tested) cliff(a, b) else NA_real_,
         p_two_sided = r$p, p_one_sided = r1$p,
         p_run_adj = mw(e1[g == v$a], e1[g == v$b])$p,
         p_run_yield_adj = mw(e2[g == v$a], e2[g == v$b])$p,
         p_baseline_adj = mw(e3[g == v$a], e3[g == v$b])$p)
  }, by = keys][, comparison := v$name][]
}
spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  t <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(t$estimate), p = t$p.value, n = sum(ok))
}
correlate <- function(d, keys) d[, {
  out <- list(n = .N)
  for (cc in COR_COLS) {
    a <- spearman(value, get(cc))
    out[[paste0("n_", cc)]] <- a$n; out[[paste0("rho_", cc)]] <- a$rho; out[[paste0("p_", cc)]] <- a$p
  }
  y <- spearman(value, log(n_reference))
  yp <- spearman(resid_rank(value, run = factor(run)), resid_rank(log(n_reference), run = factor(run)))
  c(out, list(rho_yield = y$rho, p_yield = y$p,
              rho_yield_run_partial = yp$rho, p_yield_run_partial = yp$p))
}, by = keys]

K_AB <- c("level", "parent", "unit", "metric")
E_ab <- rbindlist(lapply(COMPARISONS, test_one, d = Hl, keys = K_AB))
E_ab[, p_one_sided := NA_real_]      ## no directional hypothesis for abundances
E_ab[, significant := is.finite(p_two_sided) & p_two_sided < 0.05]
E_ab[, robust_run := significant & is.finite(p_run_adj) & p_run_adj < 0.05]
E_ab[, robust_run_yield := significant & is.finite(p_run_yield_adj) & p_run_yield_adj < 0.05]
C_ab <- correlate(Hl, K_AB)

K_EX <- c("level", "parent", "unit", "marker", "defining", "metric")
E_ex <- rbindlist(lapply(COMPARISONS, test_one, d = Xl, keys = K_EX))
E_ex[, one_sided_allowed := marker %in% DESIGN$one_sided$markers & comparison %in% DESIGN$one_sided$comparisons]
E_ex[one_sided_allowed == FALSE, p_one_sided := NA_real_]
E_ex[, significant := is.finite(p_two_sided) & p_two_sided < 0.05]
E_ex[, significant_directional := one_sided_allowed & is.finite(p_one_sided) & p_one_sided < 0.05]
E_ex[, robust_run := significant & is.finite(p_run_adj) & p_run_adj < 0.05]
E_ex[, robust_run_yield := significant & is.finite(p_run_yield_adj) & p_run_yield_adj < 0.05]
## Consistency: z median and % positive significant in the same direction?
kz <- E_ex[metric == "z_median", .(level, parent, unit, marker, comparison, sig_z = significant, d_z = difference)]
kp <- E_ex[metric == "pct_positive", .(level, parent, unit, marker, comparison, sig_p = significant, d_p = difference)]
kk <- merge(kz, kp, by = c("level", "parent", "unit", "marker", "comparison"))
kk[, consistent := sig_z & sig_p & sign(d_z) == sign(d_p)]
E_ex <- merge(E_ex, kk[, .(level, parent, unit, marker, comparison, consistent)],
              by = c("level", "parent", "unit", "marker", "comparison"), all.x = TRUE)
C_ex <- correlate(Xl[metric == "z_median"], setdiff(K_EX, "metric"))
setorder(E_ab, p_two_sided, na.last = TRUE); setorder(E_ex, p_two_sided, na.last = TRUE)

## ---- Nominal hits (p < 0.05, before any multiplicity assessment) ------------------
HITS <- rbind(
  E_ab[significant == TRUE, .(type = "abundance", level, parent, unit, marker = NA_character_, metric, comparison,
                              n_a, n_b, median_a, median_b, difference, cliffs_delta, p_two_sided, p_one_sided,
                              p_run_adj, p_run_yield_adj, robust_run, robust_run_yield, consistent = NA)],
  E_ex[(significant | significant_directional) & defining == FALSE,
       .(type = "expression", level, parent, unit, marker, metric, comparison, n_a, n_b, median_a, median_b,
         difference, cliffs_delta, p_two_sided, p_one_sided, p_run_adj, p_run_yield_adj, robust_run,
         robust_run_yield, consistent)])
setorder(HITS, p_two_sided)
message(sprintf("Abundance: %d tests, %d with p < 0.05 (%d robust to run, %d to run + yield)",
                sum(is.finite(E_ab$p_two_sided)), sum(E_ab$significant), sum(E_ab$robust_run), sum(E_ab$robust_run_yield)))
message(sprintf("Expression (non-defining, z): %d tests, %d with p < 0.05, %d consistent in %% positive, %d robust to run",
                sum(is.finite(E_ex[defining == FALSE & metric == "z_median"]$p_two_sided)),
                sum(E_ex[defining == FALSE & metric == "z_median"]$significant),
                sum(E_ex[defining == FALSE & metric == "z_median"]$consistent, na.rm = TRUE),
                sum(E_ex[defining == FALSE & metric == "z_median"]$robust_run)))
message("Expected under H0 at 5 %: abundance ", round(0.05 * sum(is.finite(E_ab$p_two_sided))),
        ", expression ", round(0.05 * sum(is.finite(E_ex[defining == FALSE & metric == "z_median"]$p_two_sided))))

settings <- data.table(parameter = c("MIN_L (reference cells per sample)", "MIN_N (per group)",
                                     "MIN_CELLS (per sample x population, expression)", "significance", "correction",
                                     "directional hypothesis", "difference / cliffs_delta", "p_run_adj", "p_run_yield_adj",
                                     "correlations", "yield", "CLR", "artefacts", "autofluorescence load", "FSC fraction",
                                     "date"),
                       value = c(MIN_L, MIN_N, MIN_CELLS, "two-sided Mann-Whitney p < 0.05", "none (see METHODS.md)",
                                 if (length(DESIGN$one_sided$markers)) paste0(paste(DESIGN$one_sided$markers, collapse = "/"),
                                   " in ", paste(DESIGN$one_sided$comparisons, collapse = "; "), ", alternative '",
                                   DESIGN$one_sided$alternative, "' (a vs b)") else "none",
                                 "b minus a; delta > 0 = b higher",
                                 "Mann-Whitney on rank residuals after run",
                                 "Mann-Whitney on rank residuals after run + log(reference cell count)",
                                 if (length(COR_COLS)) paste("Spearman with", paste(COR_COLS, collapse = ", ")) else "none",
                                 "Spearman with log(reference cell count), raw and run-partial",
                                 "log(N + 0.5) centred within the lineage composition or within the parent",
                                 sprintf("%d artefact cells removed before all denominators", n_art),
                                 sprintf("%d cells flagged; only %s set to NA for them, all other markers and the abundance unchanged",
                                         n_af, paste(AFM, collapse = "/")),
                                 if (FSC_MODE == "high") "large fraction only (sensitivity)" else "all cells",
                                 format(Sys.time(), "%Y-%m-%d %H:%M")))
## Comparison definitions for the Python modules (python/)
COMP_TAB <- rbindlist(lapply(COMPARISONS, function(v) data.table(
  comparison = v$name, var = v$var, a = v$a, b = v$b,
  stratum_var = if (is.null(v$str)) NA_character_ else names(v$str),
  stratum_value = if (is.null(v$str)) NA_character_ else unname(v$str))))
write_xlsx(list(Nominal_hits = as.data.frame(HITS), Tests_abundance = as.data.frame(E_ab),
                Tests_expression = as.data.frame(E_ex), Correlations_abundance = as.data.frame(C_ab),
                Correlations_expression = as.data.frame(C_ex), Counts = as.data.frame(H),
                Expression_per_sample = as.data.frame(X), Samples = as.data.frame(P),
                Defining_markers = as.data.frame(DEFt), Comparisons = as.data.frame(COMP_TAB),
                Settings = as.data.frame(settings)),
           file.path(OUT, paste0("20_tests", SUF, ".xlsx")))
fwrite(Hl, file.path(OUT, paste0("20_abundance_long", SUF, ".csv.gz")), compress = "gzip")
fwrite(Xl, file.path(OUT, paste0("20_expression_long", SUF, ".csv.gz")), compress = "gzip")
message("written: 20_tests.xlsx, 20_abundance_long.csv.gz, 20_expression_long.csv.gz")
