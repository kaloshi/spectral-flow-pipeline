## ============================================================================
## 21_factorial_layer.R -- two-by-two layer per population: Kruskal-Wallis
## across the four cells with Dunn post-hoc, and a composition-level PERMANOVA
## ============================================================================
## Usage:  Rscript R/21_factorial_layer.R <folder>
## Input:  <folder>/20_abundance_long.csv.gz; design file with FACTORIAL
## Output: 21_factorial_layer.xlsx
##
## The four cells are formed by an outer and an inner factor (design file).
## Only four contrasts are tested, not all pairs:
##   1, 2  inner factor within each level of the outer factor
##   3, 4  outer factor within each level of the inner factor
## The two diagonal contrasts mix both factors and are not reported.
##
## Per population x metric:
##   Kruskal-Wallis across the four cells (unpaired)
##   -> Dunn post-hoc on the global ranks, restricted to the four contrasts
##   -> Mann-Whitney per contrast for direct comparison with R/20_tests.R
##
## Multiplicity, several procedures side by side, so that the cost of each
## choice is visible: Bonferroni, Benjamini-Hochberg, Storey's adaptive q value
## (estimates the proportion of true null hypotheses) and independent hypothesis
## weighting (IHW), with the median cell count of the population as covariate
## (independent of the p value under H0, but informative about power). In
## addition a hierarchical procedure: lineages first, and sublineages only
## below a lineage whose omnibus test was rejected.
##
## Composition-level omnibus test before any single test: PERMANOVA on the
## Aitchison distance (Euclidean distance of the CLR values) of the lineage
## composition, with a test of homogeneity of dispersion (betadisper).
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(writexl); library(rstatix) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
FAC <- DESIGN$factorial
if (is.null(FAC)) stop("FACTORIAL is not set in the design file")
args <- commandArgs(TRUE); OUT <- args[1]
OV <- FAC$outer$var; OL <- FAC$outer$levels; IV <- FAC$inner$var; IL <- FAC$inner$levels

H <- fread(file.path(OUT, "20_abundance_long.csv.gz"), colClasses = c(sample_id = "character"))[included == TRUE]
## fread reads missing metadata as empty strings, not NA: a filter on !is.na()
## would let samples without metadata through. Filter positively on the
## allowed levels instead.
n_before <- uniqueN(H$file)
H <- H[get(OV) %in% OL & get(IV) %in% IL & unit != "unassigned"]
message("Samples without metadata excluded: ", n_before - uniqueN(H$file))
CELLS <- c(paste0(OL[1], ":", IL[1]), paste0(OL[1], ":", IL[2]), paste0(OL[2], ":", IL[1]), paste0(OL[2], ":", IL[2]))
H[, cell := factor(paste0(get(OV), ":", get(IV)), levels = CELLS)]
CONTRASTS <- data.table(
  group1 = CELLS[c(1, 3, 1, 2)], group2 = CELLS[c(2, 4, 3, 4)],
  comparison = c(sprintf("%s = %s: %s %s vs %s", OV, OL[1], IV, IL[2], IL[1]),
                 sprintf("%s = %s: %s %s vs %s", OV, OL[2], IV, IL[2], IL[1]),
                 sprintf("%s = %s: %s %s vs %s", IV, IL[1], OV, OL[2], OL[1]),
                 sprintf("%s = %s: %s %s vs %s", IV, IL[2], OV, OL[2], OL[1])))
message(sprintf("Samples: %d (%s)", uniqueN(H$file),
                paste(sprintf("%s %d", CELLS, vapply(CELLS, function(g) uniqueN(H[cell == g, file]), 1L)), collapse = ", ")))

## ---- Tests per population x metric ------------------------------------------------
keys <- c("metric", "level", "parent", "unit")
H[, parent := fifelse(is.na(parent), "", parent)]
combos <- unique(H[, ..keys])
res <- list(); omni <- list()
for (i in seq_len(nrow(combos))) {
  k <- combos[i]
  d <- H[metric == k$metric & level == k$level & parent == k$parent & unit == k$unit & is.finite(value)]
  n_gr <- d[, .N, by = cell]
  if (nrow(n_gr) < 4 || min(n_gr$N) < MIN_N) next
  kw <- kruskal.test(value ~ cell, data = d)
  omni[[length(omni) + 1]] <- cbind(k, data.table(
    n = nrow(d), chi2 = unname(kw$statistic), df = unname(kw$parameter), p_kruskal = kw$p.value,
    eta2 = (unname(kw$statistic) - 3) / (nrow(d) - 4)))          ## eta^2 = (H - k + 1) / (n - k) with k = 4 groups
  dn <- as.data.table(dunn_test(d, value ~ cell, p.adjust.method = "none"))
  dn <- merge(CONTRASTS, dn[, .(group1, group2, statistic, p)], by = c("group1", "group2"))
  ## Mann-Whitney on the same two cells, as in R/20_tests.R
  dn[, p_mw := vapply(seq_len(.N), function(j) {
    a <- d[cell == group1[j], value]; b <- d[cell == group2[j], value]
    suppressWarnings(wilcox.test(a, b)$p.value) }, numeric(1))]
  dn[, `:=`(n1 = vapply(group1, function(g) sum(d$cell == g), 1L),
            n2 = vapply(group2, function(g) sum(d$cell == g), 1L),
            median1 = vapply(group1, function(g) median(d[cell == g, value]), 1.0),
            median2 = vapply(group2, function(g) median(d[cell == g, value]), 1.0))]
  ## Cliff's delta (group2 vs group1) with bootstrap interval
  ci <- t(vapply(seq_len(nrow(dn)), function(j) {
    a <- d[cell == dn$group1[j], value]; b <- d[cell == dn$group2[j], value]
    cliffs_delta_ci(a, b, B = as.integer(Sys.getenv("SFP_BOOT", "2000"))) }, numeric(3)))
  dn[, `:=`(delta = ci[, "delta"], ci_lower = ci[, "lower"], ci_upper = ci[, "upper"])]
  res[[length(res) + 1]] <- cbind(k[rep(1, nrow(dn))], dn,
                                  data.table(p_kruskal = kw$p.value, n_cells_median = median(d$N)))
}
R <- rbindlist(res); O <- rbindlist(omni)
message(sprintf("Population x metric combinations: %d, single contrasts: %d", nrow(O), nrow(R)))

## ---- Multiplicity --------------------------------------------------------------------
R[, q_BH_unit   := p.adjust(p, "BH"), by = .(metric, level, parent, unit)]
R[, q_BH_metric := p.adjust(p, "BH"), by = metric]
R[, q_bonferroni := p.adjust(p, "bonferroni"), by = metric]
R[, q_BH_all    := p.adjust(p, "BH")]
## Storey: falls back to pi0 = 1 (equivalent to BH) when the p value
## distribution is too sparse for the pi0 estimate.
st <- try(qvalue::qvalue(R$p), silent = TRUE)
if (!inherits(st, "try-error")) { R[, q_storey := st$qvalues]; pi0 <- st$pi0 } else {
  st <- try(qvalue::qvalue(R$p, pi0 = 1), silent = TRUE)
  R[, q_storey := if (inherits(st, "try-error")) q_BH_all else st$qvalues]
  pi0 <- 1 }
message(sprintf("Storey pi0 (estimated proportion of true null hypotheses): %.3f", pi0))
## IHW splits the hypotheses into random folds, and its weights also depend on
## the numerical solver; without a fixed seed the q values change between runs
## and platforms. They are reported as a side-by-side column only.
set.seed(as.integer(Sys.getenv("SFP_SEED", "1234")))
ih <- try(IHW::ihw(R$p, covariates = log10(pmax(R$n_cells_median, 1)), alpha = 0.05, nbins = 4), silent = TRUE)
R[, q_ihw := if (inherits(ih, "try-error")) NA_real_ else IHW::adj_pvalues(ih)]
## Hierarchical: lineages first, sublineages only below a rejected parent
O[, q_kruskal_BH := p.adjust(p_kruskal, "BH"), by = .(metric, level)]
parents_ok <- O[level == "lineage" & q_kruskal_BH < 0.05, paste(metric, unit)]
R[, parent_key := fifelse(level == "sublineage", paste(metric, parent), paste(metric, unit))]
R[, tested_hierarchically := level == "lineage" | parent_key %in% parents_ok]
R[, q_hierarchical := NA_real_]
R[tested_hierarchically == TRUE, q_hierarchical := p.adjust(p, "BH"), by = .(metric, level)]

## ---- PERMANOVA on the Aitchison distance of the lineage composition --------------
CL <- dcast(H[metric == "clr" & level == "lineage"], as.formula(paste("file +", OV, "+", IV, "~ unit")), value.var = "value")
mat <- as.matrix(CL[, setdiff(names(CL), c("file", OV, IV)), with = FALSE])
complete <- complete.cases(mat)
set.seed(1)
pm <- vegan::adonis2(as.formula(paste("stats::dist(mat[complete, ]) ~", OV, "*", IV)),
                     data = as.data.frame(CL[complete, c(OV, IV), with = FALSE]),
                     permutations = 9999, by = "terms")
PM <- as.data.table(as.data.frame(pm), keep.rownames = "Term")
message("\nPERMANOVA on the Aitchison distance of the lineages (", sum(complete), " samples):")
print(PM)
## Outer factor alone (one-way PERMANOVA). It is computed before the dispersion
## test, which draws its permutations from the same random stream.
pm1 <- vegan::adonis2(stats::dist(mat[complete, ]) ~ group,
                      data = data.frame(group = CL[[OV]][complete]), permutations = 9999)
PM1 <- as.data.table(as.data.frame(pm1), keep.rownames = "Term")
disp <- vegan::betadisper(stats::dist(mat[complete, ]), interaction(CL[[OV]][complete], CL[[IV]][complete]))
pdisp <- vegan::permutest(disp, permutations = 999)$tab[1, "Pr(>F)"]
message("Homogeneity of dispersion (betadisper): p = ", signif(pdisp, 3))

summ <- rbindlist(lapply(c("p", "p_mw", "q_BH_unit", "q_BH_metric", "q_BH_all", "q_storey",
                           "q_ihw", "q_bonferroni", "q_hierarchical"), function(s)
  data.table(quantity = s, significant = sum(R[[s]] < 0.05, na.rm = TRUE), tested = sum(!is.na(R[[s]])))))
print(summ)

## ---- Dunn versus Mann-Whitney ----------------------------------------------------------
R[, switch := fcase(p < 0.05 & p_mw >= 0.05, "Dunn only",
                    p >= 0.05 & p_mw < 0.05, "Mann-Whitney only",
                    p < 0.05 & p_mw < 0.05, "both",
                    default = "neither")]
print(R[, .N, by = switch][order(-N)])
message(sprintf("Rank correlation of the p values: %.3f | largest difference: %.3f",
                cor(R$p, R$p_mw, method = "spearman"), max(abs(R$p - R$p_mw))))

## ---- Output ------------------------------------------------------------------------
setorder(R, metric, p)
write_xlsx(list(
  Contrasts = as.data.frame(R[, .(metric, level, parent, unit, comparison, n1, n2,
                                  median1, median2, cliffs_delta = round(delta, 3),
                                  ci_lower = round(ci_lower, 3), ci_upper = round(ci_upper, 3),
                                  p_kruskal, p_dunn = p, p_mann_whitney = p_mw,
                                  q_BH_unit, q_BH_metric, q_BH_all, q_storey, q_ihw,
                                  q_bonferroni, q_hierarchical, switch)]),
  Kruskal_Wallis = as.data.frame(O[order(metric, p_kruskal)]),
  Method_summary = as.data.frame(summ),
  PERMANOVA = as.data.frame(PM),
  PERMANOVA_outer_factor = as.data.frame(PM1),
  Dunn_vs_MW = as.data.frame(R[switch %in% c("Dunn only", "Mann-Whitney only"),
                               .(metric, level, unit, comparison, p_dunn = p, p_mann_whitney = p_mw, switch)]),
  Reading_guide = data.frame(
    item = c("cells", "contrasts", "omnibus", "post-hoc", "Storey pi0", "IHW covariate", "hierarchical",
             "composition", "effect size"),
    value = c(paste(CELLS, collapse = ", "),
              "only the four pre-specified contrasts; the two diagonal contrasts are not reported",
              "Kruskal-Wallis across the four cells, unpaired",
              "Dunn on the global ranks; for comparison the same contrast as Mann-Whitney",
              sprintf("%.3f", pi0),
              "log10 of the median cell count of the population (independent of the p value under H0)",
              "lineages first; sublineages only below a rejected lineage",
              sprintf("PERMANOVA on the Aitchison distance (sequential terms): %s p = %s, %s p = %s, interaction p = %s; dispersion p = %s",
                      OV, signif(PM$`Pr(>F)`[1], 3), IV, signif(PM$`Pr(>F)`[2], 3), signif(PM$`Pr(>F)`[3], 3), signif(pdisp, 3)),
              "Cliff's delta of group2 versus group1 (see comparison names), percentile bootstrap interval"))),
  file.path(OUT, "21_factorial_layer.xlsx"))
message("written: 21_factorial_layer.xlsx")
