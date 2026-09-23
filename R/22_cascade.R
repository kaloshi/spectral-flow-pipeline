## ============================================================================
## 22_cascade.R -- confirmatory layer: targeted instead of everything at once
## ============================================================================
## Usage:  Rscript R/22_cascade.R <folder>
## Input:  <folder>/20_abundance_long.csv.gz; design file with FACTORIAL
## Output: 22_cascade.xlsx
##
## The analysis has three layers, and only the middle one needs control of
## multiplicity:
##   primary       a pre-specified directional hypothesis (design file,
##                 ONE_SIDED; one-sided test in R/20_tests.R)
##   confirmatory  this cascade
##   descriptive   all sublineages, metrics and contrasts (R/21_factorial_layer.R),
##                 with effect size and interval, without a claim of significance
##
## The cascade, specified in advance:
##   stage 0  PERMANOVA on the Aitchison distance of the lineage composition:
##            does a factor shift the composition at all? One test per factor.
##            Only a factor with p < 0.05 is followed further; the others stay
##            descriptive.
##   stage 1  for that factor: the lineages one by one, metric CLR (the
##            compositional quantity), Mann-Whitney, Benjamini-Hochberg across
##            the lineages.
##   stage 2  only below a rejected lineage: its sublineages, metric % of the
##            parent, Benjamini-Hochberg across the children of this parent.
## This is hierarchical FDR testing (Yekutieli 2008): every family stays small,
## and testing descends only where the tree justifies it.
## Every test carries Cliff's delta with a bootstrap interval.
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(writexl) })
skriptdir <- dirname(gsub("~+~", " ", sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE))
source(file.path(skriptdir, "common.R"))
DESIGN <- load_design()
FAC <- DESIGN$factorial
if (is.null(FAC)) stop("FACTORIAL is not set in the design file")
args <- commandArgs(TRUE); OUT <- args[1]
B_BOOT <- as.integer(Sys.getenv("SFP_BOOT", "2000"))
ALPHA <- 0.05
OV <- FAC$outer$var; IV <- FAC$inner$var
CONTRAST <- setNames(list(FAC$outer$contrast %||% FAC$outer$levels, FAC$inner$contrast %||% FAC$inner$levels), c(OV, IV))

H <- fread(file.path(OUT, "20_abundance_long.csv.gz"), colClasses = c(sample_id = "character"))[included == TRUE]
H <- H[get(OV) %in% FAC$outer$levels & get(IV) %in% FAC$inner$levels & unit != "unassigned"]
H[, parent := fifelse(is.na(parent), "", parent)]

## ---- Stage 0: PERMANOVA --------------------------------------------------------------
CL <- dcast(H[metric == "clr" & level == "lineage"], as.formula(paste("file +", OV, "+", IV, "~ unit")), value.var = "value")
mat <- as.matrix(CL[, setdiff(names(CL), c("file", OV, IV)), with = FALSE])
ok <- complete.cases(mat)
set.seed(1)
pm <- vegan::adonis2(as.formula(paste("stats::dist(mat[ok, ]) ~", OV, "*", IV)),
                     data = as.data.frame(CL[ok, c(OV, IV), with = FALSE]), permutations = 9999, by = "terms")
TERMS <- c(OV, IV, paste0(OV, ":", IV))
S0 <- as.data.table(as.data.frame(pm), keep.rownames = "Term")[Term %in% TERMS]
setnames(S0, c("Term", "Df", "SumOfSqs", "R2", "F", "p"))
disp <- vegan::betadisper(stats::dist(mat[ok, ]), interaction(CL[[OV]][ok], CL[[IV]][ok]))
S0[, dispersion_p := vegan::permutest(disp, permutations = 999)$tab[1, "Pr(>F)"]]
licensed <- S0[Term %in% c(OV, IV) & p < ALPHA, Term]
message("Stage 0 -- PERMANOVA: ", paste(sprintf("%s p = %.4f (R2 %.1f %%)", S0$Term, S0$p, 100 * S0$R2), collapse = " | "))
message("followed further: ", if (length(licensed)) paste(licensed, collapse = ", ") else "nothing -- everything stays descriptive")

## ---- Stages 1 and 2 ----------------------------------------------------------------------
test_ab <- function(d, factor, a, b) {
  x <- d[get(factor) == a, value]; y <- d[get(factor) == b, value]
  if (length(x) < MIN_N || length(y) < MIN_N) return(NULL)
  ci <- cliffs_delta_ci(x, y, B = B_BOOT)
  data.table(n_a = length(x), n_b = length(y), median_a = median(x), median_b = median(y),
             delta = ci[["delta"]], ci_lower = ci[["lower"]], ci_upper = ci[["upper"]],
             p = suppressWarnings(wilcox.test(x, y)$p.value))
}
S1 <- list(); S2 <- list()
for (fk in licensed) {
  a <- CONTRAST[[fk]][1]; b <- CONTRAST[[fk]][2]
  L <- H[metric == "clr" & level == "lineage"]
  s1 <- rbindlist(lapply(sort(unique(L$unit)), function(u) {
    r <- test_ab(L[unit == u], fk, a, b); if (is.null(r)) NULL else cbind(data.table(factor = fk, unit = u), r) }))
  s1[, q_BH := p.adjust(p, "BH")][, rejected := q_BH < ALPHA]
  setorder(s1, p); S1[[fk]] <- s1
  message(sprintf("Stage 1 (%s, %d lineages, CLR): rejected %s", fk, nrow(s1),
                  if (any(s1$rejected)) paste(s1[rejected == TRUE, sprintf("%s (q %.3f)", unit, q_BH)], collapse = ", ") else "none"))
  for (mu in s1[rejected == TRUE, unit]) {
    K <- H[metric == "pct_parent" & level == "sublineage" & parent == mu]
    if (!nrow(K)) next
    s2 <- rbindlist(lapply(sort(unique(K$unit)), function(u) {
      r <- test_ab(K[unit == u], fk, a, b); if (is.null(r)) NULL else cbind(data.table(factor = fk, parent = mu, unit = u), r) }))
    s2[, q_BH := p.adjust(p, "BH")][, rejected := q_BH < ALPHA]
    setorder(s2, p); S2[[paste(fk, mu)]] <- s2
    message(sprintf("Stage 2 (%s, sublineages of %s, %% of parent): rejected %s", fk, mu,
                    if (any(s2$rejected)) paste(s2[rejected == TRUE, sprintf("%s (q %.3f)", unit, q_BH)], collapse = ", ") else "none"))
  }
}
S1 <- rbindlist(S1); S2 <- rbindlist(S2, fill = TRUE)

## ---- Counter-check: the same question without the cascade ------------------------------
## All populations in one family, as one would compute it without structure.
flat <- rbindlist(lapply(licensed, function(fk) {
  a <- CONTRAST[[fk]][1]; b <- CONTRAST[[fk]][2]
  A <- H[metric == "clr"]
  rbindlist(lapply(unique(A[, .(level, parent, unit)])[, seq_len(.N)], function(i) {
    k <- unique(A[, .(level, parent, unit)])[i]
    r <- test_ab(A[level == k$level & parent == k$parent & unit == k$unit], fk, a, b)
    if (is.null(r)) NULL else cbind(data.table(factor = fk), k, r) }))[, q_BH_flat := p.adjust(p, "BH"), by = factor][]
}))
if (nrow(flat)) message(sprintf("Without cascade (all %d populations in one family): %d rejected",
                                nrow(flat), sum(flat$q_BH_flat < ALPHA)))

## ---- Output ------------------------------------------------------------------------
os <- DESIGN$one_sided
decision <- data.table(
  layer = c("primary", "confirmatory stage 0", "confirmatory stage 1", "confirmatory stage 2", "descriptive"),
  what = c(if (length(os$markers)) sprintf("%s, one-sided in %s (R/20_tests.R)", paste(os$markers, collapse = "/"),
                                           paste(os$comparisons, collapse = "; ")) else "none specified",
           "PERMANOVA per factor on the lineage composition",
           "lineages one by one for the licensed factor, CLR, BH across the lineages",
           "sublineages only below rejected lineages, % of parent, BH per parent",
           "all sublineages, all metrics, all four contrasts (R/21_factorial_layer.R)"),
  result = c("see R/20_tests.R",
             paste(sprintf("%s p = %.4f", S0$Term, S0$p), collapse = "; "),
             if (nrow(S1)) paste(S1[rejected == TRUE, sprintf("%s q = %.3f, delta %.2f [%.2f, %.2f]", unit, q_BH, delta, ci_lower, ci_upper)], collapse = "; ") else "not licensed",
             if (nrow(S2)) { x <- S2[rejected == TRUE]; if (nrow(x)) paste(x[, sprintf("%s q = %.3f", unit, q_BH)], collapse = "; ") else
               paste("none rejected; closest:", S2[which.min(q_BH), sprintf("%s q = %.3f, delta %.2f [%.2f, %.2f]", unit, q_BH, delta, ci_lower, ci_upper)]) } else "not reached",
             "effect sizes with intervals, without a claim of significance"))
write_xlsx(list(Decision = as.data.frame(decision), Stage0_PERMANOVA = as.data.frame(S0),
                Stage1_lineages = as.data.frame(S1), Stage2_sublineages = as.data.frame(S2),
                Without_cascade = as.data.frame(flat)),
           file.path(OUT, "22_cascade.xlsx"))
message("written: 22_cascade.xlsx")
