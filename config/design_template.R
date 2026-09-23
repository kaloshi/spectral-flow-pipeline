## Study design for the statistics scripts (template).
## Copy this file outside the repository, replace the placeholders with column
## names and values of your sample sheet, and set SFP_DESIGN to the copy.
## Nothing in the pipeline depends on what the columns mean.

## Comparisons. Either list them one by one with comparison(), or generate all
## six comparisons of a two-by-two design with comparisons_2x2().
COMPARISONS <- comparisons_2x2("<column_1>", "<level_a>", "<level_b>",
                               "<column_2>", "<level_a>", "<level_b>")
# COMPARISONS <- list(comparison("<name>", var = "<column>", a = "<level_a>", b = "<level_b>"))

## Sample-sheet columns protected during batch correction (cyCombine `covar`).
## Their combination is passed as one covariate, so that run effects are
## removed without removing differences between these groups. Leave empty to
## correct without a covariate. The runs should be balanced across these groups.
PROTECTED_COVARIATES <- character(0)

## Numeric sample-sheet columns correlated (Spearman) with every readout.
CORRELATE_WITH <- character(0)

## Additional sample-level covariates for a secondary check of every nominal
## hit (R/25_extra_adjustment.R): Mann-Whitney on rank residuals after each
## covariate alone, after all of them, and after all of them plus run and log
## cell yield. Named vector: column name = "factor" or "numeric".
EXTRA_ADJUSTMENT <- character(0)
# EXTRA_ADJUSTMENT <- c("<column>" = "factor", "<column>" = "numeric")

## A directional hypothesis specified before the analysis: one-sided p values
## are reported only for these markers in these comparisons. "less" means
## that group a is expected to be lower than group b.
ONE_SIDED <- list(markers = character(0), comparisons = character(0), alternative = "less")

## Two-by-two layer (R/21_factorial_layer.R, optional): the four cells are
## formed by an outer and an inner factor. Tested are the inner factor within
## each outer level and the outer factor within each inner level (four
## contrasts); the two diagonal contrasts mix both factors and are not tested.
## `contrast` (optional) gives the direction used by the confirmatory cascade
## (R/22_cascade.R): level a versus level b; default is the order of `levels`.
## Leave NULL when the design has no two-by-two structure.
FACTORIAL <- NULL
# FACTORIAL <- list(outer = list(var = "<column_2>", levels = c("<level_1>", "<level_2>"), contrast = c("<a>", "<b>")),
#                   inner = list(var = "<column_1>", levels = c("<level_1>", "<level_2>")))
