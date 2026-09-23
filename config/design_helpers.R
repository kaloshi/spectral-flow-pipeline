## Helpers for design files (sourced before the design file).

## One comparison: `var` is a sample-sheet column, `a` and `b` are two of its
## values (b is compared against a; positive effect sizes mean b is higher).
## `strata` optionally restricts the comparison to samples with the given
## values of other columns, e.g. strata = c(other_column = "value").
comparison <- function(name, var, a, b, strata = NULL)
  list(name = name, var = var, str = strata, a = a, b = b)

## All comparisons of a two-by-two design: each factor across all samples and
## each factor within every level of the other factor (six comparisons).
comparisons_2x2 <- function(var1, first_a, first_b, var2, second_a, second_b, label = function(v, a, b) paste0(v, ": ", b, " vs ", a)) {
  s <- function(v, lev) setNames(lev, v)
  list(
    comparison(label(var1, first_a, first_b), var1, first_a, first_b),
    comparison(label(var2, second_a, second_b), var2, second_a, second_b),
    comparison(paste0(label(var2, second_a, second_b), " | ", var1, " = ", first_a), var2, second_a, second_b, s(var1, first_a)),
    comparison(paste0(label(var2, second_a, second_b), " | ", var1, " = ", first_b), var2, second_a, second_b, s(var1, first_b)),
    comparison(paste0(label(var1, first_a, first_b), " | ", var2, " = ", second_a), var1, first_a, first_b, s(var2, second_a)),
    comparison(paste0(label(var1, first_a, first_b), " | ", var2, " = ", second_b), var1, first_a, first_b, s(var2, second_b)))
}
