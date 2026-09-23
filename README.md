# spectral-flow-pipeline

Analysis pipeline for spectral flow cytometry data from two panels, an NK/T-cell panel and a myeloid/B-cell
(MYE/LYM) panel. It covers transformation, the positivity thresholds, a correction of a spectral unmixing
artefact, batch correction, clustering, filtering of debris, population annotation and the statistical
comparison of abundance and marker expression between groups of samples.

> **Status: work in progress.** This repository currently contains the shared module, the panel
> definitions and the statistics layer. The upstream steps (gating import, transformation, unmixing
> correction, thresholds, clustering, debris filter, annotation), the Python modules and the full method
> documentation will follow.

The repository contains code only. It contains no data and no results, and no information about the samples
analysed. Grouping variables and covariates are defined by the user in a local design file (see below).

## Requirements

- R **>= 4.6.0**. From R 4.6.0 on, `wilcox.test()` computes exact p values also in the presence of ties;
  older versions use the normal approximation and give different results.
- R packages: data.table, readxl, writexl, rstatix, vegan, qvalue, IHW (statistics); flowCore, CytoML,
  flowWorkspace, FlowSOM, ConsensusClusterPlus, cyCombine, harmony, bluster, igraph, uwot (upstream steps).

## Layout

| Path | Content |
|---|---|
| `R/common.R` | shared definitions: panel, sample sheet, threshold estimator, statistical helpers, autofluorescence rules |
| `R/20_tests.R` | group comparisons of abundance and expression (Mann-Whitney, Cliff's delta, run and yield adjustment, correlations) |
| `R/21_factorial_layer.R` | two-by-two layer: Kruskal-Wallis, Dunn, several multiplicity procedures, PERMANOVA |
| `R/22_cascade.R` | confirmatory hierarchical cascade (PERMANOVA, then lineages, then sublineages) |
| `R/23_robustness.R` | every expression contrast staged by cell cut-off and sample statistic |
| `R/24_cutoff_consequences.R` | what a higher cell cut-off costs in samples and testable contrasts |
| `R/25_extra_adjustment.R` | re-test of nominal hits after additional sample covariates |
| `R/26_cell_count_bias.R` | dependence of nominal hits on cell numbers and denominators |
| `config/panel_ntnk.R`, `config/panel_mye.R` | panel definitions (channel, marker, role) |
| `config/design_template.R` | template for the study design (comparisons, covariates) |
| `tools/check_public.py` | pre-publication check (no data files, identifiers, paths or study information) |

## Configuration

Set before running any step:

- `SFP_PANEL` = `ntnk` or `mye`
- `SFP_SAMPLES` = sample sheet (CSV, one row per FCS file, columns `file`, `sample_id`, `run` plus your own
  grouping and covariate columns)
- `SFP_DESIGN` = a copy of `config/design_template.R` filled in with your column names and levels

Keep the sample sheet and the design file outside the repository.

## Methods and references (statistics layer)

| Method | Used in | Reference |
|---|---|---|
| Mann-Whitney / Wilcoxon rank-sum test | 20-26 | Wilcoxon 1945, doi:10.2307/3001968; Mann & Whitney 1947, doi:10.1214/aoms/1177730491 |
| exact distribution with ties (R >= 4.6.0) | common.R | Streitberg & Röhmel 1987, EDV in Medizin und Biologie 18:12-19 |
| Cliff's delta | 20-26 | Cliff 1993, doi:10.1037/0033-2909.114.3.494 |
| percentile bootstrap | 21, 22 | Efron & Tibshirani 1993, An Introduction to the Bootstrap, Chapman & Hall |
| rank analysis of covariance (tests on rank residuals) | 20, 25, 26 | Quade 1967, doi:10.1080/01621459.1967.10500925 |
| Spearman rank correlation | 20, 25, 26 | Spearman 1904, doi:10.2307/1412159 |
| centred log-ratio, Aitchison distance | 20-22 | Aitchison 1982, doi:10.1111/j.2517-6161.1982.tb01195.x; Gloor et al. 2017, doi:10.3389/fmicb.2017.02224 |
| Kruskal-Wallis test | 21, 25 | Kruskal & Wallis 1952, doi:10.1080/01621459.1952.10483441 |
| Dunn's post-hoc test | 21 | Dunn 1964, doi:10.1080/00401706.1964.10490181 |
| Bonferroni correction | 21 | Dunn 1961, doi:10.1080/01621459.1961.10482090 |
| Benjamini-Hochberg FDR | 21, 22 | Benjamini & Hochberg 1995, doi:10.1111/j.2517-6161.1995.tb02031.x |
| Storey q value | 21 | Storey 2002, doi:10.1111/1467-9868.00346; Storey & Tibshirani 2003, doi:10.1073/pnas.1530509100 |
| independent hypothesis weighting | 21 | Ignatiadis et al. 2016, doi:10.1038/nmeth.3885 |
| hierarchical FDR | 21, 22 | Yekutieli 2008, doi:10.1198/016214507000001373 |
| PERMANOVA | 21, 22 | Anderson 2001, doi:10.1111/j.1442-9993.2001.01070.pp.x |
| homogeneity of dispersion (betadisper) | 21, 22 | Anderson 2006, doi:10.1111/j.1541-0420.2005.00440.x |
| per-sample aggregation, marker classes | 20 | Weber et al. 2019, doi:10.1038/s42003-019-0415-5; Nowicka et al., F1000Research 2017;6:748 (version 3), doi:10.12688/f1000research.11622.3 |
| negative control population | 20 | Lipsitch et al. 2010, doi:10.1097/EDE.0b013e3181d61eeb |
| detection limit instead of post-hoc power | (upcoming) | Hoenig & Heisey 2001, doi:10.1198/000313001300339897 |

Every reference above was checked against two independent sources. References for the upstream steps
(transformation, unmixing, clustering, annotation) follow with those steps.

## Third-party software

The pipeline calls the following packages; none of their code is included in this repository, and each
remains under its own license. Please cite them when you use the pipeline (`citation("<package>")` in R).

| Package | License |
|---|---|
| R (R Core Team) | GPL-2 / GPL-3 |
| data.table | MPL-2.0 |
| readxl | MIT |
| writexl | BSD-2-Clause |
| rstatix | GPL-2 |
| vegan | GPL-2 |
| qvalue | LGPL |
| IHW | Artistic-2.0 |

## Authorship

Dr. Benjamin Krämer. The code was developed with the assistance of Claude Opus 5 (Anthropic) in Visual Studio Code.

## License

MIT, see `LICENSE`.
