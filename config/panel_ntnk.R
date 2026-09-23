## NK/T panel: 33 channels, spectral cytometer.
## channel = FCS parameter name ($PnN), marker = canonical marker name,
## role = "lineage" (clustering), "state" (activation/function) or "qc".
PANEL_ID <- "ntnk"
PANEL <- data.frame(stringsAsFactors = FALSE, matrix(ncol = 3, byrow = TRUE, dimnames = list(NULL, c("channel", "marker", "role")), c(
  "[LD UV455 ]-A",      "Viability",  "qc",
  "BUV395-A",           "CD45",       "lineage",
  "BUV563-A",           "CD45RA",     "lineage",
  "SBUV665-A",          "CD25",       "lineage",
  "SBUV740-A",          "CD4",        "lineage",
  "BV421-A",            "CD69",       "state",
  "VioBlue-A",          "NKp80",      "lineage",
  "BV510-A",            "TCRgd",      "lineage",
  "BV570-A",            "HLA-DR",     "state",
  "BV605-A",            "CD56",       "lineage",
  "BV650-A",            "c-kit",      "lineage",
  "SB780-A",            "CD278",      "state",
  "[RV828]-A",          "CD8a",       "lineage",
  "FITC-A",             "Lineage",    "lineage",
  "RB705-A",            "CXCR5",      "lineage",
  "[RB670]-A",          "NKp44",      "lineage",
  "PerCP-A",            "GranzymeB",  "state",
  "BB700-A",            "FOXP3",      "lineage",
  "PerCP-eFluor710-A",  "CD49a",      "lineage",
  "RB780-A",            "TIM-3",      "state",
  "PerCP-Fire806-A",    "CD45RO",     "lineage",
  "PE-A",               "AndrogenR",  "state",
  "SparkYG-581-A",      "CD127",      "lineage",
  "PE-Dazzle594-A",     "TIGIT",      "state",
  "AF594-A",            "EstrogenR",  "state",
  "[RY655]-A",          "CRTH2",      "lineage",
  "[RY703]-A",          "CD94",       "lineage",
  "PE-Cy7-A",           "PD-1",       "state",
  "PE-Fire810-A",       "CD39",       "state",
  "[RR688]-A",          "CD223",      "state",
  "APC-A",              "CD38",       "state",
  "R718-A",             "CD103",      "lineage",
  "APC-Cy7-A",          "CD3",        "lineage")))
## The lineage dump channel ("Lineage", FITC) contains CD14, CD19, CD20, CD34,
## CD123, FceRIa, BDCA-2 (CD303) and CD66b; it contains no CD3, CD56 or CD16.

## Markers without signal above noise (see R/05_thresholds.R): NKp80 spans only
## about two spreads above its negative mode and has no positive population. In
## a percentile-scaled clustering space, noise would get the same weight as CD3.
DISCARDED <- c("NKp80")
