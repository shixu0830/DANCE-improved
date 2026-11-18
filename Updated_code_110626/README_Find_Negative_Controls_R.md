# DANCE negative control triplet discovery — R-only

`FindNegativeControls()` identifies candidate negative control triplets using marginal or baseline-conditional correlations and tetrad tests. This file documents the implementation that runs entirely in R.

## Choose one implementation

Both implementations are called from **R** and expose the same `FindNegativeControls()` arguments and output format. Choose according to your computing environment:

| Implementation | Requirements | R source file |
|---|---|---|
| **R-only** | R and `MASS`; no C++ compiler | `Find_Negative_Controls_R.R` |
| **Rcpp (R + C++)** | R, `MASS`, `Rcpp`, `RcppParallel`, a C++11-compatible compiler, and `tetrad_helpers.cpp` | `find_negative_control_conditional_110526.R` |

**Source only one implementation in an R session.** Both define `FindNegativeControls()`; sourcing the other replaces that function. To switch, start a fresh session and source the desired file. Knowledge of C++ is not required to use the Rcpp implementation.

## Setup (R-only)

```r
# Install once if needed
install.packages("MASS")

# Update this path to your local R-only file
source("path/to/Find_Negative_Controls_R.R")
```

The screening helpers are included in this file. No separate C++ file or compilation is needed. This implementation retains the supplied R helpers and the latest main-function safeguards and optional diagnostics. Interpreted triplet screening can take longer with large candidate pools.

## Input and examples

`data` must be a data frame with numeric treatment, outcome, and at least three numeric candidate NC columns. Conditional DANCE also requires numeric baseline X columns. In the examples, `df` is your data, `Z_names` lists candidate names, and `X_names` lists baseline names.

```r
# Completely marginal DANCE: supply only A, Y, and intended NC candidates
res_marginal <- FindNegativeControls(
  data = df[, c("A", "Y", Z_names), drop = FALSE],
  tname = "A", oname = "Y",
  baseline = character(0), conditional = FALSE
)

# Completely conditional DANCE: include baseline X
res_conditional <- FindNegativeControls(
  data = df[, unique(c("A", "Y", Z_names, X_names)), drop = FALSE],
  tname = "A", oname = "Y",
  baseline = X_names, conditional = TRUE,
  return_diagnostics = TRUE
)

attr(res_conditional, "status")
attr(res_conditional, "excluded_rare")
attr(res_conditional, "excluded_zero_variance")
```

## Main arguments (same in both implementations)

| Argument | Default | Description |
|---|---|---|
| `baseline` | `character(0)` | Names of baseline X variables |
| `conditional` | `TRUE` | Condition on baseline X; supply `baseline` or explicitly choose marginal DANCE |
| `vtalpha` | `0.05` | All six tetrad p-values must exceed this threshold; also used for optional pair screening |
| `independence` | `TRUE` | Screen candidates for association with both treatment and outcome |
| `alpha_corr` | `0.2` | Maximum association-screen p-value |
| `exclude_rare` | `TRUE` | Exclude binary 0/1 candidates with either level below 1% prevalence |
| `use_pair` | `FALSE` | Enable pairwise pre-screening |
| `FindBest` | `FALSE` | Return only the highest-ranked triplet when `TRUE` |
| `return_diagnostics` | `FALSE` | Attach status and candidate-exclusion attributes |

## Output and notes

- Returns a **3 × K character matrix** with one triplet per column, ranked by descending minimum tetrad p-value. Empty results have dimension 3 × 0.
- Optional `status` is `success`, `no_triplets`, or `calculation_failed`. A calculation failure is not evidence of no discoveries.
- Correlations use complete cases across treatment, outcome, current candidates, and any conditioning variables. Candidates with nonpositive variance are removed and complete cases are recalculated.
- **Simulation-specific exclusion:** candidate names beginning with `U` are excluded. For other applications, review this rule; set `Uvars <- character(0)` in the chosen R file to disable it.
- For marginal DANCE, extra numeric columns may be treated as candidates. Supply only the intended inputs.
- Both implementations follow the same screening rules, but use different numerical routines for determinants. Exact numerical equality has not been established; small differences near thresholds or ties may affect selection or ranking.
