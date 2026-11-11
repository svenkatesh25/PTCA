# PTCA — Progressive Thresholding Correlation Analysis

PTCA measures how concordant two association studies are among their strongest shared signals.

1. Each feature (e.g. a gene) is meta-analyzed across the two studies (fixed- or random-effects, via [`remaCor`](https://cran.r-project.org/package=remaCor)).
2. Features are ranked by meta-analysis p-value.
3. The correlation between the two studies' Z-scores is computed over the top *k* features, for *k* = `stepSize`, 2·`stepSize`, … up to all features.
4. Correlation is plotted against *k* (raw and log10 scale).

A curve that stays high as *k* shrinks means the two studies agree most on their top hits.

## Requirements

R ≥ 4.1 with:

```r
install.packages(c("data.table", "remaCor", "ggplot2", "scales"))
```

## Usage

```r
source("PTCA.R")

# `work` has one row per feature with columns BETA1, BETA2, SE1, SE2, Z1, Z2
result <- PTCA(data = work, method = "BETA", rankMethod = "Fixed", outDir = "output")
```

Run the bundled simulated example:

```bash
Rscript example.R
```

### Key arguments

| Argument | Default | Description |
|---|---|---|
| `method` | `"BETA"` | `"BETA"`: meta-analyze effect sizes with SEs. `"Z"`: meta-analyze Z-scores (SE = 1). |
| `rankMethod` | `"Fixed"` | `"Fixed"` (`remaCor::LS`) or `"Random"` (`remaCor::RE2C`). |
| `stepSize` | `10` | Step increment, in features (or in % if `stepPercent = TRUE`). |
| `stepPercent` | `FALSE` | Treat steps as percentages of all features. |
| `corMethod` | `"pearson"` | `"pearson"`, `"spearman"`, or `"kendall"`. |
| `minPlottedStep` | `50` | Steps ≤ this are omitted from plots (same units as `stepSize`). |
| `restrictY` | `TRUE` | Fix y-axis to [0, 1]. Negative correlations are clipped. |
| `BETA1`…`Z2` | same names | Column names in `data`. `Z1`/`Z2` are always required. |
| `outDir` | `"output"` | Output directory. |

### Outputs (in `outDir`)

- `rankedMetaAnalysisTable.tsv`: input data plus `MA_beta`, `MA_se`, `MA_Z`, `MA_p`, sorted by `MA_p`
- `PTCATable.tsv`: correlation at every step (`step`, `nTop`, `corVal`, `logStep`)
- `PTCA_rawStep.pdf/.png`: correlation vs. number of top features
- `PTCA_logStep.pdf/.png`: correlation vs. log10 number of top features

### Notes

- Features with data in only one study are kept and use that study's own statistics.
- Correlations use pairwise-complete observations, so missing Z values are skipped.
- For `rankMethod = "Random"`, `MA_Z` is signed using the fixed-effect estimate (RE2C itself is unsigned), and `MA_se` is `NA`.
