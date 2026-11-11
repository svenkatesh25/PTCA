#' PTCA: Progressive Thresholding Correlation Analysis
#'
#' Meta-analyzes two sets of per-feature (e.g. per-gene) association statistics,
#' ranks features by meta-analysis p-value, and then computes the correlation
#' between the two input Z-scores across progressively larger sets of the
#' top-ranked features. The resulting curve shows how concordant the two
#' studies are among their strongest shared signals.
#'
#' Outputs written to `outDir`:
#'   - rankedMetaAnalysisTable.tsv : input data + meta-analysis results, ranked
#'   - PTCATable.tsv               : correlation at every step (unfiltered)
#'   - PTCA_rawStep.pdf/.png       : correlation vs. number of top features
#'   - PTCA_logStep.pdf/.png       : correlation vs. log10(number of top features)
#'
#' @param data          data.frame / data.table with one row per feature.
#' @param method        "BETA" (meta-analyze BETA1/BETA2 with SE1/SE2) or
#'                      "Z" (meta-analyze Z1/Z2 with unit standard errors).
#'                      Z scores may be pre-scaled by the user if desired.
#' @param rankMethod    "Fixed" (fixed-effect, remaCor::LS) or
#'                      "Random" (random-effects, remaCor::RE2C).
#' @param stepSize      Step increment. Number of features if stepPercent = FALSE,
#'                      percentage of features if stepPercent = TRUE.
#' @param corMethod     Passed to cor(method = ...): "pearson", "spearman", "kendall".
#' @param minPlottedStep Steps <= this value are dropped from the plots (small
#'                      top sets give noisy correlations). Same units as stepSize.
#' @param restrictY     If TRUE, fix the y axis to [0, 1]. Note: negative
#'                      correlations will then be clipped from the plot.
#' @param stepPercent   If TRUE, steps are percentages of the total feature count.
#' @param BETA1,BETA2,SE1,SE2,Z1,Z2 Column names in `data`. Z1/Z2 are always
#'                      required (they are what gets correlated); BETA/SE columns
#'                      are required only when method = "BETA".
#' @param outDir        Output directory (created if needed).
#' @param plotHeight,plotWidth,plotDPI Passed to ggsave (units = "px").
#'
#' @return data.table of the plotted steps (step, nTop, corVal, logStep).
#'         All steps, including those below minPlottedStep, are in PTCATable.tsv.

library(data.table)
library(remaCor)
library(ggplot2)
library(scales)

PTCA <- function(data,
		method = c("BETA", "Z"), rankMethod = c("Fixed", "Random"),
		stepSize = 10, corMethod = "pearson", minPlottedStep = 50,
		restrictY = TRUE, stepPercent = FALSE,
		BETA1 = "BETA1", BETA2 = "BETA2", SE1 = "SE1", SE2 = "SE2", Z1 = "Z1", Z2 = "Z2",
		outDir = "output", plotHeight = 3360, plotWidth = 3360, plotDPI = 600
	) {

	### 1) Input checks / preprocessing
	method <- match.arg(method)
	rankMethod <- match.arg(rankMethod)
	corMethod <- match.arg(corMethod, c("pearson", "spearman", "kendall"))

	data <- as.data.table(data)
	requiredCols <- c(Z1, Z2)
	if (method == "BETA") requiredCols <- c(requiredCols, BETA1, BETA2, SE1, SE2)
	missingCols <- setdiff(requiredCols, colnames(data))
	if (length(missingCols) > 0) {
		stop("Columns not found in data: ", paste(missingCols, collapse = ", "))
	}
	if (stepPercent && (stepSize <= 0 || stepSize > 100)) {
		stop("With stepPercent = TRUE, stepSize must be in (0, 100].")
	}
	if (stepSize <= 0) stop("stepSize must be positive.")

	dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
	corLabel <- paste0(toupper(substr(corMethod, 1, 1)), substr(corMethod, 2, nchar(corMethod)))

	### 2) Per-feature meta-analysis
	n <- nrow(data)
	MA_p <- rep(NA_real_, n)
	MA_Z <- rep(NA_real_, n)
	MA_beta <- rep(NA_real_, n)
	MA_se <- rep(NA_real_, n)

	for (a in seq_len(n)) {
		if (method == "BETA") {
			BETAV <- c(data[[BETA1]][a], data[[BETA2]][a])
			SEV <- c(data[[SE1]][a], data[[SE2]][a])
		} else {
			BETAV <- c(data[[Z1]][a], data[[Z2]][a])
			SEV <- c(1, 1)
		}

		# Drop studies with a missing estimate OR SE (keeps BETA/SE pairs aligned)
		keep <- !is.na(BETAV) & !is.na(SEV)
		BETAV <- BETAV[keep]
		SEV <- SEV[keep]

		if (length(BETAV) == 0) {
			next
		} else if (length(BETAV) == 1) {
			# Only one study available: use its own statistics
			MA_beta[a] <- BETAV
			MA_se[a] <- SEV
			MA_Z[a] <- BETAV / SEV
			MA_p[a] <- 2 * pnorm(-abs(MA_Z[a]))
		} else if (rankMethod == "Fixed") {
			meta <- LS(BETAV, SEV)
			MA_beta[a] <- meta$beta
			MA_se[a] <- meta$se
			MA_Z[a] <- meta$beta / meta$se
			MA_p[a] <- meta$p
		} else {
			# RE2Cp is the random-effects association p-value
			# (QEp is Cochran's Q heterogeneity p-value and must not be used for ranking).
			# RE2C is not a signed test; the sign of the fixed-effect estimate is used for MA_Z.
			meta <- RE2C(BETAV, SEV)
			fe <- LS(BETAV, SEV)
			MA_beta[a] <- fe$beta
			MA_p[a] <- meta$RE2Cp
			MA_Z[a] <- sign(fe$beta) * qnorm(meta$RE2Cp / 2, lower.tail = FALSE)
		}
	}
	data[, `:=`(MA_beta = MA_beta, MA_se = MA_se, MA_Z = MA_Z, MA_p = MA_p)]

	# Rank by p-value; break ties (e.g. p underflowing to 0) by |Z|. NAs go last.
	data <- data[order(MA_p, -abs(MA_Z), na.last = TRUE)]
	fwrite(data, file = file.path(outDir, "rankedMetaAnalysisTable.tsv"), sep = "\t", quote = FALSE)

	### 3) Correlation across progressively larger top-ranked sets
	if (stepPercent) {
		steps <- seq(stepSize, 100, stepSize)
		nTop <- pmax(1, round(n * steps / 100))
	} else {
		steps <- seq(stepSize, n, stepSize)
		nTop <- steps
	}
	if (length(steps) == 0) stop("stepSize is larger than the number of rows in data.")

	corVal <- vapply(nTop, function(k) {
		suppressWarnings(cor(data[[Z1]][seq_len(k)], data[[Z2]][seq_len(k)],
			method = corMethod, use = "pairwise.complete.obs"))
	}, numeric(1))

	stepTable <- data.table(step = steps, nTop = nTop, corVal = corVal, logStep = log10(steps))
	fwrite(stepTable, file = file.path(outDir, "PTCATable.tsv"), sep = "\t", quote = FALSE)

	### 4) Plots
	plotTable <- stepTable[step > minPlottedStep]
	if (nrow(plotTable) == 0) {
		warning("No steps above minPlottedStep; skipping plots.")
		return(stepTable)
	}

	yScale <- if (restrictY) {
		scale_y_continuous(expand = c(0, 0), limits = c(0, 1))
	} else {
		scale_y_continuous(expand = c(0, 0))
	}
	xUnit <- if (stepPercent) "% Top Features" else "N(Top Features)"

	savePlot <- function(graph, stem) {
		for (dev in c("pdf", "png")) {
			ggsave(filename = file.path(outDir, paste0(stem, ".", dev)), plot = graph, device = dev,
				width = plotWidth, height = plotHeight, units = "px", dpi = plotDPI)
		}
	}

	rawPlot <- ggplot(plotTable, aes(x = step, y = corVal)) +
		geom_line() +
		theme_classic() +
		ylab(paste0(corLabel, " Correlation")) +
		yScale +
		scale_x_reverse(labels = comma) +
		xlab(xUnit)
	savePlot(rawPlot, "PTCA_rawStep")

	logPlot <- ggplot(plotTable, aes(x = logStep, y = corVal)) +
		geom_line() +
		theme_classic() +
		ylab(paste0(corLabel, " Correlation")) +
		yScale +
		scale_x_reverse() +
		xlab(paste0("log10 ", xUnit))
	savePlot(logPlot, "PTCA_logStep")

	return(plotTable)
}
