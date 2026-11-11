# Minimal example: run PTCA on simulated, partially concordant data.
# Run from the repository root:  Rscript example.R

source("PTCA.R")
set.seed(1)

n <- 1000
trueEffect <- c(rnorm(100, sd = 0.15), rep(0, n - 100))  # 100 shared true signals
SE1 <- runif(n, 0.01, 0.2)
SE2 <- runif(n, 0.01, 0.2)
BETA1 <- trueEffect + rnorm(n, sd = SE1)
BETA2 <- trueEffect + rnorm(n, sd = SE2)

work <- data.table(BETA1, BETA2, SE1, SE2)
work[, Z1 := BETA1 / SE1]
work[, Z2 := BETA2 / SE2]

result <- PTCA(data = work, method = "BETA", rankMethod = "Fixed", outDir = "example_output")
print(head(result))
