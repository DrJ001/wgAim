source("scratch/aom_stopping/aom_null_diag.R")
set.seed(1)
for (L in list(c(.5, .3, .2), c(.9, .05, .05), c(1/3, 1/3, 1/3)))
    for (x in c(2, 5, 10, 15)) {
        s <- colSums(L * matrix(rnorm(3 * 4e5), 3)^2)
        cat(sprintf("lam=%-16s x=%4.1f imhof=%.6f  mc=%.6f\n",
                    paste(round(L, 3), collapse = ","), x, .imhof_upper(x, L), mean(s > x)))
    }
cat("1-df check:", .imhof_upper(4, 1), pchisq(4, 1, lower.tail = FALSE), "\n")
cat("3 equal = chi2_3/3 check:", .imhof_upper(5, rep(1/3, 3)), pchisq(15, 3, lower.tail = FALSE), "\n")
