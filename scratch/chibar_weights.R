# Sensitivity of the 3-trait boundary LRT critical value to the mixture weights.
# Binomial weights (current pchisq.mixture) assume the three tested variance
# estimators are uncorrelated. With correlation rho between them, the exact
# chi-bar-squared weights for k = 3 follow Kudo (1963):
#   w3 = P(all 3 > 0)              = 1/8 + sum(asin(rho_ij)) / (4 pi)
#   w0 = same with partial correlations (from the inverse correlation matrix)
#   w1 = 1/2 - w3,  w2 = 1/2 - w0
cb <- function(w, a = 0.95)
    uniroot(function(c) sum(w * pchisq(c, 0:3)) - a, c(1e-8, 30))$root
wb <- c(1, 3, 3, 1) / 8
cat(sprintf("binomial weights        : crit = %.3f\n", cb(wb)))
cat(sprintf("bound (1/2,1/2,0,0)     : crit = %.3f\n", cb(c(.5, .5, 0, 0))))
cat(sprintf("bound (0,0,1/2,1/2)     : crit = %.3f\n", cb(c(0, 0, .5, .5))))
for (r in c(-0.3, 0.3, 0.6, 0.9)) {
    R  <- matrix(r, 3, 3); diag(R) <- 1
    w3 <- 1/8 + 3 * asin(r) / (4 * pi)
    P  <- -cov2cor(solve(R))
    w0 <- 1/8 + 3 * asin(P[1, 2]) / (4 * pi)
    w  <- c(w0, 0.5 - w3, 0.5 - w0, w3)
    cat(sprintf("rho = %4.1f weights = %s  crit = %.3f  actual size at binomial crit = %.4f\n",
                r, paste(sprintf("%.3f", w), collapse = "/"), cb(w),
                1 - sum(w * pchisq(cb(wb), 0:3))))
}
cat("pchisq(0, 0) =", pchisq(0, 0), "  pchisq(-1, 0) =", pchisq(-1, 0), "\n")
