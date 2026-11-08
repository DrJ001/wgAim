# Summarise the null simulations (scratch/aom_stopping/aom_null_k*.csv).
# Weak background: seeds 1005, 1022, 1023, 1050 never reached selection because
# the initial LRT statistic was <= 0 (p = 1 > TypeI = 0.999) -- gate stops,
# i.e. no false positive.  They are added here as non-rejections.
for (k in c("1", "0.05")) {
    R <- read.csv(sprintf("scratch/aom_stopping/results/aom_null_k%s.csv", k))
    R <- R[!duplicated(R$seed), ]
    n.gate <- if (k == "0.05") length(setdiff(c(1005, 1022, 1023, 1050), R$seed)) else 0L
    n  <- nrow(R) + n.gate
    pf <- c(R$p.fwe, rep(1, n.gate)); pb <- c(R$p.bonf, rep(1, n.gate))
    pm <- c(R$p.maha.bonf, rep(1, n.gate))
    gate <- c(R$p.lrt0 <= 0.05, rep(FALSE, n.gate))
    ci <- function(x) sprintf("%d/%d = %.3f  (95%% CI %.3f-%.3f)", sum(x), length(x),
                             mean(x), binom.test(sum(x), length(x))$conf.int[1],
                             binom.test(sum(x), length(x))$conf.int[2])
    cat(sprintf("\n==== Ga_K x %s : %d replicates (%d gate stops) ====\n", k, n, n.gate))
    cat("Reject at 0.05, AOM simulated max-stat :", ci(pf <= 0.05), "\n")
    cat("Reject at 0.05, AOM Bonferroni         :", ci(pb <= 0.05), "\n")
    cat("Reject at 0.05, Mahalanobis Bonferroni :", ci(pm <= 0.05), "\n")
    cat("Initial LRT significant (gate passed)  :", ci(gate), "\n")
    cat("Full procedure (gate AND AOM max-stat) :", ci(gate & pf <= 0.05), "\n")
    cat("p.fwe quantiles (excl. gate stops):\n")
    print(round(quantile(R$p.fwe, c(0, .05, .1, .25, .5, .75, .9, 1)), 3))
    cat(sprintf("mean p.fwe = %.3f (0.5 if uniform);  KS vs U(0,1): p = %.3f\n",
                mean(R$p.fwe), suppressWarnings(ks.test(R$p.fwe, "punif")$p.value)))
    cat(sprintf("null mean of oint (should be 1): %.4f;  diag reproduces package oint: %s\n",
                mean(R$null.mean), all(R$chk)))
    cat(sprintf("Meff median %.0f (range %.0f-%.0f) of 1501 markers\n",
                median(R$Meff), min(R$Meff), max(R$Meff)))
}
