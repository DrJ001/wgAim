# Summarise Stage 1 null calibration (scratch/stage1/*.csv).
files <- list.files("C:/jtaylor/Claude/wgAim/scratch/stage1", "\\.csv$", full.names = TRUE)
alpha <- 0.05
ci <- function(x) {
    x <- x[!is.na(x)]; k <- sum(x); n <- length(x)
    if (!n) return("  -  ")
    b <- binom.test(k, n)$conf.int
    sprintf("%5.3f (%3d/%3d; %.3f-%.3f)", k / n, k, n, b[1], b[2])
}
rows <- list()
for (f in files) {
    R  <- read.csv(f, stringsAsFactors = FALSE)
    R  <- R[!duplicated(R$seed), ]
    err  <- grepl("^error", R$status)
    gate <- R$status == "gate-stop"
    ok   <- R$status == "ok"
    U    <- R[!err, ]                                  # usable replicates
    lab  <- sprintf("%s  t=%d  Ga_K x %s", U$pop[1], U$ntrait[1], U$kscale[1])
    cat("\n==========================================================================\n")
    cat(lab, sprintf("   replicates %d  (ok %d, gate-stop %d, errors %d)\n",
                     nrow(R), sum(ok), sum(gate), sum(err)))
    if (any(err)) print(table(substr(R$status[err], 1, 70)))
    cat("--------------------------------------------------------------------------\n")
    cat("AOM rule, false positives at 0.05 (gate stops = no rejection):\n")
    cat("  simulated max-stat    :", ci(U$p.fwe <= alpha), "\n")
    cat("  Bonferroni            :", ci(U$p.bonf <= alpha), "\n")
    cat("  standardised S_j Bonf :", ci(U$p.maha.bonf <= alpha), "\n")
    if (U$pop[1] == "dh")
        cat("  chromosome statistic  :", ci(U$p.fwe.chr <= alpha), "\n")
    cat("Initial LRT significant :", ci(U$p.lrt0 <= alpha), "\n")
    cat("Full procedure (LRT gate at 0.05 AND AOM sim):",
        ci(U$p.lrt0 <= alpha & U$p.fwe <= alpha), "\n")
    O <- U[U$status == "ok", ]
    if (nrow(O) > 5) {
        q <- quantile(O$p.fwe, c(.05, .1, .25, .5, .75, .9))
        cat(sprintf("p.fwe (ok reps): mean %.3f  quantiles 5/10/25/50/75/90%%: %s\n",
                    mean(O$p.fwe), paste(sprintf("%.3f", q), collapse = " ")))
        cat(sprintf("  P(p.fwe <= 0.10) = %.3f   KS vs U(0,1): p = %.3f\n",
                    mean(O$p.fwe <= 0.10), suppressWarnings(ks.test(O$p.fwe, "punif")$p.value)))
        cat(sprintf("null mean of oint: median %.4f (range %.4f-%.4f);  max rel diff vs package: %.2g\n",
                    median(O$null.mean), min(O$null.mean), max(O$null.mean), max(O$chk)))
        cat(sprintf("Meff median %.0f (IQR %.0f-%.0f);  secs/rep median %.1f\n",
                    median(O$Meff, na.rm = TRUE), quantile(O$Meff, .25, na.rm = TRUE),
                    quantile(O$Meff, .75, na.rm = TRUE), median(R$secs)))
    }
    rows[[lab]] <- data.frame(cell = lab, n = nrow(U), gate = sum(gate),
        fwer.sim = mean(U$p.fwe <= alpha), fwer.bonf = mean(U$p.bonf <= alpha),
        fwer.maha = mean(U$p.maha.bonf <= alpha),
        fwer.chr = if (U$pop[1] == "dh") mean(U$p.fwe.chr <= alpha) else NA,
        lrt.sig = mean(U$p.lrt0 <= alpha),
        full = mean(U$p.lrt0 <= alpha & U$p.fwe <= alpha),
        mean.p = if (nrow(O)) mean(O$p.fwe) else NA,
        Meff = if (nrow(O)) median(O$Meff, na.rm = TRUE) else NA)
}
cat("\n\n==== STAGE 1 OVERVIEW ====\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 3)
