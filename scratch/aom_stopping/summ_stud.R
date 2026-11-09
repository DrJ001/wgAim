# Compare the original simulated null with the studentised nulls
# (scratch/aom_stopping/results/stage1_stud/*.csv; same seeds as stage1/).
files <- list.files("scratch/aom_stopping/results/stage1_stud", "\\.csv$", full.names = TRUE)
rate <- function(p) { p <- p[!is.na(p)]; b <- binom.test(sum(p <= 0.05), length(p))$conf.int
    sprintf("%.3f (%d/%d; %.3f-%.3f)", mean(p <= 0.05), sum(p <= 0.05), length(p), b[1], b[2]) }
for (f in files) {
    R <- read.csv(f, stringsAsFactors = FALSE)
    U <- R[!grepl("^error", R$status), ]
    O <- U[U$status == "ok", ]
    cat(sprintf("\n==== %s  (usable %d, ok %d, gate-stop %d) ====\n", basename(f),
                nrow(U), nrow(O), sum(U$status == "gate-stop")))
    cat("FWER at 0.05   original sim :", rate(U$p.fwe), "\n")
    cat("               studentised  :", rate(U$p.fwe.stud), "\n")
    if (U$ntrait[1] > 1)
        cat("               stud. matrix :", rate(U$p.fwe.studM), "\n")
    cat("               Bonferroni   :", rate(U$p.bonf), "\n")
    for (v in c("p.fwe", "p.fwe.stud", if (U$ntrait[1] > 1) "p.fwe.studM")) {
        p <- O[[v]]
        cat(sprintf("  %-12s mean %.3f  P(p<=.05) %.3f  P(p<=.10) %.3f  quantiles 5/25/50/75/95: %s  KS p %.3f\n",
                    v, mean(p), mean(p <= .05), mean(p <= .10),
                    paste(sprintf("%.3f", quantile(p, c(.05, .25, .5, .75, .95))), collapse = " "),
                    suppressWarnings(ks.test(p, "punif")$p.value)))
    }
    cat(sprintf("  threshold: original median %.3f, studentised median %.3f (ratio %.3f)\n",
                median(O$fwe.thr), median(O$fwe.thr.stud), median(O$fwe.thr.stud / O$fwe.thr)))
    cat(sprintf("  observed sum(num)/sum(den): median %.5f  range %.5f-%.5f\n",
                median(O$obs.ratio), min(O$obs.ratio), max(O$obs.ratio)))
    cat(sprintf("  simulated null mean of oint: median %.4f\n", median(O$null.mean)))
}
