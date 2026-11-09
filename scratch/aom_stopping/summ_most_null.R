# Summarise the standalone-MOST null grid (results/most_null/*.csv).
files <- list.files("scratch/aom_stopping/results/most_null", "\\.csv$", full.names = TRUE)
# cells rerun with run_most_null_v2.R (base fallback on non-convergence) replace v1
v2 <- list.files("scratch/aom_stopping/results/most_null_v2", "\\.csv$", full.names = TRUE)
files <- c(files[!basename(files) %in% basename(v2)], v2)
ci <- function(x) { x <- x[!is.na(x)]; b <- binom.test(sum(x), length(x))$conf.int
    sprintf("%.3f (%.3f-%.3f)", mean(x), b[1], b[2]) }
rows <- list()
for (f in files) {
    R <- read.csv(f, stringsAsFactors = FALSE); R <- R[!duplicated(R$seed), ]
    err <- grepl("^error", R$status); U <- R[!err, ]; O <- U[U$status == "ok", ]
    lab <- sprintf("%-4s t=%d  K x %-4s L x %s", U$pop[1], U$ntrait[1], U$kscale[1], U$lscale[1])
    cat("\n=====================================================================\n", lab,
        sprintf("  reps %d: ok %d, no-signal %d, errors %d | base us %d / diag %d | vm conv %.3f | trials bounded (ok reps): %s\n",
                nrow(R), nrow(O), sum(U$status == "no-signal"), sum(err),
                sum(U$base.str == "us", na.rm = TRUE), sum(U$base.str == "diag", na.rm = TRUE),
                mean(U$vm.conv, na.rm = TRUE),
                paste(names(table(O$n.bounded)), table(O$n.bounded), sep = ":", collapse = " ")), sep = "")
    if (any(err)) print(table(substr(R$status[err], 1, 70)))
    cat("  MOST (diag-studentised) reject @0.01/0.05/0.10:",
        sprintf("%.3f / %s / %.3f", mean(U$p.most <= .01), ci(U$p.most <= .05), mean(U$p.most <= .10)), "\n")
    cat("  comparison @0.05:  unstudentised", sprintf("%.3f", mean(U$p.orig <= .05)),
        "  matrix-stud", sprintf("%.3f", mean(U$p.studM <= .05)),
        "  Bonferroni", sprintf("%.3f", mean(U$p.bonf <= .05)),
        "  | legacy LRT reject", ci(U$p.lrt <= .05), "\n")
    if (nrow(O) > 5)
        cat(sprintf("  ok reps: MOST p mean %.3f, P(p<=.10) %.3f, KS p %.3f | pin.dev median %.5f max %.4f | Meff median %.0f\n",
                    mean(O$p.most), mean(O$p.most <= .10),
                    suppressWarnings(ks.test(O$p.most, "punif")$p.value),
                    median(O$diag.pin.dev), max(O$diag.pin.dev), median(O$Meff, na.rm = TRUE)))
    rows[[lab]] <- data.frame(cell = lab, n = nrow(U), no.signal = sum(U$status == "no-signal"),
        MOST.01 = mean(U$p.most <= .01), MOST.05 = mean(U$p.most <= .05),
        MOST.10 = mean(U$p.most <= .10), unstud.05 = mean(U$p.orig <= .05),
        bonf.05 = mean(U$p.bonf <= .05), LRT.05 = mean(U$p.lrt <= .05, na.rm = TRUE))
}
cat("\n\n==== STANDALONE MOST NULL GRID ====\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 3)
