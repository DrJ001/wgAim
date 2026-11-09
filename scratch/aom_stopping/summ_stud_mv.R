# Scalar vs matrix studentisation for MV cells (3 trials, us).
# Combines results/stage1_stud, stage1_stud_gA, stage1_stud_gB by cell.
dirs  <- file.path("scratch/aom_stopping/results",
                   c("stage1_stud", "stage1_stud_gA", "stage1_stud_gB"))
files <- unlist(lapply(dirs, list.files, pattern = "_t3_.*\\.csv$", full.names = TRUE))
cells <- split(files, basename(files))
alphas <- c(.01, .05, .10, .20)
methods <- c(original = "p.fwe", scalar = "p.fwe.stud", matrix = "p.fwe.studM",
             Bonferroni = "p.bonf")
ci <- function(x) { b <- binom.test(sum(x), length(x))$conf.int
    sprintf("%.3f (%.3f-%.3f)", mean(x), b[1], b[2]) }
pool <- list()
for (cl in names(cells)) {
    R <- do.call(rbind, lapply(cells[[cl]], read.csv, stringsAsFactors = FALSE))
    R <- R[!duplicated(R$seed), ]
    err <- grepl("^error", R$status)
    U <- R[!err, ]; O <- U[U$status == "ok", ]
    cat(sprintf("\n==== %s : %d replicates (ok %d, gate-stop %d, errors %d) ====\n",
                sub("\\.csv$", "", cl), nrow(R), nrow(O), sum(U$status == "gate-stop"), sum(err)))
    tab <- sapply(methods, function(v) sapply(alphas, function(a) mean(U[[v]] <= a)))
    rownames(tab) <- paste0("alpha=", alphas)
    print(round(tab, 3))
    cat("at 0.05 with 95% CI:\n")
    for (nm in names(methods)) cat(sprintf("  %-10s %s\n", nm, ci(U[[methods[nm]]] <= .05)))
    cat("KS uniformity (ok reps):",
        paste(sprintf("%s p=%.3f", names(methods)[1:3],
              sapply(methods[1:3], function(v) suppressWarnings(ks.test(O[[v]], "punif")$p.value))),
              collapse = "  "), "\n")
    pool[[cl]] <- U
}
P <- do.call(rbind, pool)
cat(sprintf("\n==== POOLED over %d MV cells: %d replicates ====\n", length(pool), nrow(P)))
tab <- sapply(methods, function(v) sapply(alphas, function(a) mean(P[[v]] <= a)))
rownames(tab) <- paste0("alpha=", alphas); print(round(tab, 3))
for (nm in names(methods)) cat(sprintf("  %-10s at 0.05: %s\n", nm, ci(P[[methods[nm]]] <= .05)))
