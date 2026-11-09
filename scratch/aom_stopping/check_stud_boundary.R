# Does miscalibration in the weak MV cells concentrate in replicates where the
# REML pin fails (a per-trial genomic variance at the boundary, so
# sum num != sum den and S.obs != S.target)?
dirs  <- file.path("scratch/aom_stopping/results",
                   c("stage1_stud", "stage1_stud_gA", "stage1_stud_gB"))
files <- unlist(lapply(dirs, list.files, pattern = "_t3_.*\\.csv$", full.names = TRUE))
for (cl in split(files, basename(files))) {
    R <- do.call(rbind, lapply(cl, read.csv, stringsAsFactors = FALSE))
    O <- R[R$status == "ok", ]
    O$pinned <- abs(O$obs.ratio - 1) < 0.005 & O$obs.Sdev < 0.005
    cat(sprintf("\n==== %s : ok %d, pin holds %d, pin fails %d ====\n",
                sub("\\.csv$", "", basename(cl[1])), nrow(O), sum(O$pinned), sum(!O$pinned)))
    for (g in c(TRUE, FALSE)) {
        S <- O[O$pinned == g, ]
        if (!nrow(S)) next
        cat(sprintf("  pin %-5s n=%3d  rate@0.05: orig %.3f  scalar %.3f  matrix %.3f   rate@0.10: orig %.3f  scalar %.3f  matrix %.3f\n",
                    if (g) "holds" else "FAILS", nrow(S),
                    mean(S$p.fwe <= .05), mean(S$p.fwe.stud <= .05), mean(S$p.fwe.studM <= .05),
                    mean(S$p.fwe <= .10), mean(S$p.fwe.stud <= .10), mean(S$p.fwe.studM <= .10)))
    }
    if (any(!O$pinned))
        cat(sprintf("  pin fails: obs.ratio median %.3f (range %.3f-%.3f)\n",
                    median(O$obs.ratio[!O$pinned]), min(O$obs.ratio[!O$pinned]),
                    max(O$obs.ratio[!O$pinned])))
}
