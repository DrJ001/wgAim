# Compare gpAim() vm vs mbf fits for us / fa(Trial,1) on the reduced panel,
# with and without ungenotyped lines (Gsave), and trace the fa mbf error.
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)
for (ch in names(panel$geno)) {
    g <- panel$geno[[ch]]; keep <- seq(1L, ncol(g$imputed.data), by = 6L)
    g$imputed.data <- g$imputed.data[, keep, drop = FALSE]
    if (!is.null(g$data)) g$data <- g$data[, keep, drop = FALSE]
    if (!is.null(g$map))  g$map  <- g$map[keep]
    panel$geno[[ch]] <- g
}
quiet <- function(expr) { out <- NULL; invisible(capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }
fx <- yield ~ Trial + Trial:Rep; rs <- ~ dsum(~units | Trial)
for (dat in c("all lines (Gsave)", "genotyped only")) {
    ph <- if (dat == "genotyped only") droplevels(pheno[pheno$id %in% geno.ids, ]) else pheno
    for (s in c("us", "fa(Trial,1)")) {
        rfs <- as.formula(if (grepl("\\(", s)) sprintf("~ %s:id", s) else sprintf("~ %s(Trial):id", s))
        b <- asreml(fx, random = rfs, residual = rs, na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE)
        for (frc in c(FALSE, TRUE)) {
            g <- tryCatch(quiet(gpAim(b, panel, "id", Trait = "Trial", trace = FALSE, force = frc, maxit = 50)),
                          error = function(e) { tb <- sys.calls(); paste("ERROR:", conditionMessage(e)) })
            lab <- sprintf("[%s] base %-11s %s", dat, s, if (frc) "mbf" else "vm ")
            if (is.character(g)) { cat(lab, g, "\n"); next }
            Gc <- g$GP$Gcor
            cat(sprintf("%s converge=%-5s logL=%9.3f  var.g=%s  Gcor=%s\n", lab, g$converge, g$loglik,
                        paste(round(g$GP$var.genetic, 4), collapse = "/"),
                        paste(round(Gc[lower.tri(Gc)], 3), collapse = "/")))
        }
    }
}
cat("\n--- traceback for fa(Trial,1) mbf with Gsave ---\n")
b <- asreml(fx, random = ~ fa(Trial, 1):id, residual = rs, na.action = na.method(x = "include"),
            data = pheno, workspace = "1gb", trace = FALSE)
withCallingHandlers(
    tryCatch(quiet(gpAim(b, panel, "id", Trait = "Trial", trace = FALSE, force = TRUE, maxit = 50)),
             error = function(e) cat("final error:", conditionMessage(e), "\n")),
    error = function(e) { cs <- sys.calls(); cat("call stack (innermost 8):\n");
        for (cl in tail(head(cs, -2), 8)) cat("  ", substr(paste(deparse(cl, width.cutoff = 200L)[1]), 1, 150), "\n") })
