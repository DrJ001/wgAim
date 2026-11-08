# Inspect .qtlSelect() internals on the MV mbf path (reduced panel).
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
ns <- asNamespace("wgAim")
trace(".qtlSelect", where = ns, print = FALSE, exit = quote({
    cat("\n--- .qtlSelect exit ---\n")
    cat("Ga:\n"); print(round(Ga, 6)); cat("Ginv diag:", round(diag(Ginv), 2), "\n")
    if (exists("atilde")) { cat("atilde dim:", dim(atilde), " sd by col:", round(apply(atilde, 2, sd), 6), "\n") }
    if (exists("pevar.v")) cat("pevar.v summary:", round(summary(pevar.v), 6), "\n")
    cat("qtilde summary:", signif(summary(qtilde), 4), "\n")
    cat("vqtilde summary:", signif(summary(vqtilde), 4), "  #>0:", sum(vqtilde > 0), "\n")
    cat("oint summary:", signif(summary(oint[as.logical(oint != 0) | TRUE]), 4), "\n")
}))
b <- asreml(yield ~ Trial + Trial:Rep, random = ~ us(Trial):id, residual = ~ dsum(~units | Trial),
            na.action = na.method(x = "include"), data = pheno, workspace = "1gb", trace = FALSE)
g <- suppressWarnings(suppressMessages(
    gwasAim(b, panel, "id", Trait = "Trial", trace = FALSE, force = TRUE, breakout = 1, maxit = 50)))
