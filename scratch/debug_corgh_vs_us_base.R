# corgh(Trial):id and us(Trial):id are the same covariance model, so a
# properly converged fit must give the same REML log-likelihood.
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
for (nt in c("3", "4")) {
    commandArgs <- function(...) c("corgh", "NULL", nt)
    invisible(capture.output(eval(parse(text = src[1:stop.at]))))
    rm(commandArgs)
    covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
    ph <- pheno; whg <- levels(ph$id) %in% geno.ids
    ph$Gomit <- ph$Gsave <- ph$id
    levels(ph$Gsave)[!whg] <- NA; levels(ph$Gomit)[whg] <- "GEN"
    cat(sprintf("\n=== %s trials ===\n", nt))
    for (rf in c("~ corgh(Trial):Gsave", "~ us(Trial):Gsave",
                 "~ diag(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave",
                 "~ diag(Trial):vm(Gsave, covObj) + us(Trial):Gsave")) {
        m <- tryCatch(suppressWarnings(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                 random = as.formula(rf), residual = ~ dsum(~units | Trial),
                 na.action = na.method(x = "include"), data = ph,
                 workspace = "1gb", trace = FALSE, maxit = 100)),
                 error = function(e) NULL)
        if (is.null(m)) { cat(sprintf("%-55s ABORT\n", rf)); next }
        m2 <- tryCatch(suppressWarnings(update(m)), error = function(e) m)
        cat(sprintf("%-55s converge=%-5s logL=%.3f  (after update: %.3f)\n",
                    rf, m$converge, m$loglik, m2$loglik))
    }
}
