# How does ASReml name vm terms for each variance structure, and which
# string does predict(only=) honour?
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
phg <- droplevels(pheno[pheno$id %in% geno.ids, ])

for (s in c("diag(Trial)", "corgh(Trial)", "us(Trial)", "fa(Trial,1)")) {
    rf <- as.formula(sprintf("~ %s:vm(id, covObj) + diag(Trial):id", s))
    m <- tryCatch(suppressWarnings(asreml(yield ~ Trial + Trial:Rep, random = rf,
                    residual = ~ dsum(~units | Trial), data = phg,
                    workspace = "1gb", trace = FALSE, maxit = 30)),
                  error = function(e) { cat(s, "fit ERROR\n"); NULL })
    if (is.null(m)) next
    raw   <- attr(terms.formula(rf), "term.labels")[1]
    gkey  <- names(m$G.param)[grep("vm.*covObj", names(m$G.param))]
    cat(sprintf("\n[%s] raw='%s'  G.param key='%s'\n", s, raw, gkey))
    for (o in unique(c(raw, gkey))) {
        pv <- tryCatch(suppressWarnings(predict(m, classify = "Trial:id", only = o,
                              vcov = TRUE, maxit = 1)),
                       error = function(e) { cat("   only='", o, "' ERROR: ", conditionMessage(e), "\n", sep = ""); NULL })
        if (!is.null(pv))
            cat(sprintf("   only='%s': mean pred=%.4f  sd=%.4f  vcov NULL=%s\n", o,
                        mean(pv$pvals$predicted.value, na.rm = TRUE),
                        sd(pv$pvals$predicted.value, na.rm = TRUE), is.null(pv$vcov)))
    }
}
