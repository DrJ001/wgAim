# vparameter names/order for each vm structure (3 trials)
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
phg <- droplevels(pheno[pheno$id %in% geno.ids, ])
for (s in c("diag(Trial)", "corh(Trial)", "corgh(Trial)", "us(Trial)", "fa(Trial,1)", "fa(Trial,2)")) {
    rf <- as.formula(sprintf("~ %s:vm(id, covObj) + diag(Trial):id", s))
    m <- tryCatch(suppressWarnings(asreml(yield ~ Trial + Trial:Rep, random = rf,
                    residual = ~ dsum(~units | Trial), data = phg,
                    workspace = "1gb", trace = FALSE, maxit = 20)),
                  error = function(e) { cat(s, "fit ERROR\n"); NULL })
    if (is.null(m)) next
    vp <- m$vparameters; sel <- grep("vm", names(vp))
    cat("\n", s, "  G.param key:", names(m$G.param)[1], "\n")
    print(data.frame(name = names(vp)[sel], value = round(vp[sel], 4),
                     con = m$vparameters.con[sel], row.names = NULL))
    cat("sigma2 =", m$sigma2, "\n")
}
