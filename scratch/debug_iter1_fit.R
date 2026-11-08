# Is the collapse of the T1 vm variance after adding the Chr8 QTL a property
# of the model, or of the warm-started update() path?  Fit the iteration-1
# model directly (fresh) and via update() from the iteration-0 fit.
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "4")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)

covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
ph <- pheno; whg <- levels(ph$id) %in% geno.ids
ph$Gomit <- ph$Gsave <- ph$id
levels(ph$Gsave)[!whg] <- NA; levels(ph$Gomit)[whg] <- "GEN"
x8 <- M[, qtl.cols[3]]                       # planted Chr8 marker (-1/0/1)
ph$X.Chr8.73 <- x8[match(as.character(ph$id), geno.ids)]
ph$X.Chr8.73[is.na(ph$X.Chr8.73)] <- 0

ctl <- list(na.action = na.method(x = "include"), workspace = "1gb", trace = FALSE, maxit = 60)
vmvar <- function(m) round(m$vparameters[grep("vm.*!Trial_T[0-9]$", names(m$vparameters))], 4)
rep_fit <- function(m, lab) {
    if (is.null(m)) return(invisible())
    cat(sprintf("%-40s converge=%s logL=%.3f\n   vm vars: %s\n", lab, m$converge, m$loglik,
                paste(vmvar(m), collapse = " ")))
}
f0 <- do.call(asreml, c(list(fixed = yield ~ Trial:Gomit + Trial + Trial:Rep,
          random = ~ corgh(Trial):vm(Gsave, covObj) + us(Trial):Gsave,
          residual = ~ dsum(~units | Trial), data = ph), ctl))
rep_fit(f0, "iter0 (no QTL)")
f1 <- do.call(asreml, c(list(fixed = yield ~ Trial:Gomit + Trial + Trial:Rep + X.Chr8.73 + Trial:X.Chr8.73,
          random = ~ corgh(Trial):vm(Gsave, covObj) + us(Trial):Gsave,
          residual = ~ dsum(~units | Trial), data = ph), ctl))
rep_fit(f1, "iter1 fresh")
u1 <- suppressWarnings(update(f0, fixed. = . ~ . + X.Chr8.73 + Trial:X.Chr8.73))
rep_fit(u1, "iter1 update(iter0)")
# alternative residual genetic structure
f1d <- do.call(asreml, c(list(fixed = yield ~ Trial:Gomit + Trial + Trial:Rep + X.Chr8.73 + Trial:X.Chr8.73,
          random = ~ corgh(Trial):vm(Gsave, covObj) + diag(Trial):Gsave,
          residual = ~ dsum(~units | Trial), data = ph), ctl))
rep_fit(f1d, "iter1 fresh, diag(Trial):Gsave")
cat("\nus(Trial):Gsave variances iter0 / iter1:\n")
g <- function(m) round(m$vparameters[grep("^Trial:Gsave!Trial_T([0-9]):T\\1$", names(m$vparameters))], 4)
print(rbind(iter0 = g(f0), iter1 = g(f1)))
