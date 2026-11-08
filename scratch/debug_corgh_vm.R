# Diagnose convergence of diag(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })

src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^gw <- gwasAim", src) - 1L
invisible(capture.output(eval(parse(text = src[1:stop.at]))))

covObj <- K
dimnames(covObj) <- list(geno.ids, geno.ids)
ph <- pheno
whg <- levels(ph$id) %in% geno.ids
ph$Gomit <- ph$Gsave <- ph$id
levels(ph$Gsave)[!whg] <- NA
levels(ph$Gomit)[whg]  <- "GEN"

rpt <- function(m, lab) {
    if (is.null(m)) return(invisible())
    cat(sprintf("%-50s converge=%s logL=%.3f\n", lab, m$converge, m$loglik))
}
tryfit <- function(expr, lab) {
    m <- tryCatch(suppressWarnings(expr),
                  error = function(e) { cat(sprintf("%-50s ERROR: %s\n", lab, conditionMessage(e))); NULL })
    rpt(m, lab); invisible(m)
}

b0 <- tryfit(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep, random = ~ corgh(Trial):Gsave,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE), "base corgh(Gsave)")
f1 <- tryfit(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                    random = ~ diag(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE, maxit = 50),
             "fresh diag:vm + corgh(Gsave)")
u1 <- tryfit(update(b0, random. = ~ diag(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave, maxit = 50),
             "update(base) -> diag:vm + corgh(Gsave)")
f2 <- tryfit(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                    random = ~ diag(Trial):vm(Gsave, covObj) + diag(Trial):Gsave,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE, maxit = 50),
             "fresh diag:vm + diag(Gsave)")
f3 <- tryfit(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                    random = ~ corgh(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE, maxit = 50),
             "fresh corgh:vm + corgh(Gsave)")
phg <- droplevels(pheno[pheno$id %in% geno.ids, ])
g1 <- tryfit(asreml(yield ~ Trial + Trial:Rep,
                    random = ~ diag(Trial):vm(id, covObj) + corgh(Trial):id,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = phg, workspace = "1gb", trace = FALSE, maxit = 50),
             "genotyped-only diag:vm + corgh(id)")
g2 <- tryfit(asreml(yield ~ Trial + Trial:Rep,
                    random = ~ diag(Trial):vm(id, covObj) + fa(Trial,1):id,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = phg, workspace = "1gb", trace = FALSE, maxit = 50),
             "genotyped-only diag:vm + fa1(id)")
g3 <- tryfit(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                    random = ~ diag(Trial):vm(Gsave, covObj) + fa(Trial,1):Gsave,
                    residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                    data = ph, workspace = "1gb", trace = FALSE, maxit = 50),
             "Gsave diag:vm + fa1(Gsave)")
if (!is.null(f1)) print(summary(f1)$varcomp)
if (!is.null(u1)) print(summary(u1)$varcomp)
