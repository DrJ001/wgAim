# Can ASReml predict line-level values of an MV mbf term (classify = Trial:Gmbf,
# only = mbf G.param key) with a full vcov?  If so, check that
#   q = (M'M)^-1 M' u   reproduces the mbf marker BLUPs exactly.
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)
keep <- unlist(lapply(panel$geno, function(g) colnames(g$imputed.data)[seq(1L, ncol(g$imputed.data), by = 6L)]))
Ms <- M[, keep]
ph <- pheno; whg <- levels(ph$id) %in% geno.ids
ph$Gomit <- ph$Gsave <- ph$id
levels(ph$Gsave)[!whg] <- NA; levels(ph$Gomit)[whg] <- "GEN"
ph$Gmbf <- factor(ifelse(is.na(ph$Gsave), ".NoGeno", as.character(ph$Gsave)))
covObj <- wgAim:::.mbfCovObj(Ms, "Gmbf")
m <- suppressWarnings(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
        random = ~ us(Trial):mbf('ints') + us(Trial):Gsave,
        mbf = list(ints = list(key = c("Gmbf", "Gmbf"), cov = "covObj")),
        residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
        data = ph, workspace = "1gb", trace = FALSE, maxit = 50))
cat("converge:", m$converge, " logL:", m$loglik, "\n")
key <- names(m$G.param)[grep("mbf", names(m$G.param))]
cat("G.param key:", key, "\n")
for (cl in c("Trial:Gmbf", "Trial:Gsave")) {
    pv <- tryCatch(suppressWarnings(predict(m, classify = cl, only = key, vcov = TRUE, maxit = 1)),
                   error = function(e) conditionMessage(e))
    if (is.character(pv)) { cat(cl, "PREDICT ERROR:", pv, "\n"); next }
    cat(sprintf("%s: pvals %d rows, vcov %s, mean pred %.5f\n", cl, nrow(pv$pvals),
                paste(dim(pv$vcov), collapse = "x"), mean(pv$pvals$predicted.value, na.rm = TRUE)))
    # Reconstruct marker effects from line values for trial T1
    lc <- names(pv$pvals)[2L]
    p1 <- pv$pvals[pv$pvals$Trial == "T1" & as.character(pv$pvals[[lc]]) %in% geno.ids, ]
    u  <- p1$predicted.value[match(geno.ids, as.character(p1[[lc]]))]
    q.rec <- MASS::ginv(Ms) %*% u
    rn <- rownames(m$coefficients$random)
    er <- wgAim:::.mbfEffectRows(rn, "Trial", levels(ph$Trial), colnames(Ms))
    q.blup <- m$coefficients$random[er[, "T1"], 1L]
    cat(sprintf("   T1: max|q.rec - q.blup| = %.3g   cor = %.6f\n", max(abs(q.rec - q.blup)), cor(q.rec, q.blup)))
}
