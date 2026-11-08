# (a) fa(Trial,1):mbf coefficient names
# (b) us:mbf / fa1:mbf FRESH fits on data with ungenotyped lines (Gsave NA),
#     vs the vm equivalents -- path / start-value issue or model difference?
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)
keep <- unlist(lapply(panel$geno, function(g) colnames(g$imputed.data)[seq(1L, ncol(g$imputed.data), by = 6L)]))
Ms <- M[, keep]
relm <- tcrossprod(Ms); sc <- mean(diag(relm)); relm <- relm / sc
ev <- eigen(relm, symmetric = TRUE); ev$values[ev$values < 0] <- 0
relm <- ev$vectors %*% diag(ev$values + max(1e-6, 1e-10 * max(ev$values))) %*% t(ev$vectors)

ph <- pheno; whg <- levels(ph$id) %in% geno.ids
ph$Gomit <- ph$Gsave <- ph$id
levels(ph$Gsave)[!whg] <- NA; levels(ph$Gomit)[whg] <- "GEN"
covObj <- relm; dimnames(covObj) <- list(geno.ids, geno.ids)
ints <- data.frame(Gsave = geno.ids, Ms, check.names = FALSE)
mbfl <- list(ints = list(key = c("Gsave", "Gsave"), cov = "ints"))
ctl  <- list(residual = ~ dsum(~units | Trial), workspace = "1gb", trace = FALSE, maxit = 100,
             na.action = na.method(x = "include"))
fx   <- yield ~ Trial:Gomit + Trial + Trial:Rep
fit  <- function(rf, mbf = FALSE, d = ph) tryCatch(suppressWarnings(do.call(asreml,
            c(list(fixed = fx, random = rf, data = d), if (mbf) list(mbf = mbfl), ctl))),
            error = function(e) conditionMessage(e))
rep  <- function(m, lab) {
    if (is.character(m)) { cat(sprintf("%-40s ERROR: %s\n", lab, gsub("\\s+", " ", m))); return(invisible()) }
    cat(sprintf("%-40s converge=%-5s logL=%9.3f\n", lab, m$converge, m$loglik)); invisible(m)
}
cat("(b) with ungenotyped lines (Gsave NA)\n")
rep(fit(~ us(Trial):vm(Gsave, covObj) + us(Trial):Gsave), "us:vm + us(Gsave)")
m.us <- rep(fit(~ us(Trial):mbf('ints') + us(Trial):Gsave, TRUE), "us:mbf + us(Gsave) [fresh]")
rep(fit(~ diag(Trial):vm(Gsave, covObj) + us(Trial):Gsave), "diag:vm + us(Gsave)")
rep(fit(~ diag(Trial):mbf('ints') + us(Trial):Gsave, TRUE), "diag:mbf + us(Gsave)")
rep(fit(~ fa(Trial, 1):vm(Gsave, covObj) + fa(Trial, 1):Gsave), "fa1:vm + fa1(Gsave)")
rep(fit(~ fa(Trial, 1):mbf('ints') + fa(Trial, 1):Gsave, TRUE), "fa1:mbf + fa1(Gsave) [fresh]")
cat("\nmbf random coefficient count (us, Gsave data):",
    if (!is.character(m.us)) length(grep("mbf", rownames(m.us$coefficients$random))) else NA, "\n")

cat("\n(a) fa(Trial,1):mbf coefficient names (genotyped only)\n")
phg <- droplevels(ph[ph$id %in% geno.ids, ]); phg$Gsave <- droplevels(phg$Gsave)
fx <- yield ~ Trial + Trial:Rep
mf <- fit(~ fa(Trial, 1):mbf('ints') + fa(Trial, 1):Gsave, TRUE, phg)
if (!is.character(mf)) {
    rn <- rownames(mf$coefficients$random); mb <- rn[grep("mbf", rn)]
    cat("n =", length(mb), " (3 trials x", ncol(Ms), "markers =", 3 * ncol(Ms), ")\n")
    pre <- sub(":mbf.*$", "", mb); print(table(pre))
    print(head(mb, 2))
} else cat(mf, "\n")
