# Probe (1) model$formulae when random= is passed as a variable, and whether
# update() works after substituting the formula object into model$call;
# (2) naming / ordering of MV mbf random coefficients.
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)

cat("==== (1) formulae ====\n")
rf <- ~ us(Trial):id
fx <- yield ~ Trial + Trial:Rep
m <- asreml(fx, random = rf, residual = ~ dsum(~units | Trial),
            na.action = na.method(x = "include"), data = pheno, trace = FALSE)
cat("call$fixed class:", class(m$call$fixed), "  call$random class:", class(m$call$random), "\n")
cat("names(formulae):", names(m$formulae), "\n")
for (nm in names(m$formulae)) {
    f <- m$formulae[[nm]]
    cat(sprintf("  %-9s class=%s : %s\n", nm, paste(class(f), collapse = "/"),
                paste(deparse(f, width.cutoff = 500L), collapse = "")))
}
m2 <- m
m2$call$random <- m$formulae$random
m2$call$fixed  <- m$formulae$fixed
u <- tryCatch(update(m2, trace = FALSE), error = function(e) conditionMessage(e))
cat("update() after substitution:", if (is.character(u)) u else sprintf("OK logL=%.3f (orig %.3f)", u$loglik, m$loglik), "\n")
cat("deparse(call$random[[2]]):", deparse(m2$call$random[[2]]), "\n")
u2 <- tryCatch(update(m2, random. = ~ . + diag(Trial):Rep, trace = FALSE), error = function(e) conditionMessage(e))
cat("update(random. = ~ . + x):", if (is.character(u2)) u2 else deparse(u2$call$random), "\n")

cat("\n==== (2) MV mbf coefficients ====\n")
sub.cols <- round(seq(1, ncol(M), length.out = 150))
ints <- data.frame(id = geno.ids, M[, sub.cols], check.names = TRUE)
phg <- droplevels(pheno[pheno$id %in% geno.ids, ])
mm <- suppressWarnings(asreml(yield ~ Trial + Trial:Rep,
        random = ~ us(Trial):mbf('ints') + us(Trial):id,
        residual = ~ dsum(~units | Trial), mbf = list(ints = list(key = c("id", "id"), cov = "ints")),
        data = phg, workspace = "1gb", trace = FALSE, maxit = 30))
rn <- rownames(mm$coefficients$random)
mb <- grep("mbf", rn)
cat("n mbf coefs:", length(mb), " (ntrait x nmarkers =", 3 * length(sub.cols), ")\n")
print(head(rn[mb], 4)); print(rn[mb][c(149:152)]); print(tail(rn[mb], 2))
cat("vcoeff$random length:", length(mm$vcoeff$random), " coef length:", length(rn), "\n")
