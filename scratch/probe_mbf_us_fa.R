# (1) Do us vm and us mbf fits reach the same optimum?
# (2) Does fa(Trial,1):mbf fail fresh, or only via update() from diag?
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)

keep <- unlist(lapply(panel$geno, function(g) colnames(g$imputed.data)[seq(1L, ncol(g$imputed.data), by = 6L)]))
Ms   <- M[, keep]
cat("markers:", ncol(Ms), "\n")
phg <- droplevels(pheno[pheno$id %in% geno.ids, ])
# vm covObj exactly as .constructCM()
relm <- tcrossprod(Ms); sc <- mean(diag(relm)); relm <- relm / sc
ev <- eigen(relm, symmetric = TRUE); ev$values[ev$values < 0] <- 0
relm <- ev$vectors %*% diag(ev$values + max(1e-6, 1e-10 * max(ev$values))) %*% t(ev$vectors)
covObj <- relm; dimnames(covObj) <- list(geno.ids, geno.ids)
ints <- data.frame(id = geno.ids, Ms, check.names = FALSE)
mbfl <- list(ints = list(key = c("id", "id"), cov = "ints"))
ctl  <- list(residual = ~ dsum(~units | Trial), workspace = "1gb", trace = FALSE, maxit = 100)
fx   <- yield ~ Trial + Trial:Rep
fit  <- function(rf, mbf = FALSE) tryCatch(suppressWarnings(do.call(asreml,
            c(list(fixed = fx, random = rf, data = phg), if (mbf) list(mbf = mbfl), ctl))),
            error = function(e) conditionMessage(e))
rep  <- function(m, lab) {
    if (is.character(m)) { cat(sprintf("%-42s ERROR: %s\n", lab, gsub("\\s+", " ", m))); return() }
    vp <- m$vparameters[grep("vm.*covObj|mbf", names(m$vparameters))]
    sc.f <- if (grepl("mbf", lab)) sc else 1    # mbf per-marker var x scale = vm var
    Ga <- wgAim:::.extractVmGa(m, "Trial", levels(phg$Trial), 1) * sc.f
    d  <- 1 / sqrt(diag(Ga)); R <- Ga * tcrossprod(d)
    cat(sprintf("%-42s converge=%-5s logL=%9.3f  diag(Ga)=%s  cor=%s\n", lab, m$converge, m$loglik,
                paste(round(diag(Ga), 4), collapse = "/"), paste(round(R[lower.tri(R)], 3), collapse = "/")))
}
cat("\n(1) us: vm vs mbf (Ga of mbf rescaled to the vm scale)\n")
rep(fit(~ diag(Trial):vm(id, covObj) + us(Trial):id), "diag:vm + us(id)")
rep(fit(~ us(Trial):vm(id, covObj) + us(Trial):id), "us:vm + us(id)")
rep(fit(~ diag(Trial):mbf('ints') + us(Trial):id, TRUE), "diag:mbf + us(id)")
m.us.mbf <- fit(~ us(Trial):mbf('ints') + us(Trial):id, TRUE)
rep(m.us.mbf, "us:mbf + us(id)")
if (!is.character(m.us.mbf)) rep(tryCatch(suppressWarnings(update(m.us.mbf)), error = function(e) conditionMessage(e)),
                                 "us:mbf + us(id)  [after extra update]")

cat("\n(2) fa(Trial,1) on mbf\n")
d.mbf <- fit(~ diag(Trial):mbf('ints') + fa(Trial, 1):id, TRUE)
rep(d.mbf, "diag:mbf + fa1(id)")
rep(fit(~ fa(Trial, 1):mbf('ints') + fa(Trial, 1):id, TRUE), "fa1:mbf + fa1(id)  [fresh]")
if (!is.character(d.mbf)) {
    rep(tryCatch(suppressWarnings(update(d.mbf, random. = ~ fa(Trial, 1):mbf('ints') + fa(Trial, 1):id)),
                 error = function(e) conditionMessage(e)), "fa1:mbf  [update(diag)]")
    rep(tryCatch(suppressWarnings(update(d.mbf, random. = ~ fa(Trial, 1):mbf('ints') + fa(Trial, 1):id,
                                         G.param = NULL, R.param = NULL)),
                 error = function(e) conditionMessage(e)), "fa1:mbf  [update(diag), G.param=NULL]")
}
cat("\n(2b) fa(Trial,1) on vm for reference\n")
rep(fit(~ fa(Trial, 1):vm(id, covObj) + fa(Trial, 1):id), "fa1:vm + fa1(id)  [fresh]")
