# =============================================================================
# scratch/probe_us_genomic.R
#
# Probe: is us(Trait) supported on the genomic term in ASReml-R 4.2?
# (engine_model.R comment claims "us() is not supported as a variance
#  structure on vm()/mbf() terms" and mirrors a us base model to corgh.)
#
# For vm (markers > lines) and mbf (lines >= markers) paths, 3 and 4 trials:
#   - fresh fit of us(Trial):<genomic> + diag(Trial):id
#   - update() from diag(Trial):<genomic> fit (wgAim's upgrade route)
#   - same with corgh for comparison
#   - G.param key, vparameter names, predict(only = key) behaviour
#   - .extractVmGa() result (PD?), logLik equality us vs corgh
# Run: Rscript scratch/probe_us_genomic.R
# =============================================================================
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
cat("asreml version:", as.character(packageVersion("asreml")), "\n")

src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L

ctl <- list(residual = ~ dsum(~units | Trial), workspace = "1gb", trace = FALSE, maxit = 60)
tryfit <- function(expr) tryCatch(suppressWarnings(expr),
                                  error = function(e) structure(list(msg = conditionMessage(e)), class = "fitfail"))
ok <- function(m) !inherits(m, "fitfail")

probe <- function(m, lab, classify, tl) {
    if (!ok(m)) { cat(sprintf("  %-28s FIT FAILED: %s\n", lab, gsub("\\s+", " ", m$msg))); return(invisible(NULL)) }
    key <- names(m$G.param)[grep("vm.*covObj|mbf.*ints", names(m$G.param))]
    vpn <- names(m$vparameters)[grep("vm.*covObj|mbf.*ints", names(m$vparameters))]
    pv  <- tryfit(predict(m, classify = classify, only = key, vcov = TRUE, maxit = 1))
    pr  <- if (ok(pv)) sprintf("mean=%.4f sd=%.4f vcov=%s", mean(pv$pvals$predicted.value, na.rm = TRUE),
                               sd(pv$pvals$predicted.value, na.rm = TRUE),
                               if (is.null(pv$vcov)) "NULL" else paste(dim(pv$vcov), collapse = "x"))
           else paste("PREDICT FAILED:", pv$msg)
    Ga  <- tryfit(wgAim:::.extractVmGa(m, "Trial", tl, 1))
    ga  <- if (ok(Ga)) sprintf("min eig=%.4f", min(eigen(Ga, only.values = TRUE)$values))
           else paste("Ga FAILED:", Ga$msg)
    cat(sprintf("  %-28s converge=%-5s logL=%9.3f  key='%s'\n", lab, m$converge, m$loglik, key))
    cat(sprintf("  %-28s vparams: %s%s\n", "", paste(head(sub("^[^!]*", "", vpn), 4), collapse = " "),
                if (length(vpn) > 4) sprintf(" ... (%d)", length(vpn)) else ""))
    cat(sprintf("  %-28s predict(only=key): %s | Ga %s\n", "", pr, ga))
    invisible(m)
}

for (nt in c("3", "4")) {
    commandArgs <- function(...) c("corgh", "NULL", nt)
    invisible(capture.output(eval(parse(text = src[1:stop.at]))))
    rm(commandArgs)
    tl  <- levels(pheno$Trial)
    phg <- droplevels(pheno[pheno$id %in% geno.ids, ])

    # ---------------- vm path ----------------
    covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
    cat(sprintf("\n================ %s trials : vm(id, covObj) ================\n", nt))
    fx <- yield ~ Trial + Trial:Rep
    d0 <- tryfit(do.call(asreml, c(list(fixed = fx, random = ~ diag(Trial):vm(id, covObj) + diag(Trial):id, data = phg), ctl)))
    probe(d0, "diag (start)", "Trial:id", tl)
    for (s in c("us", "corgh")) {
        rf <- as.formula(sprintf("~ %s(Trial):vm(id, covObj) + diag(Trial):id", s))
        probe(tryfit(do.call(asreml, c(list(fixed = fx, random = rf, data = phg), ctl))),
              paste(s, "fresh"), "Trial:id", tl)
        if (ok(d0)) probe(tryfit(update(d0, random. = rf)), paste(s, "update(diag)"), "Trial:id", tl)
    }

    # ---------------- mbf path ----------------
    # lines >= markers: use a 150-marker subset so wgAim would choose mbf
    sub.cols <- round(seq(1, ncol(M), length.out = 150))
    ints <- data.frame(id = geno.ids, M[, sub.cols], check.names = TRUE)
    cat(sprintf("\n================ %s trials : mbf('ints') (%d markers) ================\n", nt, length(sub.cols)))
    mbfl <- list(ints = list(key = c("id", "id"), cov = "ints"))
    d1 <- tryfit(do.call(asreml, c(list(fixed = fx, random = ~ diag(Trial):mbf('ints') + diag(Trial):id,
                                        mbf = mbfl, data = phg), ctl)))
    probe(d1, "diag (start)", "Trial:id", tl)
    for (s in c("us", "corgh")) {
        rf <- as.formula(sprintf("~ %s(Trial):mbf('ints') + diag(Trial):id", s))
        probe(tryfit(do.call(asreml, c(list(fixed = fx, random = rf, mbf = mbfl, data = phg), ctl))),
              paste(s, "fresh"), "Trial:id", tl)
        if (ok(d1)) probe(tryfit(update(d1, random. = rf)), paste(s, "update(diag)"), "Trial:id", tl)
    }
}
