# =============================================================================
# scratch/probe_fixes_symbol_mbf.R
# (A) base model fitted with formulae held in variables (no update()) now
#     works in gwasAim() and gpAim().
# (B) gpAim() MV mbf path returns per-trial GEBVs, Ga/Gcor, and its GEBVs
#     agree with the vm path (same G-BLUP model in two parameterisations).
# Run: Rscript scratch/probe_fixes_symbol_mbf.R
# =============================================================================
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L
commandArgs <- function(...) c("us", "NULL", "3")
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
rm(commandArgs)
quiet <- function(expr) { out <- NULL; invisible(capture.output(out <- suppressWarnings(suppressMessages(expr)))); out }

cat("==== (A) formulae passed as variables, no update() ====\n")
fx <- yield ~ Trial + Trial:Rep
rf <- ~ us(Trial):id
rs <- ~ dsum(~units | Trial)
base.sym <- asreml(fx, random = rf, residual = rs, na.action = na.method(x = "include"),
                   data = pheno, workspace = "1gb", trace = FALSE)
cat("call$fixed:", class(base.sym$call$fixed), " call$random:", class(base.sym$call$random), "\n")
gw <- tryCatch(quiet(gwasAim(base.sym, panel, "id", Trait = "Trial", trace = FALSE,
                             breakout = 3, maxit = 50)),
               error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(gw)) cat("gwasAim:", gw, "\n") else {
    cat("gwasAim OK. QTL:", paste(gw$QTL$qtl, collapse = " "), "\n")
    print(summary(gw, genObj = panel)[, 1:6])
}
gp <- tryCatch(quiet(gpAim(base.sym, panel, "id", Trait = "Trial", trace = FALSE, maxit = 50)),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("gpAim:", if (is.character(gp)) gp else
    paste("OK, random =", paste(deparse(gp$call$random, width.cutoff = 500L), collapse = "")), "\n")

cat("\n==== (B) MV mbf path vs vm path ====\n")
# Reduced panel (~250 markers, every 6th): still > 191 lines so the default
# path is vm, while force = TRUE gives mbf on the same markers -- keeps the
# mbf fits fast.
panel.small <- panel
for (ch in names(panel.small$geno)) {
    g <- panel.small$geno[[ch]]
    keep <- seq(1L, ncol(g$imputed.data), by = 6L)
    g$imputed.data <- g$imputed.data[, keep, drop = FALSE]
    if (!is.null(g$data)) g$data <- g$data[, keep, drop = FALSE]
    if (!is.null(g$map))  g$map  <- g$map[keep]
    panel.small$geno[[ch]] <- g
}
cat("reduced panel markers:", sum(vapply(panel.small$geno, function(g) ncol(g$imputed.data), 1L)),
    "  lines:", nrow(panel.small$pheno), "\n")
panel <- panel.small
for (s in c("diag", "us", "fa(Trial,1)")) {
    rfs <- as.formula(if (grepl("\\(", s)) sprintf("~ %s:id", s) else sprintf("~ %s(Trial):id", s))
    b <- asreml(fx, random = rfs, residual = rs, na.action = na.method(x = "include"),
                data = pheno, workspace = "1gb", trace = FALSE)
    g.vm  <- tryCatch(quiet(gpAim(b, panel, "id", Trait = "Trial", trace = FALSE, maxit = 50)),
                      error = function(e) paste("ERROR:", conditionMessage(e)))
    g.mbf <- tryCatch(quiet(gpAim(b, panel, "id", Trait = "Trial", trace = FALSE,
                                  force = TRUE, maxit = 50)),
                      error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("\n-- base %s --\n", s))
    if (is.character(g.mbf)) { cat("mbf:", g.mbf, "\n"); next }
    if (is.character(g.vm))  { cat("vm:", g.vm, "\n"); next }
    G1 <- g.vm$GP; G2 <- g.mbf$GP
    cat("mbf random:", paste(deparse(g.mbf$call$random, width.cutoff = 500L), collapse = ""), "\n")
    cat("path vm/mbf:", G1$path, "/", G2$path, "  gebv rows:", nrow(G1$gebv), "/", nrow(G2$gebv),
        "  cols:", paste(names(G2$gebv), collapse = ","), "\n")
    k  <- paste(G1$gebv[[1]], G1$gebv$Trial); k2 <- paste(G2$gebv[[1]], G2$gebv$Trial)
    m  <- match(k, k2)
    cat("per-trial cor(GEBV vm, GEBV mbf):",
        paste(sprintf("%s=%.4f", G1$trait.levels, vapply(G1$trait.levels, function(t) {
            i <- G1$gebv$Trial == t; cor(G1$gebv$GEBV[i], G2$gebv$GEBV[m[i]]) }, 0)), collapse = " "), "\n")
    cat("Gcor vm (lower tri): ", paste(round(G1$Gcor[lower.tri(G1$Gcor)], 3), collapse = " "), "\n")
    cat("Gcor mbf (lower tri):", paste(round(G2$Gcor[lower.tri(G2$Gcor)], 3), collapse = " "), "\n")
    cat("var.genetic vm: ", paste(round(G1$var.genetic, 4), collapse = " "), "\n")
    cat("var.genetic mbf:", paste(round(G2$var.genetic, 4), collapse = " "), "\n")
    cat("mean accuracy vm/mbf:", round(mean(G1$gebv$Accuracy), 3), "/", round(mean(G2$gebv$Accuracy), 3), "\n")
    pr <- tryCatch({ capture.output(print(g.mbf)); capture.output(summary(g.mbf)); "OK" },
                   error = function(e) conditionMessage(e))
    pl <- tryCatch({ p <- plot(g.mbf); "OK" }, error = function(e) conditionMessage(e))
    si <- tryCatch({ x <- selIndex(g.mbf, weights = setNames(rep(1, 3), G2$trait.levels)); "OK" },
                   error = function(e) conditionMessage(e))
    cat("print/summary:", pr, "  plot:", pl, "  selIndex:", si, "\n")
}

cat("\n==== (C) gwasAim MV on the mbf path (force = TRUE), ungenotyped lines present ====\n")
for (s in c("us", "fa(Trial,1)")) {
    rfs <- as.formula(if (grepl("\\(", s)) sprintf("~ %s:id", s) else sprintf("~ %s(Trial):id", s))
    b <- asreml(fx, random = rfs, residual = rs, na.action = na.method(x = "include"),
                data = pheno, workspace = "1gb", trace = FALSE)
    for (frc in c(FALSE, TRUE)) {
        w <- NULL
        g <- tryCatch(withCallingHandlers(quiet(gwasAim(b, panel, "id", Trait = "Trial", trace = FALSE,
                                     force = frc, breakout = 3, maxit = 50)),
                      warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }),
                      error = function(e) paste("ERROR:", conditionMessage(e)))
        lab <- sprintf("base %-11s %s", s, if (frc) "mbf" else "vm ")
        if (is.character(g)) { cat(lab, g, "\n"); next }
        lk <- vapply(g$QTL$diag$lik, function(x) x[3L], 0)
        cat(sprintf("%s QTL: %-40s LRTs: %s\n", lab, paste(g$QTL$qtl, collapse = " "),
                    paste(round(lk, 1), collapse = " ")))
    }
}
