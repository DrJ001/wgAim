# =============================================================================
# scratch/run_stage0.R
#
# Stage 0 verification: does the standalone AOM diagnostic
# (scratch/aom_null_diag.R) reproduce the package's outlier statistics on
# every engine path?  One selection step per case (breakout = 1,
# TypeI = 0.999).  Diagnostic attached to wgAim:::.qtlSelect for this R
# session only; package code unchanged.
#
# Cases: gwasAim  UV/MV x vm/mbf(force)
#        qtlAim   interval UV/MV x vm/mbf(force), selection = "chromosome"
#                 UV/MV, gen.type = "marker" UV
#
# Run: Rscript scratch/run_stage0.R [B]
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
B    <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L

.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("C:/jtaylor/Claude/wgAim/scratch/aom_null_diag.R")

# ---- attach diagnostic ------------------------------------------------------
ns   <- asNamespace("wgAim")
orig <- get(".qtlSelect", envir = ns)
cap  <- NULL
wrap <- function(asm, phenoData, intervalObj, gen.type, selection,
                 exclusion.window, state, verboseLev, merge.by = NULL,
                 Trait = NULL, ntrait = 1L) {
    t0 <- proc.time()[["elapsed"]]
    d <- tryCatch(aomNullDiag(asm, phenoData, intervalObj, state, merge.by,
                              Trait, ntrait, selection = selection, B = B),
                  error = function(e) { cat("  aomNullDiag ERROR:", conditionMessage(e), "\n"); NULL })
    t.diag <- proc.time()[["elapsed"]] - t0
    out <- orig(asm, phenoData, intervalObj, gen.type, selection,
                exclusion.window, state, verboseLev, merge.by, Trait, ntrait)
    if (is.null(cap)) cap <<- list(d = d, out = out, selection = selection,
                                   t.diag = t.diag, asm = asm)
    out
}
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", wrap, envir = ns)

# Thin a wgPanel / wgCross to every `by`-th column so that lines >= markers
# and the engine takes the mbf path NATURALLY (no force = TRUE).
thin_obj <- function(obj, by) {
    for (ch in names(obj$geno)) {
        g <- obj$geno[[ch]]
        for (slot in c("imputed.data", "data", "interval.data")) {
            if (is.null(g[[slot]])) next
            keep <- seq(1L, ncol(g[[slot]]), by = by)
            g[[slot]] <- g[[slot]][, keep, drop = FALSE]
            if (slot == "interval.data" && !is.null(g$inferred.map))
                g$inferred.map <- g$inferred.map[keep]
            if (slot %in% c("imputed.data", "data") && !is.null(g$map) &&
                length(g$map) > length(keep)) g$map <- g$map[keep]
        }
        obj$geno[[ch]] <- g
    }
    obj
}

results <- list()
run_case <- function(label, fun) {
    cat(sprintf("\n---- %s ----\n", label))
    cap <<- NULL
    t0 <- proc.time()[["elapsed"]]
    # evaluate in the global env: gwasAim/qtlAim assign covObj/phenoData into
    # the caller's frame, which ASReml must be able to see
    ok <- tryCatch({ invisible(capture.output(suppressMessages(suppressWarnings(
              invisible(eval(fun, envir = globalenv())))))); TRUE },
                   error = function(e) { cat("  RUN ERROR:", conditionMessage(e), "\n"); FALSE })
    if (is.null(cap) || is.null(cap$d)) {
        results[[label]] <<- data.frame(case = label, status = "NO DIAG")
        cat("  no diagnostic captured\n"); return(invisible())
    }
    d <- cap$d; out <- cap$out
    pk <- out$oint[names(d$oint)]
    rel <- max(abs(d$oint - pk) / pmax(abs(pk), 1e-8))
    if (rel > 1e-6) {   # detail of the discrepancy
        dd <- data.frame(diag = d$oint, pkg = pk, absdiff = abs(d$oint - pk),
                         ratio = d$oint / pk)
        cat(sprintf("  DISCREPANCY: max abs diff %.4g; ratio diag/pkg: min %.4f median %.4f max %.4f; cor %.6f\n",
                    max(dd$absdiff), min(dd$ratio, na.rm = TRUE), median(dd$ratio, na.rm = TRUE),
                    max(dd$ratio, na.rm = TRUE), cor(dd$diag, dd$pkg)))
        cat(sprintf("  argmax diag %s, argmax pkg %s\n", names(d$oint)[which.max(d$oint)],
                    names(pk)[which.max(pk)]))
        print(head(dd[order(-dd$absdiff), ], 5))
        if (d$summary$path == "mbf" && d$summary$ntrait == 1L) {
            # package UV mbf uses the fitted coefficients / vcoeff directly
            asm <- cap$asm
            s2  <- asm$sigma2
            if (asm$vparameters.con[length(asm$vparameters.con)] == 4) s2 <- 1
            rr  <- grep("mbf", rownames(asm$coefficients$random))
            qpk <- asm$coefficients$random[rr, 1L]
            vpk <- asm$vparameters[grep("mbf.*ints", names(asm$vparameters))] * s2 -
                   s2 * asm$vcoeff$random[rr]
            qdg <- d$Q[, 1]; vdg <- d$den * d$Ga[1, 1]      # den = v / Ga for UV
            cat(sprintf("  UV mbf numerator q: cor %.8f, max rel diff %.3g\n",
                        cor(qdg, qpk), max(abs(qdg - qpk) / max(abs(qpk)))))
            cat(sprintf("  UV mbf variance  v: diag/pkg ratio min %.4f median %.4f max %.4f\n",
                        min(vdg / vpk), median(vdg / vpk), max(vdg / vpk)))
            cat(sprintf("  sigma2 = %.5g, vcoeff range %.4g-%.4g, prior %.4g\n",
                        s2, min(asm$vcoeff$random[rr]), max(asm$vcoeff$random[rr]),
                        asm$vparameters[grep("mbf.*ints", names(asm$vparameters))] * s2))
        }
        assign(paste0("disc_", gsub("[^A-Za-z]", "", label)), list(d = d, out = out),
               envir = globalenv())
    }
    chr.rel <- NA_real_
    if (identical(cap$selection, "chromosome") && !is.null(out$ochr)) {
        pc <- out$ochr[names(d$ochr)]
        chr.rel <- max(abs(d$ochr - pc) / pmax(abs(pc), 1e-8))
    }
    s <- d$summary
    r <- data.frame(case = label, status = if (ok) "ok" else "run error",
                    path = s$path, t = s$ntrait, m = s$n.markers, n = s$n.lines,
                    max.rel.diff.oint = signif(rel, 3),
                    max.rel.diff.ochr = signif(chr.rel, 3),
                    max.oint = round(s$max.oint, 3), p.fwe = round(s$p.fwe, 4),
                    p.bonf = signif(s$p.bonf, 3), Meff = round(s$Meff),
                    p.fwe.chr = round(s$p.fwe.chr, 4),
                    null.mean = round(s$null.mean.oint, 4),
                    t.diag = round(cap$t.diag, 1),
                    t.factor = round(d$timing[["factor"]], 2),
                    t.sim = round(d$timing[["simulate"]], 2),
                    t.total = round(proc.time()[["elapsed"]] - t0, 1))
    results[[label]] <<- r
    print(r, row.names = FALSE)
}

# =============================================================================
# GWAS data: shipped panel + simulated 3-trial phenotypes (planted QTL, strong
# background, 9 ungenotyped lines) -- as in scratch/run_aom_diag.R
# =============================================================================
data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos"))
panel <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
set.seed(7)
M <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
geno.ids <- rownames(M); n.g <- nrow(M); p <- ncol(M); n.t <- 3L
Mc <- sweep(M, 2L, colMeans(M))
K  <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
trials <- paste0("T", 1:n.t)
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1; D <- diag(sqrt(v)); D %*% S %*% D }
Ga_K <- mkcov(c(0.12, 0.08, 0.15), 0.6); Ga_L <- mkcov(0.05, 0.6)
U <- t(chol(K)) %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)
qc <- c(round(p * 0.15), round(p * 0.45), round(p * 0.75))
Q  <- M[, qc] %*% rbind(rep(0.20, n.t), rep(-0.18, n.t), seq(0.25, -0.25, length.out = n.t))
G  <- U + L + Q; rownames(G) <- geno.ids
G.extra <- matrix(rnorm(9 * n.t), 9) %*% chol(Ga_K + Ga_L)
rownames(G.extra) <- sprintf("NoGeno%02d", 1:9)
G <- rbind(G, G.extra); ids <- rownames(G)
gph <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
gph$yield <- c(T1 = 4.5, T2 = 4.8, T3 = 4.2)[gph$Trial] +
    G[cbind(match(gph$id, ids), match(gph$Trial, trials))] +
    ifelse(gph$Rep == "2", 0.1, 0) + rnorm(nrow(gph), sd = sqrt(0.15))
gph$id <- factor(gph$id); gph$Rep <- factor(gph$Rep); gph$Trial <- factor(gph$Trial)
gph1 <- droplevels(subset(gph, Trial == "T1"))

g_uv <- asreml(yield ~ Rep, random = ~ id, na.action = na.method(x = "include"),
               data = gph1, trace = FALSE)
g_mv <- asreml(yield ~ Trial + Trial:Rep, random = ~ us(Trial):id,
               residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
               data = gph, workspace = "1gb", trace = FALSE)
g_mv <- update(g_mv, trace = FALSE)

run_case("gwas UV vm", quote(
    gwasAim(g_uv, panel, merge.by = "id", TypeI = 0.49, breakout = 1, trace = FALSE)))
run_case("gwas MV vm (us)", quote(
    gwasAim(g_mv, panel, merge.by = "id", Trait = "Trial", TypeI = 0.999,
            breakout = 1, trace = FALSE, maxit = 50)))
panel.small <- thin_obj(panel, 10L)
cat("thinned panel markers:", sum(sapply(panel.small$geno, function(g) ncol(g$imputed.data))), "\n")
run_case("gwas UV mbf (thinned)", quote(
    gwasAim(g_uv, panel.small, merge.by = "id", TypeI = 0.49,
            breakout = 1, trace = FALSE)))
run_case("gwas MV mbf (thinned, us)", quote(
    gwasAim(g_mv, panel.small, merge.by = "id", Trait = "Trial",
            TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)))

# =============================================================================
# qtlAim data: genoMVsim / phenoMVsim (120 DH lines, 600 markers, 8 trials)
# =============================================================================
data(genoMVsim, package = "wgAim"); data(phenoMVsim, package = "wgAim")
qint <- primeCross(genoMVsim, type = "interval", id = "id")
qmrk <- primeCross(genoMVsim, type = "marker",   id = "id")
phenoMVsim$id <- factor(phenoMVsim$id); phenoMVsim$Rep <- factor(phenoMVsim$Rep)
phenoMVsim$Trial <- factor(phenoMVsim$Trial, levels = paste0("Trial", 1:8))
qph1 <- droplevels(subset(phenoMVsim, Trial == "Trial1"))
qph3 <- droplevels(subset(phenoMVsim, Trial %in% paste0("Trial", 1:3)))

q_uv <- asreml(yield ~ Rep, random = ~ id, na.action = na.method(x = "include"),
               data = qph1, trace = FALSE)
q_mv <- asreml(yield ~ Trial + Trial:Rep, random = ~ fa(Trial, 1):id,
               residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
               data = qph3, maxit = 80, trace = FALSE)
q_mv <- update(q_mv, trace = FALSE)
q_dg <- asreml(yield ~ Trial + Trial:Rep, random = ~ diag(Trial):id,
               residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
               data = qph3, maxit = 80, trace = FALSE)

run_case("qtl interval UV vm", quote(
    qtlAim(q_uv, qint, merge.by = "id", gen.type = "interval", TypeI = 0.49,
           breakout = 1, trace = FALSE)))
run_case("qtl interval MV vm (fa1)", quote(
    qtlAim(q_mv, qint, merge.by = "id", gen.type = "interval", Trait = "Trial",
           TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)))
qint.small <- thin_obj(qint, 6L)
cat("thinned intervals:", sum(sapply(qint.small$geno, function(g) ncol(g$interval.data))), "\n")
run_case("qtl interval UV mbf (thinned)", quote(
    qtlAim(q_uv, qint.small, merge.by = "id", gen.type = "interval",
           TypeI = 0.49, breakout = 1, trace = FALSE)))
run_case("qtl interval MV mbf (thinned, diag)", quote(
    qtlAim(q_dg, qint.small, merge.by = "id", gen.type = "interval", Trait = "Trial",
           TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)))
run_case("qtl chromosome UV vm", quote(
    qtlAim(q_uv, qint, merge.by = "id", gen.type = "interval", selection = "chromosome",
           TypeI = 0.49, breakout = 1, trace = FALSE)))
run_case("qtl chromosome MV vm (fa1)", quote(
    qtlAim(q_mv, qint, merge.by = "id", gen.type = "interval", selection = "chromosome",
           Trait = "Trial", TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)))
run_case("qtl marker UV vm", quote(
    qtlAim(q_uv, qmrk, merge.by = "id", gen.type = "marker", TypeI = 0.49,
           breakout = 1, trace = FALSE)))

# =============================================================================
# Summary
# =============================================================================
R <- do.call(rbind, lapply(results, function(r) {
    full <- c("case", "status", "path", "t", "m", "n", "max.rel.diff.oint",
              "max.rel.diff.ochr", "max.oint", "p.fwe", "p.bonf", "Meff",
              "p.fwe.chr", "null.mean", "t.diag", "t.factor", "t.sim", "t.total")
    r[setdiff(full, names(r))] <- NA; r[full] }))
cat("\n\n==== STAGE 0 SUMMARY (B =", B, ") ====\n")
print(R, row.names = FALSE)
saveRDS(R, "C:/jtaylor/Claude/wgAim/scratch/stage0_results.rds")
