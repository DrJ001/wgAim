# =============================================================================
# scratch/aom_stopping/run_stage1.R
#
# STAGE 1: null calibration of the AOM stopping rule.  NO planted QTL.
# Phenotypes = K-structured polygenic background (Ga_K x kscale, K from the
# population's own markers) + IID line effect + residual.
#
# Under this global null the family-wise error rate of the proposed rule is
# P(step 1 rejects), so each replicate runs ONE selection step (breakout = 1)
# with the standalone diagnostic (scratch/aom_stopping/aom_null_diag.R) attached to
# wgAim:::.qtlSelect for this session only.  Package code is not modified.
#
#   pop = "gwas" : shipped GWAS panel (191 lines, 1501 markers) + 9
#                  ungenotyped lines; gwasAim(); vm path
#   pop = "dh"   : genoMVsim DH cross (120 lines, 600 markers, 594
#                  intervals); qtlAim(gen.type = "interval"); vm path
#
# The entry LRT is forced open (UV TypeI = 0.49, MV TypeI = 0.999) so that the
# AOM rule is evaluated in every replicate; the initial LRT p-value is
# recorded so the full procedure (LRT gate at 0.05 AND AOM) can be assessed.
#
# Run : Rscript scratch/aom_stopping/run_stage1.R <pop> <ntrait> <kscale> <nrep> <seed0> [struct] [B]
# Out : scratch/aom_stopping/stage1/<pop>_t<ntrait>_k<kscale>.csv   (appended per replicate)
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
pop    <- args[1L]
n.t    <- as.integer(args[2L])
kscale <- as.numeric(args[3L])
nrep   <- as.integer(args[4L])
seed0  <- as.integer(args[5L])
struct <- if (length(args) >= 6L) args[6L] else "us"
B      <- if (length(args) >= 7L) as.integer(args[7L]) else 1000L
out.sub <- if (length(args) >= 8L) args[8L] else "stage1"
stopifnot(pop %in% c("gwas", "dh"), n.t %in% c(1L, 3L))

out.dir <- file.path("scratch/aom_stopping/results", out.sub)
dir.create(out.dir, showWarnings = FALSE)
out.csv <- file.path(out.dir, sprintf("%s_t%d_k%s.csv", pop, n.t, kscale))

.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("scratch/aom_stopping/aom_null_diag.R")

# ---- genotypes ---------------------------------------------------------------
if (pop == "gwas") {
    data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
    flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                         map.id = "marker", map.chr = "chr", map.pos = "pos"))
    genObj <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                        map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
    M <- do.call(cbind, lapply(genObj$geno, function(el) el$imputed.data))
    n.extra <- 9L
} else {
    data(genoMVsim, package = "wgAim")
    genObj <- primeCross(genoMVsim, type = "interval", id = "id")
    gm     <- primeCross(genoMVsim, type = "marker",   id = "id")
    M <- do.call(cbind, lapply(gm$geno, function(el) el$imputed.data))
    n.extra <- 0L
}
geno.ids <- rownames(M); n.g <- nrow(M); p <- ncol(M)
Mc <- sweep(M, 2L, colMeans(M))
K  <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
cK <- t(chol(K))
trials <- paste0("T", seq_len(n.t))
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1
                          D <- diag(sqrt(v), n.t); D %*% S %*% D }
Ga_K <- kscale * mkcov(c(0.12, 0.08, 0.15), 0.6)
Ga_L <- mkcov(0.05, 0.6)

sim_pheno <- function(seed) {
    set.seed(seed)
    U <- cK %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
    L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)
    G <- U + L; rownames(G) <- geno.ids                       # NO QTL
    if (n.extra > 0L) {
        Ge <- matrix(rnorm(n.extra * n.t), n.extra) %*% chol(Ga_K + Ga_L)
        rownames(Ge) <- sprintf("NoGeno%02d", seq_len(n.extra))
        G <- rbind(G, Ge)
    }
    ids <- rownames(G)
    ph <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
    mu <- setNames(c(4.5, 4.8, 4.2)[seq_len(n.t)], trials)
    ph$yield <- mu[ph$Trial] + G[cbind(match(ph$id, ids), match(ph$Trial, trials))] +
        ifelse(ph$Rep == "2", 0.1, 0) + rnorm(nrow(ph), sd = sqrt(0.15))
    ph$id <- factor(ph$id); ph$Rep <- factor(ph$Rep); ph$Trial <- factor(ph$Trial)
    ph
}

# ---- attach diagnostic (session only) --------------------------------------
ns   <- asNamespace("wgAim")
orig <- get(".qtlSelect", envir = ns)
last <- NULL
wrap <- function(asm, phenoData, intervalObj, gen.type, selection,
                 exclusion.window, state, verboseLev, merge.by = NULL,
                 Trait = NULL, ntrait = 1L) {
    d <- tryCatch(aomNullDiag(asm, phenoData, intervalObj, state, merge.by,
                              Trait, ntrait, selection = selection, B = B),
                  error = function(e) { cat("  aomNullDiag ERROR:", conditionMessage(e), "\n"); NULL })
    out <- orig(asm, phenoData, intervalObj, gen.type, selection,
                exclusion.window, state, verboseLev, merge.by, Trait, ntrait)
    if (!is.null(d) && is.null(last)) {
        s <- d$summary
        s$chk <- max(abs(d$oint - out$oint[names(d$oint)]) / pmax(abs(out$oint[names(d$oint)]), 1e-8))
        last <<- s
    }
    out
}
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", wrap, envir = ns)

run_one <- function(ph) {
    # ASReml resolves data/formula objects in the global env: fit everything there
    assign("pheno", ph, envir = globalenv())
    assign("rf.", as.formula(sprintf("~ %s(Trial):id", struct)), envir = globalenv())
    fit <- if (n.t == 1L)
        quote(asreml(yield ~ Rep, random = ~ id, na.action = na.method(x = "include"),
                     data = pheno, trace = FALSE))
    else
        quote(asreml(yield ~ Trial + Trial:Rep, random = rf.,
                     residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                     data = pheno, workspace = "1gb", maxit = 50, trace = FALSE))
    TI <- if (n.t == 1L) 0.49 else 0.999
    base <- suppressWarnings(eval(fit, envir = globalenv()))
    assign("base", base, envir = globalenv())
    base <- suppressWarnings(eval(quote(update(base, trace = FALSE)), envir = globalenv()))
    assign("base", base, envir = globalenv())
    assign("genObj", genObj, envir = globalenv())
    call <- if (pop == "gwas") {
        if (n.t == 1L) quote(gwasAim(base, genObj, merge.by = "id", TypeI = TI.,
                                     breakout = 1, trace = FALSE))
        else           quote(gwasAim(base, genObj, merge.by = "id", Trait = "Trial",
                                     TypeI = TI., breakout = 1, trace = FALSE, maxit = 50))
    } else {
        if (n.t == 1L) quote(qtlAim(base, genObj, merge.by = "id", gen.type = "interval",
                                    TypeI = TI., breakout = 1, trace = FALSE))
        else           quote(qtlAim(base, genObj, merge.by = "id", gen.type = "interval",
                                    Trait = "Trial", TypeI = TI., breakout = 1,
                                    trace = FALSE, maxit = 50))
    }
    assign("TI.", TI, envir = globalenv())
    res <- NULL
    invisible(capture.output(suppressMessages(suppressWarnings(
        res <- eval(call, envir = globalenv())))))
    res
}

cat(sprintf("STAGE 1 NULL: pop=%s ntrait=%d kscale=%s nrep=%d seeds %d..%d struct=%s B=%d\n",
            pop, n.t, kscale, nrep, seed0, seed0 + nrep - 1L, struct, B))
for (r in seq_len(nrep)) {
    seed <- seed0 + r - 1L
    t0 <- Sys.time()
    last <- NULL
    row <- tryCatch({
        gw  <- run_one(sim_pheno(seed))
        lik <- gw$QTL$diag$lik.mat      # NULL from qtlAim() when the gate stops
        base.row <- data.frame(pop = pop, ntrait = n.t, kscale = kscale, seed = seed,
                               lrt0   = if (is.null(lik)) NA_real_ else lik[1, "Statistic"],
                               p.lrt0 = if (is.null(lik)) 1 else lik[1, "Pvalue"])
        if (is.null(last)) {
            cbind(base.row, status = "gate-stop", marker = NA, max.oint = NA, p.marker = NA,
                  p.bonf = 1, fwe.thr = NA, p.fwe = 1, Meff = NA, p.maha.bonf = 1,
                  p.fwe.chr = 1, null.mean = NA, chk = NA,
                  p.fwe.stud = 1, p.fwe.studM = 1, fwe.thr.stud = NA, fwe.thr.studM = NA,
                  obs.ratio = NA, obs.Sdev = NA, p.fwe.studO = 1, p.fwe.studP = 1,
                  n.unpinned = NA, w.min = NA, w.max = NA)
        } else {
            cbind(base.row, status = "ok", marker = last$marker, max.oint = last$max.oint,
                  p.marker = last$p.marker, p.bonf = last$p.bonf, fwe.thr = last$fwe.thr,
                  p.fwe = last$p.fwe, Meff = last$Meff, p.maha.bonf = last$p.maha.bonf,
                  p.fwe.chr = last$p.fwe.chr, null.mean = last$null.mean.oint, chk = last$chk,
                  p.fwe.stud = last$p.fwe.stud, p.fwe.studM = last$p.fwe.studM,
                  fwe.thr.stud = last$fwe.thr.stud, fwe.thr.studM = last$fwe.thr.studM,
                  obs.ratio = last$obs.ratio, obs.Sdev = last$obs.Sdev,
                  p.fwe.studO = last$p.fwe.studO, p.fwe.studP = last$p.fwe.studP,
                  n.unpinned = last$n.unpinned, w.min = last$w.min, w.max = last$w.max)
        }
    }, error = function(e) {
        data.frame(pop = pop, ntrait = n.t, kscale = kscale, seed = seed, lrt0 = NA,
                   p.lrt0 = NA, status = paste("error:", substr(conditionMessage(e), 1, 80)),
                   marker = NA, max.oint = NA, p.marker = NA, p.bonf = NA, fwe.thr = NA,
                   p.fwe = NA, Meff = NA, p.maha.bonf = NA, p.fwe.chr = NA,
                   null.mean = NA, chk = NA, p.fwe.stud = NA, p.fwe.studM = NA,
                   fwe.thr.stud = NA, fwe.thr.studM = NA, obs.ratio = NA, obs.Sdev = NA,
                   p.fwe.studO = NA, p.fwe.studP = NA, n.unpinned = NA, w.min = NA, w.max = NA)
    })
    row$secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    write.table(row, out.csv, sep = ",", row.names = FALSE,
                col.names = !file.exists(out.csv), append = file.exists(out.csv))
    cat(sprintf("rep %3d seed %d  %-10s LRT0 %7.2f  p.fwe %.4f  studM %.4f  studO %.4f  studP %.4f  unpinned %s  w %.3f-%.3f  [%.0fs]\n",
                r, seed, substr(row$status, 1, 10), row$lrt0, row$p.fwe,
                row$p.fwe.studM, row$p.fwe.studO, row$p.fwe.studP,
                as.character(row$n.unpinned), row$w.min, row$w.max, row$secs))
}
cat("DONE\n")
