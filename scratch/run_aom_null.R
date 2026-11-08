# =============================================================================
# scratch/run_aom_null.R
#
# Null simulation for the AOM stopping rule: NO planted QTL.
# Phenotypes = K-structured polygenic background (Ga_K x kscale) + IID line
# effects + residual.  Under this global null the family-wise error rate of
# the proposed procedure is P(step 1 rejects), so each replicate runs one
# iteration of gwasAim() (breakout = 1) with the standalone diagnostic
# (scratch/aom_null_diag.R) attached to .qtlSelect for this session only.
#
# Run: Rscript scratch/run_aom_null.R [kscale] [nrep] [first.seed] [struct]
# Output: scratch/aom_null_k<kscale>.csv (appended per replicate)
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
kscale <- if (length(args) >= 1L) as.numeric(args[1L]) else 1
nrep   <- if (length(args) >= 2L) as.integer(args[2L]) else 20L
seed0  <- if (length(args) >= 3L) as.integer(args[3L]) else 1001L
struct <- if (length(args) >= 4L) args[4L] else "us"
n.t    <- 3L
out.csv <- sprintf("C:/jtaylor/Claude/wgAim/scratch/aom_null_k%s.csv", kscale)

.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("C:/jtaylor/Claude/wgAim/scratch/aom_null_diag.R")

data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos"))
panel <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
M <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
geno.ids <- rownames(M); n.g <- nrow(M); p <- ncol(M)
Mc <- sweep(M, 2L, colMeans(M))
K  <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
cK <- t(chol(K))
trials <- paste0("T", seq_len(n.t))
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1; D <- diag(sqrt(v)); D %*% S %*% D }
Ga_K <- kscale * mkcov(c(0.12, 0.08, 0.15), 0.6)
Ga_L <- mkcov(c(0.05, 0.05, 0.05), 0.6)

sim_pheno <- function(seed) {
    set.seed(seed)
    U <- cK %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
    L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)
    G <- U + L; rownames(G) <- geno.ids                 # NO QTL
    extra <- sprintf("NoGeno%02d", 1:9)
    G.extra <- matrix(rnorm(9 * n.t), 9) %*% chol(Ga_K + Ga_L); rownames(G.extra) <- extra
    G <- rbind(G, G.extra); ids <- rownames(G)
    ph <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
    mu <- setNames(c(4.5, 4.8, 4.2), trials)
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
                              Trait, ntrait, B = 2000L),
                  error = function(e) { cat("aomNullDiag failed:", conditionMessage(e), "\n"); NULL })
    out <- orig(asm, phenoData, intervalObj, gen.type, selection,
                exclusion.window, state, verboseLev, merge.by, Trait, ntrait)
    if (!is.null(d)) { last <<- d$summary; last$oint.pkg <<- max(out$oint) }
    out
}
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", wrap, envir = ns)

cat(sprintf("NULL simulation: kscale = %s, nrep = %d, seeds %d..%d, struct = %s\n",
            kscale, nrep, seed0, seed0 + nrep - 1L, struct))
for (r in seq_len(nrep)) {
    seed <- seed0 + r - 1L
    t0 <- Sys.time()
    pheno <- sim_pheno(seed)
    res <- tryCatch({
        base <- asreml(yield ~ Trial + Trial:Rep,
                       random = as.formula(sprintf("~ %s(Trial):id", struct)),
                       residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                       data = pheno, workspace = "1gb", trace = FALSE)
        base <- update(base, trace = FALSE)
        last <- NULL
        gw <- suppressMessages(gwasAim(baseModel = base, genObj = panel, merge.by = "id",
                      Trait = "Trial", str = NULL, TypeI = 0.999,
                      trace = FALSE, breakout = 1, maxit = 50))
        lik <- gw$QTL$diag$lik.mat
        if (is.null(last)) {  # initial LRT failed even at TypeI = 0.999 (stat <= 0)
            data.frame(seed = seed, kscale = kscale,
                   lrt0 = lik[1, "Statistic"], p.lrt0 = lik[1, "Pvalue"],
                   marker = "GATE-STOP", max.oint = NA, chk = TRUE,
                   p.marker = NA, p.bonf = 1, fwe.thr = NA, p.fwe = 1, Meff = NA,
                   p.maha.bonf = 1, null.mean = NA,
                   secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
        } else data.frame(seed = seed, kscale = kscale,
                   lrt0 = lik[1, "Statistic"], p.lrt0 = lik[1, "Pvalue"],
                   marker = last$marker, max.oint = last$max.oint,
                   chk = abs(last$max.oint - last$oint.pkg) < 1e-6,
                   p.marker = last$p.marker, p.bonf = last$p.bonf,
                   fwe.thr = last$fwe.thr, p.fwe = last$p.fwe, Meff = last$Meff,
                   p.maha.bonf = last$p.maha.bonf, null.mean = last$null.mean.oint,
                   secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    }, error = function(e) {
        message("seed ", seed, " failed: ", conditionMessage(e)); NULL })
    if (!is.null(res)) {
        write.table(res, out.csv, sep = ",", row.names = FALSE,
                    col.names = !file.exists(out.csv), append = file.exists(out.csv))
        cat(sprintf("rep %2d seed %d  LRT0 %6.2f (p %.3g)  max.oint %5.2f  p.fwe %.4f  p.bonf %.4f  [%.0fs]%s\n",
                    r, seed, res$lrt0, res$p.lrt0, res$max.oint, res$p.fwe, res$p.bonf, res$secs,
                    if (res$marker == "GATE-STOP") "  (initial LRT <= 0: gate stop)" else ""))
    }
}

R <- read.csv(out.csv)
cat("\n==== NULL summary (kscale =", kscale, ") ====\n")
cat(sprintf("replicates: %d   diagnostic matches package oint: %s\n", nrow(R), all(R$chk)))
cat(sprintf("rejections at 0.05:  FWE-sim %d   Bonferroni %d   Mahalanobis-Bonf %d   initial LRT %d\n",
            sum(R$p.fwe <= 0.05), sum(R$p.bonf <= 0.05), sum(R$p.maha.bonf <= 0.05),
            sum(R$p.lrt0 <= 0.05)))
cat("p.fwe quantiles:\n"); print(round(quantile(R$p.fwe, c(0, .1, .25, .5, .75, .9, 1)), 3))
cat(sprintf("KS test p.fwe ~ U(0,1): p = %.3f\n", suppressWarnings(ks.test(R$p.fwe, "punif")$p.value)))
