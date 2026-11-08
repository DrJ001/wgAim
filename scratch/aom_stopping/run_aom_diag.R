# =============================================================================
# scratch/aom_stopping/run_aom_diag.R
#
# Runs gwasAim() on the simulated 3-trial data with the standalone AOM-null
# diagnostic (scratch/aom_stopping/aom_null_diag.R) attached to wgAim:::.qtlSelect for this
# R session only.  Package code is NOT modified; selection is unchanged.
#
# TypeI = 0.999 so the LRT never stops the loop; breakout = 8.  For each
# iteration we report the LRT and the AOM-based rules side by side, so we can
# see where each WOULD have stopped.
#
# Run: Rscript scratch/aom_stopping/run_aom_diag.R [us|corgh|...] [n.trials] [Ga_K scale] [breakout]
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
struct <- if (length(args) >= 1L) args[1L] else "us"
n.t    <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
kscale <- if (length(args) >= 3L) as.numeric(args[3L]) else 1
brk    <- if (length(args) >= 4L) as.integer(args[4L]) else 8L
cat("struct:", struct, " trials:", n.t, " Ga_K scale:", kscale, " breakout:", brk, "\n")

.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("scratch/aom_stopping/aom_null_diag.R")

# ---- data: identical to repro_corgh_3trial.R / diag_lrt_rise.R -------------
data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos"))
panel <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
set.seed(7)
M <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
geno.ids <- rownames(M); n.g <- nrow(M); p <- ncol(M)
Mc <- sweep(M, 2L, colMeans(M))
K  <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
trials <- paste0("T", seq_len(n.t))
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1; D <- diag(sqrt(v)); D %*% S %*% D }
Ga_K <- kscale * mkcov(c(0.12, 0.08, 0.15, 0.10, 0.14, 0.09), 0.6)
Ga_L <- mkcov(c(0.05, 0.05, 0.05), 0.6)
U <- t(chol(K)) %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)
qtl.cols <- c(round(p * 0.15), round(p * 0.45), round(p * 0.75))
beta <- rbind(rep(0.20, n.t), rep(-0.18, n.t), seq(0.25, -0.25, length.out = n.t))
Q <- M[, qtl.cols] %*% beta
cat("Planted QTL markers:", colnames(M)[qtl.cols], "\n")
G <- U + L + Q; rownames(G) <- geno.ids
extra <- sprintf("NoGeno%02d", 1:9)
G.extra <- matrix(rnorm(9 * n.t), 9) %*% chol(Ga_K + Ga_L); rownames(G.extra) <- extra
G <- rbind(G, G.extra); ids <- rownames(G)
pheno <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
mu <- setNames(rep_len(c(4.5, 4.8, 4.2, 4.6, 4.9, 4.3), n.t), trials)
pheno$yield <- mu[pheno$Trial] + G[cbind(match(pheno$id, ids), match(pheno$Trial, trials))] +
    ifelse(pheno$Rep == "2", 0.1, 0) + rnorm(nrow(pheno), sd = sqrt(0.15))
pheno$id <- factor(pheno$id); pheno$Rep <- factor(pheno$Rep); pheno$Trial <- factor(pheno$Trial)

base <- asreml(yield ~ Trial + Trial:Rep, random = as.formula(sprintf("~ %s(Trial):id", struct)),
               residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
               data = pheno, workspace = "1gb", trace = FALSE)
base <- update(base, trace = FALSE)

# ---- attach diagnostic to .qtlSelect (session only) ------------------------
ns   <- asNamespace("wgAim")
orig <- get(".qtlSelect", envir = ns)
aom  <- list()
wrap <- function(asm, phenoData, intervalObj, gen.type, selection,
                 exclusion.window, state, verboseLev, merge.by = NULL,
                 Trait = NULL, ntrait = 1L) {
    d <- tryCatch(aomNullDiag(asm, phenoData, intervalObj, state, merge.by,
                              Trait, ntrait, B = 2000L),
                  error = function(e) { message("aomNullDiag failed: ", conditionMessage(e)); NULL })
    out <- orig(asm, phenoData, intervalObj, gen.type, selection,
                exclusion.window, state, verboseLev, merge.by, Trait, ntrait)
    if (!is.null(d)) {
        s <- d$summary
        s$selected  <- out$qtl
        s$oint.pkg  <- max(out$oint)            # package's own max oint (sanity)
        aom[[length(aom) + 1L]] <<- s
    }
    out
}
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", wrap, envir = ns)

gw <- gwasAim(baseModel = base, genObj = panel, merge.by = "id",
              Trait = "Trial", str = NULL, TypeI = 0.999,
              trace = FALSE, breakout = brk, maxit = 50)

# ---- report -----------------------------------------------------------------
lik <- gw$QTL$diag$lik.mat
A   <- do.call(rbind, aom)
planted <- colnames(M)[qtl.cols]
# engine names "Chr.<chr>.<idx>" follow the column order of M chromosome by chromosome
eng.names <- unlist(lapply(names(panel$geno), function(ch)
    paste("Chr", ch, seq_len(ncol(panel$geno[[ch]]$imputed.data)), sep = ".")))
snp.of <- setNames(colnames(M), eng.names)
mk <- A$marker
A$snp <- snp.of[mk]
A$planted <- ifelse(A$snp %in% planted, "PLANTED", "")
rep <- data.frame(
    iter       = seq_len(nrow(A)),
    selected   = A$selected,
    snp        = A$snp,
    planted    = A$planted,
    LRT.before = round(lik[seq_len(nrow(A)), "Statistic"], 2),
    p.LRT      = signif(lik[seq_len(nrow(A)), "Pvalue"], 3),
    max.oint   = round(A$max.oint, 2),
    chk.pkg    = round(A$oint.pkg, 2),
    p.marker   = signif(A$p.marker, 3),
    p.bonf     = signif(A$p.bonf, 3),
    fwe.thr    = round(A$fwe.thr, 2),
    p.fwe      = signif(A$p.fwe, 3),
    Meff       = round(A$Meff),
    p.maha.bonf = signif(A$p.maha.bonf, 3),
    null.mean  = round(A$null.mean.oint, 3))
cat("\n==== AOM stopping-rule diagnostic (", struct, ",", n.t, "trials, Ga_K x", kscale, ") ====\n")
cat("markers:", A$n.markers[1], "  genotyped lines:", A$n.lines[1], "\n")
cat("Planted:", planted, "\n\n")
print(rep, row.names = FALSE)
cat("\nLRT after final iteration:", round(lik[nrow(lik), "Statistic"], 2),
    " p =", signif(lik[nrow(lik), "Pvalue"], 3), "\n")
saveRDS(list(rep = rep, A = A, lik = lik),
        sprintf("scratch/aom_stopping/results/aom_diag_%s_%dt_k%s.rds", struct, n.t, kscale))
