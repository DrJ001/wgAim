# =============================================================================
# research/aom_stopping/diag_lrt_rise.R
#
# Why does the MV LRT statistic RISE after the true QTL are fitted?
# Wraps wgAim:::.buildLRTModels to record, at every LRT step:
#   L0        : baseForLRT logL (warm-started chain via .addEffect)
#   L1        : qtlForLRT logL as used by gwasAim (cold-start refit)
#   L1.conv   : its converge flag
#   L1.more   : qtlForLRT continued with extra iterations (maxit = 100)
#   L0.more   : baseForLRT continued with extra iterations
#
# Run: Rscript research/aom_stopping/diag_lrt_rise.R [corgh|us] [n.trials]
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
struct <- if (length(args) >= 1L) args[1L] else "corgh"
n.t    <- if (length(args) >= 2L) as.integer(args[2L]) else 3L
kscale <- if (length(args) >= 3L) as.numeric(args[3L]) else 1     # scales Ga_K
cat("struct:", struct, " trials:", n.t, " Ga_K scale:", kscale, "\n")

.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })

# ---- same data as repro_corgh_3trial.R -------------------------------------
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

# ---- instrument .buildLRTModels --------------------------------------------
ns   <- asNamespace("wgAim")
orig <- get(".buildLRTModels", envir = ns)
rec  <- list(); null.rec <- list()
wrap <- function(qtlModel, baseModel, ntrait, Trait, merge.by, n.fa, vmterms, phenoData, ...) {
    out <- orig(qtlModel, baseModel, ntrait, Trait, merge.by, n.fa, vmterms, phenoData, ...)
    q1  <- out$qtlForLRT; b0 <- out$baseForLRT
    q2  <- suppressWarnings(update(q1, maxit = 100, trace = FALSE))
    q3  <- suppressWarnings(update(q2, maxit = 100, trace = FALSE))
    b2  <- suppressWarnings(update(b0, maxit = 100, trace = FALSE))
    rec[[length(rec) + 1L]] <<- data.frame(
        step = length(rec), L0 = b0$loglik, L0.conv = b0$converge,
        L0.more = b2$loglik,
        L1 = q1$loglik, L1.conv = q1$converge,
        L1.more = q2$loglik, L1.more2 = q3$loglik, L1.more.conv = q3$converge,
        full.struct.logL = qtlModel$loglik)

    # Null-marker experiment (at step 2): add a marker that is NOT a QTL and
    # far from every planted QTL as X + Trial:X to BOTH LRT models and see
    # how the LRT statistic moves.
    if (length(rec) == 3L) {
        nulls <- setdiff(round(p * c(0.05, 0.30, 0.60, 0.90, 0.97)), qtl.cols)
        for (j in nulls) {
            mk  <- colnames(M)[j]
            xv  <- M[, j][match(as.character(phenoData[[merge.by]]), geno.ids)]
            xv[is.na(xv)] <- 0
            phenoData$XNULL <- xv
            ff  <- as.formula(paste(". ~ . + XNULL +", Trait, ":XNULL"))
            qn  <- suppressWarnings(update(q3, fixed. = ff, maxit = 100, trace = FALSE))
            bn  <- suppressWarnings(update(b2, fixed. = ff, maxit = 100, trace = FALSE))
            null.rec[[length(null.rec) + 1L]] <<- data.frame(
                marker = mk, dL0 = bn$loglik - b2$loglik, dL1 = qn$loglik - q3$loglik,
                stat.before = 2 * (q3$loglik - b2$loglik),
                stat.after  = 2 * (qn$loglik - bn$loglik))
        }
    }
    out
}
unlockBinding(".buildLRTModels", ns); assign(".buildLRTModels", wrap, envir = ns)

gw <- gwasAim(baseModel = base, genObj = panel, merge.by = "id",
              Trait = "Trial", str = NULL, TypeI = 0.05,
              trace = FALSE, breakout = 4, maxit = 50)

res <- do.call(rbind, rec)
res$stat.used  <- 2 * (res$L1 - res$L0)
res$stat.conv  <- 2 * (res$L1.more2 - res$L0.more)
cat("\n==== LRT diagnostics (", struct, ",", n.t, "trials ) ====\n")
print(res, digits = 7, row.names = FALSE)
cat("\nQTL:", gw$QTL$qtl, "\n")
if (length(null.rec)) {
    cat("\n==== Null-marker experiment (added at LRT step 2) ====\n")
    print(do.call(rbind, null.rec), digits = 5, row.names = FALSE)
}
