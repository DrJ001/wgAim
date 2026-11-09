# Why is the REML pin (S.obs = S.target) far from holding in weak MV fits?
# One DH 3-trial weak-background replicate (as run_stage1.R, seed given):
# inspect the us(Trial):vm fit used by .qtlSelect -- convergence, the vm
# variance parameters and their constraint codes, eigenvalues of Ga, and
# whether continuing the fit (more iterations) restores the pin.
args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args)) as.integer(args[1]) else 2002L
.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("scratch/aom_stopping/aom_null_diag.R")

data(genoMVsim, package = "wgAim")
genObj <- primeCross(genoMVsim, type = "interval", id = "id")
gm <- primeCross(genoMVsim, type = "marker", id = "id")
M <- do.call(cbind, lapply(gm$geno, function(el) el$imputed.data))
n.g <- nrow(M); p <- ncol(M); n.t <- 3L; trials <- paste0("T", 1:3)
Mc <- sweep(M, 2L, colMeans(M)); K <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1
                          D <- diag(sqrt(v), n.t); D %*% S %*% D }
Ga_K <- 0.05 * mkcov(c(0.12, 0.08, 0.15), 0.6); Ga_L <- mkcov(0.05, 0.6)
set.seed(seed)                                          # identical to run_stage1.R sim_pheno
U <- t(chol(K)) %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)
G <- U + L; rownames(G) <- rownames(M); ids <- rownames(G)
pheno <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
pheno$yield <- setNames(c(4.5, 4.8, 4.2), trials)[pheno$Trial] +
    G[cbind(match(pheno$id, ids), match(pheno$Trial, trials))] +
    ifelse(pheno$Rep == "2", 0.1, 0) + rnorm(nrow(pheno), sd = sqrt(0.15))
pheno$id <- factor(pheno$id); pheno$Rep <- factor(pheno$Rep); pheno$Trial <- factor(pheno$Trial)

base <- asreml(yield ~ Trial + Trial:Rep, random = ~ us(Trial):id,
               residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
               data = pheno, workspace = "1gb", maxit = 50, trace = FALSE)
base <- update(base, trace = FALSE)

ns <- asNamespace("wgAim"); orig <- get(".qtlSelect", envir = ns); cap <- NULL
wrap <- function(asm, phenoData, intervalObj, gen.type, selection, exclusion.window,
                 state, verboseLev, merge.by = NULL, Trait = NULL, ntrait = 1L) {
    if (is.null(cap)) cap <<- list(asm = asm, phenoData = phenoData, intervalObj = intervalObj,
                                   state = state, merge.by = merge.by)
    orig(asm, phenoData, intervalObj, gen.type, selection, exclusion.window, state,
         verboseLev, merge.by, Trait, ntrait)
}
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", wrap, envir = ns)
invisible(capture.output(suppressWarnings(suppressMessages(
    gw <- qtlAim(base, genObj, merge.by = "id", gen.type = "interval", Trait = "Trial",
                 TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)))))

asm <- cap$asm
cat("seed", seed, "\nvm model converged:", asm$converge, "  logL", asm$loglik, "\n")
vp <- asm$vparameters; sel <- grep("vm", names(vp))
print(data.frame(parameter = names(vp)[sel], value = signif(vp[sel], 4),
                 con = as.character(asm$vparameters.con[sel])), row.names = FALSE)
cat("\nvarcomp (bound column):\n"); vc <- summary(asm)$varcomp
print(vc[grep("vm", rownames(vc)), c("component", "bound", "%ch")])
Ga <- get(".extractVmGa", envir = ns)(asm, "Trial", levels(pheno$Trial), 1)
cat("\neigenvalues of Ga:", signif(eigen(Ga, only.values = TRUE)$values, 4), "\n")

pin <- function(model, lab) {
    d <- aomNullDiag(model, cap$phenoData, cap$intervalObj, cap$state, cap$merge.by,
                     "Trial", 3L, B = 50L)$summary
    cat(sprintf("%-28s converge %-5s logL %.4f  ratio %.4f  W eig %.3f-%.3f  unpinned %d\n",
                lab, model$converge, model$loglik, d$obs.ratio, d$w.min, d$w.max, d$n.unpinned))
}
pin(asm, "as used by .qtlSelect")
assign("phenoData", cap$phenoData, envir = globalenv())
asm2 <- suppressWarnings(update(asm, maxit = 200, trace = FALSE)); pin(asm2, "continued (maxit 200)")
asm3 <- suppressWarnings(update(asm2, maxit = 200, trace = FALSE)); pin(asm3, "continued again")
vp3 <- asm3$vparameters; s3 <- grep("vm", names(vp3))
print(data.frame(parameter = names(vp3)[s3], value = signif(vp3[s3], 4),
                 con = as.character(asm3$vparameters.con[s3])), row.names = FALSE)
