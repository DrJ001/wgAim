.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
data(genoMVsim, package = "wgAim")
genObj <- primeCross(genoMVsim, type = "interval", id = "id")
gm <- primeCross(genoMVsim, type = "marker", id = "id")
M <- do.call(cbind, lapply(gm$geno, function(el) el$imputed.data))
n.g <- nrow(M); p <- ncol(M); n.t <- 3
Mc <- sweep(M, 2L, colMeans(M)); K <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1; D <- diag(sqrt(v), n.t); D %*% S %*% D }
set.seed(9)
U <- t(chol(K)) %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(0.05 * mkcov(c(.12, .08, .15), .6))
L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(mkcov(.05, .6))
G <- U + L; rownames(G) <- rownames(M)
pheno <- expand.grid(id = rownames(G), Rep = c("1", "2"), Trial = paste0("T", 1:3), stringsAsFactors = FALSE)
pheno$yield <- 4.5 + G[cbind(match(pheno$id, rownames(G)), match(pheno$Trial, paste0("T", 1:3)))] +
    rnorm(nrow(pheno), sd = sqrt(.15))
pheno$id <- factor(pheno$id); pheno$Rep <- factor(pheno$Rep); pheno$Trial <- factor(pheno$Trial)
cat("ids in pheno:", nlevels(pheno$id), " in genObj:", nrow(genObj$pheno), "\n")
cat("head ids pheno:", head(levels(pheno$id)), "\n")
cat("head ids genObj:", head(as.character(genObj$pheno$id)), "\n")
base <- asreml(yield ~ Trial + Trial:Rep, random = ~ us(Trial):id,
               residual = ~ dsum(~units | Trial), data = pheno, maxit = 50, trace = FALSE)
base <- update(base, trace = FALSE)
gw <- qtlAim(base, genObj, merge.by = "id", gen.type = "interval", Trait = "Trial",
             TypeI = 0.999, breakout = 1, trace = FALSE, maxit = 50)
str(gw$QTL$diag$lik.mat)
print(gw$QTL$diag$lik.mat)
