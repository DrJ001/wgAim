# =============================================================================
# research/aom_stopping/repro_corgh_3trial.R
#
# Reproduce colleague's failure: gwasAim() with a 3-trial data set where the
# genotype x trial term in the base model is corgh(Trial):id or us(Trial):id.
#
# Reported error (after "Predict step for outlier statistics"):
#   Error in array(x, c(length(x), 1L), ...) :
#     'data' must be of a vector type, was 'NULL'
#
# Data: shipped GWAS panel genotypes (filtered + KNN imputed) with NEW
# 3-trial phenotypes simulated here so the K-structured genetic signal is
# strong (initial LRT large, as in the colleague's run).  9 phenotyped lines
# have no genotypes, so the Gsave / Gomit path in .fixLines() is exercised.
#
# Run:  Rscript research/aom_stopping/repro_corgh_3trial.R [corgh|us] [str|NULL] [n.trials]
# =============================================================================

args   <- commandArgs(trailingOnly = TRUE)
struct <- if (length(args) >= 1L) args[1L] else "corgh"
str.arg <- if (length(args) >= 2L && args[2L] != "NULL") args[2L] else NULL
n.t     <- if (length(args) >= 3L) as.integer(args[3L]) else 3L

.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({
    library(asreml)
    library(wgAim)
})
cat("wgAim loaded from:", find.package("wgAim"), "\n")
cat("base structure:", struct, "  str:", if (is.null(str.arg)) "NULL" else str.arg, "\n")

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")

flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos"))
panel <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos",
                    impute = "knn", knn.k = 5L))

# ---------------------------------------------------------------------------
# Simulate 3-trial phenotypes
# ---------------------------------------------------------------------------
set.seed(7)   # seed 2026 gave a draw where corgh + vm aborts in ASReml
M   <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))  # -1/0/1
geno.ids <- rownames(M)
n.g <- nrow(M); p <- ncol(M)
Mc  <- sweep(M, 2L, colMeans(M))
K   <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)

trials <- paste0("T", seq_len(n.t))
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1; D <- diag(sqrt(v)); D %*% S %*% D }
Ga_K <- mkcov(c(0.12, 0.08, 0.15, 0.10, 0.14, 0.09), 0.6)   # K-structured polygenic
Ga_L <- mkcov(c(0.05, 0.05, 0.05), 0.6)                     # IID line polygenic
U <- t(chol(K)) %*% matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_K)
L <- matrix(rnorm(n.g * n.t), n.g) %*% chol(Ga_L)

# Planted QTL: two main, one crossover
qtl.cols <- c(round(p * 0.15), round(p * 0.45), round(p * 0.75))
beta <- rbind(rep( 0.20, n.t),
              rep(-0.18, n.t),
              seq(0.25, -0.25, length.out = n.t))   # crossover
Q <- M[, qtl.cols] %*% beta
cat("Planted QTL markers:", colnames(M)[qtl.cols], "\n")

G <- U + L + Q
rownames(G) <- geno.ids

# 9 extra phenotyped lines with no genotypes (IID genetic only)
extra <- sprintf("NoGeno%02d", 1:9)
G.extra <- matrix(rnorm(9 * n.t), 9) %*% chol(Ga_K + Ga_L)
rownames(G.extra) <- extra
G <- rbind(G, G.extra)

ids <- rownames(G)
pheno <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials,
                     stringsAsFactors = FALSE)
mu <- setNames(rep_len(c(4.5, 4.8, 4.2, 4.6, 4.9, 4.3), n.t), trials)
pheno$yield <- mu[pheno$Trial] +
    G[cbind(match(pheno$id, ids), match(pheno$Trial, trials))] +
    ifelse(pheno$Rep == "2", 0.1, 0) +
    rnorm(nrow(pheno), sd = sqrt(0.15))
pheno$id    <- factor(pheno$id)
pheno$Rep   <- factor(pheno$Rep)
pheno$Trial <- factor(pheno$Trial)
cat(sprintf("pheno rows: %d  (trials=%d, lines=%d, genotyped=%d)\n",
            nrow(pheno), nlevels(pheno$Trial), nlevels(pheno$id), n.g))

# ---------------------------------------------------------------------------
# Base model + gwasAim
# ---------------------------------------------------------------------------
rand.form <- as.formula(if (grepl("\\(", struct)) sprintf("~ %s:id", struct)
                        else sprintf("~ %s(Trial):id", struct))
base <- asreml(yield ~ Trial + Trial:Rep,
               random   = rand.form,
               residual = ~ dsum(~units | Trial),
               na.action = na.method(x = "include"),
               data = pheno, workspace = "1gb", trace = FALSE)
base <- update(base, trace = FALSE)
cat("base converge:", base$converge, "\n")
print(summary(base)$varcomp)

gw <- gwasAim(baseModel = base, genObj = panel, merge.by = "id",
              Trait = "Trial", str = str.arg, TypeI = 0.05,
              trace = FALSE, breakout = 4, maxit = 50)

cat("\nDetected QTL:\n"); print(gw$QTL$qtl)
print(gw$QTL$is.interaction)
print(summary(gw, genObj = panel))
