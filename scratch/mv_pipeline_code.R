## ---- setup ----
asreml_ok <- requireNamespace("asreml", quietly = TRUE)
if (requireNamespace("wgAim", quietly = TRUE)) {
  library(wgAim)
} else {
  devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                     quiet = TRUE)
}
library(ggplot2)
if (asreml_ok) {
  library(asreml)
  asreml.options(trace = FALSE)
}

## ---- load_data ----
data(genoMVsim,  package = "wgAim")
data(phenoMVsim, package = "wgAim")

## ---- primeCross ----
genObj <- primeCross(genoMVsim, type = "interval", id = "id")

## ---- map_summary ----
n_chr <- length(genObj$geno)
n_int <- sum(sapply(genObj$geno, function(ch) ncol(ch$interval.data)))
n_mk  <- sum(sapply(genObj$geno, function(ch) ncol(ch$data)))
cat(sprintf("Chromosomes : %d\n", n_chr))
cat(sprintf("Markers     : %d\n", n_mk))
cat(sprintf("Intervals   : %d\n", n_int))
cat(sprintf("Lines       : %d\n", nrow(genObj$pheno)))

## ---- pheno_setup ----
phenoMVsim$id    <- factor(phenoMVsim$id)
phenoMVsim$Trial <- factor(phenoMVsim$Trial,
                           levels = paste0("Trial", 1:8))
phenoMVsim$Rep   <- factor(phenoMVsim$Rep)

## ---- base_diag ----
base_diag <- asreml(
  fixed     = yield ~ Trial + Trial:Rep,
  random    = ~ diag(Trial):id,
  residual  = ~ dsum(~units | Trial),
  data      = phenoMVsim,
  na.action = na.method(x = "include"),
  maxit     = 80
)
base_diag <- update(base_diag)

## ---- base_fa1 ----
base_fa1 <- asreml(
  fixed     = yield ~ Trial + Trial:Rep,
  random    = ~ fa(Trial, 1):id,
  residual  = ~ dsum(~units | Trial),
  data      = phenoMVsim,
  na.action = na.method(x = "include"),
  maxit     = 80
)
base_fa1 <- update(base_fa1)

## ---- base_fa2 ----
base_fa2 <- asreml(
  fixed     = yield ~ Trial + Trial:Rep,
  random    = ~ fa(Trial, 2):id,
  residual  = ~ dsum(~units | Trial),
  data      = phenoMVsim,
  na.action = na.method(x = "include"),
  maxit     = 80
)
base_fa2 <- update(base_fa2)

fasum <- fa.asreml(base_fa2, trunc.char = NULL)
fasum$gammas[[1]]

## ---- model_comparison ----
aic <- function(mod, k) -2 * mod$loglik + 2 * k
# Parameters: diag = 8 var; fa1 = 8 psi + 8 load = 16; fa2 = 8 psi + 16 load = 24
mod_comp <- data.frame(
  Model      = c("diag(Trial)", "fa(Trial, 1)", "fa(Trial, 2)"),
  Parameters = c(8L, 16L, 24L),
  LogLik     = round(c(base_diag$loglik, base_fa1$loglik, base_fa2$loglik), 2),
  AIC        = round(c(aic(base_diag, 8), aic(base_fa1, 16), aic(base_fa2, 24)), 2)
)
print(mod_comp, row.names = FALSE)

## ---- qtlAim ----
qtl_mv <- qtlAim(
  baseModel        = base_fa2,
  genObj           = genObj,
  merge.by         = "id",
  Trait            = "Trial",
  exclusion.window = 20,
  TypeI            = 0.01,
  trace            = FALSE,
  na.action        = na.method(x = "include"),
  maxit            = 80
)

summary.asreml(qtl_mv, vparameters = TRUE)$vparameters

## ---- qtlAim_class ----
cat("Class        :", paste(class(qtl_mv), collapse = ", "), "\n")
cat("Engine       :", qtl_mv$QTL$path, "\n")
cat("QTL detected :", length(qtl_mv$QTL$qtl), "\n")
cat("QTL keys     :", paste(qtl_mv$QTL$qtl, collapse = ", "), "\n")

## ---- summary_qtl ----
summ_mv <- summary(qtl_mv, genObj = genObj)
print(summ_mv)

## ---- plot_outlier ----
plot(qtl_mv, genObj = genObj, type = "outlier", sig.col = "firebrick")

## ---- plot_blups ----
plot(qtl_mv, genObj = genObj, type = "blups", sig.col = "firebrick")

## ---- plot_effects ----
plot(qtl_mv, genObj = genObj, type = "effects")

## ---- plot_contrast ----
plot(qtl_mv, genObj = genObj, type = "contrast", data = phenoData, ncol = 3)

## ---- aimtrace_lrt ----
aimTrace(qtl_mv, plot = "lrt", print.out = FALSE)

## ---- aimtrace_stability ----
aimTrace(qtl_mv, plot = "stability", print.out = FALSE)

