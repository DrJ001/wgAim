# Third smoke test -- try mbf path (force=TRUE) plus simple base model.
suppressPackageStartupMessages({
    devtools::load_all()
    library(asreml)
})
asreml.options(trace = FALSE)

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

chk   <- checkPanel(genoGWASraw, mapGWASsim, id = "id")
flt   <- filterPanel(chk)
panel <- primePanel(flt, impute = "knn", knn.k = 5)

env1 <- droplevels(phenoGWASsim[phenoGWASsim$Site == "Env01", ])
env1$Row    <- factor(env1$Row)
env1$Column <- factor(env1$Column)
env1$Rep    <- factor(env1$Rep)

base <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ id + Row + Column,
    data      = env1,
    na.action = na.method(x = "include"),
    maxit     = 60
)
base <- update(base)

cat("\n--- Try force = TRUE (mbf path) ---\n")
gwas.mbf <- gwasAim(
    baseModel        = base,
    genObj           = panel,
    merge.by         = "id",
    force            = TRUE,
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = FALSE,
    na.action        = na.method(x = "include"),
    maxit            = 80
)
cat("\nmbf path detected:", length(gwas.mbf$QTL$qtl), "QTL\n")
if (!is.null(gwas.mbf$QTL$qtl)) {
    print(gwas.mbf, genObj = panel)
    print(summary(gwas.mbf, genObj = panel))
}

truth <- c("SNP_C02_0095", "SNP_C05_0183", "SNP_C07_0060",
            "SNP_C03_0214", "SNP_C09_0123")
det <- if (!is.null(gwas.mbf$QTL$qtl))
    summary(gwas.mbf, genObj = panel)$Marker else character(0)
cat("\nRecovered planted:\n")
for (m in truth) cat(sprintf("  %s: %s\n", m, m %in% det))
