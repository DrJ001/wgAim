# =============================================================================
# scratch/smoke_gwasAim_uv2.R
# Second smoke test -- upgrade base model to AR1(Row):AR1(Column) spatial to
# absorb the biomAid spatial trend and boost detection power.
# =============================================================================

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

# ---- Inspect Env01 layout ---------------------------------------------------
cat(sprintf("Row range:    %d - %d  (unique: %d)\n",
            min(env1$Row), max(env1$Row), length(unique(env1$Row))))
cat(sprintf("Column range: %d - %d  (unique: %d)\n",
            min(env1$Column), max(env1$Column), length(unique(env1$Column))))
cat(sprintf("Total plots:  %d  (should be %d)\n",
            nrow(env1), length(unique(env1$Row)) * length(unique(env1$Column))))

# ---- Prep for AR1 (factors + proper sort order) -----------------------------
env1$Rowf <- factor(env1$Row)
env1$Colf <- factor(env1$Column)
env1$Rep  <- factor(env1$Rep)
env1      <- env1[order(env1$Rowf, env1$Colf), ]

# ---- AR1(Row):AR1(Col) base model ------------------------------------------
cat("\nFitting AR1 spatial base model ...\n")
base <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ id,
    residual  = ~ ar1v(Rowf):ar1(Colf),
    data      = env1,
    na.action = na.method(x = "include"),
    maxit     = 80
)
base <- update(base)
base <- update(base)
cat(sprintf("Base converged: %s   logL = %.3f\n", base$converge, base$loglik))
print(summary(base)$varcomp)

# ---- gwasAim ----------------------------------------------------------------
cat("\nRunning gwasAim ...\n")
gwas <- gwasAim(
    baseModel        = base,
    genObj           = panel,
    merge.by         = "id",
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = FALSE,
    na.action        = na.method(x = "include"),
    maxit            = 80
)

cat("\nnQTL:", length(gwas$QTL$qtl), "\n")
print(gwas, genObj = panel)
cat("\nSummary:\n"); print(summary(gwas, genObj = panel))

# Ground-truth recovery check
truth <- data.frame(
    marker = c("SNP_C02_0095", "SNP_C05_0183", "SNP_C07_0060",
                "SNP_C03_0214", "SNP_C09_0123"),
    type   = c("MAIN-strong (+100)", "MAIN-neg (-80)", "MAIN-weak (+50)",
                "GxE-fade (Env01: +150)", "GxE-cross (Env01: -120)"),
    stringsAsFactors = FALSE
)
detected <- if (!is.null(gwas$QTL$qtl))
    summary(gwas, genObj = panel)$Marker else character(0)
truth$detected <- truth$marker %in% detected
cat("\nRecovery:\n"); print(truth)
