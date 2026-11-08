# =============================================================================
# scratch/smoke_gwasAim_uv.R
# End-to-end smoke test of the single-environment GWAS pipeline on Env01.
# Confirms the pipeline runs and reports which planted QTL are recovered.
# =============================================================================

suppressPackageStartupMessages({
    devtools::load_all()
    library(asreml)
})
asreml.options(trace = FALSE)

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

# ---- 1.  Panel prep ---------------------------------------------------------
cat("\n--------- 1. checkPanel ------------------------------------------------\n")
chk <- checkPanel(genoGWASraw, mapGWASsim, id = "id")
print(chk)

cat("\n--------- 2. filterPanel -----------------------------------------------\n")
flt <- filterPanel(chk)
print(flt)

cat("\n--------- 3. primePanel (knn) ------------------------------------------\n")
panel <- primePanel(flt, impute = "knn", knn.k = 5)

# ---- 2.  Env01 phenotype subset ---------------------------------------------
env1 <- droplevels(phenoGWASsim[phenoGWASsim$Site == "Env01", ])
env1$Row    <- factor(env1$Row)
env1$Column <- factor(env1$Column)
env1$Rep    <- factor(env1$Rep)
cat(sprintf("\nEnv01 subset: %d rows, %d unique lines, %d reps\n",
            nrow(env1), length(unique(env1$id)), length(unique(env1$Rep))))

# ---- 3.  Base model ---------------------------------------------------------
cat("\n--------- 4. Base ASReml model -----------------------------------------\n")
base <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ id + Row + Column,
    data      = env1,
    na.action = na.method(x = "include"),
    maxit     = 60
)
base <- update(base)
cat(sprintf("Base model converged: %s   logL = %.3f\n",
            base$converge, base$loglik))
print(summary(base)$varcomp)

# ---- 4.  gwasAim ------------------------------------------------------------
cat("\n--------- 5. gwasAim ---------------------------------------------------\n")
gwas <- gwasAim(
    baseModel        = base,
    genObj           = panel,
    merge.by         = "id",
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = FALSE,
    na.action        = na.method(x = "include"),
    maxit            = 60
)

# ---- 5.  Results ------------------------------------------------------------
cat("\n--------- 6. Results ---------------------------------------------------\n")
cat("class :", paste(class(gwas), collapse = ", "), "\n")
cat("nQTL  :", length(gwas$QTL$qtl), "\n")
cat("keys  :\n"); print(gwas$QTL$qtl)

cat("\n--------- 7. Print / Summary -------------------------------------------\n")
print(gwas, genObj = panel)
cat("\nSummary table:\n")
print(summary(gwas, genObj = panel))

# ---- 6.  Ground-truth check -------------------------------------------------
cat("\n--------- 8. Planted-QTL recovery check --------------------------------\n")
truth <- data.frame(
    marker = c("SNP_C02_0095", "SNP_C05_0183", "SNP_C07_0060",
                "SNP_C03_0214", "SNP_C09_0123"),
    type   = c("MAIN-strong", "MAIN-neg", "MAIN-weak",
                "GxE-fade (env1: +150)", "GxE-cross (env1: -120)"),
    stringsAsFactors = FALSE
)
detected.markers <- summary(gwas, genObj = panel)$Marker
truth$detected <- truth$marker %in% detected.markers
print(truth)
