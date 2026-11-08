# =============================================================================
# scratch/smoke_gwas_mv_new_sim.R
#
# MV smoke test on the 200-line AlphaSimR panel (post-S36 sim rewrite).
#
# What the sim now provides (see data-raw/simulate_GWASsim.R):
#   * Panel  : AlphaSimR-generated, 200 lines x 1500 markers, 10 chr x 150 cM.
#   * K      : realistic full-rank kinship (off-diag SD ~0.086, top-5
#              eigenvalues ~5-6, min ~0).
#   * Polyg. : TWO Kronecker components matching the wgAim MV analysis's
#              fa(Site,k):vm(id, K) + fa(Site,k):id decomposition --
#                U ~ MVN(0, Ga_K (x) K)   K-weighted polygenic
#                L ~ MVN(0, Ga_L (x) I)   IID-line polygenic
#              Ga_K and Ga_L are 6 x 6 Wishart-style with mean diag ~0.07
#              each (total polygenic variance ~0.14, on par with residual),
#              moderate positive cross-site correlations.  Real K-shape in
#              U means vm(id, K) has genome-wide signal at iteration 1 so
#              the initial LRT is expected to fire.
#   * QTL    : 5 planted at fixed cM positions (effects in t/ha) --
#                Chr2 @  45 cM  MAIN         +0.080 all sites
#                Chr5 @  90 cM  MAIN         -0.070 all sites
#                Chr7 @  30 cM  MAIN-weak    +0.055 all sites
#                Chr3 @ 110 cM  GxE-fade     +0.135 -> 0 across sites
#                Chr9 @  60 cM  GxE-cross    -/+ 0.110 crossover
#            Marker names drift across sim runs (coalescent picks
#            different segregating sites); positions are stable.
#   * Yield  : tonnes per hectare (site means 4.5-4.9), residual variance
#              0.142 (SD ~0.377 t/ha, ~8% of the ~4.7 t/ha trait mean).
#
# Notes on `str` for MV GWAS (see session_summary S36):
#   * str = NULL   : mirrors the base model -- with fa(Site,2) here, this
#                    gives per-site psi variances that CAN boundary at 0.
#                    If Ginv[j,j] blows up for any j, all oint go to 0
#                    and forward selection stalls on Chr1.1 spuriously.
#   * str = "corh" : pools to a single variance + single correlation,
#                    no per-site boundary risk.  Fallback if str = NULL
#                    misbehaves.
# =============================================================================

suppressPackageStartupMessages({
    library(asreml)
    asreml.options(trace = TRUE)
    devtools::load_all()
})

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

flt   <- filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos")
panel <- primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos",
                    impute = "knn", knn.k = 5L)
panel.ids <- as.character(panel$pheno$id)

pheno       <- phenoGWASsim
pheno       <- pheno[as.character(pheno$id) %in% panel.ids, ]
pheno       <- droplevels(pheno)
pheno$Site  <- factor(pheno$Site)
pheno$Rep   <- factor(pheno$Rep)
cat(sprintf("MV pheno rows: %d  (sites=%d, lines=%d, reps=%d)\n",
            nrow(pheno), nlevels(pheno$Site),
            nlevels(pheno$id), nlevels(pheno$Rep)))

# --- base model progression: diag -> fa1 -> fa2
n.sites <- nlevels(pheno$Site)

cat("\n--- base_diag ---\n")
base_diag <- asreml(yield ~ Site + Site:Rep,
                    random   = ~ diag(Site):id,
                    residual = ~ dsum(~units | Site),
                    data = pheno, workspace = "1gb")
summary(base_diag)$varcomp

cat("converge base_diag:", base_diag$converge,
    "  logLik:", summary(base_diag)$logl, "\n")

cat("\n--- base_fa1 ---\n")
base_fa1 <- asreml(yield ~ Site + Site:Rep,
                   random   = ~ fa(Site, 1):id,
                   residual = ~ dsum(~units | Site),
                   data = pheno, workspace = "1gb")
cat("converge base_fa1 :", base_fa1$converge,
    "  logLik:", summary(base_fa1)$logl, "\n")

cat("\n--- base_fa2 ---\n")
base_fa2 <- asreml(yield ~ Site + Site:Rep,
                   random   = ~ fa(Site, 2):id,
                   residual = ~ dsum(~units | Site),
                   data = pheno, workspace = "1gb")
cat("converge base_fa2 :", base_fa2$converge,
    "  logLik:", summary(base_fa2)$logl, "\n")

# Parameter counts: n.sites psi + k*n.sites loadings + n.sites resid
np   <- c(diag = 2L*n.sites, fa1 = 3L*n.sites, fa2 = 4L*n.sites)
ll   <- c(summary(base_diag)$logl, summary(base_fa1)$logl,
          summary(base_fa2)$logl)
comp <- data.frame(model = c("diag(Site):id", "fa(Site,1):id", "fa(Site,2):id"),
                   params = np, logLik = round(ll, 2),
                   AIC    = round(-2 * ll + 2 * np, 2))
cat("\n--- base-model comparison ---\n")
print(comp, row.names = FALSE)

# ----------------------------------------------------------------------------
# Diagnostic: fit the initial diag(Site):vm(id, K) + fa(Site,2):id model
# that gwasAim() builds internally in the pre-loop LRT step, and report
# the per-site vm variance components with their ASReml constraint code.
# Constraint codes (asreml-R): 1 = P (positive, ok), 2 = U (unrestricted),
# 3 = B (BOUNDARY -- variance at 0), 4 = F (fixed), 5 = C (constrained).
# If any vm site variance is at the boundary (code 3), MASS::ginv(Ga) will
# inflate the corresponding Ginv[j,j], drive vqtilde <= 0 for every marker,
# and .qtlSelect() will return all-zero oint  -->  which.max picks Chr.Chr1.1.
# This diagnostic makes the failure mode visible BEFORE calling gwasAim().
# ----------------------------------------------------------------------------
cat("\n--- Diagnostic: pre-loop diag(Site):vm(id,K) fit ---\n")

# Build covObj = kinship K exactly as .constructCM() does internally.
M     <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
Mc    <- sweep(M, 2L, colMeans(M), "-")
covObj <- tcrossprod(Mc) / ncol(Mc)
covObj <- covObj / mean(diag(covObj))
ev    <- eigen(covObj, symmetric = TRUE)
ev$values[ev$values < 0] <- 0
eps   <- max(1e-6, 1e-10 * max(ev$values))
covObj <- ev$vectors %*% diag(ev$values + eps) %*% t(ev$vectors)
dimnames(covObj) <- list(rownames(M), rownames(M))
cat(sprintf("covObj: %d x %d, mean diag = %.3f, off-diag SD = %.3f\n",
            nrow(covObj), ncol(covObj),
            mean(diag(covObj)), sd(covObj[lower.tri(covObj)])))

diag_vm <- asreml(fixed    = yield ~ Site + Site:Rep,
                  random   = ~ fa(Site, 1):vm(id, covObj) + fa(Site, 2):id,
                  residual = ~ dsum(~units | Site),
                  data     = pheno, workspace = "1gb")
diag_vm <- update(diag_vm)
summary(diag_vm)$varcomp

vp   <- diag_vm$vparameters
vpc  <- diag_vm$vparameters.con
mask <- grepl("vm\\(.*covObj\\)", names(vp))
if (!any(mask)) mask <- grepl("vm.*covObj", names(vp))
diag_tab <- data.frame(component = names(vp)[mask],
                       estimate  = signif(vp[mask], 4),
                       con_code  = vpc[mask],
                       flag      = ifelse(vpc[mask] == 3L, "BOUNDARY",
                                   ifelse(vp[mask] < 1e-4, "near-zero", "ok")),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
cat("\ndiag(Site):vm(id, K) variance components:\n")
print(diag_tab, row.names = FALSE)
cat(sprintf("\nBoundary flag summary: any-boundary = %s, any-near-zero = %s\n",
            any(diag_tab$flag == "BOUNDARY"),
            any(diag_tab$flag == "near-zero")))
cat("(If TRUE with str='diag' -> Chr1.1 will be spuriously picked.\n",
    " Use str=NULL to upgrade vm to fa(Site,2) [MVsim recipe],\n",
    " or str='corh' as a robust fallback.)\n", sep = "")

cat("\n--- gwasAim(Trait='Site') ---\n")
t0 <- proc.time()
gwas_mv <- gwasAim(
    baseModel        = base_fa2,
    genObj           = panel,
    merge.by         = "id",
    Trait            = "Site",
    str              = NULL,
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = TRUE,
    verboseLev       = 0
)
t1 <- proc.time()
cat(sprintf("\ngwasAim() MV elapsed: %.1f s\n", (t1 - t0)[3]))

# --- Planted-QTL reference (positions are stable; marker names drift) ---
truth <- data.frame(
    chr    = c("Chr2", "Chr5", "Chr7",       "Chr3",     "Chr9"),
    pos    = c(   45,     90,     30,          110,         60),
    type   = c("MAIN", "MAIN", "MAIN-weak", "GxE-fade", "GxE-cross"),
    expect = c("MAIN", "MAIN", "MAIN",      "INT",      "INT"),
    stringsAsFactors = FALSE
)
truth$nearest_marker <- mapply(function(ch, ps) {
    sub <- mapGWASsim[mapGWASsim$chr == ch, , drop = FALSE]
    sub$marker[which.min(abs(sub$pos - ps))]
}, truth$chr, truth$pos, USE.NAMES = FALSE)
cat("\n--- Planted QTL (truth) ---\n")
print(truth, row.names = FALSE)

cat("\n--- Detection ---\n")
cat(sprintf("Detected QTL: %d\n", length(gwas_mv$QTL$qtl)))
if (length(gwas_mv$QTL$qtl) > 0L) {
    print(gwas_mv$QTL$qtl)
    cat("\nis.interaction:\n"); print(gwas_mv$QTL$is.interaction)
    cat("\nSummary:\n")
    print(summary(gwas_mv, genObj = panel))
}
