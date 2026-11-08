# =============================================================================
# scratch/smoke_gwas_new_sim.R
#
# UV smoke test on the Env1 subset of the S37-calibrated 200-line panel.
# The MV counterpart is scratch/smoke_gwas_mv_new_sim.R.
#
# What the sim provides at Env1 (see data-raw/simulate_GWASsim.R):
#   * Panel  : 200 lines x ~1691 markers (all F1 MAF > 0.15, median ~0.32).
#   * Base   : yield ~ Rep + id + Row + Column (id random, Row/Column
#              random for good measure though the sim has no row/column
#              trend so their variance components should sit near zero).
#   * QTL (effects at Env1 -- ALL five have signal at this site) :
#              Chr2 @  45 cM  MAIN         +0.145  (constant across sites)
#              Chr5 @  90 cM  MAIN         -0.145  (constant across sites)
#              Chr7 @  30 cM  MAIN-weak    +0.125  (constant across sites)
#              Chr3 @ 110 cM  GxE-fade     +0.230  (PEAK at Env1, fading)
#              Chr9 @  60 cM  GxE-cross    -0.180  (Env1..3 -ve; Env4..6 +ve)
#            Marker names drift across sim regenerations; positions stable.
#   * Yield  : tonnes per hectare, Env1 mean ~4.6, residual SD ~0.377.
#
# Expected univariate detection at Env1: all 5 QTL should land (Env1 is
# actually the STRONGEST environment for both G x E QTL, since it's the
# peak of the fade and one end of the crossover).
#
# Prints:
#   - initial LRT statistic (should be well above 3.84 for chi-sq(1))
#   - detected QTL count + summary table
#   - truth / detection comparison at the end
suppressPackageStartupMessages({
    library(asreml)
    asreml.options(trace = TRUE)
    devtools::load_all()
})

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

cat("\n----- Panel prep -----\n")
flt   <- filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos")
print(flt)
panel <- primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos",
                    impute = "knn", knn.k = 5L)
panel.ids <- as.character(panel$pheno$id)
cat(sprintf("Panel: %d lines x %d markers after primePanel()\n",
            length(panel.ids),
            sum(vapply(panel$geno, function(g) ncol(g$imputed.data), integer(1)))))

cat("\n----- Env1 base model -----\n")
env1 <- droplevels(phenoGWASsim[phenoGWASsim$Site == "Env1", ])
# For the smoke test, restrict pheno to genotyped lines to sidestep the
# .fixLines() path (whose interaction with `random = ~ id + Row + Column`
# may or may not be an existing bug -- investigating separately).
env1 <- env1[as.character(env1$id) %in% panel.ids, ]
env1 <- droplevels(env1)
env1$Row    <- factor(env1$Row)
env1$Column <- factor(env1$Column)
env1$Rep    <- factor(env1$Rep)
cat(sprintf("env1 rows after restricting to genotyped lines: %d\n", nrow(env1)))

base_env1 <- asreml(
    fixed    = yield ~ Rep,
    random   = ~ id + Row + Column,
    data     = env1,
    workspace = "512mb"
)
base_env1 <- update(base_env1)
cat("Base varcomp:\n")
print(summary(base_env1)$varcomp)

cat("\n----- gwasAim() -----\n")
t0 <- proc.time()
gwas_env1 <- gwasAim(
    baseModel        = base_env1,
    genObj           = panel,
    merge.by         = "id",
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = FALSE,
    verboseLev       = 0
)
t1 <- proc.time()
cat(sprintf("\ngwasAim() elapsed: %.1f s\n", (t1 - t0)[3]))

cat("\n----- Detection summary -----\n")
cat(sprintf("Detected QTL: %d\n", length(gwas_env1$QTL$qtl)))
if (length(gwas_env1$QTL$qtl) > 0L) {
    cat("QTL keys:\n")
    print(gwas_env1$QTL$qtl)
    cat("\nSummary table:\n")
    print(summary(gwas_env1, genObj = panel))
}

# --- Planted-QTL reference (positions stable; marker names drift; Env1 effs) ---
truth <- data.frame(
    chr           = c("Chr2", "Chr5", "Chr7",       "Chr3",     "Chr9"),
    pos           = c(   45,     90,     30,          110,         60),
    type          = c("MAIN", "MAIN", "MAIN-weak", "GxE-fade", "GxE-cross"),
    env1_effect   = c( 0.145, -0.145,  0.125,        0.230,     -0.180),
    stringsAsFactors = FALSE
)
truth$nearest_marker <- mapply(function(ch, ps) {
    sub <- mapGWASsim[mapGWASsim$chr == ch, , drop = FALSE]
    sub$marker[which.min(abs(sub$pos - ps))]
}, truth$chr, truth$pos, USE.NAMES = FALSE)
cat("\n----- Planted QTL (truth, Env1 effects) -----\n")
print(truth, row.names = FALSE)

cat("\n----- LRT diagnostics -----\n")
if (!is.null(gwas_env1$QTL$diag$lik.list)) {
    cat("LRT trace across iterations:\n")
    print(gwas_env1$QTL$diag$lik.list)
}
