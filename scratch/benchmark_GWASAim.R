# =============================================================================
# benchmark_gwasAim.R
#
# Benchmarks gwasAim.asreml() against wgaim.asreml(gen.type = "marker")
# using the RxK wheat dataset. Both use the same marker data (imputed ±1
# scores from the biparental cross) and should produce identical QTL
# selections, effect sizes, log-likelihoods, and LRT values.
#
# Key mapping between old and new:
#   wgaim(gen.type = "marker")  <->  gwasAim()
#   intervalObj$imputed.data    <->  wgPanel$imputed.data (encoding = "pm1")
#   method  = "fixed"  [default in wgaim, hardcoded in gwasAim]
#   selection = "interval" [default in wgaim, hardcoded in gwasAim]
#
# Run from the wgaim project root:
#   Rscript scratch/benchmark_gwasAim.R
# =============================================================================

library(asreml)
library(wgaim)

cat("\n============================================================\n")
cat(" wgAim Benchmark: gwasAim.asreml() vs wgaim.asreml(marker)\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------
data(phenoRxK, package = "wgaim")
data(genoRxK,  package = "wgaim")

genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))

# Old wgaim reference object (cross2int keeps imputed.data for marker mode)
genoRxK_old <- cross2int(genoRxK, impute = "Martinez", id = "Genotype")

cat("Data loaded:\n")
cat("  Phenotypic lines :", nrow(phenoRxK), "\n")
cat("  Chromosomes      :", length(genoRxK$geno), "\n")
cat("  Total markers    :",
    sum(sapply(genoRxK_old$geno, function(el) ncol(el$imputed.data))), "\n\n")

# -----------------------------------------------------------------------------
# 2. Fit base model
# -----------------------------------------------------------------------------
cat("Fitting base model...\n")
base.mod <- asreml(yld ~ Type + lrow,
                   random   = ~ Genotype + Range,
                   residual = ~ ar1(Range):ar1(Row),
                   data     = phenoRxK)

# -----------------------------------------------------------------------------
# 3. Run wgaim.asreml(gen.type = "marker") — reference
# -----------------------------------------------------------------------------
cat("\n--- Running wgaim.asreml(gen.type='marker') [reference] ---\n")
t_wgaim <- system.time({
    rkyld.wgaim <- wgaim(base.mod,
                          intervalObj = genoRxK_old,
                          merge.by    = "Genotype",
                          gen.type    = "marker",
                          na.action   = na.method(x = "include"))
})

# -----------------------------------------------------------------------------
# 4. Build wgPanel from the same marker data and run gwasAim.asreml()
#
# The old wgaim uses imputed.data (±1 encoding) directly. We pass the same
# matrix to primePanel with encoding = "pm1" so no transformation is applied.
# maf = 0 ensures no markers are dropped that wgaim would have kept.
# -----------------------------------------------------------------------------
devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                   quiet = TRUE)

data(phenoRxK, package = "wgaim")   # reload clean copy

cat("\n--- Building wgPanel from marker data ---\n")

# Extract imputed marker matrix (lines x markers, ±1 encoding) and map
gdat <- lapply(genoRxK_old$geno, function(el) el$imputed.data)
marker.mat <- do.call(cbind, gdat)

map.list <- lapply(names(genoRxK_old$geno), function(chr) {
    m <- genoRxK_old$geno[[chr]]$map
    data.frame(marker = names(m),
               chr    = chr,
               pos    = as.numeric(m),
               stringsAsFactors = FALSE)
})
map.df <- do.call(rbind, map.list)

# primePanel with encoding = "pm1": no transformation, no MAF filter
panel <- primePanel(geno     = marker.mat,
                    map      = map.df,
                    id       = "Genotype",
                    map.id   = "marker",
                    map.chr  = "chr",
                    map.pos  = "pos",
                    encoding = "pm1",
                    maf      = 0)

cat("\n--- Running gwasAim.asreml() [new] ---\n")
base.mod2 <- asreml(yld ~ Type + lrow,
                    random   = ~ Genotype + Range,
                    residual = ~ ar1(Range):ar1(Row),
                    data     = phenoRxK)

t_gwasaim <- system.time({
    rkyld.gwasaim <- gwasAim(base.mod2,
                              genObj   = panel,
                              merge.by = "Genotype",
                              na.action = na.method(x = "include"))
})

# -----------------------------------------------------------------------------
# 5. Compare results
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat(" Comparison Results\n")
cat("============================================================\n\n")

# --- Timing ---
cat("Timing:\n")
cat(sprintf("  wgaim.asreml  : %.2f sec\n", t_wgaim["elapsed"]))
cat(sprintf("  gwasAim.asreml: %.2f sec\n", t_gwasaim["elapsed"]))

# --- QTL/markers detected ---
cat("\nMarkers detected:\n")
cat("  wgaim   :", paste(rkyld.wgaim$QTL$qtl,    collapse = ", "), "\n")
cat("  gwasAim :", paste(rkyld.gwasaim$QTL$qtl,  collapse = ", "), "\n")
qtl_match <- identical(rkyld.wgaim$QTL$qtl, rkyld.gwasaim$QTL$qtl)
cat("  Match   :", ifelse(qtl_match, "YES", "NO"), "\n")

# --- Effect sizes ---
cat("\nEffect sizes:\n")
cat("  wgaim   :", round(rkyld.wgaim$QTL$effects,   4), "\n")
cat("  gwasAim :", round(rkyld.gwasaim$QTL$effects, 4), "\n")
eff_match <- isTRUE(all.equal(unname(rkyld.wgaim$QTL$effects),
                               unname(rkyld.gwasaim$QTL$effects),
                               tolerance = 1e-3))
cat("  Match   :", ifelse(eff_match, "YES", "NO"), "\n")

# --- Final log-likelihood ---
cat("\nFinal model log-likelihood:\n")
cat(sprintf("  wgaim   : %.6f\n", rkyld.wgaim$loglik))
cat(sprintf("  gwasAim : %.6f\n", rkyld.gwasaim$loglik))
ll_match <- isTRUE(all.equal(rkyld.wgaim$loglik, rkyld.gwasaim$loglik,
                              tolerance = 1e-4))
cat("  Match   :", ifelse(ll_match, "YES", "NO"), "\n")

# --- LRT values per iteration ---
cat("\nLRT values per iteration:\n")
lik_w <- do.call(rbind, rkyld.wgaim$QTL$diag$lik)
lik_q <- do.call(rbind, rkyld.gwasaim$QTL$diag$lik)
colnames(lik_w) <- colnames(lik_q) <- c("baseLogL", "qtlLogL", "stat", "pvalue")
cat("  wgaim:\n");   print(round(lik_w, 4))
cat("  gwasAim:\n"); print(round(lik_q, 4))
lrt_match <- isTRUE(all.equal(lik_w, lik_q, tolerance = 1e-3,
                               check.attributes = FALSE))
cat("  Match   :", ifelse(lrt_match, "YES", "NO"), "\n")

# --- Iterations ---
cat("\nIterations:\n")
cat("  wgaim   :", rkyld.wgaim$QTL$iterations,   "\n")
cat("  gwasAim :", rkyld.gwasaim$QTL$iterations, "\n")

# --- Summary tables ---
cat("\n--- summary(wgaim, gen.type='marker') ---\n")
print(summary(rkyld.wgaim, intervalObj = genoRxK_old))

cat("\n--- summary(gwasAim) ---\n")
print(summary(rkyld.gwasaim, genObj = panel))

# --- Overall verdict ---
cat("\n============================================================\n")
all_match <- qtl_match & eff_match & ll_match & lrt_match
cat(sprintf(" Overall: %s\n", ifelse(all_match,
    "ALL CHECKS PASSED - gwasAim matches wgaim(marker)",
    "DIFFERENCES DETECTED - check output above")))
cat("============================================================\n\n")
