# =============================================================================
# benchmark_QTLAim.R
#
# Benchmarks qtlAim.asreml() against wgaim.asreml() using the RxK wheat
# dataset from the package man page example. Both functions should produce
# identical QTL selections, effect sizes, log-likelihoods, and diagnostics.
#
# Run from the wgaim project root:
#   Rscript scratch/benchmark_QTLaim.R
# =============================================================================

library(asreml)
library(wgaim)

cat("\n============================================================\n")
cat(" wgAim Benchmark: qtlAim.asreml() vs wgaim.asreml()\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
data(phenoRxK, package = "wgaim")
data(genoRxK,  package = "wgaim")

# Subset to four chromosomes (as per man page example)
genoRxK <- subset(genoRxK, chr = c("1A", "2D1", "2D2", "3B"))

# Old wgaim uses cross2int(); new wgAim uses primeCross()
genoRxK_old <- cross2int(genoRxK, impute = "Martinez", id = "Genotype")

cat("Data loaded:\n")
cat("  Phenotypic lines :", nrow(phenoRxK), "\n")
cat("  Chromosomes      :", length(genoRxK$geno), "\n")
cat("  Total intervals  :",
    sum(sapply(genoRxK_old$geno, function(el) ncol(el$interval.data))), "\n\n")

# -----------------------------------------------------------------------------
# 2. Fit base model (shared by both analyses)
# -----------------------------------------------------------------------------
cat("Fitting base model...\n")
base.mod <- asreml(yld ~ Type + lrow,
                   random   = ~ Genotype + Range,
                   residual = ~ ar1(Range):ar1(Row),
                   data     = phenoRxK)

# -----------------------------------------------------------------------------
# 3. Run wgaim.asreml() — the reference
# -----------------------------------------------------------------------------
cat("\n--- Running wgaim.asreml() [reference] ---\n")
t_wgaim <- system.time({
    rkyld.wgaim <- wgaim(base.mod,
                          intervalObj = genoRxK_old,
                          merge.by    = "Genotype",
                          na.action   = na.method(x = "include"))
})

# -----------------------------------------------------------------------------
# 4. Run qtlAim.asreml() — the new function
# Reload package and data to ensure a clean environment.
# -----------------------------------------------------------------------------
devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                   quiet = TRUE)

data(phenoRxK, package = "wgaim")   # reload clean copy

cat("\n--- Running qtlAim.asreml() [new] ---\n")
base.mod2 <- asreml(yld ~ Type + lrow,
                    random   = ~ Genotype + Range,
                    residual = ~ ar1(Range):ar1(Row),
                    data     = phenoRxK)

genoRxKa <- primeCross(genoRxK, impute = "Martinez", id = "Genotype")

t_qtlaim <- system.time({
    rkyld.qtlaim <- qtlAim(base.mod2,
                            genObj   = genoRxKa,
                            merge.by = "Genotype",
                            gen.type = "interval",
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
cat(sprintf("  qtlAim.asreml : %.2f sec\n", t_qtlaim["elapsed"]))

# --- QTL detected ---
cat("\nQTL detected:\n")
cat("  wgaim  :", paste(rkyld.wgaim$QTL$qtl,  collapse = ", "), "\n")
cat("  qtlAim :", paste(rkyld.qtlaim$QTL$qtl, collapse = ", "), "\n")
qtl_match <- identical(rkyld.wgaim$QTL$qtl, rkyld.qtlaim$QTL$qtl)
cat("  Match  :", ifelse(qtl_match, "YES", "NO"), "\n")

# --- Effect sizes ---
cat("\nEffect sizes:\n")
cat("  wgaim  :", round(rkyld.wgaim$QTL$effects,  4), "\n")
cat("  qtlAim :", round(rkyld.qtlaim$QTL$effects, 4), "\n")
eff_match <- isTRUE(all.equal(unname(rkyld.wgaim$QTL$effects),
                               unname(rkyld.qtlaim$QTL$effects),
                               tolerance = 1e-3))
cat("  Match  :", ifelse(eff_match, "YES", "NO"), "\n")

# --- Log-likelihoods (final model) ---
cat("\nFinal model log-likelihood:\n")
cat(sprintf("  wgaim  : %.6f\n", rkyld.wgaim$loglik))
cat(sprintf("  qtlAim : %.6f\n", rkyld.qtlaim$loglik))
ll_match <- isTRUE(all.equal(rkyld.wgaim$loglik, rkyld.qtlaim$loglik,
                              tolerance = 1e-4))
cat("  Match  :", ifelse(ll_match, "YES", "NO"), "\n")

# --- LRT matrix ---
cat("\nLRT values per iteration:\n")
lik_w <- do.call(rbind, rkyld.wgaim$QTL$diag$lik)
lik_q <- do.call(rbind, rkyld.qtlaim$QTL$diag$lik)
colnames(lik_w) <- colnames(lik_q) <- c("baseLogL", "qtlLogL", "stat", "pvalue")
cat("  wgaim:\n");  print(round(lik_w, 4))
cat("  qtlAim:\n"); print(round(lik_q, 4))
lrt_match <- isTRUE(all.equal(lik_w, lik_q, tolerance = 1e-3,
                               check.attributes = FALSE))
cat("  Match  :", ifelse(lrt_match, "YES", "NO"), "\n")

# --- Iterations ---
cat("\nIterations:\n")
cat("  wgaim  :", rkyld.wgaim$QTL$iterations,  "\n")
cat("  qtlAim :", rkyld.qtlaim$QTL$iterations, "\n")

# --- Summary output ---
cat("\n--- summary(wgaim) ---\n")
print(summary(rkyld.wgaim, intervalObj = genoRxK_old))

cat("\n--- summary(qtlAim) ---\n")
print(summary(rkyld.qtlaim, genObj = genoRxKa))

# --- Overall verdict ---
cat("\n============================================================\n")
all_match <- qtl_match & eff_match & ll_match & lrt_match
cat(sprintf(" Overall: %s\n", ifelse(all_match,
    "ALL CHECKS PASSED - qtlAim matches wgaim",
    "DIFFERENCES DETECTED - check output above")))
cat("============================================================\n\n")
