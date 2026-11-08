# =============================================================================
# demo_selIndex.R
# Selection Index demo for selIndex() — works directly from the 8-trial
# gpAim object produced by demo_gpAim_mv.R (Section 3).
#
# Run demo_gpAim_mv.R first (or at least Sections 0 and 3) so that gp8
# is in your workspace.
#
# Sections:
#   1 — Weighted index        :  direct economic weights
#   2 — Smith-Hazel index     :  optimal combining weights via P^{-1} Ga w
#   3 — Desired-gains index   :  Pesek-Baker index via Ga^{-1} d
#   4 — Univariate edge case  :  selIndex on univariate gp.uni
#   5 — Comparison table      :  rankings from all three MV indices
#   6 — Plots                 :  all three plot types for each index
# =============================================================================

devtools::load_all()

outdir <- "scratch/plots_demo_gp"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 8, h = 6) {
    ggplot2::ggsave(file.path(outdir, paste0(name, ".png")),
                    plot = p, width = w, height = h, dpi = 150)
    cat("  saved:", file.path(outdir, paste0(name, ".png")), "\n")
    invisible(p)
}

# Confirm gp8 is available
if (!exists("gp8"))
    stop("gp8 not found. Please run demo_gpAim_mv.R (Sections 0 & 3) first.")

trials <- gp8$GP$trait.levels
cat(sprintf("Using gp8: %d lines x %d trials\n",
            length(unique(gp8$GP$gebv$id)), length(trials)))
cat("Trials:", paste(trials, collapse = ", "), "\n\n")

# =============================================================================
# Section 1 — Weighted index: equal weights & custom economic weights
# =============================================================================
cat("=== Section 1: Weighted index ===\n\n")

# 1a: Equal weights (default -- NULL weights)
si.eq <- selIndex(gp8)
cat("--- 1a: Equal weights (default) ---\n")
print(si.eq)

cat("\n--- summary (top 10%) ---\n")
print(summary(si.eq, prop.select = 0.10))

# 1b: Custom economic weights
#     Suppose Trial1 and Trial2 are high-value commercial environments (w=3),
#     Trials 3-6 are secondary (w=1.5), Trials 7-8 are exploratory (w=0.5)
w.econ <- setNames(
    c(3, 3, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5),
    trials
)
cat("\n--- 1b: Custom economic weights ---\n")
cat("Weights:", paste(sprintf("%s=%.1f", names(w.econ), w.econ), collapse = "  "), "\n\n")

si.wt <- selIndex(gp8, weights = w.econ)
print(si.wt)

cat("\n--- summary (top 10%) ---\n")
print(summary(si.wt, prop.select = 0.10))

# =============================================================================
# Section 2 — Smith-Hazel index
#   Same economic weights as 1b but optimal combining weights derived via
#   b = P^{-1} Ga w.  Compare b to the raw economic weights.
# =============================================================================
cat("\n=== Section 2: Smith-Hazel index ===\n\n")

si.sh <- selIndex(gp8, weights = w.econ, type = "smith-hazel")
print(si.sh)

cat("\n--- summary (top 10%) ---\n")
print(summary(si.sh, prop.select = 0.10))

# Compare economic weights vs Smith-Hazel derived weights
cat("\nWeight comparison: economic vs Smith-Hazel\n")
w.comp <- data.frame(
    Trial    = trials,
    Economic = round(w.econ, 3),
    SH       = round(si.sh$weights, 3),
    row.names = NULL
)
print(w.comp)

# =============================================================================
# Section 3 — Desired-gains index (Pesek-Baker)
#   Breeder specifies proportional gains wanted per trial.
#   e.g. "I want the biggest gain in Trials 1 & 2; modest in 3-6; Trial 7 & 8
#         are OK as-is."
# =============================================================================
cat("\n=== Section 3: Desired-gains index ===\n\n")

d.gains <- setNames(
    c(0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.1, 0.1),
    trials
)
cat("Desired gains:", paste(sprintf("%s=%.1f", names(d.gains), d.gains),
                            collapse = "  "), "\n\n")

si.dg <- selIndex(gp8, desired = d.gains, type = "desired-gains")
print(si.dg)

cat("\n--- summary (top 10%) ---\n")
print(summary(si.dg, prop.select = 0.10))

# =============================================================================
# Section 4 — Univariate edge case
# =============================================================================
if (exists("gp.uni")) {
    cat("\n=== Section 4: Univariate edge case ===\n\n")
    si.uni <- selIndex(gp.uni)
    cat("--- print ---\n")
    print(si.uni)
    cat("\n--- summary ---\n")
    print(summary(si.uni, prop.select = 0.10))
} else {
    cat("\n=== Section 4: Skipped (gp.uni not in workspace) ===\n")
}

# =============================================================================
# Section 5 — Ranking comparison: all three MV indices
#   Show how changing the index type re-orders the top candidates.
# =============================================================================
cat("\n=== Section 5: Ranking comparison ===\n\n")

gterm <- gp8$GP$genetic.term
top.n <- 15

# Extract top-15 line IDs from each index
top.eq <- head(si.eq$index[[gterm]], top.n)
top.wt <- head(si.wt$index[[gterm]], top.n)
top.sh <- head(si.sh$index[[gterm]], top.n)
top.dg <- head(si.dg$index[[gterm]], top.n)

all.lines <- unique(c(top.eq, top.wt, top.sh, top.dg))

rank.comp <- data.frame(
    Line    = all.lines,
    Equal   = match(all.lines, si.eq$index[[gterm]]),
    Wtd     = match(all.lines, si.wt$index[[gterm]]),
    SmithH  = match(all.lines, si.sh$index[[gterm]]),
    DesGain = match(all.lines, si.dg$index[[gterm]]),
    stringsAsFactors = FALSE
)
# Sort by equal-weight rank (primary ranking)
rank.comp <- rank.comp[order(rank.comp$Equal), ]
rownames(rank.comp) <- NULL

cat("Top", top.n, "lines across all indices (rank, lower = better):\n")
print(rank.comp)

# Pairwise Spearman correlations between full index rankings
all.ranks <- data.frame(
    Equal   = si.eq$index$Rank,
    Wtd     = match(si.eq$index[[gterm]], si.wt$index[[gterm]]),
    SmithH  = match(si.eq$index[[gterm]], si.sh$index[[gterm]]),
    DesGain = match(si.eq$index[[gterm]], si.dg$index[[gterm]])
)
cat("\nSpearman rank correlations between indices:\n")
print(round(cor(all.ranks, method = "spearman"), 3))

# =============================================================================
# Section 6 — Plots: all three types for each index
# =============================================================================
cat("\n=== Section 6: Plots ===\n")

# ---- Equal weights ----
cat("\n--- Equal weights ---\n")
save_plot(
    plot(si.eq, type = "index", prop.select = 0.1),
    "si_eq_index", w = 8, h = 6)

save_plot(
    plot(si.eq, type = "heatmap", prop.select = 0.1),
    "si_eq_heatmap", w = 9, h = 8)

save_plot(
    plot(si.eq, type = "weights"),
    "si_eq_weights", w = 6, h = 5)

# ---- Custom economic weights ----
cat("\n--- Custom economic weights ---\n")
save_plot(
    plot(si.wt, type = "index", prop.select = 0.1),
    "si_wt_index", w = 8, h = 6)

save_plot(
    plot(si.wt, type = "heatmap", prop.select = 0.1),
    "si_wt_heatmap", w = 9, h = 8)

save_plot(
    plot(si.wt, type = "weights"),
    "si_wt_weights", w = 6, h = 5)

# ---- Smith-Hazel ----
cat("\n--- Smith-Hazel ---\n")
save_plot(
    plot(si.sh, type = "index", prop.select = 0.1),
    "si_sh_index", w = 8, h = 6)

save_plot(
    plot(si.sh, type = "heatmap", prop.select = 0.1),
    "si_sh_heatmap", w = 9, h = 8)

save_plot(
    plot(si.sh, type = "weights"),
    "si_sh_weights", w = 6, h = 5)

# ---- Desired gains ----
cat("\n--- Desired gains ---\n")
save_plot(
    plot(si.dg, type = "index", prop.select = 0.1),
    "si_dg_index", w = 8, h = 6)

save_plot(
    plot(si.dg, type = "heatmap", prop.select = 0.1),
    "si_dg_heatmap", w = 9, h = 8)

save_plot(
    plot(si.dg, type = "weights"),
    "si_dg_weights", w = 6, h = 5)

cat("\n=== demo_selIndex.R complete. Plots saved to:", outdir, "===\n")
