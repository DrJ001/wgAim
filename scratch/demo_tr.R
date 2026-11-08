# =============================================================================
# demo_tr.R
# Demonstration of aimTrace() with diagnostic plots.
#
# Loads pre-fitted qtl_yld (11 QTL, RxK wheat) and gwas_fit (3 markers,
# simulated panel) from the saved workspace to avoid re-running models.
#
# Plots saved to scratch/plots/:
#   TR01_qtlAim_lrt.png         — LRT statistic trace
#   TR02_qtlAim_stability.png   — QTL effect stability
#   TR03_gwasAim_lrt.png        — GWAS LRT statistic trace
#   TR04_gwasAim_stability.png  — GWAS marker effect stability
# =============================================================================

devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                   quiet = TRUE)
library(asreml)
library(ggplot2)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

# Load pre-fitted objects
cat("Loading saved workspace...\n")
load("scratch/.RData")
cat(sprintf("  qtl_yld : %d QTL detected, %d iterations stored\n",
            length(qtl_yld$QTL$qtl),
            length(qtl_yld$QTL$diag$oint)))
cat(sprintf("  gwas_fit: %d markers detected, %d iterations stored\n",
            length(gwas_fit$QTL$qtl),
            length(gwas_fit$QTL$diag$oint)))

# =============================================================================
# PART 1 — aimTrace.qtlAim(): console output + plots
# =============================================================================
cat("\n--- Part 1: aimTrace.qtlAim() ---\n")

cat("\n1a. Console output (p-value matrix + LRT table)\n")
aimTrace(qtl_yld)

cat("\n1b. LRT trace plot\n")
p_lrt_qtl <- aimTrace(qtl_yld, plot = "lrt", sig.col = "firebrick")
save_plot(p_lrt_qtl, "TR01_qtlAim_lrt", width = 10, height = 5)

cat("\n1c. Effect stability plot\n")
n_qtl <- length(qtl_yld$QTL$qtl)
p_stab_qtl <- aimTrace(qtl_yld, plot = "stability", sig.col = "firebrick")
save_plot(p_stab_qtl, "TR02_qtlAim_stability",
          width = 12, height = 4 * ceiling(n_qtl / 3))

cat("\n1d. Both plots returned as a list\n")
plist <- aimTrace(qtl_yld, plot = "both", sig.col = "firebrick")
cat("  List names:", paste(names(plist), collapse = ", "), "\n")

# =============================================================================
# PART 2 — aimTrace.gwasAim(): console output + plots
# =============================================================================
cat("\n--- Part 2: aimTrace.gwasAim() ---\n")

cat("\n2a. Console output (p-value matrix + LRT table)\n")
aimTrace(gwas_fit)

cat("\n2b. LRT trace plot\n")
p_lrt_gwas <- aimTrace(gwas_fit, plot = "lrt", sig.col = "red")
save_plot(p_lrt_gwas, "TR03_gwasAim_lrt", width = 8, height = 5)

cat("\n2c. Effect stability plot\n")
n_sig <- length(gwas_fit$QTL$qtl)
p_stab_gwas <- aimTrace(gwas_fit, plot = "stability", sig.col = "red")
save_plot(p_stab_gwas, "TR04_gwasAim_stability",
          width = max(8, n_sig * 3.5), height = 5)

cat("\nDone. Plots saved to", out_dir, "\n")
