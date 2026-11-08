# =============================================================================
# demo_fineMap_top4.R
# Fine-map the four largest QTL from the RxK wheat qtlAim analysis.
#
# "Largest" = four QTL with the greatest absolute fixed effect estimate.
# Window = 20 cM each side, inferred grid step = 0.5 cM.
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

set.seed(42)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = NULL) {
    if (is.null(height)) height <- attr(p, "suggested_height") %||% 4
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

# %||% helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

cat("\n")
cat("================================================================\n")
cat("  RxK wheat — fineMap() for four largest QTL\n")
cat("================================================================\n\n")

# =============================================================================
# 1. Data preparation
# =============================================================================
cat("1. Loading RxK data and running primeCross()\n")
data(genoRxK,  package = "wgAim")
data(phenoRxK, package = "wgAim")

pheno          <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)

rxk <- primeCross(genoRxK, type = "interval", id = "Genotype")
cat("    Chromosomes:", length(rxk$geno),
    "  Intervals:", sum(sapply(rxk$geno, function(e) ncol(e$interval.data))), "\n\n")

# =============================================================================
# 2. qtlAim
# =============================================================================
cat("2. Fitting base ASReml model and running qtlAim()\n")
base_yld <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)

qtl_yld <- qtlAim(
    baseModel = base_yld,
    genObj    = rxk,
    merge.by  = "Genotype",
    TypeI     = 0.05,
    trace     = FALSE
)

n_qtl <- length(qtl_yld$QTL$qtl)
cat("    QTL detected:", n_qtl, "\n")

# =============================================================================
# 3. Select four largest QTL by absolute effect
# =============================================================================
cat("\n3. Selecting four largest QTL by LOD score (from summary)\n")
qtl_summ  <- summary(qtl_yld, genObj = rxk)
print(qtl_summ)

# Map summary LOD scores back to qtl_keys via the Infer. Marker column.
# e.g. "3B.midm30" -> "Chr.3B.30", "4D.midm2" -> "Chr.4D.2"
infer_col <- if ("Infer. Marker" %in% names(qtl_summ)) "Infer. Marker" else "Marker"
summ_keys <- sub("^(.+)\\.midm([0-9]+)$", "Chr.\\1.\\2",
                 qtl_summ[[infer_col]])
lod_qtl   <- setNames(qtl_summ$LOD, summ_keys)

top4_idx  <- order(lod_qtl, decreasing = TRUE)[1:4]
top4_keys <- summ_keys[top4_idx]

cat("\n    LOD scores per QTL (descending):\n")
for (i in order(lod_qtl, decreasing = TRUE)) {
    cat(sprintf("      %s  LOD = %.2f%s\n",
                summ_keys[i], lod_qtl[i],
                if (i %in% top4_idx) "  *" else ""))
}
cat("\n    Top 4 selected:", paste(top4_keys, collapse = ", "), "\n\n")

# =============================================================================
# 4. fineMap() — 20 cM window, 0.5 cM step
# =============================================================================
cat("4. fineMap() — window = 20 cM, step = 0.5 cM\n")
fm_top4 <- fineMap(
    object           = qtl_yld,
    genObj           = rxk,
    qtl              = top4_keys,
    window           = 20,
    step             = 0.5,
    exclusion.window = 20,
    trace            = FALSE
)

cat("\n    Results summary:\n")
for (nm in names(fm_top4)) {
    df    <- fm_top4[[nm]]
    peak  <- df[which.max(df$LOD), ]
    orig  <- attr(df, "qtl_pos")
    shift <- round(abs(peak$dist - orig), 2)
    cat(sprintf("    %s  %d positions  peak=%.1f cM  qtl=%.1f cM  shift=%.1f cM\n",
                nm, nrow(df), peak$dist, orig, shift))
}

# =============================================================================
# 5. print.fineMap()
# =============================================================================
cat("\n5. print.fineMap()\n")
print(fm_top4)

# =============================================================================
# 6. plot.fineMap()
# =============================================================================
cat("\n6. plot.fineMap() — LOD profiles\n")
p_top4 <- plot(fm_top4, col = "steelblue", peak.col = "firebrick",
               qtl.col = "navy")
save_plot(p_top4, "FMtop4_LOD", width = 14)

cat("\n")
cat("================================================================\n")
cat("  Done.\n")
cat("================================================================\n\n")
