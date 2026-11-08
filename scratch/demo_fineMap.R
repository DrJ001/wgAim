# =============================================================================
# demo_fineMap.R
# Demonstration of fineMap() for qtlAim and gwasAim models.
#
# Structure:
#   PART 1 — qtlAim + fineMap()
#     1a.  primeCross() with type = "interval"  (infer = "mid")
#     1b.  qtlAim() on RxK wheat data
#     1c.  fineMap() — all QTL, default window/step
#     1d.  fineMap() — single nominated QTL
#     1e.  print.fineMap() + plot.fineMap() (LOD and p-value profiles)
#     1f.  fineMap() on a marker-type wgCross — shows the grid is always built
#
#   PART 2 — gwasAim + fineMap()
#     2a.  Simulate panel data with 3 planted effects
#     2b.  primePanel() + gwasAim()
#     2c.  fineMap() — all significant markers
#     2d.  print.fineMap() + plot.fineMap()
#
# Run interactively after  devtools::load_all()  or via source().
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

set.seed(42)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 12, height = 5) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

cat("\n")
cat("================================================================\n")
cat("  wgAim fineMap() demonstration\n")
cat("================================================================\n\n")

# =============================================================================
# PART 1 — qtlAim: fine mapping of QTL from a wheat cross
# =============================================================================
cat("------ PART 1: fineMap() for qtlAim (RxK wheat) ------\n\n")

data(genoRxK,  package = "wgAim")
data(phenoRxK, package = "wgAim")

# ---- 1a. primeCross() — interval type with midpoint inference ---------------
cat("1a. primeCross() — type = 'interval', infer = 'mid'\n")
rxk_int <- primeCross(genoRxK, type = "interval", id = "Genotype")
cat("    Type:", rxk_int$type, "\n")
cat("    Chromosomes:", length(rxk_int$geno), "\n\n")

# ---- 1b. qtlAim() -----------------------------------------------------------
cat("1b. qtlAim() — forward-selection QTL detection (yld)\n")
pheno       <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)

base_yld <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)

qtl_yld <- qtlAim(
    baseModel = base_yld,
    genObj    = rxk_int,
    merge.by  = "Genotype",
    TypeI     = 0.05,
    trace     = FALSE
)

n_qtl <- length(qtl_yld$QTL$qtl)
cat("    QTL detected:", n_qtl, "\n")
cat("    QTL keys:", paste(qtl_yld$QTL$qtl, collapse = ", "), "\n\n")

# ---- 1c. fineMap() — first four QTL ----------------------------------------
cat("1c. fineMap() — fine-map first four QTL\n")
cat("    window = 50 cM, step = 2 cM, exclusion.window = 20 cM\n")
fm_all <- fineMap(
    object           = qtl_yld,
    genObj           = rxk_int,
    qtl              = qtl_yld$QTL$qtl[1:4],
    window           = 50,
    step             = 2,
    exclusion.window = 20,
    trace            = FALSE
)
cat("    Class:", class(fm_all), "\n")
cat("    QTL fine-mapped:", length(fm_all), "\n")
for (nm in names(fm_all)) {
    df   <- fm_all[[nm]]
    peak <- df[which.min(df$pvalue), ]
    cat(sprintf("    %s  — %d positions scanned, peak at %.1f cM (p = %.4f, LOD = %.4f)\n",
                nm, nrow(df), peak$dist, peak$pvalue, peak$LOD))
}

# ---- 1d. fineMap() — single nominated QTL -----------------------------------
cat("\n1d. fineMap() — fine-map a single nominated QTL\n")
target_qtl <- qtl_yld$QTL$qtl[1]
cat("    Target:", target_qtl, "\n")
fm_one <- fineMap(
    object           = qtl_yld,
    genObj           = rxk_int,
    qtl              = target_qtl,
    window           = 40,
    step             = 1,
    exclusion.window = 20,
    trace            = FALSE
)
cat("    Positions scanned:", nrow(fm_one[[1]]), "\n\n")

# ---- 1e. print and plot the fine-map results --------------------------------
cat("1e. print.fineMap()\n")
print(fm_all)

cat("\n1f. plot.fineMap() — LOD profiles (first four QTL)\n")
p_fm_lod <- plot(fm_all, col = "steelblue", peak.col = "firebrick")
save_plot(p_fm_lod, "FM01_qtlAim_finemap_LOD",
          width = 16, height = attr(p_fm_lod, "suggested_height"))

cat("\n1g. plot.fineMap() — single QTL (step = 1 cM)\n")
p_fm_one <- plot(fm_one, col = "steelblue", peak.col = "darkred")
save_plot(p_fm_one, "FM02_qtlAim_finemap_single_QTL",
          width = 7, height = attr(p_fm_one, "suggested_height"))

# ---- 1i. fineMap() on a marker-type primeCross output -----------------------
cat("\n1h. fineMap() on a marker-type wgCross\n")
cat("    The inferred grid is always built from lambdaf(), regardless of\n")
cat("    whether the original primeCross() used type = 'marker'.\n")
rxk_mkr <- primeCross(genoRxK, type = "marker", id = "Genotype")
base_mkr <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)
qtl_mkr <- qtlAim(
    baseModel = base_mkr,
    genObj    = rxk_mkr,
    merge.by  = "Genotype",
    TypeI     = 0.05,
    trace     = FALSE
)
cat("    QTL detected (marker type):", length(qtl_mkr$QTL$qtl), "\n")

fm_mkr <- fineMap(
    object           = qtl_mkr,
    genObj           = rxk_mkr,
    qtl              = qtl_mkr$QTL$qtl[1],
    window           = 40,
    step             = 2,
    exclusion.window = 20,
    trace            = FALSE
)
cat("    Positions scanned (inferred grid from imputed markers):",
    nrow(fm_mkr[[1]]), "\n")
df_mkr <- fm_mkr[[1]]
peak_m <- df_mkr[which.min(df_mkr$pvalue), ]
cat(sprintf("    Peak: %.1f cM  p = %.4f  LOD = %.4f\n\n",
            peak_m$dist, peak_m$pvalue, peak_m$LOD))

p_fm_mkr <- plot(fm_mkr, col = "purple", peak.col = "darkred")
save_plot(p_fm_mkr, "FM03_qtlAim_finemap_markercross",
          width = 7, height = attr(p_fm_mkr, "suggested_height"))

# =============================================================================
# PART 2 — gwasAim: fine mapping of significant markers
# =============================================================================
cat("------ PART 2: fineMap() for gwasAim (simulated panel) ------\n\n")

# Simulate 300 lines × 5 chromosomes × 100 markers, 2 reps
# Plant 2 effects large enough to be reliably detected
n.lines   <- 300
n.chr     <- 5
n.per.chr <- 100
n.markers <- n.chr * n.per.chr
line.ids  <- paste0("Line_", sprintf("%03d", seq_len(n.lines)))
mk.ids    <- paste0("M", sprintf("%04d", seq_len(n.markers)))

mafs     <- runif(n.markers, 0.2, 0.5)
geno.mat <- matrix(NA_integer_, nrow = n.lines, ncol = n.markers,
                   dimnames = list(line.ids, mk.ids))
for (j in seq_len(n.markers))
    geno.mat[, j] <- rbinom(n.lines, 2, mafs[j])

# Plant two strong effects
planted.marks   <- c("M0050", "M0175")
planted.effects <- c(1.2, -0.9)
geno.pm1        <- geno.mat - 1L

pheno.gwas <- do.call(rbind, lapply(1:2, function(r) {
    g <- rowSums(sweep(geno.pm1[, planted.marks], 2, planted.effects, "*"))
    data.frame(line  = factor(line.ids),
               rep   = factor(r),
               yield = 10 + g + rnorm(n.lines, 0, 0.7),
               stringsAsFactors = FALSE)
}))

map.df <- data.frame(
    marker = mk.ids,
    chr    = rep(paste0("Chr", seq_len(n.chr)), each = n.per.chr),
    pos    = rep(seq_len(n.per.chr), times = n.chr),
    stringsAsFactors = FALSE
)

# ---- 2a. primePanel() -------------------------------------------------------
cat("2a. primePanel() — preparing wgPanel object\n")
panel <- primePanel(
    geno     = geno.mat,
    map      = map.df,
    id       = "line",
    map.id   = "marker",
    map.chr  = "chr",
    map.pos  = "pos",
    encoding = "012",
    maf      = 0.05,
    impute   = "none"
)
panel$pheno$line <- factor(line.ids)
cat("    Class:", paste(class(panel), collapse = ", "), "\n")
cat("    Planted effects at:", paste(planted.marks, collapse = ", "), "\n\n")

# ---- 2b. gwasAim() ----------------------------------------------------------
cat("2b. gwasAim() — forward-selection GWAS\n")
base.gwas <- asreml::asreml(
    fixed     = yield ~ rep,
    random    = ~ line,
    data      = pheno.gwas,
    trace     = FALSE
)
gwas_fit <- gwasAim(
    baseModel        = base.gwas,
    genObj           = panel,
    merge.by         = "line",
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE,
    na.action        = asreml::na.method(x = "include")
)
n_sig <- length(gwas_fit$QTL$qtl)
cat("    Significant markers detected:", n_sig, "\n")
cat("    Keys:", paste(gwas_fit$QTL$qtl, collapse = ", "), "\n\n")

if (n_sig == 0L) {
    cat("    No significant markers — skipping fineMap demo for gwasAim.\n\n")
} else {

    # ---- 2c. fineMap() — all significant markers ----------------------------
    cat("2c. fineMap() — fine-map all significant markers\n")
    cat("    window = 30 cM (= ±30 actual markers here), exclusion.window = 15\n")
    fm_gwas <- fineMap(
        object           = gwas_fit,
        genObj           = panel,
        window           = 30,
        exclusion.window = 15,
        trace            = FALSE,
        na.action        = asreml::na.method(x = "include")
    )
    cat("    Class:", class(fm_gwas), "\n")
    cat("    Type attribute:", attr(fm_gwas, "type"), "\n")
    for (nm in names(fm_gwas)) {
        df   <- fm_gwas[[nm]]
        peak <- df[which.min(df$pvalue), ]
        cat(sprintf("    %s  — %d markers scanned, peak: %s at %.1f cM (p = %.4f, LOD = %.4f)\n",
                    nm, nrow(df), peak$mark, peak$dist, peak$pvalue, peak$LOD))
    }

    # ---- 2d. print and plot --------------------------------------------------
    cat("\n2d. print.fineMap()\n")
    print(fm_gwas)

    cat("\n2e. plot.fineMap() — LOD profiles (gwasAim)\n")
    p_gfm_lod <- plot(fm_gwas, col = "darkgreen", peak.col = "firebrick")
    save_plot(p_gfm_lod, "FM04_gwasAim_finemap_LOD",
              width = 12, height = attr(p_gfm_lod, "suggested_height"))

    # ---- 2g. Single marker fine-map -----------------------------------------
    if (n_sig >= 1L) {
        cat("\n2g. fineMap() — single nominated marker\n")
        tgt_mk <- gwas_fit$QTL$qtl[1]
        cat("    Target:", tgt_mk, "\n")
        fm_gwas_one <- fineMap(
            object           = gwas_fit,
            genObj           = panel,
            qtl              = tgt_mk,
            window           = 20,
            exclusion.window = 15,
            trace            = FALSE,
            na.action        = asreml::na.method(x = "include")
        )
        cat("    Markers scanned:", nrow(fm_gwas_one[[1]]), "\n")
        p_gfm_one <- plot(fm_gwas_one, col = "darkgreen", peak.col = "darkred")
        save_plot(p_gfm_one, "FM05_gwasAim_finemap_single_marker",
                  width = 7, height = attr(p_gfm_one, "suggested_height"))
    }
}

# =============================================================================
# Summary
# =============================================================================
fm_plots <- list.files(out_dir, pattern = "^FM.*\\.png$")
cat("\n")
cat("================================================================\n")
cat("  fineMap demo complete —", length(fm_plots), "plots saved to", out_dir, "\n")
cat("================================================================\n\n")
for (f in sort(fm_plots))
    cat("  ", f, "\n")
