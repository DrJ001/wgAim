# =============================================================================
# demo_wgAim.R
# Comprehensive demonstration of the wgAim package.
#
# Covers (in order):
#   1.  primeCross()          — prepare a cross object
#   2.  qtlAim()              — QTL detection
#   3.  print / summary       — QTL results
#   4.  plot.qtlAim()         — outlier, blups, chr, heatmap plots
#   5.  linkMap (cross/qtl)   — base map, single-trait, multi-trait overlays
#   6.  primePanel()          — prepare a marker panel
#   7.  gwasAim()             — forward-selection GWAS
#   8.  print / summary       — GWAS results
#   9.  plot.gwasAim()        — manhattan, outlier, blups, heatmap
#  10.  manhattan()           — classic Manhattan plot
#  11.  linkMap.gwasAim()     — significant marker overlay
#  12.  gpAim()               — genomic prediction (vm + mbf paths)
#  13.  print / summary       — GP results
#  14.  plot.gpAim()          — caterpillar, heatmap, density, blups
#
# Run interactively after  devtools::load_all()  or via:
#   Rscript scratch/demo_wgAim.R
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

set.seed(42)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

cat("\n")
cat("================================================================\n")
cat("  wgAim comprehensive demonstration\n")
cat("================================================================\n\n")

# =============================================================================
# PART 1 — qtlAim: QTL detection on RxK wheat dataset
# =============================================================================
cat("------ PART 1: primeCross + qtlAim (RxK wheat) ------\n\n")

data(genoRxK,  package = "wgAim")
data(phenoRxK, package = "wgAim")

# ---- 1a. primeCross() -------------------------------------------------------
cat("1a. primeCross() — preparing wgCross object (interval type)\n")
rxk_cross <- primeCross(genoRxK, type = "interval", id = "Genotype")
cat("    Class:", paste(class(rxk_cross), collapse = ", "), "\n")
cat("    Stored type:", rxk_cross$type, "\n")
cat("    Chromosomes:", length(rxk_cross$geno), "\n")
cat("    Interval data present on chr 1A:",
    !is.null(rxk_cross$geno[["1A"]]$interval.data), "\n\n")

# Also demo marker-type for contrast
cat("1b. primeCross() — marker type (no interval data generated)\n")
rxk_marker <- primeCross(genoRxK, type = "marker", id = "Genotype")
cat("    Stored type:", rxk_marker$type, "\n")
cat("    Interval data present:", !is.null(rxk_marker$geno[["1A"]]$interval.data), "\n\n")

# ---- 1c. Base ASReml model --------------------------------------------------
cat("1c. Fitting base ASReml model (yld ~ Range + Genotype + Row)\n")
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
cat("    Converged:", base_yld$converge, "\n\n")

# ---- 1d. qtlAim() -----------------------------------------------------------
cat("1d. qtlAim() — forward-selection QTL detection\n")
qtl_yld <- qtlAim(
    baseModel = base_yld,
    genObj    = rxk_cross,
    exclusion.window = 50,
    merge.by  = "Genotype",
    TypeI     = 0.05,
    trace     = FALSE
)

qtl_yld$QTL$diag$oint[[2]]
rxk_cross$geno[["3B"]]$map
cat("    Class:", paste(class(qtl_yld), collapse = ", "), "\n")
cat("    gen.type auto-resolved from genObj:", qtl_yld$QTL$type, "\n\n")

# ---- 1e. print / summary ----------------------------------------------------
cat("1e. print.qtlAim()\n")
print(qtl_yld, genObj = rxk_cross)

cat("\n1f. summary.qtlAim()\n")
qtl_summ <- summary(qtl_yld, genObj = rxk_cross)
print(qtl_summ)

# ---- 1g. plot.qtlAim() ------------------------------------------------------
cat("\n1g. plot.qtlAim() — three plot types\n")

cat("    type = 'outlier'\n")
p_oint <- plot(qtl_yld, genObj = rxk_cross, type = "outlier",
               sig.col = "firebrick")
save_plot(p_oint, "Q01_qtlAim_outlier", width = 14,
          height = 3 * length(qtl_yld$QTL$diag$oint))

cat("    type = 'blups'\n")
p_blup <- plot(qtl_yld, genObj = rxk_cross, type = "blups",
               sig.col = "firebrick")
save_plot(p_blup, "Q02_qtlAim_blups", width = 14,
          height = 3 * length(qtl_yld$QTL$diag$blups))

# chr type only works with selection = "chromosome" — skip if not applicable
if (!is.null(qtl_yld$QTL$diag$ochr[[1]])) {
    cat("    type = 'chr'\n")
    p_chr <- plot(qtl_yld, genObj = rxk_cross, type = "chr")
    save_plot(p_chr, "Q03_qtlAim_chr", width = 10,
              height = 3 * length(qtl_yld$QTL$diag$oint))
} else {
    cat("    type = 'chr' skipped (selection = 'interval' used)\n")
}

# ---- 1h. plot.qtlAim() — heatmap --------------------------------------------
cat("\n1h. plot.qtlAim() — type = 'heatmap'\n")

# Full genome: every iteration, all chromosomes
p_heat <- plot(qtl_yld, genObj = rxk_cross, type = "heatmap",
               sig.col = "firebrick")
save_plot(p_heat, "Q03_qtlAim_heatmap_full",
          width = 18, height = 3 + length(qtl_yld$QTL$diag$oint) * 0.6)

# Subset to the two chromosomes where most signal lives — handy for focusing
summ_yld <- summary(qtl_yld, genObj = rxk_cross)
sig_chrs <- unique(as.character(summ_yld$Chromosome))
cat("    Significant chromosomes:", paste(sig_chrs, collapse = ", "), "\n")
p_heat_sig <- plot(qtl_yld, genObj = rxk_cross, type = "heatmap",
                   chr     = sig_chrs,
                   sig.col = "firebrick")
save_plot(p_heat_sig, "Q03b_qtlAim_heatmap_sigchrs",
          width = max(8, length(sig_chrs) * 3),
          height = 3 + length(qtl_yld$QTL$diag$oint) * 0.6)

# ---- 1j. linkMap — base, single-trait, multi-trait -------------------------
cat("\n1h. linkMap.cross() — base genetic map\n")
p_map <- linkMap(rxk_cross, marker.names = "none") +
    ggtitle("RAC875 x Kukri — Genetic Map")
save_plot(p_map, "Q04_linkMap_base", width = 20, height = 6)

cat("    linkMap.qtlAim() — single-trait QTL overlay\n")
p_qtl <- linkMap(qtl_yld, genObj = rxk_cross,
                 marker.names = "markers", flanking = TRUE,
                 qtl.colour = "steelblue", m.cex = 6) +
    ggtitle("RAC875 x Kukri — QTL Map (yld)")
save_plot(p_qtl, "Q05_linkMap_qtl_single", width = 20, height = 6)

# Fit tgw model for multi-trait demo
cat("    Fitting tgw base model for multi-trait demo ...\n")
base_tgw <- asreml::asreml(
    fixed     = tgw ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)
qtl_tgw <- qtlAim(
    baseModel = base_tgw,
    genObj    = rxk_cross,
    merge.by  = "Genotype",
    TypeI     = 0.05,
    trace     = FALSE
)
cat("    linkMap.qtlAim() — multi-trait QTL overlay\n")
p_multi <- linkMap(qtl_yld, qtl_tgw,
                   genObj       = rxk_cross,
                   marker.names = "none") +
    ggtitle("RAC875 x Kukri — Multi-trait QTL (yld + tgw)")
save_plot(p_multi, "Q06_linkMap_qtl_multi", width = 20, height = 6)

# ---- 1k. class guard demo ---------------------------------------------------
cat("\n1i. Class guard — passing marker-type wgCross with gen.type='interval'\n")
tryCatch(
    qtlAim(base_yld, genObj = rxk_marker, merge.by = "Genotype",
           gen.type = "interval", TypeI = 0.05, trace = FALSE),
    error = function(e) cat("    Caught expected error:", conditionMessage(e), "\n")
)

# =============================================================================
# PART 2 — gwasAim: GWAS on simulated panel data
# =============================================================================
cat("\n")
cat("------ PART 2: primePanel + gwasAim (simulated panel) ------\n\n")

# Simulate 300 lines × 5 chromosomes × 100 markers, 2 reps
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

map.df <- data.frame(
    marker = mk.ids,
    chr    = rep(paste0("Chr", seq_len(n.chr)), each = n.per.chr),
    pos    = rep(seq_len(n.per.chr), times = n.chr),
    stringsAsFactors = FALSE
)

planted.marks   <- c("M0050", "M0175", "M0320")
planted.effects <- c(0.8, -0.6, 0.5)
geno.pm1        <- geno.mat - 1L

pheno.gwas <- do.call(rbind, lapply(1:2, function(r) {
    g <- rowSums(sweep(geno.pm1[, planted.marks], 2, planted.effects, "*"))
    data.frame(line  = factor(line.ids),
               rep   = factor(r),
               yield = 10 + g + rnorm(n.lines, 0, 0.9),
               stringsAsFactors = FALSE)
}))

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
cat("    Inherits 'cross':", inherits(panel, "cross"), "\n\n")

# ---- 2b. gwasAim() ----------------------------------------------------------
cat("2b. gwasAim() — forward-selection GWAS\n")
base.gwas <- asreml::asreml(
    fixed     = yield ~ rep,
    random    = ~ line,
    data      = pheno.gwas,
    trace     = TRUE
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
cat("    Class:", paste(class(gwas_fit), collapse = ", "), "\n\n")

# ---- 2c. print / summary ----------------------------------------------------
cat("2c. print.gwasAim()\n")
print(gwas_fit, genObj = panel)

cat("\n2d. summary.gwasAim()\n")
print(summary(gwas_fit, genObj = panel))

# ---- 2e. plot.gwasAim() -----------------------------------------------------
cat("\n2e. plot.gwasAim() — three plot types\n")

cat("    type = 'manhattan'\n")
p_mh <- plot(gwas_fit, genObj = panel, type = "manhattan",
             sig.col = "red", pt.col = c("steelblue", "grey60"))
save_plot(p_mh, "G01_gwasAim_manhattan", width = 14,
          height = 3 * length(gwas_fit$QTL$diag$oint))

cat("    type = 'outlier'\n")
p_go <- plot(gwas_fit, genObj = panel, type = "outlier",
             sig.col = "firebrick")
save_plot(p_go, "G02_gwasAim_outlier", width = 14,
          height = 3 * length(gwas_fit$QTL$diag$oint))

cat("    type = 'blups'\n")
p_gb <- plot(gwas_fit, genObj = panel, type = "blups",
             sig.col = "firebrick")
save_plot(p_gb, "G03_gwasAim_blups", width = 14,
          height = 3 * length(gwas_fit$QTL$diag$blups))

# ---- 2f. plot.gwasAim() — heatmap -------------------------------------------
cat("\n2f. plot.gwasAim() — type = 'heatmap'\n")

p_gheat <- plot(gwas_fit, genObj = panel, type = "heatmap",
                sig.col = "red")
save_plot(p_gheat, "G04_gwasAim_heatmap",
          width = 14, height = 3 + length(gwas_fit$QTL$diag$oint) * 0.6)

# ---- 2g. linkMap.gwasAim() --------------------------------------------------
cat("\n2g. linkMap.gwasAim() — significant marker map overlay\n")
p_gmap <- linkMap(gwas_fit, genObj = panel, marker.names = "markers", m.cex = 7) +
    ggtitle("GWAS — Significant Marker Map")
save_plot(p_gmap, "G05_linkMap_gwas", width = 14, height = 6)

sig_chrs <- unique(summary(gwas_fit, genObj = panel)$Chromosome)
p_gmap2  <- linkMap(gwas_fit, genObj = panel,
                    chr = sig_chrs, marker.names = "dist", m.cex = 7) +
    ggtitle("GWAS — Significant Chromosomes Only")
save_plot(p_gmap2, "G06_linkMap_gwas_sig_chrs",
          width = max(4, length(sig_chrs) * 3), height = 6)

# ---- 2h. class guard --------------------------------------------------------
cat("\n2h. Class guard — passing wgCross to gwasAim\n")
tryCatch(
    gwasAim(base.gwas, genObj = rxk_cross, merge.by = "line",
            TypeI = 0.05, trace = FALSE),
    error = function(e) cat("    Caught expected error:", conditionMessage(e), "\n")
)

# =============================================================================
# PART 3 — gpAim: Genomic Prediction
# =============================================================================
cat("\n")
cat("------ PART 3: gpAim (simulated data, vm + mbf paths) ------\n\n")

# Simulate 200 lines × 300 markers, h²=0.5, 30 causal markers, 2 reps
n.gp   <- 200
n.mk   <- 300
n.caus <- 30
gp.ids <- paste0("G_", sprintf("%03d", seq_len(n.gp)))
mk.gp  <- paste0("MG", sprintf("%04d", seq_len(n.mk)))

mafs.gp <- runif(n.mk, 0.2, 0.5)
M <- matrix(NA_real_, nrow = n.gp, ncol = n.mk,
            dimnames = list(gp.ids, mk.gp))
for (j in seq_len(n.mk))
    M[, j] <- rbinom(n.gp, 2, mafs.gp[j]) - 1L

causal.idx     <- sample(n.mk, n.caus)
marker.effects <- rnorm(n.caus, 0, sqrt(0.5 / n.caus))
g.true         <- as.numeric(M[, causal.idx] %*% marker.effects)
e.var          <- var(g.true) * (1 - 0.5) / 0.5

map.gp <- data.frame(
    marker = mk.gp,
    chr    = rep(paste0("C", 1:5), each = n.mk / 5),
    pos    = rep(seq(1, n.mk / 5), times = 5),
    stringsAsFactors = FALSE
)

pheno.gp <- do.call(rbind, lapply(1:2, function(r) {
    data.frame(
        line  = factor(gp.ids),
        rep   = factor(r),
        trait = g.true + rnorm(n.gp, 0, sqrt(e.var)),
        stringsAsFactors = FALSE
    )
}))

# ---- 3a. primePanel() for GP ------------------------------------------------
cat("3a. primePanel() for GP panel (200 lines, 300 markers)\n")
panel.gp <- primePanel(
    geno     = M + 1L,         # back to 0/1/2 for encoding="012"
    map      = map.gp,
    id       = "line",
    map.id   = "marker",
    map.chr  = "chr",
    map.pos  = "pos",
    encoding = "012",
    maf      = 0.02,
    impute   = "none"
)
panel.gp$pheno$line <- factor(gp.ids)

# ---- 3b. vm path (markers > lines) ------------------------------------------
cat("\n3b. gpAim() — vm path (300 markers > 200 lines)\n")
base.gp <- asreml::asreml(
    fixed     = trait ~ rep,
    random    = ~ line,
    data      = pheno.gp,
    trace     = FALSE
)
gp_vm <- gpAim(
    baseModel = base.gp,
    genObj    = panel.gp,
    merge.by  = "line",
    gen.type  = "marker",
    trace     = FALSE,
    na.action = asreml::na.method(x = "include")
)
cat("    Path used:", gp_vm$GP$path, "\n")
cat("    Estimated h:", round(gp_vm$GP$heritability, 3), "\n")
cat("    marker.effects stored:",
    !is.null(gp_vm$GP$marker.effects), "\n")

cat("\n3c. print.gpAim()\n")
print(gp_vm)

cat("\n3d. summary.gpAim()\n")
print(summary(gp_vm))

cat("\n3e. plot.gpAim() — four plot types (vm path)\n")

cat("    type = 'caterpillar'\n")
p_cat <- plot(gp_vm, type = "caterpillar", top.n = 15)
save_plot(p_cat, "P01_gpAim_caterpillar_vm", width = 10, height = 6)

cat("    type = 'density'\n")
p_den <- plot(gp_vm, type = "density", prop.select = 0.10, max.labels = 5)
save_plot(p_den, "P02_gpAim_density_vm", width = 8, height = 5)

cat("    type = 'heatmap'\n")
p_hm <- plot(gp_vm, type = "heatmap")
save_plot(p_hm, "P03_gpAim_heatmap_vm", width = 8, height = 7)

cat("    type = 'blups' (marker effects across genome)\n")
p_bp <- plot(gp_vm, genObj = panel.gp, type = "blups",
             pt.col = c("steelblue", "grey60"), pt.cex = 0.8)
save_plot(p_bp, "P04_gpAim_blups_vm", width = 12, height = 5)

# ---- 3f. mbf path (force lines >= markers scenario) ------------------------
cat("\n3f. gpAim() — mbf path (forced, to demonstrate both paths)\n")
gp_mbf <- gpAim(
    baseModel = base.gp,
    genObj    = panel.gp,
    merge.by  = "line",
    gen.type  = "marker",
    force     = TRUE,           # force mbf path
    trace     = FALSE,
    na.action = asreml::na.method(x = "include")
)
cat("    Path used:", gp_mbf$GP$path, "\n")
cat("    Estimated h²:", round(gp_mbf$GP$heritability, 3), "\n")

cat("\n    plot.gpAim() type = 'caterpillar' (mbf path)\n")
p_cat2 <- plot(gp_mbf, type = "caterpillar", top.n = 15)
save_plot(p_cat2, "P05_gpAim_caterpillar_mbf", width = 10, height = 6)

cat("    plot.gpAim() type = 'blups' (mbf path, direct marker BLUPs)\n")
p_bp2 <- plot(gp_mbf, genObj = panel.gp, type = "blups",
              pt.col = c("steelblue", "grey60"), pt.cex = 0.8)
save_plot(p_bp2, "P06_gpAim_blups_mbf", width = 12, height = 5)

# ---- 3g. class guard --------------------------------------------------------
cat("\n3g. Class guard — genObj = NULL for blups type\n")
tryCatch(
    plot(gp_vm, genObj = NULL, type = "blups"),
    error = function(e) cat("    Caught expected error:", conditionMessage(e), "\n")
)

cat("\n3h. Class guard — interval request on marker-type wgCross\n")
tryCatch(
    qtlAim(base_yld, genObj = rxk_marker, merge.by = "Genotype",
           gen.type = "interval", TypeI = 0.05, trace = FALSE),
    error = function(e) cat("    Caught expected error:", conditionMessage(e), "\n")
)

# =============================================================================
# Summary
# =============================================================================
plots_made <- list.files(out_dir, pattern = "\\.png$")
cat("\n")
cat("================================================================\n")
cat("  Demo complete —", length(plots_made), "plots saved to", out_dir, "\n")
cat("================================================================\n\n")
for (f in sort(plots_made))
    cat("  ", f, "\n")
