# =============================================================================
# demo_linkMap.R
# Demonstration of the new ggplot2-based linkMap.cross() and
# linkMap.QTLAim() functions.
#
# Run from the package root with:
#   Rscript scratch/demo_linkMap.R
# or source() interactively after devtools::load_all().
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

## Output directory — saves plots as PNGs into scratch/
out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("Saved: ", path)
}

# =============================================================================
# 1.  linkMap.cross() — SxT dataset (287 markers, 21 chromosomes)
#     Shows: full map, all marker names via ggrepel
# =============================================================================

data(genoSxT, package = "wgAim")

cat("\n--- Plot 1: SxT full map, all marker labels ---\n")
p1 <- linkMap(genoSxT, marker.names = "markers", m.cex = 6) +
    ggplot2::ggtitle("Sunco x Tasman — Full Genetic Map (287 markers)")
save_plot(p1, "01_SxT_full_markers", width = 27, height = 7)

# =============================================================================
# 2.  linkMap.cross() — SxT subset: chromosomes 1A, 1B, 1D, 2A, 2B
#     Shows: chr subsetting, distances on y-axis labels only
# =============================================================================

cat("--- Plot 2: SxT subset of 5 chromosomes, cM distances ---\n")
p2 <- linkMap(genoSxT,
              chr          = c("1A", "1B", "1D", "2A", "2B"),
              marker.names = "dist",
              m.cex        = 6) +
    ggplot2::ggtitle("Sunco x Tasman — Chromosomes 1A–2B (cM distances)")
save_plot(p2, "02_SxT_subset_dist", width = 10, height = 7)

# =============================================================================
# 3.  linkMap.cross() — RxK dataset (500 markers, 21 chromosomes)
#     Shows: dense map with marker.names = "none" (clean, no labels)
# =============================================================================

data(genoRxK, package = "wgAim")

cat("--- Plot 3: RxK dense map, no marker labels ---\n")
p3 <- linkMap(genoRxK, marker.names = "none") +
    ggplot2::ggtitle("RAC875 x Kukri — Full Genetic Map (500 markers, no labels)")
save_plot(p3, "03_RxK_full_nolabels", width = 18, height = 7)

# =============================================================================
# 4.  linkMap.QTLAim() — RxK QTL analysis then overlay
#     Fits a QTLAim model using the RxK phenotypic/genotypic data,
#     then plots the genetic map with QTL overlaid.
# =============================================================================

cat("--- Running QTLAim on RxK dataset ---\n")
data(phenoRxK, package = "wgAim")

## Build interval object
genoRxK_int <- cross2int(genoRxK, id = "Genotype")

## Base ASReml model (row-column field design)
pheno <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)

base_mod <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)

## QTL detection
qtl_rxk <- QTLAim(
    baseModel   = base_mod,
    intervalObj = genoRxK_int,
    merge.by    = "Genotype",
    TypeI       = 0.05,
    trace       = FALSE
)

cat("--- Plot 4: RxK QTL map — flanking markers only ---\n")
p4 <- linkMap(qtl_rxk,
              intervalObj   = genoRxK_int,
              marker.names  = "markers",
              flanking      = TRUE,
              qtl.colour    = "steelblue",
              marker.colour = "firebrick",
              m.cex         = 6) +
    ggplot2::ggtitle("RAC875 x Kukri — QTL Map (yield, TypeI = 0.05)")
save_plot(p4, "04_RxK_QTL_flanking", width = 18, height = 7)

# =============================================================================
# 5.  linkMap.QTLAim() — same result but no marker labels (clean version)
# =============================================================================

cat("--- Plot 5: RxK QTL map — no marker labels ---\n")
p5 <- linkMap(qtl_rxk,
              intervalObj  = genoRxK_int,
              marker.names = "none",
              qtl.colour   = "steelblue",
              m.cex        = 6) +
    ggplot2::ggtitle("RAC875 x Kukri — QTL Map (yield, no labels)")
save_plot(p5, "05_RxK_QTL_nolabels", width = 18, height = 7)

# =============================================================================
# 6.  linkMap.QTLAim() — chromosomes carrying QTL only
#     Shows: chr subsetting to focus on QTL-carrying chromosomes
# =============================================================================

qtl_chrs <- unique(summary(qtl_rxk, genoRxK_int)$Chromosome)

cat("--- Plot 6: RxK QTL chromosomes only ---\n")
p6 <- linkMap(qtl_rxk,
              intervalObj   = genoRxK_int,
              chr           = qtl_chrs,
              marker.names  = "markers",
              flanking      = TRUE,
              qtl.colour    = "steelblue",
              marker.colour = "firebrick",
              m.cex         = 7) +
    ggplot2::ggtitle("RAC875 x Kukri — QTL Chromosomes Only")
save_plot(p6, "06_RxK_QTL_chrs_only",
          width  = max(4, length(qtl_chrs) * 2.5),
          height = 7)

# =============================================================================
# 7.  linkMap.default() — multi-trait overlay (yld + tgw)
#     Run QTLAim on a second trait then call linkMap() with a list of models.
# =============================================================================

cat("--- Running QTLAim on tgw (thousand grain weight) ---\n")
base_mod_tgw <- asreml::asreml(
    fixed     = tgw ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)

qtl_tgw <- QTLAim(
    baseModel   = base_mod_tgw,
    intervalObj = genoRxK_int,
    merge.by    = "Genotype",
    TypeI       = 0.05,
    trace       = FALSE
)

cat("--- Plot 7: Multi-trait (yld + tgw) — all QTL chromosomes ---\n")
qtl_chrs_multi <- unique(c(
    as.character(summary(qtl_rxk, genoRxK_int)$Chromosome),
    as.character(summary(qtl_tgw, genoRxK_int)$Chromosome)
))

p7 <- linkMap(
    list(yld = qtl_rxk, tgw = qtl_tgw),
    intervalObj = genoRxK_int,
    chr         = qtl_chrs_multi,
    marker.names = "markers",
    flanking    = TRUE,
    m.cex       = 7
) + ggplot2::ggtitle("RAC875 x Kukri — Multi-trait QTL Map (yld + tgw)")

save_plot(p7, "07_RxK_multitrait_QTL_chrs",
          width  = max(4, length(qtl_chrs_multi) * 2.5),
          height = 7)

cat("--- Plot 8: Multi-trait (yld + tgw) — full map, no marker labels ---\n")
p8 <- linkMap(
    list(yld = qtl_rxk, tgw = qtl_tgw),
    intervalObj  = genoRxK_int,
    marker.names = "none",
    m.cex        = 7
) + ggplot2::ggtitle("RAC875 x Kukri — Multi-trait QTL Map (yld + tgw, all chromosomes)")

save_plot(p8, "08_RxK_multitrait_QTL_full", width = 27, height = 7)

# =============================================================================
# 9.  nrow — SxT full map split across 3 rows automatically
# =============================================================================
cat("--- Plot 9: SxT full map, nrow = 3 ---\n")
p9 <- linkMap(genoSxT, marker.names = "markers", m.cex = 6, nrow = 3) +
    ggplot2::ggtitle("Sunco x Tasman — Full Genetic Map (3 rows)")
save_plot(p9, "09_SxT_nrow3", width = 10, height = 18)

# =============================================================================
# 10.  nrow = "auto" — RxK QTL chromosomes, automatic row count
# =============================================================================
cat("--- Plot 10: RxK QTL map, nrow = 'auto' ---\n")
p10 <- linkMap(qtl_rxk,
               intervalObj  = genoRxK_int,
               marker.names = "markers",
               flanking     = TRUE,
               nrow         = "auto",
               m.cex        = 7) +
    ggplot2::ggtitle("RAC875 x Kukri — QTL Map (auto rows)")
save_plot(p10, "10_RxK_QTL_auto_nrow", width = 14, height = 14)

# =============================================================================
# 11.  row.chr — RxK full map grouped by homoeologous chromosome number
#      Builds row groups dynamically: all chrs starting with "1" → row 1, etc.
# =============================================================================
cat("--- Plot 11: RxK full map, row.chr by homoeologous group ---\n")
all_rxk_chrs <- names(qtl::pull.map(genoRxK))
chr_prefix   <- sub("[A-Za-z].*$", "", all_rxk_chrs)   # leading digits
chr_groups   <- split(all_rxk_chrs, chr_prefix)
row_chr_rxk  <- unname(chr_groups[order(as.integer(names(chr_groups)))])

p11 <- linkMap(genoRxK,
               marker.names = "none",
               row.chr      = row_chr_rxk) +
    ggplot2::ggtitle("RAC875 x Kukri — Full Map by Homoeologous Group")
save_plot(p11, "11_RxK_row_chr_groups",
          width  = max(lengths(row_chr_rxk)) * 1.5,
          height = length(row_chr_rxk) * 3)

# =============================================================================
# 12.  row.chr + multi-trait — RxK yld + tgw QTL, grouped by chromosome number
# =============================================================================
cat("--- Plot 12: RxK multi-trait QTL, row.chr by homoeologous group ---\n")

## Only include groups that contain at least one QTL chromosome
qtl_chrs_all <- unique(c(
    as.character(summary(qtl_rxk, genoRxK_int)$Chromosome),
    as.character(summary(qtl_tgw, genoRxK_int)$Chromosome)
))
row_chr_qtl <- lapply(row_chr_rxk, function(r) r[r %in% qtl_chrs_all])
row_chr_qtl <- Filter(function(r) length(r) > 0, row_chr_qtl)

p12 <- linkMap(
    list(yld = qtl_rxk, tgw = qtl_tgw),
    intervalObj  = genoRxK_int,
    marker.names = "markers",
    flanking     = TRUE,
    row.chr      = row_chr_qtl,
    m.cex        = 7
) + ggplot2::ggtitle("RAC875 x Kukri — Multi-trait QTL by Homoeologous Group")
save_plot(p12, "12_RxK_multitrait_row_chr",
          width  = max(lengths(row_chr_qtl)) * 2.5,
          height = length(row_chr_qtl) * 3.5)

cat("\nAll plots saved to", out_dir, "\n")
