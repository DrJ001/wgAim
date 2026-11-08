# Minimal script to generate heatmap plots only.
# Run from the wgaim project root.

devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim",
                   quiet = TRUE)
library(asreml)
library(ggplot2)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

# =============================================================================
# PART 1 — qtlAim heatmap: RxK wheat (many QTL, real data)
# =============================================================================
cat("\n--- Part 1: qtlAim heatmap (RxK wheat) ---\n")

data(genoRxK,  package = "wgAim")
data(phenoRxK, package = "wgAim")

rxk_cross <- primeCross(genoRxK, type = "interval", id = "Genotype")

pheno        <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)

base_yld <- asreml::asreml(
    fixed     = yld ~ Type + lrow,
    random    = ~ Genotype + Range,
    residual  = ~ ar1(Range):ar1(Row),
    data      = pheno,
    trace     = FALSE
)

cat("  Running qtlAim()...\n")
qtl_yld <- qtlAim(
    baseModel        = base_yld,
    genObj           = rxk_cross,
    merge.by         = "Genotype",
    TypeI            = 0.05,
    trace            = FALSE,
    na.action        = asreml::na.method(x = "include")
)

n_qtl  <- length(qtl_yld$QTL$qtl)
n_iter <- length(qtl_yld$QTL$diag$oint)
cat(sprintf("  QTL detected: %d  |  iterations stored: %d\n", n_qtl, n_iter))
print(summary(qtl_yld, genObj = rxk_cross))

# Full genome heatmap — cap at 3 so gradient spreads across the weaker signals
cat("  Generating full-genome heatmap...\n")
p_heat_full <- plot(qtl_yld, genObj = rxk_cross, type = "heatmap",
                    sig.col = "firebrick", cap = 10)
save_plot(p_heat_full, "HM01_qtlAim_heatmap_full",
          width = 18, height = 3 + n_iter * 0.6)

# Significant chromosomes only
sig_chrs <- unique(as.character(summary(qtl_yld, genObj = rxk_cross)$Chromosome))
cat(sprintf("  Significant chromosomes: %s\n", paste(sig_chrs, collapse = ", ")))
p_heat_sig <- plot(qtl_yld, genObj = rxk_cross, type = "heatmap",
                   chr     = sig_chrs,
                   sig.col = "firebrick", cap = 3)
save_plot(p_heat_sig, "HM02_qtlAim_heatmap_sigchrs",
          width = max(8, length(sig_chrs) * 3),
          height = 3 + n_iter * 0.6)

# =============================================================================
# PART 2 — gwasAim heatmap: simulated panel
# =============================================================================
cat("\n--- Part 2: gwasAim heatmap (simulated panel) ---\n")

set.seed(42)
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

base.gwas <- asreml::asreml(
    fixed  = yield ~ rep,
    random = ~ line,
    data   = pheno.gwas,
    trace  = FALSE
)

cat("  Running gwasAim()...\n")
gwas_fit <- gwasAim(
    baseModel        = base.gwas,
    genObj           = panel,
    merge.by         = "line",
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE,
    na.action        = asreml::na.method(x = "include")
)

n_sig  <- length(gwas_fit$QTL$qtl)
n_giter <- length(gwas_fit$QTL$diag$oint)
cat(sprintf("  Markers detected: %d  |  iterations stored: %d\n", n_sig, n_giter))
print(summary(gwas_fit, genObj = panel))

cat("  Generating GWAS heatmap...\n")
p_gheat <- plot(gwas_fit, genObj = panel, type = "heatmap",
                sig.col = "red")
save_plot(p_gheat, "HM03_gwasAim_heatmap",
          width = 14, height = 3 + n_giter * 0.6)

cat("\nDone. Plots saved to", out_dir, "\n")
