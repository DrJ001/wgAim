# =============================================================================
# demo_linkMap_GWAS.R
# Demonstration of linkMap.GWASAim().
#
# Uses the same simulation setup as benchmark_GWASAim.R:
#   300 lines × 500 markers across 5 chromosomes, 2 crossed reps,
#   3 planted SNP effects.
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

set.seed(2025)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("Saved: ", path)
}

## ---------------------------------------------------------------------------
## 1.  Simulate marker genotypes — 5 chromosomes × 100 markers (1 cM spacing)
## ---------------------------------------------------------------------------
n.lines   <- 300
n.chr     <- 5
n.per.chr <- 100
n.markers <- n.chr * n.per.chr

line.ids   <- paste0("Line_", sprintf("%03d", seq_len(n.lines)))
marker.ids <- paste0("M", sprintf("%04d", seq_len(n.markers)))

mafs     <- runif(n.markers, 0.2, 0.5)
geno.mat <- matrix(NA_integer_, nrow = n.lines, ncol = n.markers,
                   dimnames = list(line.ids, marker.ids))
for (j in seq_len(n.markers))
    geno.mat[, j] <- rbinom(n.lines, 2, mafs[j])

map.df <- data.frame(
    marker = marker.ids,
    chr    = rep(paste0("Chr", seq_len(n.chr)), each = n.per.chr),
    pos    = rep(seq_len(n.per.chr), times = n.chr),
    stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------------
## 2.  Plant 3 marker effects (different chromosomes)
## ---------------------------------------------------------------------------
planted.markers <- c("M0050", "M0175", "M0320")   # Chr1-50, Chr2-75, Chr4-20
planted.effects <- c(0.8, -0.6, 0.5)
geno.pm1        <- geno.mat - 1L

## ---------------------------------------------------------------------------
## 3.  Simulate 2 crossed reps — each line appears once in each rep
## ---------------------------------------------------------------------------
pheno.df <- do.call(rbind, lapply(1:2, function(r) {
    genetic  <- rowSums(sweep(geno.pm1[, planted.markers], 2,
                              planted.effects, "*"))
    data.frame(
        line  = factor(line.ids),
        rep   = factor(r),
        yield = 10 + genetic + rnorm(n.lines, 0, 0.8),
        stringsAsFactors = FALSE
    )
}))

## ---------------------------------------------------------------------------
## 4.  Build panel object
## ---------------------------------------------------------------------------
cat("--- Building panel object ---\n")
panel <- makePanel(
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

## ---------------------------------------------------------------------------
## 5.  Fit GWASAim
## ---------------------------------------------------------------------------
cat("--- Fitting base model ---\n")
base.mod <- asreml::asreml(
    fixed     = yield ~ rep,
    random    = ~ line,
    data      = pheno.df,
    trace     = FALSE
)

cat("--- Running GWASAim ---\n")
gwas.fit <- GWASAim(
    baseModel        = base.mod,
    panelObj         = panel,
    merge.by         = "line",
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE,
    na.action        = asreml::na.method(x = "include")
)

cat("\nSignificant markers:\n")
print(gwas.fit, panelObj = panel)

## ---------------------------------------------------------------------------
## 6.  Plot G01 — all chromosomes, marker names labelled
## ---------------------------------------------------------------------------
cat("\n--- Plot G01: all chromosomes ---\n")
pg1 <- linkMap(gwas.fit,
               panelObj     = panel,
               marker.names = "markers",
               m.cex        = 7) +
    ggplot2::ggtitle("GWAS Linkage Map — Significant Markers")
save_plot(pg1, "G01_GWAS_all_chrs", width = 14, height = 6)

## ---------------------------------------------------------------------------
## 7.  Plot G02 — significant-marker chromosomes only, cM distances
## ---------------------------------------------------------------------------
cat("--- Plot G02: significant-marker chromosomes only ---\n")
sig_chrs <- unique(summary(gwas.fit, panelObj = panel)$Chromosome)

pg2 <- linkMap(gwas.fit,
               panelObj     = panel,
               chr          = sig_chrs,
               marker.names = "dist",
               m.cex        = 7) +
    ggplot2::ggtitle("GWAS Linkage Map — Significant Marker Chromosomes (cM)")
save_plot(pg2, "G02_GWAS_sig_chrs",
          width = max(3, length(sig_chrs) * 2.5), height = 6)

## ---------------------------------------------------------------------------
## 8.  Plot G03 — nrow = "auto" (all chromosomes in multiple rows)
## ---------------------------------------------------------------------------
cat("--- Plot G03: nrow = 'auto' ---\n")
pg3 <- linkMap(gwas.fit,
               panelObj     = panel,
               marker.names = "markers",
               nrow         = "auto",
               m.cex        = 7) +
    ggplot2::ggtitle("GWAS Linkage Map — Auto Row Layout")
save_plot(pg3, "G03_GWAS_auto_nrow", width = 10, height = 10)

cat("\nAll GWAS plots saved to", out_dir, "\n")
