# =============================================================================
# demo_effects_contrast.R
# Demonstration of the two new result-level plot types:
#
#   plot(object, genObj, type = "effects")   — effect direction lollipop
#   plot(object, genObj, type = "contrast",
#        data = <pheno.data>)                — allele contrast plot
#
# Both types work for qtlAim and gwasAim objects.
#
# Structure:
#   PART 1 — qtlAim on RxK wheat (real data)
#     1a.  primeCross + qtlAim
#     1b.  type = "effects"   (interval-type wgCross, method = "fixed")
#     1c.  type = "contrast"  (interval-type wgCross)
#     1d.  type = "effects"   (marker-type wgCross — same helpers, different input)
#     1e.  type = "contrast"  (marker-type wgCross)
#
#   PART 2 — gwasAim on simulated diversity panel
#     2a.  primePanel + gwasAim
#     2b.  type = "effects"   (GWAS panel, 0/1/2 encoding)
#     2c.  type = "contrast"  (GWAS panel)
#
#   PART 3 — edge-case guards
#     3a.  Missing data argument for contrast
#     3b.  No QTL found
#
# Run interactively after  devtools::load_all()  or via source().
# =============================================================================

devtools::load_all()
library(asreml)
library(ggplot2)

set.seed(42)

out_dir <- file.path("scratch", "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p, width = width, height = height,
                    dpi = 150, bg = "white")
    message("  Saved: ", path)
}

cat("\n")
cat("================================================================\n")
cat("  wgAim — effects & contrast plot demonstration\n")
cat("================================================================\n\n")


# =============================================================================
# PART 1 — qtlAim: RxK wheat cross
# =============================================================================
cat("------ PART 1: qtlAim effects & contrast (RxK wheat) ------\n\n")

data(genoRxK,  package = "wgAim")
data(phenoRxK, package = "wgAim")

pheno          <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)

# ---- 1a. primeCross (interval type) + qtlAim --------------------------------
cat("1a. primeCross() — interval type\n")
rxk_int <- primeCross(genoRxK, type = "interval", id = "Genotype")
cat("    Type:", rxk_int$type,
    "| Chromosomes:", length(rxk_int$geno), "\n\n")

cat("    Fitting base ASReml model (yld ~ Range + Genotype + Row)\n")
base_yld <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)

cat("    qtlAim() — forward selection (method = 'fixed')\n")
qtl_yld <- qtlAim(
    baseModel        = base_yld,
    genObj           = rxk_int,
    merge.by         = "Genotype",
    exclusion.window = 50,
    TypeI            = 0.05,
    trace            = FALSE
)
n_qtl <- length(qtl_yld$QTL$qtl)
cat("    QTL detected:", n_qtl, "\n")
cat("    QTL keys:", paste(qtl_yld$QTL$qtl, collapse = ", "), "\n\n")

# phenoData written to caller env as yld.data
cat("    phenoData variable in environment: yld.data\n\n")

# ---- 1b. type = "effects" (interval qtlAim) ---------------------------------
cat("1b. plot(qtl_yld, genObj = rxk_int, type = 'effects')\n")
cat("    Lollipop of additive QTL effects — firebrick = A favoured,\n")
cat("    steelblue = B favoured.  ±1 SE bars.  % var on heads.\n\n")

p_eff_int <- plot(qtl_yld, genObj = rxk_int, type = "effects")
print(p_eff_int)
save_plot(p_eff_int, "EC01_qtlAim_interval_effects",
          width = max(7, n_qtl * 1.4), height = 5)

# ---- 1c. type = "contrast" (interval qtlAim) --------------------------------
cat("1c. plot(..., type = 'contrast') — ncol = 1 (default, vertical)\n")
cat("    Allele class from interval.data (what the model used): B(−1) / A(+1)\n\n")

p_con_int <- plot(qtl_yld, genObj = rxk_int, type = "contrast",
                  data = yld.data, ncol = 1L)
print(p_con_int)
save_plot(p_con_int, "EC02_qtlAim_interval_contrast_vertical",
          width = 8, height = n_qtl * 2.8)

cat("    ncol = 3 — grid layout\n\n")
p_con_int3 <- plot(qtl_yld, genObj = rxk_int, type = "contrast",
                   data = yld.data, ncol = 3L)
print(p_con_int3)
save_plot(p_con_int3, "EC02b_qtlAim_interval_contrast_ncol3",
          width = 14, height = ceiling(n_qtl / 3) * 3.2)

# ---- 1d. primeCross (marker type) + qtlAim -----------------------------------
cat("1d. primeCross() — marker type (no interval midpoints)\n")
rxk_mkr <- primeCross(genoRxK, type = "marker", id = "Genotype")
cat("    Type:", rxk_mkr$type, "\n\n")

cat("    qtlAim() — marker-type analysis\n")
base_yld2 <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)
qtl_mkr <- qtlAim(
    baseModel        = base_yld2,
    genObj           = rxk_mkr,
    merge.by         = "Genotype",
    exclusion.window = 50,
    TypeI            = 0.05,
    trace            = FALSE
)
n_mkr <- length(qtl_mkr$QTL$qtl)
cat("    QTL detected:", n_mkr, "\n\n")

cat("1e. type = 'effects' — marker-type qtlAim\n")
p_eff_mkr <- plot(qtl_mkr, genObj = rxk_mkr, type = "effects")
print(p_eff_mkr)
save_plot(p_eff_mkr, "EC03_qtlAim_marker_effects",
          width = max(7, n_mkr * 1.4), height = 5)

cat("\n1f. type = 'contrast' — marker-type qtlAim (ncol = 1)\n")
p_con_mkr <- plot(qtl_mkr, genObj = rxk_mkr, type = "contrast",
                  data = yld.data, ncol = 1L)
print(p_con_mkr)
save_plot(p_con_mkr, "EC04_qtlAim_marker_contrast",
          width = 8, height = n_mkr * 2.8)


# =============================================================================
# PART 2 — gwasAim: simulated diversity panel
# =============================================================================
cat("\n")
cat("------ PART 2: gwasAim effects & contrast (simulated panel) ------\n\n")

# Simulate 300 lines × 5 chromosomes × 100 markers, 2 reps.
# Three planted additive effects of varying magnitude and direction.
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

# Three planted effects: two large (reliably detected), one moderate
planted.marks   <- c("M0050", "M0175", "M0320")
planted.effects <- c(1.2, -0.9, 0.6)     # +, -, +  → tests both directions
geno.pm1        <- geno.mat - 1L          # centre to ±1 for effect computation

pheno.gwas <- do.call(rbind, lapply(seq_len(2), function(r) {
    g <- rowSums(sweep(geno.pm1[, planted.marks], 2, planted.effects, "*"))
    data.frame(
        line  = factor(line.ids),
        rep   = factor(r),
        yield = 10 + g + rnorm(n.lines, 0, 0.7),
        stringsAsFactors = FALSE
    )
}))

map.df <- data.frame(
    marker = mk.ids,
    chr    = rep(paste0("Chr", seq_len(n.chr)), each = n.per.chr),
    pos    = rep(seq_len(n.per.chr), times = n.chr),
    stringsAsFactors = FALSE
)

cat("2a. primePanel() — 300 lines × 500 markers, 3 planted effects\n")
cat("    Planted markers:", paste(planted.marks, collapse = ", "), "\n")
cat("    Planted effects:", paste(planted.effects, collapse = ", "),
    " (units: yield)\n\n")

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
if (n_sig > 0)
    cat("    Keys:", paste(gwas_fit$QTL$qtl, collapse = ", "), "\n")
cat("\n")

if (n_sig == 0L) {
    cat("    No significant markers — skipping effects/contrast for gwasAim.\n\n")
} else {

    # ---- 2c. type = "effects" (gwasAim) -------------------------------------
    cat("2c. plot(gwas_fit, genObj = panel, type = 'effects')\n")
    cat("    For GWAS: positive effect = alt allele (dosage 2) favoured,\n")
    cat("    negative = ref allele (dosage 0) favoured.\n\n")

    p_geff <- plot(gwas_fit, genObj = panel, type = "effects")
    print(p_geff)
    save_plot(p_geff, "EC05_gwasAim_effects",
              width = max(6, n_sig * 1.6), height = 5)

    # ---- 2d. type = "contrast" (gwasAim) ------------------------------------
    cat("2d. plot(gwas_fit, genObj = panel, type = 'contrast', ncol = 1)\n")
    cat("    Allele classes: 0 (ref hom) / 1 (het) / 2 (alt hom).\n")
    cat("    y-axis = total genetic value.  Group mean ± SE connected by line.\n\n")

    p_gcon <- plot(gwas_fit, genObj = panel, type = "contrast",
                   data = yield.data, ncol = 1L)
    print(p_gcon)
    save_plot(p_gcon, "EC06_gwasAim_contrast",
              width = 8, height = n_sig * 3.0)

    # ---- 2e. Side-by-side comparison: planted vs recovered ------------------
    cat("2e. Planted vs recovered marker comparison\n")
    summ_gwas <- summary(gwas_fit, genObj = panel)
    cat("    Summary table:\n")
    print(summ_gwas)
    cat("\n")

    recovered <- summ_gwas$Marker
    for (pm in planted.marks) {
        hit <- any(sapply(recovered, function(r) grepl(pm, r, fixed = TRUE)))
        cat(sprintf("    %-8s  (planted effect %+.1f) — %s\n",
                    pm, planted.effects[match(pm, planted.marks)],
                    if (hit) "RECOVERED" else "not recovered"))
    }
    cat("\n")
}


# =============================================================================
# PART 3 — Edge-case guards
# =============================================================================
cat("------ PART 3: Edge-case guards ------\n\n")

cat("3a. Missing 'data' argument for type = 'contrast'\n")
tryCatch(
    plot(qtl_yld, genObj = rxk_int, type = "contrast"),
    error = function(e) cat("    Caught:", conditionMessage(e), "\n\n")
)

cat("3b. Object with no significant QTL\n")
base_tiny <- asreml::asreml(
    fixed     = yld ~ Range,
    random    = ~ Genotype + Row,
    data      = pheno,
    na.action = asreml::na.method(x = "include"),
    trace     = FALSE
)
qtl_none <- qtlAim(
    baseModel        = base_tiny,
    genObj           = rxk_int,
    merge.by         = "Genotype",
    exclusion.window = 50,
    TypeI            = 1e-10,   # impossibly strict — no QTL will pass
    trace            = FALSE
)
cat("    QTL detected at TypeI=1e-10:", length(qtl_none$QTL$qtl), "\n")
tryCatch(
    plot(qtl_none, genObj = rxk_int, type = "effects"),
    error = function(e) cat("    Caught:", conditionMessage(e), "\n\n")
)


# =============================================================================
# Summary
# =============================================================================
ec_plots <- list.files(out_dir, pattern = "^EC.*\\.png$")
cat("================================================================\n")
cat("  Demo complete —", length(ec_plots),
    "effects/contrast plots saved to", out_dir, "\n")
cat("================================================================\n\n")
for (f in sort(ec_plots))
    cat("  ", f, "\n")
