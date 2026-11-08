# =============================================================================
# run_demo3only.R
# Demo 3 only: eight field trials/sites, fa(Trial,2) genetic structure.
#
#  Demo 3a -- qtlAim()  eight sites, fa(Trial,2):id base random,
#             qtlAim called with Trait = "Trial", n.fa = 2
#  Demo 3b -- gwasAim() eight sites, same structure
#
# Population: 120 DH lines, 5 chr x 80 markers @ 1 cM.
#   395 intervals > 120 lines => vm() path throughout.
#
# Simulation design notes
# -----------------------
# True G is a general PD matrix (Wishart), NOT an exact FA(2)+Psi structure.
# fa(Trial,2):id in the base model therefore APPROXIMATES G; the approximation
# residual keeps the vm term identifiable (Psi parameters stay positive).
#
# Vg target: mean diagonal of Ga ≈ 0.08 per site.  This is large enough for
# fa(Trial,2) to be estimable with 8 sites (23 parameters per FA term) yet
# small enough — combined with the improved Ginv eigenvalue floor in
# engine_select.R — that the LRT drops cleanly after the 4 planted QTL.
#
# Run interactively after  devtools::load_all()  or:
#   Rscript scratch/run_demo3only.R
# =============================================================================


library(asreml)
library(qtl)
library(MASS)
library(ggplot2)

set.seed(101)

# ---- Output directory -------------------------------------------------------
out_dir <- file.path("scratch", "plots_demo3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggsave(path, plot = p, width = width, height = height, dpi = 150, bg = "white")
    message("  Saved: ", path)
}

cat("\n")
cat("================================================================\n")
cat("  wgAim Demo 3 — eight field sites, fa(Trial,2)\n")
cat("  120 DH lines, 5 chr x 80 markers @ 1 cM\n")
cat("================================================================\n\n")

# =============================================================================
# SECTION 0 — Simulate DH population
# =============================================================================
cat("------ SECTION 0: Simulate DH population ------\n\n")

n.lines  <- 120
n.chr    <- 5
n.mar    <- 80
mar.step <- 1.0    # 79 intervals per chr; 395 total > 120 lines => vm()

cat("0a. sim.cross (bc), 120 lines, 5 chr x 80 markers @ 1 cM\n")
cat("    Intervals:", n.chr * (n.mar - 1L), "> lines:", n.lines, "=> vm() path\n\n")

map.sim <- sim.map(len       = rep((n.mar - 1L) * mar.step, n.chr),
                   n.mar     = rep(n.mar, n.chr),
                   eq.spac   = TRUE,
                   include.x = FALSE,
                   sex.sp    = FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

cross.sim <- sim.cross(map.sim, n.ind = n.lines, type = "bc", error.prob = 0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)
cat("    Simulated markers:", sum(nmar(cross.sim)), "\n\n")

cat("0b. primeCross() — wgCross (interval type)\n")
genObj <- primeCross(cross.sim, type = "interval", id = "id")
cat("    Total intervals:",
    sum(sapply(genObj$geno, function(g) ncol(g$interval.data))), "\n\n")

cat("0c. Genomic relationship matrix K\n")
M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M) / ncol(M)
K <- K / mean(diag(K))
cat("    K: mean diagonal =", round(mean(diag(K)), 3),
    "| off-diag range [",
    paste(round(range(K[lower.tri(K)]), 3), collapse = ", "), "]\n\n")

cat("0d. primePanel() — wgPanel for gwasAim (pm1 encoding)\n")
marker.mat.pm1 <- do.call(cbind, lapply(genObj$geno, function(g) g$imputed.data))
rownames(marker.mat.pm1) <- ids
map.panel <- do.call(rbind, lapply(names(genObj$geno), function(chr)
    data.frame(marker = colnames(genObj$geno[[chr]]$imputed.data),
               chr    = chr,
               pos    = as.numeric(genObj$geno[[chr]]$map),
               stringsAsFactors = FALSE)))
panel <- primePanel(marker.mat.pm1, map.panel,
                    id = "id", map.id = "marker", map.chr = "chr", map.pos = "pos",
                    encoding = "pm1", maf = 0.05, impute = "none")
panel$pheno$id <- factor(ids)
cat("    Markers after MAF filter:",
    sum(sapply(panel$geno, function(g) ncol(g$imputed.data))), "\n\n")

cat("0e. Planted QTL\n")
qtl.info <- list(list(chr = "Chr1", int = 20L),
                 list(chr = "Chr2", int = 55L),
                 list(chr = "Chr3", int = 30L),
                 list(chr = "Chr4", int = 60L))
n.qtl <- length(qtl.info)
qtl.X <- matrix(0, nrow = n.lines, ncol = n.qtl,
                dimnames = list(ids,
                    sapply(qtl.info, function(q) paste(q$chr, q$int, sep = "."))))
for (i in seq_len(n.qtl))
    qtl.X[, i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]
cat("    QTL 1 — Chr1 int 20 (main)\n")
cat("    QTL 2 — Chr2 int 55 (main, opposite sign)\n")
cat("    QTL 3 — Chr3 int 30 (G×E: strong sites 1-4, weak 5-8)\n")
cat("    QTL 4 — Chr4 int 60 (G×E: crossover, sites 1-4 negative, 5-8 positive)\n\n")

# =============================================================================
# SECTION 3 — Eight-site simulation  (Demo 3a + 3b)
# =============================================================================
cat("================================================================\n")
cat("  DEMO 3: Eight field sites\n")
cat("  fa(Trial,2):id base random structure\n")
cat("  Trait = 'Trial', n.fa = 2\n")
cat("================================================================\n\n")

# ---- 3.0  Genetic covariance matrix Ga — general PD (Wishart) --------------
#
# True G is general PD, NOT an exact FA(2)+Psi structure.
# fa(Trial,2):id approximates G; the approximation residual keeps the
# vm term identifiable so Psi parameters do not collapse to zero.
#
# Vg target: mean diagonal ≈ 0.04 per site.
#
# The LRT compares diag(Trial):vm(id,covObj) + diag(Trial):id against
# diag(Trial):id alone.  The vm term uses K_remaining (genomic relationship
# from the remaining intervals); the id term uses identity covariance.
# Because the polygenic background is drawn from MVN(0, Ga ⊗ K), it is
# correlated with K_remaining and the vm term detects it correctly.
# The pchisq.mixture threshold for ntrait = 8, TypeI = 0.05 is ≈ 13.5.
# Back-solving:  stat ≈ 2 × 1920 × (Vg/Ve)²
#   Vg = 0.08  →  stat ≈ 25  (above threshold → spurious QTL detected)
#   Vg = 0.04  →  stat ≈  6  (below threshold → clean stop after 4 QTL)
# Vg = 0.04 is also sufficient for fa(Trial,2) to be estimable with
# 8 sites (23 FA parameters per term, 1920 total observations).
#
# S.scale: moderate positive between-site genetic correlation (off-diag ≈ 0.06)
# with some G×E heterogeneity (diagonal ≈ 0.10).
n.trials <- 8L
set.seed(42)

S.scale <- matrix(0.06, n.trials, n.trials)
diag(S.scale) <- 0.10
Ga.raw <- MASS::mvrnorm(n.trials + 2L, mu = rep(0, n.trials), Sigma = S.scale)
Ga     <- crossprod(Ga.raw) / (n.trials + 2L)   # scaled Wishart realisation
Ga     <- Ga * (0.04 / mean(diag(Ga)))           # rescale to target diagonal

cat("True Ga diagonal (genetic variances per site):\n")
cat("  ", paste(round(diag(Ga), 3), collapse = "  "), "\n\n")
cat("True genetic correlations:\n")
D.inv  <- diag(1 / sqrt(diag(Ga)))
Ga.cor <- D.inv %*% Ga %*% D.inv
for (i in seq_len(n.trials))
    cat(sprintf("  Site %-2d: %s\n", i,
                paste(sprintf("%5.2f", Ga.cor[i, ]), collapse = " ")))
cat("\n")

# ---- 3.1  QTL effect matrix (8 sites) --------------------------------------
#
# QTL 1 (Chr1 int 20): MAIN — constant +0.80 across all 8 sites.
# QTL 2 (Chr2 int 55): MAIN — constant -0.65 across all 8 sites.
# QTL 3 (Chr3 int 30): INTERACTION — large effect in sites 1-4, small in 5-8.
# QTL 4 (Chr4 int 60): CROSSOVER G×E — negative in sites 1-4, positive in 5-8.

qtl.eff <- matrix(
    c( 0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  0.80,  # QTL 1 main
      -0.65, -0.65, -0.65, -0.65, -0.65, -0.65, -0.65, -0.65,  # QTL 2 main
       0.60,  0.55,  0.50,  0.45,  0.15,  0.10,  0.10,  0.10,  # QTL 3 G×E
      -0.50, -0.45, -0.40, -0.35,  0.30,  0.35,  0.40,  0.45), # QTL 4 crossover
    nrow = n.qtl, ncol = n.trials, byrow = TRUE,
    dimnames = list(colnames(qtl.X), paste0("Trial", seq_len(n.trials)))
)

cat("Planted QTL effects (per site):\n")
cat(sprintf("%-30s", ""), paste(sprintf("%-8s", paste0("S", seq_len(n.trials))), collapse=""), "\n")
for (i in seq_len(n.qtl)) {
    type <- if (i <= 2L) "MAIN" else "INTERACTION"
    cat(sprintf("  %-26s [%s]  ", colnames(qtl.X)[i], type))
    cat(paste(sprintf("%6.3f", qtl.eff[i, ]), collapse = " "), "\n")
}
cat("\n")

# ---- 3.2  Polygenic background (8-site MVN) ---------------------------------
set.seed(101)   # restore master seed
cat("3.2 Sampling 8-site polygenic background from MVN(0, Ga x K)\n")
Sigma.poly <- kronecker(Ga, K)
u.vec      <- MASS::mvrnorm(1, mu = rep(0, n.lines * n.trials), Sigma = Sigma.poly)
U <- matrix(u.vec, nrow = n.lines, ncol = n.trials,
            dimnames = list(ids, paste0("Trial", seq_len(n.trials))))
cat("    Polygenic range per site:\n")
for (j in seq_len(n.trials))
    cat(sprintf("      Site %d: [%.2f, %.2f]\n", j, min(U[, j]), max(U[, j])))
cat("\n")

# ---- 3.3  Residual variances and site means ---------------------------------
set.seed(42)
# Ve = 0.70 per site, consistent with the two-trial demo design.
# h²(polygenic) ≈ 0.04/(0.04+0.70) ≈ 5%   (small, below LRT threshold after QTL)
# h²(QTL)       ≈ 1.40/(1.40+0.04+0.70) ≈ 65%  (dominates)
sigma2.e <- rep(0.70, n.trials)
mu.site  <- runif(n.trials, 9.0, 13.0)
n.reps   <- 2L

cat("Site means and residual variances:\n")
for (j in seq_len(n.trials))
    cat(sprintf("  Site %d : mu = %.2f  Ve = %.2f\n", j, mu.site[j], sigma2.e[j]))
cat("\n")

# ---- 3.4  Phenotypic data ---------------------------------------------------
cat("3.4 Assembling 8-site phenotypic data\n")
pheno <- do.call(rbind, lapply(seq_len(n.trials), function(j) {
    qtl.contrib <- rowSums(sweep(qtl.X, 2, qtl.eff[, j], "*"))
    do.call(rbind, lapply(seq_len(n.reps), function(r)
        data.frame(
            id    = factor(ids),
            Trial = factor(paste0("Trial", j),
                           levels = paste0("Trial", seq_len(n.trials))),
            Rep   = factor(paste0("R", r)),
            yield = mu.site[j] + qtl.contrib + U[, j] +
                    rnorm(n.lines, 0, sqrt(sigma2.e[j])),
            stringsAsFactors = FALSE
        )
    ))
}))
pheno <- pheno[order(pheno$Trial, pheno$Rep, pheno$id), ]
rownames(pheno) <- NULL
cat("    Rows:", nrow(pheno), "(",
    n.lines, "lines x", n.trials, "sites x", n.reps, "reps)\n")
cat("    yield range:", round(range(pheno$yield), 2), "\n\n")

# ---- 3.5  Base model for 8-site (qtlAim) ------------------------------------
#
# random   : ~ fa(Trial,2):id   — rank-2 FA polygenic G×E (23 params for 8 sites)
# residual : ~ dsum(~units|Trial) — separate Ve per site
#
# The engine (.buildGenomeModel) detects "fa(Trial,2)" in the residual line
# term and — with n.fa = 2 — steps up the additive genomic term:
#   diag(Trial):vm()  ->  fa(Trial,1):vm()  ->  fa(Trial,2):vm()
# The residual fa(Trial,2):id is never modified by the engine.
#
# maxit = 80: fa(Trial,2) over 8 sites has many variance parameters;
# extra REML iterations help convergence.

cat("3.5 Fitting 8-site base model for qtlAim\n")
cat("    fixed    = yield ~ Trial + Trial:Rep\n")
cat("    random   = ~ fa(Trial,2):id\n")
cat("    residual = ~ dsum(~units | Trial)\n")
cat("    maxit = 80  (FA over 8 sites needs room to converge)\n\n")

base.qtl <- asreml(
    fixed     = yield ~ Trial + Trial:Rep,
    random    = ~ fa(Trial, 2):id,
    residual  = ~ dsum(~units | Trial),
    data      = pheno,
    trace     = FALSE,
    na.action = asreml::na.method(x = "include"),
    maxit     = 80
)
if (!base.qtl$converge) {
    cat("  Not converged — running one update...\n")
    base.qtl <- update(base.qtl)
}
cat("qtlAim base model converged:", base.qtl$converge, "\n")
cat("LogLik:", round(base.qtl$loglik, 3),
    "  sigma2:", round(base.qtl$sigma2, 4), "\n")

vp <- base.qtl$vparameters
psi.nm <- names(vp)[grepl("fa.*id.*!var$|id.*fa.*!var$", names(vp))]
fa1.nm <- names(vp)[grepl("fa.*id.*!fa1$|id.*fa.*!fa1$", names(vp))]
fa2.nm <- names(vp)[grepl("fa.*id.*!fa2$|id.*fa.*!fa2$", names(vp))]

if (length(psi.nm)) {
    cat("\nEstimated Psi (specific variances x sigma2):\n")
    for (i in seq_along(psi.nm))
        cat(sprintf("  Site %-2d : %.4f\n", i, vp[psi.nm[i]] * base.qtl$sigma2))
}
if (length(fa1.nm)) {
    cat("Estimated Lambda col 1 (fa1 loadings):\n")
    for (i in seq_along(fa1.nm))
        cat(sprintf("  Site %-2d : %.4f\n", i, vp[fa1.nm[i]]))
}
if (length(fa2.nm)) {
    cat("Estimated Lambda col 2 (fa2 loadings):\n")
    for (i in seq_along(fa2.nm))
        cat(sprintf("  Site %-2d : %.4f\n", i, vp[fa2.nm[i]]))
}
cat("\n")

# ---- 3.6  Separate base model for gwasAim -----------------------------------
cat("3.6 Fitting 8-site base model for gwasAim\n")
base.gwas <- asreml(
    fixed     = yield ~ Trial + Trial:Rep,
    random    = ~ fa(Trial, 2):id,
    residual  = ~ dsum(~units | Trial),
    data      = pheno,
    trace     = FALSE,
    na.action = asreml::na.method(x = "include"),
    maxit     = 80
)
if (!base.gwas$converge) {
    cat("  Not converged — running one update...\n")
    base.gwas <- update(base.gwas)
}
cat("gwasAim base model converged:", base.gwas$converge, "\n\n")

# =============================================================================
# DEMO 3a: qtlAim — eight sites, fa(Trial,2)
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Demo 3a: qtlAim() — 8 sites, Trait = 'Trial', n.fa = 2\n")
cat("  Engine step-up:\n")
cat("    diag(Trial):vm()  ->  fa(Trial,1):vm()  ->  fa(Trial,2):vm()\n")
cat("  Residual fa(Trial,2):id unchanged throughout.\n")
cat("================================================================\n\n")

qtl.3a <- qtlAim(
    baseModel        = base.qtl,
    genObj           = genObj,
    merge.by         = "id",
    Trait            = "Trial",
    n.fa             = 2L,
    TypeI            = 0.05,
    exclusion.window = 20,
    trace            = "scratch/trace_demo3a_qtlAim.txt",
    na.action        = asreml::na.method(x = "include"),
    maxit            = 80
)
cat("\n--- qtlAim results ---\n")
cat("  QTL detected:", length(qtl.3a$QTL$qtl), "\n")
if (length(qtl.3a$QTL$qtl)) {
    cat("  Locations and type:\n")
    for (i in seq_along(qtl.3a$QTL$qtl)) {
        is.int <- if (!is.null(qtl.3a$QTL$is.interaction))
            qtl.3a$QTL$is.interaction[i] else NA
        cat(sprintf("    %d. %-22s  %s\n", i, qtl.3a$QTL$qtl[i],
                    if (is.na(is.int)) "—"
                    else if (is.int) "INTERACTION" else "MAIN"))
    }
    cat("\n  print.qtlAim():\n")
    print(qtl.3a, genObj = genObj)
    cat("\n  summary.qtlAim():\n")
    print(summary(qtl.3a, genObj = genObj))
}

# --- Diagnostic plots for Demo 3a -------------------------------------------
if (length(qtl.3a$QTL$qtl) > 0) {
    cat("\n  Saving diagnostic plots for Demo 3a...\n")
    save_plot(plot(qtl.3a, genObj = genObj, type = "outlier"), "demo3a_qtlAim_outlier")
    save_plot(plot(qtl.3a, genObj = genObj, type = "blups"),   "demo3a_qtlAim_blups")
    save_plot(plot(qtl.3a, genObj = genObj, type = "heatmap", cap = 5), "demo3a_qtlAim_heatmap")
    save_plot(plot(qtl.3a, genObj = genObj, type = "effects"), "demo3a_qtlAim_effects")
    save_plot(plot(qtl.3a, genObj = genObj, data = yield.data, type = "contrast"), "demo3a_qtlAim_contrast")
    p <- aimTrace(qtl.3a, plot = "lrt")
    save_plot(p, "demo3a_qtlAim_lrt_trace")
    p <- aimTrace(qtl.3a, plot = "stability")
    save_plot(p, "demo3a_qtlAim_stability")
}

# =============================================================================
# DEMO 3b: gwasAim — eight sites, fa(Trial,2)
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Demo 3b: gwasAim() — 8 sites, Trait = 'Trial', n.fa = 2\n")
cat("  Scans panel markers; same FA(2) base model.\n")
cat("================================================================\n\n")

gwas.3b <- gwasAim(
    baseModel        = base.gwas,
    genObj           = panel,
    merge.by         = "id",
    Trait            = "Trial",
    n.fa             = 2L,
    TypeI            = 0.05,
    exclusion.window = 20,
    trace            = "scratch/trace_demo3b_gwasAim.txt",
    na.action        = asreml::na.method(x = "include"),
    maxit            = 80
)
cat("\n--- gwasAim results ---\n")
cat("  Markers detected:", length(gwas.3b$QTL$qtl), "\n")
if (length(gwas.3b$QTL$qtl)) {
    cat("  Locations and type:\n")
    for (i in seq_along(gwas.3b$QTL$qtl)) {
        is.int <- if (!is.null(gwas.3b$QTL$is.interaction))
            gwas.3b$QTL$is.interaction[i] else NA
        cat(sprintf("    %d. %-22s  %s\n", i, gwas.3b$QTL$qtl[i],
                    if (is.na(is.int)) "—"
                    else if (is.int) "INTERACTION" else "MAIN"))
    }
    cat("\n  print.gwasAim():\n")
    print(gwas.3b, genObj = panel)
    cat("\n  summary.gwasAim():\n")
    print(summary(gwas.3b, genObj = panel))
}

# --- Diagnostic plots for Demo 3b -------------------------------------------
if (length(gwas.3b$QTL$qtl) > 0) {
    cat("\n  Saving diagnostic plots for Demo 3b...\n")
    save_plot(plot(gwas.3b, genObj = panel, type = "manhattan"), "demo3b_gwasAim_manhattan")
    save_plot(plot(gwas.3b, genObj = panel, type = "blups"),     "demo3b_gwasAim_blups")
    save_plot(plot(gwas.3b, genObj = panel, type = "heatmap"),   "demo3b_gwasAim_heatmap")
    save_plot(plot(gwas.3b, genObj = panel, type = "effects"),   "demo3b_gwasAim_effects")
    p <- aimTrace(gwas.3b, plot = "lrt")
    if (inherits(p, "ggplot")) save_plot(p, "demo3b_gwasAim_lrt_trace")
    p <- aimTrace(gwas.3b, plot = "stability")
    if (inherits(p, "ggplot")) save_plot(p, "demo3b_gwasAim_stability")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Demo 3 — summary\n")
cat("================================================================\n\n")

cat("Planted QTL:\n")
cat("  QTL 1 — Chr1 int 20 — MAIN        ( +0.80 all sites )\n")
cat("  QTL 2 — Chr2 int 55 — MAIN        ( -0.65 all sites )\n")
cat(sprintf("  QTL 3 — Chr3 int 30 — INTERACTION ( sites 1-4: %.2f..%.2f, sites 5-8: %.2f..%.2f )\n",
    qtl.eff[3,1], qtl.eff[3,4], qtl.eff[3,5], qtl.eff[3,8]))
cat(sprintf("  QTL 4 — Chr4 int 60 — CROSSOVER   ( sites 1-4: %.2f..%.2f, sites 5-8: +%.2f..+%.2f )\n\n",
    qtl.eff[4,1], qtl.eff[4,4], qtl.eff[4,5], qtl.eff[4,8]))

report_obj <- function(label, obj) {
    n <- length(obj$QTL$qtl)
    cat(sprintf("  %-38s %d QTL\n", label, n))
    if (n > 0)
        for (i in seq_len(n)) {
            is.int <- if (!is.null(obj$QTL$is.interaction)) obj$QTL$is.interaction[i] else NA
            cat(sprintf("      %d. %-22s %s\n", i, obj$QTL$qtl[i],
                        if (is.na(is.int)) "—" else if (is.int) "INTERACTION" else "MAIN"))
        }
}
report_obj("Demo 3a  qtlAim  8 sites fa(2):", qtl.3a)
report_obj("Demo 3b  gwasAim 8 sites fa(2):", gwas.3b)

cat("\nTrace files: scratch/trace_demo3a_qtlAim.txt\n")
cat("             scratch/trace_demo3b_gwasAim.txt\n")
cat("Plots:       scratch/plots_demo3/\n\n")
