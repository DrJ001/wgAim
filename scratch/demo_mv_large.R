# =============================================================================
# demo_mv_large.R
# Large-scale wgAim demonstration: 200 lines, ~1000 markers.
#
# Six analyses are demonstrated in order:
#
#  Demo 1a -- qtlAim()  single field trial, 200 lines, 2 reps
#  Demo 1b -- gwasAim() single field trial, 200 lines, 2 reps
#
#  Demo 2a -- qtlAim()  two field trials (same 200 lines, 2 reps each)
#             base model uses corh(Trial):Line residual genetic structure
#             qtlAim called with Trait = "Trial"
#  Demo 2b -- gwasAim() two field trials, same structure
#
#  Demo 3a -- qtlAim()  eight field trials (same 200 lines, 2 reps each)
#             base model uses fa(Trial,2):Line residual genetic structure
#             qtlAim called with Trait = "Trial"
#  Demo 3b -- gwasAim() eight field trials, same structure
#
# Simulation notes
# ----------------
#  * 200 DH lines, 5 chromosomes x 200 markers at 0.5 cM spacing (99.5 cM/chr).
#  * 5 x 199 = 995 intervals > 200 lines  =>  vm() path throughout.
#  * Same 200 lines appear in every trial.
#  * Plants 4 QTL/markers with known effects (two main, two with G x E).
#  * The SAME genomic relationship matrix K used in simulation is the one
#    .constructCM() will build internally, so models are correctly specified.
#
# Run interactively after  devtools::load_all()  or:
#   Rscript scratch/demo_mv_large.R
# =============================================================================

devtools::load_all()
library(asreml)
library(qtl)
library(MASS)
library(ggplot2)

set.seed(101)

out_dir <- file.path("scratch", "plots_mv_large")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggsave(path, plot = p, width = width, height = height, dpi = 150, bg = "white")
    message("  Saved: ", path)
}

# Convenience wrapper: run expr, report errors but continue.
# Warnings are printed but execution continues (withCallingHandlers allows
# resumption after a warning; tryCatch would abort).
try_run <- function(label, expr) {
    cat("\n  [", label, "]\n", sep = "")
    withCallingHandlers(
        tryCatch(
            expr,
            error = function(e) cat("  ERROR in", label, ":", conditionMessage(e), "\n")
        ),
        warning = function(w) {
            cat("  WARNING in", label, ":", conditionMessage(w), "\n")
            invokeRestart("muffleWarning")
        }
    )
}

cat("\n")
cat("================================================================\n")
cat("  wgAim large-scale demonstration (200 lines, ~1000 markers)\n")
cat("================================================================\n\n")

# =============================================================================
# SECTION 0 — Simulate population and phenotypic data
# =============================================================================
cat("------ SECTION 0: Simulate DH population and phenotypes ------\n\n")

# ---- 0a. Genetic map and cross ---------------------------------------------
#
# 5 chromosomes x 200 markers at 0.5 cM spacing = 99.5 cM per chromosome.
# Intervals per chr = 199;  total intervals = 5 x 199 = 995.
# 995 intervals > 200 lines  =>  vm() path will be used throughout.
n.lines  <- 200
n.chr    <- 5
n.mar    <- 200          # markers per chromosome
mar.step <- 0.5          # cM spacing  -> 99.5 cM chromosomes

cat("0a. Simulating DH cross\n")
cat("    Lines:", n.lines, " | Chromosomes:", n.chr,
    " | Markers per chr:", n.mar, " | Spacing:", mar.step, "cM\n")
cat("    Total intervals:", n.chr * (n.mar - 1L),
    " > ", n.lines, " lines  =>  vm() path\n\n")

map.sim <- sim.map(
    len       = rep((n.mar - 1L) * mar.step, n.chr),
    n.mar     = rep(n.mar, n.chr),
    eq.spac   = TRUE,
    include.x = FALSE,
    sex.sp    = FALSE
)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

cross.sim <- sim.cross(
    map        = map.sim,
    n.ind      = n.lines,
    type       = "bc",
    error.prob = 0.005
)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)
cat("    Simulated markers:", sum(nmar(cross.sim)), "\n\n")

# ---- 0b. primeCross() — wgCross for qtlAim ---------------------------------
cat("0b. primeCross() — building wgCross object\n")
genObj <- primeCross(cross.sim, type = "interval", id = "id")
cat("    Class:", paste(class(genObj), collapse = ", "), "\n")
cat("    Type:", genObj$type, "\n")
cat("    Chromosomes:", paste(names(genObj$geno), collapse = ", "), "\n\n")

# ---- 0c. Interval score matrix and relationship matrix ---------------------
cat("0c. Building genomic relationship matrix K\n")
gdat <- lapply(genObj$geno, function(el) el$interval.data)
M    <- do.call("cbind", gdat)          # 200 lines x 995 intervals
rownames(M) <- ids
n.int <- ncol(M)

K      <- tcrossprod(M) / n.int
K      <- K / mean(diag(K))            # mean diagonal = 1 (same as .constructCM)
cat("    Interval matrix:", nrow(M), "x", ncol(M), "\n")
cat("    K: mean diagonal =", round(mean(diag(K)), 3),
    " | off-diag range:", round(range(K[lower.tri(K)]), 3), "\n\n")

# ---- 0d. primePanel() — wgPanel for gwasAim ---------------------------------
#
# Build the panel from the cross marker data.
# primeCross() has already recoded genotypes to additive +/-1 in $imputed.data.
# Use encoding = "pm1" so primePanel() applies no further transformation.
# This avoids the singularity that arises when bc cross codes 1/2 are naively
# shifted to 0/1 (producing markers with no +1 allele class).
cat("0d. primePanel() — building wgPanel for gwasAim\n")
marker.mat.pm1 <- do.call(cbind, lapply(genObj$geno, function(g) g$imputed.data))
rownames(marker.mat.pm1) <- ids

map.panel <- do.call(rbind, lapply(names(genObj$geno), function(chr) {
    data.frame(
        marker = colnames(genObj$geno[[chr]]$imputed.data),
        chr    = chr,
        pos    = as.numeric(genObj$geno[[chr]]$map),
        stringsAsFactors = FALSE
    )
}))

panel <- primePanel(
    geno     = marker.mat.pm1,
    map      = map.panel,
    id       = "id",
    map.id   = "marker",
    map.chr  = "chr",
    map.pos  = "pos",
    encoding = "pm1",
    maf      = 0.05,
    impute   = "none"
)
panel$pheno$id <- factor(ids)

n.panel.markers <- sum(sapply(panel$geno, function(g) ncol(g$imputed.data)))
cat("    Markers after MAF filter:", n.panel.markers, "\n\n")

# ---- 0e. Plant QTL -----------------------------------------------------------
#
# Four QTL on four different chromosomes; effects given as per-allele shifts.
# In a DH population the two allele classes are -1 (B) and +1 (A).
#
#  QTL 1 (Chr1, int 50): MAIN EFFECT    -- same in all trials
#  QTL 2 (Chr2, int 80): MAIN EFFECT    -- same in all trials (opposite sign)
#  QTL 3 (Chr3, int 30): INTERACTION    -- large in trials 1-4, small in 5-8
#  QTL 4 (Chr4, int 120):INTERACTION    -- sign reversal (crossover G x E)
#
# For the single-trial demos (1a/1b) all four QTL contribute.
# For multi-trial demos (2x, 3x) the interaction pattern varies by trial.

cat("0e. Planted QTL positions and main-effect sizes\n")
qtl.info <- list(
    list(chr = "Chr1", int = 50L,  effect.main =  0.80),
    list(chr = "Chr2", int = 80L,  effect.main = -0.65),
    list(chr = "Chr3", int = 30L,  effect.main =  0.55),
    list(chr = "Chr4", int = 120L, effect.main = -0.50)
)
n.qtl <- length(qtl.info)

# Extract interval genotype scores at each QTL position
qtl.X <- matrix(0, nrow = n.lines, ncol = n.qtl,
                dimnames = list(ids, sapply(qtl.info, function(q)
                    paste("Chr", q$chr, q$int, sep = "."))))
for (i in seq_len(n.qtl))
    qtl.X[, i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]

for (i in seq_len(n.qtl))
    cat("    QTL", i, "-", qtl.info[[i]]$chr, "int", qtl.info[[i]]$int,
        ": main effect =", qtl.info[[i]]$effect.main, "\n")
cat("\n")

# Convenience: total QTL contribution per line (main effects only, used in 1x demos)
qtl.main.contrib <- rowSums(
    sweep(qtl.X, 2,
          sapply(qtl.info, function(q) q$effect.main), "*"))

# ---- 0f. Residual polygenic effects ----------------------------------------
cat("0f. Sampling residual polygenic effects (single-trial)\n")
sigma2.g.uni <- 0.60     # polygenic genetic variance (single-trial)
sigma2.e.uni <- 0.40     # residual error variance (single-trial)

# u ~ MVN(0, sigma2.g * K), the polygenic background
u.uni <- MASS::mvrnorm(1,
                       mu    = rep(0, n.lines),
                       Sigma = sigma2.g.uni * K)
names(u.uni) <- ids
cat("    sigma2.g =", sigma2.g.uni, " | sigma2.e =", sigma2.e.uni,
    " | h2 ~=", round(sigma2.g.uni / (sigma2.g.uni + sigma2.e.uni), 2), "\n\n")

# =============================================================================
# DEMO 1: Single field trial
# =============================================================================
cat("================================================================\n")
cat("  DEMO 1: Single field trial — qtlAim and gwasAim\n")
cat("  200 lines, 2 complete reps (RCBD)\n")
cat("================================================================\n\n")

# ---- Build single-trial phenotypic data ------------------------------------
#
# Model: y_ir = mu + rep_r + QTL_i + poly_i + e_ir
#   rep_r   : fixed block effect (2 levels)
#   QTL_i   : sum of four QTL contributions (fixed signal)
#   poly_i  : polygenic BLUP u.uni[i]  (random background)
#   e_ir    : iid N(0, sigma2.e.uni)
mu.uni   <- 10.0
rep.effs.uni <- c(0.25, -0.25)     # two block effects that sum to zero

pheno.uni <- do.call(rbind, lapply(seq_len(2L), function(r) {
    data.frame(
        id    = factor(ids),
        Rep   = factor(r),
        yield = mu.uni + rep.effs.uni[r] +
                qtl.main.contrib +
                u.uni +
                rnorm(n.lines, 0, sqrt(sigma2.e.uni)),
        stringsAsFactors = FALSE
    )
}))
# Sort Line-then-Rep for ASReml
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep), ]
rownames(pheno.uni) <- NULL
cat("Single-trial data: ", nrow(pheno.uni), "rows (", n.lines, "lines x 2 reps)\n")
cat("  yield range:", round(range(pheno.uni$yield), 2), "\n\n")

# ---- Demo 1a: qtlAim single trial ------------------------------------------
cat("------ Demo 1a: qtlAim() — single field trial ------\n\n")
cat("  Base model:\n")
cat("    fixed  = yield ~ Rep\n")
cat("    random = ~ Line\n\n")

base.uni <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ id,
    data      = pheno.uni,
    trace     = FALSE
)
cat("  Base model converged:", base.uni$converge, "\n")
cat("  Vg:", round(base.uni$vparameters["id"] * base.uni$sigma2, 3),
    " | Ve:", round(base.uni$sigma2, 3), "\n\n")

qtl.1a <- try_run("Demo 1a — qtlAim", {
    fit <- qtlAim(
        baseModel        = base.uni,
        genObj           = genObj,
        merge.by         = "id",
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo1a_qtlAim.txt"
    )
    cat("  QTL detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n\n")
        cat("  print.qtlAim():\n")
        print(fit, genObj = genObj)
        cat("\n  summary.qtlAim():\n")
        print(summary(fit, genObj = genObj))
    }
    fit
})

# ---- Demo 1b: gwasAim single trial -----------------------------------------
cat("\n------ Demo 1b: gwasAim() — single field trial ------\n\n")
cat("  Same phenotypic data. primePanel wgPanel used.\n")
cat("  Base model:\n")
cat("    fixed  = yield ~ Rep\n")
cat("    random = ~ id    (matching panel$pheno$id column)\n\n")

base.gwas.uni <- asreml(
    fixed  = yield ~ Rep,
    random = ~ id,
    data   = pheno.uni,
    trace  = FALSE
)
cat("  Base model converged:", base.gwas.uni$converge, "\n\n")

gwas.1b <- try_run("Demo 1b — gwasAim", {
    fit <- gwasAim(
        baseModel        = base.gwas.uni,
        genObj           = panel,
        merge.by         = "id",
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo1b_gwasAim.txt"
    )
    cat("  Markers detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n\n")
        cat("  print.gwasAim():\n")
        print(fit, genObj = panel)
        cat("\n  summary.gwasAim():\n")
        print(summary(fit, genObj = panel))
    }
    fit
})

# =============================================================================
# DEMO 2: Two field trials with corh() genetic structure
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  DEMO 2: Two field trials — qtlAim and gwasAim\n")
cat("  200 lines, 2 reps per trial, corh(Trial) residual genetic structure\n")
cat("  Trait = 'Trial' in qtlAim / gwasAim\n")
cat("================================================================\n\n")

# ---- Simulate two-trial phenotypic data ------------------------------------
#
# Two trials.  The residual genetic structure in the BASE model is
# corh(Trial):Line — this estimates separate genetic variances per trial
# AND the between-trial genetic correlation.
#
# Engine behaviour when Trait = "Trial":
#   .buildGenomeModel() reads 'corh(Trial):Line' from the base random formula.
#   struct.prefix = "corh(Trial)"
#   upgrade.to.corgh = TRUE   (corh triggers upgrade)
#   Step 1: fits diag(Trial):vm(Line, covObj) + corh(Trial):Line
#   Step 2: upgrades additive term to corgh(Trial):vm(Line, covObj)
#
# For gwasAim the logic is identical.

n.trials.2 <- 2L

# Genetic covariance matrix for the two-trial case
sigma2.g2 <- c(0.55, 0.50)      # genetic variances per trial
r.g2      <- 0.75               # between-trial genetic correlation
Ga2 <- matrix(c(sigma2.g2[1],
                r.g2 * sqrt(prod(sigma2.g2)),
                r.g2 * sqrt(prod(sigma2.g2)),
                sigma2.g2[2]),
              2, 2)

# QTL effects per trial
#   QTL 1, 2: main effects (same across trials)
#   QTL 3: present in trial 1, reduced in trial 2
#   QTL 4: sign reversal (crossover G x E)
qtl.eff.2 <- matrix(
    c( 0.80,  0.80,    # QTL 1 main
      -0.65, -0.65,    # QTL 2 main
       0.55,  0.15,    # QTL 3 interaction
      -0.50,  0.45),   # QTL 4 crossover
    nrow = n.qtl, ncol = n.trials.2, byrow = TRUE,
    dimnames = list(colnames(qtl.X), paste0("Trial", seq_len(n.trials.2)))
)

cat("2-trial QTL effects:\n")
cat("                   Trial1  Trial2  Type\n")
cat("  QTL 1 (Chr1-50)  ", sprintf("%.3f", qtl.eff.2[1,1]),
    " ", sprintf("%.3f", qtl.eff.2[1,2]), "  Main\n")
cat("  QTL 2 (Chr2-80) ", sprintf("%.3f", qtl.eff.2[2,1]),
    " ", sprintf("%.3f", qtl.eff.2[2,2]), "  Main\n")
cat("  QTL 3 (Chr3-30)  ", sprintf("%.3f", qtl.eff.2[3,1]),
    " ", sprintf("%.3f", qtl.eff.2[3,2]), "  Interaction\n")
cat("  QTL 4 (Chr4-120)", sprintf("%.3f", qtl.eff.2[4,1]),
    " ", sprintf("%.3f", qtl.eff.2[4,2]), "  Crossover GxE\n\n")

# Polygenic G x E from MVN(0, Ga2 x K)
Sigma.poly2 <- kronecker(Ga2, K)
u2.vec      <- MASS::mvrnorm(1, mu = rep(0, n.lines * n.trials.2),
                              Sigma = Sigma.poly2)
U2 <- matrix(u2.vec, nrow = n.lines, ncol = n.trials.2,
             dimnames = list(ids, paste0("Trial", seq_len(n.trials.2))))

sigma2.e2 <- c(0.30, 0.35)
mu2       <- c(10.0, 11.5)
n.reps.2  <- 2L

pheno.2 <- do.call(rbind, lapply(seq_len(n.trials.2), function(j) {
    qtl.contrib <- rowSums(sweep(qtl.X, 2, qtl.eff.2[, j], "*"))
    do.call(rbind, lapply(seq_len(n.reps.2), function(r) {
        data.frame(
            id    = factor(ids),
            Trial = factor(paste0("Trial", j)),
            Rep   = factor(paste0("R", r)),
            yield = mu2[j] + qtl.contrib + U2[, j] +
                    rnorm(n.lines, 0, sqrt(sigma2.e2[j])),
            stringsAsFactors = FALSE
        )
    }))
}))
# Data MUST be sorted Trial-first for dsum(~units | Trial)
pheno.2 <- pheno.2[order(pheno.2$Trial, pheno.2$Rep, pheno.2$id), ]
rownames(pheno.2) <- NULL

cat("2-trial phenotypic data: ", nrow(pheno.2), "rows (",
    n.lines, "lines x", n.trials.2, "trials x", n.reps.2, "reps)\n\n")

# ---- Base model for 2-trial ------------------------------------------------
#
# The critical structure is corh(Trial):Line in the random formula.
# This tells .buildGenomeModel() to:
#   (a) start with diag(Trial):vm(Line, covObj) as the genomic term
#   (b) upgrade the genomic term to corgh(Trial):vm(Line, covObj)
#       (corh in input -> corgh on the vm term; both are hetero-corr structures)
# The residual line term corh(Trial):Line is NEVER changed — it stays as-is.
#
# fixed    : yield ~ Trial + Trial:Rep  (trial means + block effects per trial)
# random   : ~ corh(Trial):Line         (polygenic G x E with correlation)
# residual : ~ dsum(~units | Trial)     (separate Ve per trial)
cat("2-trial base model:\n")
cat("  fixed    = yield ~ Trial + Trial:Rep\n")
cat("  random   = ~ corh(Trial):id\n")
cat("  residual = ~ dsum(~units | Trial)\n\n")

base.2 <- asreml(
    fixed    = yield ~ Trial + Trial:Rep,
    random   = ~ corh(Trial):id,
    residual = ~ dsum(~units | Trial),
    data     = pheno.2,
    trace    = FALSE
)
cat("Base model converged:", base.2$converge, "\n")

# Report variance components
vp2 <- base.2$vparameters * base.2$sigma2
cat("Variance components (corh(Trial):id):\n")
for (nm in names(vp2)[grep("id", names(vp2))])
    cat(" ", nm, ":", round(vp2[nm], 4), "\n")
cat(" sigma2 (Ve):", round(base.2$sigma2, 4), "\n\n")

# ---- Demo 2a: qtlAim two trials --------------------------------------------
cat("------ Demo 2a: qtlAim() — two trials, Trait = 'Trial' ------\n\n")
cat("  Engine will upgrade additive genomic term:\n")
cat("  diag(Trial):vm()  ->  corgh(Trial):vm()\n\n")

qtl.2a <- try_run("Demo 2a — qtlAim two trials", {
    fit <- qtlAim(
        baseModel        = base.2,
        genObj           = genObj,
        merge.by         = "id",
        Trait            = "Trial",
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo2a_qtlAim.txt"
    )
    cat("  QTL detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n")
        if (!is.null(fit$QTL$is.interaction)) {
            for (i in seq_along(fit$QTL$qtl))
                cat("    ", fit$QTL$qtl[i], ":",
                    ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN"), "\n")
        }
        cat("\n  print.qtlAim():\n")
        print(fit, genObj = genObj)
        cat("\n  summary.qtlAim():\n")
        print(summary(fit, genObj = genObj))
    }
    fit
})

# ---- Demo 2b: gwasAim two trials -------------------------------------------
cat("\n------ Demo 2b: gwasAim() — two trials, Trait = 'Trial' ------\n\n")

cat("  Base model:\n")
cat("  fixed    = yield ~ Trial + Trial:Rep\n")
cat("  random   = ~ corh(Trial):id\n")
cat("  residual = ~ dsum(~units | Trial)\n\n")

base.2.gwas <- asreml(
    fixed    = yield ~ Trial + Trial:Rep,
    random   = ~ corh(Trial):id,
    residual = ~ dsum(~units | Trial),
    data     = pheno.2,
    trace    = FALSE
)
cat("Base model converged:", base.2.gwas$converge, "\n\n")

gwas.2b <- try_run("Demo 2b — gwasAim two trials", {
    fit <- gwasAim(
        baseModel        = base.2.gwas,
        genObj           = panel,
        merge.by         = "id",
        Trait            = "Trial",
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo2b_gwasAim.txt"
    )
    cat("  Markers detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n")
        if (!is.null(fit$QTL$is.interaction)) {
            for (i in seq_along(fit$QTL$qtl))
                cat("    ", fit$QTL$qtl[i], ":",
                    ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN"), "\n")
        }
        cat("\n  print.gwasAim():\n")
        print(fit, genObj = panel)
        cat("\n  summary.gwasAim():\n")
        print(summary(fit, genObj = panel))
    }
    fit
})

# =============================================================================
# DEMO 3: Eight field trials with fa(Trial,2) genetic structure
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  DEMO 3: Eight field trials — qtlAim and gwasAim\n")
cat("  200 lines, 2 reps per trial, fa(Trial,2) residual genetic structure\n")
cat("  Trait = 'Trial' with n.fa = 2 in qtlAim / gwasAim\n")
cat("================================================================\n\n")

# ---- Simulate eight-trial phenotypic data ----------------------------------
#
# Eight trials sharing the same 200 DH lines.
# The residual line structure in the BASE model is fa(Trial,2):Line.
# Factor-analytic structures are appropriate for large METs where the
# between-trial genetic covariance is complex (multiple axes of G x E).
#
# Engine behaviour when Trait = "Trial", n.fa = 2:
#   .buildGenomeModel() reads 'fa(Trial,2):Line' from the base random formula.
#   struct.prefix = "fa(Trial,2)"
#   upgrade.to.fa = TRUE
#   Steps:
#     1. Fit diag(Trial):vm() + fa(Trial,2):Line   (initial)
#     2. Upgrade: diag(Trial):vm() -> fa(Trial,1):vm()
#     3. Upgrade: fa(Trial,1):vm() -> fa(Trial,2):vm()
#   The residual fa(Trial,2):Line term is NEVER changed.
#
# Ga8 is an 8x8 covariance matrix constructed from a factor-analytic
# structure to keep simulation tractable.

n.trials.8 <- 8L

# Build an 8x8 genetic covariance matrix via FA(2) structure
# Lambda: loadings on 2 latent factors (each trial's sensitivity)
set.seed(42)
Lambda8 <- matrix(c(
    0.70, 0.15,    # Trial 1
    0.65, 0.20,    # Trial 2
    0.60, 0.30,    # Trial 3
    0.55, 0.35,    # Trial 4
    0.50, 0.40,    # Trial 5
    0.70,-0.10,    # Trial 6
    0.60,-0.20,    # Trial 7
    0.55,-0.30     # Trial 8
), nrow = n.trials.8, ncol = 2L, byrow = TRUE)
Psi8  <- diag(runif(n.trials.8, 0.05, 0.15))  # specific variances
Ga8   <- Lambda8 %*% t(Lambda8) + Psi8
sigma2.g8 <- diag(Ga8)

# QTL effect matrix for 8 trials
# QTL 1 and 2: main effects (constant across all 8 trials)
# QTL 3: interaction — scales with first FA loading
# QTL 4: interaction — scales with second FA loading (sign reversal in trials 6-8)
qtl.eff.8 <- matrix(0, nrow = n.qtl, ncol = n.trials.8,
                    dimnames = list(colnames(qtl.X),
                                    paste0("Trial", seq_len(n.trials.8))))
qtl.eff.8[1, ] <-  0.80               # QTL 1: main
qtl.eff.8[2, ] <- -0.65               # QTL 2: main
qtl.eff.8[3, ] <- Lambda8[, 1] * 1.0  # QTL 3: scales with factor 1
qtl.eff.8[4, ] <- Lambda8[, 2] * 1.2  # QTL 4: scales with factor 2 (has negative loading)

cat("8-trial QTL effects (mean across trials):\n")
for (i in seq_len(n.qtl))
    cat("  QTL", i, "(", colnames(qtl.X)[i], "):",
        "mean =", round(mean(qtl.eff.8[i,]), 3),
        " range [", round(min(qtl.eff.8[i,]), 3),
        ",", round(max(qtl.eff.8[i,]), 3), "]\n")
cat("\n")

# Polygenic G x E from MVN(0, Ga8 x K)
Sigma.poly8 <- kronecker(Ga8, K)
u8.vec      <- MASS::mvrnorm(1, mu = rep(0, n.lines * n.trials.8),
                              Sigma = Sigma.poly8)
U8 <- matrix(u8.vec, nrow = n.lines, ncol = n.trials.8,
             dimnames = list(ids, paste0("Trial", seq_len(n.trials.8))))

sigma2.e8 <- runif(n.trials.8, 0.25, 0.45)
mu8       <- runif(n.trials.8, 9.0, 13.0)
n.reps.8  <- 2L

pheno.8 <- do.call(rbind, lapply(seq_len(n.trials.8), function(j) {
    qtl.contrib <- rowSums(sweep(qtl.X, 2, qtl.eff.8[, j], "*"))
    do.call(rbind, lapply(seq_len(n.reps.8), function(r) {
        data.frame(
            id    = factor(ids),
            Trial = factor(paste0("Trial", j)),
            Rep   = factor(paste0("R", r)),
            yield = mu8[j] + qtl.contrib + U8[, j] +
                    rnorm(n.lines, 0, sqrt(sigma2.e8[j])),
            stringsAsFactors = FALSE
        )
    }))
}))
# Sort Trial-first for dsum(~units | Trial)
pheno.8 <- pheno.8[order(pheno.8$Trial, pheno.8$Rep, pheno.8$id), ]
rownames(pheno.8) <- NULL

cat("8-trial phenotypic data:", nrow(pheno.8), "rows (",
    n.lines, "lines x", n.trials.8, "trials x", n.reps.8, "reps)\n\n")

# ---- Base model for 8-trial ------------------------------------------------
#
# fa(Trial,2):Line is the residual genetic structure.
# This is a rank-2 factor-analytic structure over 8 trials.
# The additive genomic term will be stepped up:
#   diag(Trial):vm() -> fa(Trial,1):vm() -> fa(Trial,2):vm()
#
# fixed    : yield ~ Trial + Trial:Rep
# random   : ~ fa(Trial,2):Line
# residual : ~ dsum(~units | Trial)
cat("8-trial base model:\n")
cat("  fixed    = yield ~ Trial + Trial:Rep\n")
cat("  random   = ~ fa(Trial,2):id\n")
cat("  residual = ~ dsum(~units | Trial)\n\n")
cat("  (Fitting fa(Trial,2) with 8 trials may take a moment...)\n\n")

base.8 <- asreml(
    fixed    = yield ~ Trial + Trial:Rep,
    random   = ~ fa(Trial, 2):id,
    residual = ~ dsum(~units | Trial),
    data     = pheno.8,
    trace    = FALSE
)
# If not converged, try one more update
if (!base.8$converge) {
    cat("  Base model not converged on first fit; running one update...\n")
    base.8 <- update(base.8)
}
cat("Base model converged:", base.8$converge, "\n")

# Report the specific variances (Psi) and factor loadings (Lambda)
vp8 <- base.8$vparameters
psi.nms <- names(vp8)[grep("fa.*Line.*!var|var.*Line.*fa", names(vp8))]
lam.nms <- names(vp8)[grep("fa.*Line.*fa[0-9]|fa[0-9].*Line.*fa", names(vp8))]
if (length(psi.nms)) {
    cat("  Specific variances (Psi):", paste(round(vp8[psi.nms] * base.8$sigma2, 3), collapse=", "), "\n")
}
if (length(lam.nms)) {
    cat("  Factor loadings (Lambda):", paste(round(vp8[lam.nms], 3), collapse=", "), "\n")
}
cat("\n")

# ---- Demo 3a: qtlAim eight trials ------------------------------------------
cat("------ Demo 3a: qtlAim() — eight trials, Trait = 'Trial', n.fa = 2 ------\n\n")
cat("  Engine will step the additive genomic term up:\n")
cat("  diag(Trial):vm()  ->  fa(Trial,1):vm()  ->  fa(Trial,2):vm()\n")
cat("  Residual fa(Trial,2):Line is unchanged throughout.\n\n")

qtl.3a <- try_run("Demo 3a — qtlAim eight trials", {
    fit <- qtlAim(
        baseModel        = base.8,
        genObj           = genObj,
        merge.by         = "id",
        Trait            = "Trial",
        n.fa             = 2L,         # step up to fa(Trial,2) on the vm term
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo3a_qtlAim.txt"
    )
    cat("  QTL detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n")
        if (!is.null(fit$QTL$is.interaction)) {
            for (i in seq_along(fit$QTL$qtl))
                cat("    ", fit$QTL$qtl[i], ":",
                    ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN"), "\n")
        }
        cat("\n  print.qtlAim():\n")
        print(fit, genObj = genObj)
        cat("\n  summary.qtlAim():\n")
        print(summary(fit, genObj = genObj))
    }
    fit
})

# ---- Demo 3b: gwasAim eight trials -----------------------------------------
cat("\n------ Demo 3b: gwasAim() — eight trials, Trait = 'Trial', n.fa = 2 ------\n\n")

cat("  Base model:\n")
cat("  fixed    = yield ~ Trial + Trial:Rep\n")
cat("  random   = ~ fa(Trial,2):id\n")
cat("  residual = ~ dsum(~units | Trial)\n\n")

base.8.gwas <- asreml(
    fixed    = yield ~ Trial + Trial:Rep,
    random   = ~ fa(Trial, 2):id,
    residual = ~ dsum(~units | Trial),
    data     = pheno.8,
    trace    = FALSE
)
if (!base.8.gwas$converge) {
    cat("  Base model not converged on first fit; running one update...\n")
    base.8.gwas <- update(base.8.gwas)
}
cat("Base model converged:", base.8.gwas$converge, "\n\n")

gwas.3b <- try_run("Demo 3b — gwasAim eight trials", {
    fit <- gwasAim(
        baseModel        = base.8.gwas,
        genObj           = panel,
        merge.by         = "id",
        Trait            = "Trial",
        n.fa             = 2L,
        TypeI            = 0.05,
        exclusion.window = 20,
        trace            = "scratch/trace_demo3b_gwasAim.txt"
    )
    cat("  Markers detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        cat("  Locations:", paste(fit$QTL$qtl, collapse = ", "), "\n")
        if (!is.null(fit$QTL$is.interaction)) {
            for (i in seq_along(fit$QTL$qtl))
                cat("    ", fit$QTL$qtl[i], ":",
                    ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN"), "\n")
        }
        cat("\n  print.gwasAim():\n")
        print(fit, genObj = panel)
        cat("\n  summary.gwasAim():\n")
        print(summary(fit, genObj = panel))
    }
    fit
})

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Demo summary\n")
cat("================================================================\n\n")
cat("  True QTL planted:\n")
cat("    QTL 1 (Chr1, int 50): main effect  +0.80 (all trials)\n")
cat("    QTL 2 (Chr2, int 80): main effect  -0.65 (all trials)\n")
cat("    QTL 3 (Chr3, int 30): interaction  (varies by trial)\n")
cat("    QTL 4 (Chr4, int120): interaction  (sign reversal in Demo 3)\n\n")

report <- function(label, obj, genObj.arg) {
    n <- if (!is.null(obj) && !inherits(obj, "try-error") && length(obj$QTL$qtl))
             length(obj$QTL$qtl)
         else 0L
    cat("  ", label, ": ", n, " QTL/marker(s) detected", sep = "")
    if (n > 0) cat(" —", paste(obj$QTL$qtl, collapse = ", "))
    cat("\n")
}
report("Demo 1a  qtlAim  single trial", qtl.1a,  genObj)
report("Demo 1b  gwasAim single trial", gwas.1b, panel)
report("Demo 2a  qtlAim  2 trials (corh)", qtl.2a, genObj)
report("Demo 2b  gwasAim 2 trials (corh)", gwas.2b, panel)
report("Demo 3a  qtlAim  8 trials (fa,2)", qtl.3a, genObj)
report("Demo 3b  gwasAim 8 trials (fa,2)", gwas.3b, panel)

cat("\n  Notes:\n")
cat("  - All traces written to scratch/trace_demo*.txt\n")
cat("  - Demos 1x: single trial, univariate (Trait = NULL)\n")
cat("  - Demos 2x: corh() in base model -> engine upgrades vm term to corgh()\n")
cat("  - Demos 3x: fa(,2) in base model -> engine steps vm term: diag->fa(1)->fa(2)\n")
cat("  - Any errors above are caught and reported; remaining demos still run\n\n")
