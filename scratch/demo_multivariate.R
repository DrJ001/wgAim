# =============================================================================
# demo_multivariate.R
# Demonstration of univariate and multivariate QTL analysis using wgAim.
#
# The multivariate (multi-environment) case requires a phenotypic dataset
# whose genetic structure is consistent with the models we will fit.
# Specifically we need:
#
#   1. A biparental DH population with an interval-based wgCross object.
#   2. Polygenic G x E effects sampled from MVN(0, Ga x K), where
#      K = MM'/p is the genomic relationship matrix and Ga is a 2x2
#      genetic covariance matrix.  This ensures the base ASReml model
#      (diag(Site):Line) and the qtlAim genome-wide term (vm():diag(Site))
#      are appropriate for the data.
#   3. QTL effects of three kinds:
#        (a) Main-effect QTL : same effect in both sites
#        (b) Interaction QTL : large effect in Site 1, negligible in Site 2
#        (c) A third QTL with intermediate G x E behaviour
#   4. Independent residual error per site.
#
# The demo then:
#   Part 1 -- fits qtlAim() univariately for each site separately
#   Part 2 -- fits qtlAim() multivariately with Trait = "Site"
#   Part 3 -- compares the two approaches
# =============================================================================

devtools::load_all()
library(asreml)
library(qtl)       # sim.cross()
library(MASS)      # mvrnorm()
library(ggplot2)

set.seed(2025)

out_dir <- file.path("scratch", "plots_mv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 14, height = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggsave(path, plot = p, width = width, height = height, dpi = 150, bg = "white")
    message("  Saved: ", path)
}

cat("\n")
cat("================================================================\n")
cat("  wgAim — Univariate vs Multivariate QTL Demo\n")
cat("================================================================\n\n")

# =============================================================================
# SECTION 0 — Simulate a DH population and phenotypic data
# =============================================================================
cat("------ SECTION 0: Simulate DH population + multi-site phenotypes ------\n\n")

# ---- 0a. Simulate a biparental DH cross via qtl::sim.cross() ---------------
#
# 100 DH lines, 4 chromosomes, 30 markers per chromosome at 3 cM spacing.
# 30 markers gives 29 intervals per chromosome = 4 x 29 = 116 intervals.
# 116 intervals > 100 lines  =>  qtlAim() will use the vm() path.
#
# For the vm() path (genome-wide relationship matrix), the number of intervals
# (columns of M) must exceed the number of lines (rows of M). If intervals <
# lines the engine automatically switches to the mbf() path instead.
cat("0a. Simulating DH cross (100 lines, 4 chr x 30 markers @ 3 cM)\n")
cat("    116 intervals > 100 lines  =>  vm() path will be used\n\n")

n.lines  <- 100
n.chr    <- 4
n.mar    <- 30          # markers per chromosome
mar.step <- 3           # cM between markers (87 cM per chromosome)

map.sim <- sim.map(
    len     = rep((n.mar - 1) * mar.step, n.chr),
    n.mar   = rep(n.mar, n.chr),
    eq.spac = TRUE,
    include.x = FALSE,
    sex.sp  = FALSE
)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

cross.sim <- sim.cross(
    map      = map.sim,
    n.ind    = n.lines,
    type     = "dh",
    error.prob = 0.01
)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
n.int.expected <- n.chr * (n.mar - 1L)
cat("    Lines:", n.lines, " | Chromosomes:", n.chr,
    " | Markers per chr:", n.mar, "\n")
cat("    Expected intervals:", n.int.expected, "(>", n.lines, "lines => vm path)\n")
cat("    Total markers:", sum(sapply(cross.sim$geno, function(g) ncol(g$data))), "\n\n")

# ---- 0b. primeCross() — build the wgCross object ---------------------------
cat("0b. primeCross() — building wgCross object (interval type)\n")
genObj <- primeCross(cross.sim, type = "interval", id = "id")
cat("    Class:", paste(class(genObj), collapse = ", "), "\n")
cat("    Type:", genObj$type, "\n")
cat("    Chromosomes:", names(genObj$geno), "\n\n")

# ---- 0c. Extract the interval-score matrix ----------------------------------
#
# We need this now to build the genomic relationship matrix K.
# The interval score matrix M has rows = lines, columns = intervals.
# Scores are in {-1, +1} for DH.
gdat <- lapply(genObj$geno, function(el) el$interval.data)
M    <- do.call("cbind", gdat)          # 200 lines x 56 intervals
n.int <- ncol(M)
ids   <- as.character(genObj$pheno$id)
rownames(M) <- ids

cat("0c. Interval score matrix: ", nrow(M), "lines x", ncol(M), "intervals\n")
cat("    Score range:", range(M), "\n\n")

# ---- 0d. Build genomic relationship matrix K = MM'/p -----------------------
#
# This is the SAME matrix the vm() path will build internally via .constructCM().
# By using it here for simulation we ensure the polygenic G x E component is
# exactly what the model will try to fit.
p <- ncol(M)
K <- tcrossprod(M) / p            # 200 x 200 additive kinship matrix
scale.K <- mean(diag(K))
K <- K / scale.K                  # standardise so mean diagonal = 1

cat("0d. Genomic relationship matrix K (", nrow(K), "x", ncol(K), ")\n")
cat("    Mean diagonal:", round(mean(diag(K)), 3),
    " | Range off-diagonal:", round(range(K[lower.tri(K)]), 3), "\n\n")

# ---- 0e. Define QTL effects -----------------------------------------------
#
# We choose 3 QTL on 3 different chromosomes.
# Each QTL is at a known interval-column index.
#
# QTL 1 (Chr1, interval 7):  MAIN EFFECT -- same positive effect in both sites
# QTL 2 (Chr2, interval 10): INTERACTION -- large in Site 1, near-zero in Site 2
# QTL 3 (Chr3, interval 4):  INTERACTION -- negative in Site 1, positive in Site 2
#                             (sign reversal = classic crossover G x E)
#
# Effect sizes are expressed as per-allele effects (A vs B DH class, -1 vs +1).
# A value of 0.6 means one allele substitution shifts the phenotype by 0.6 sigma.

qtl.chrs  <- c("Chr1", "Chr2", "Chr3")
qtl.ints  <- c(7L,     10L,    4L)       # interval index within chromosome
qtl.names <- paste("Chr", qtl.chrs, qtl.ints, sep = ".")

# Effect matrix: rows = QTL, cols = sites.  effect[i,j] = effect of QTL i in site j
qtl.effects <- matrix(
    c( 0.7,  0.7,    # QTL 1: main effect (same both sites)
       0.9,  0.05,   # QTL 2: interaction (large site1, absent site2)
      -0.6,  0.55),  # QTL 3: crossover G x E
    nrow = 3, ncol = 2, byrow = TRUE,
    dimnames = list(qtl.names, c("Site1", "Site2"))
)

cat("0e. QTL effects:\n")
cat("    QTL             Site1   Site2   Type\n")
cat("    ----------      -----   -----   ----\n")
cat("   ", qtl.names[1], "   ", formatC(qtl.effects[1,1], 3, format="f"),
    " ", formatC(qtl.effects[1,2], 3, format="f"), "  Main effect\n")
cat("   ", qtl.names[2], "   ", formatC(qtl.effects[2,1], 3, format="f"),
    " ", formatC(qtl.effects[2,2], 3, format="f"), "  Interaction (absent Site2)\n")
cat("   ", qtl.names[3], "   ", formatC(qtl.effects[3,1], 3, format="f"),
    " ", formatC(qtl.effects[3,2], 3, format="f"), "  Crossover GxE\n\n")

# Extract the genotype scores at QTL positions
n.sites <- 2
qtl.X <- matrix(0, nrow = n.lines, ncol = 3)
colnames(qtl.X) <- qtl.names
rownames(qtl.X) <- ids
for (i in 1:3) {
    chr.dat    <- genObj$geno[[qtl.chrs[i]]]$interval.data
    qtl.X[,i]  <- chr.dat[, qtl.ints[i]]
}

# ---- 0f. Polygenic G x E effects -------------------------------------------
#
# Genetic covariance matrix for residual polygenic G x E:
#   Ga[1,1] = variance of polygenic effects in site 1
#   Ga[2,2] = variance of polygenic effects in site 2
#   Ga[1,2] = genetic covariance (correlation ~ 0.7, typical of wheat METs)
#
# We sample u = vec(U) ~ MVN(0, Ga x K)
# where U is n.lines x n.sites and u is stacked column-by-column.
# This is what the vm():diag(Site) random term models.

sigma2.g <- c(0.50, 0.45)   # genetic variances per site
r.g      <- 0.70             # genetic correlation between sites (on same scale)
Ga <- matrix(c(sigma2.g[1],
               r.g * sqrt(prod(sigma2.g)),
               r.g * sqrt(prod(sigma2.g)),
               sigma2.g[2]),
             2, 2)

cat("0f. Polygenic G x E covariance matrix Ga:\n")
cat("        Site1   Site2\n")
cat("  Site1", formatC(Ga[1,1], 4, format="f"), " ", formatC(Ga[1,2], 4, format="f"), "\n")
cat("  Site2", formatC(Ga[2,1], 4, format="f"), " ", formatC(Ga[2,2], 4, format="f"), "\n")
cat("  Genetic correlation:", round(cov2cor(Ga)[1,2], 3), "\n\n")

# Sample polygenic effects as vec(U) ~ MVN(0, Ga x K).
#
# U is n.lines x n.sites.  The stacked vector vec(U) has covariance
# kronecker(Ga, K) -- a (n.lines*n.sites) x (n.lines*n.sites) matrix.
# For 100 lines x 2 sites this is 200x200, which is trivial.
# MASS::mvrnorm() handles near-singular K gracefully via its own tolerance.
Sigma.poly <- kronecker(Ga, K)
u.vec      <- MASS::mvrnorm(1, mu = rep(0, n.lines * n.sites),
                             Sigma = Sigma.poly)
U.poly <- matrix(u.vec, nrow = n.lines, ncol = n.sites,
                  dimnames = list(ids, c("Site1", "Site2")))

# ---- 0g. Replication and residual errors ------------------------------------
#
# Without replication within each site ASReml cannot separate genetic variance
# from residual variance -- the model will be singular or fail to converge.
# We simulate n.reps complete blocks (RCBD) per site.
#
# The model for observation (i, j, r) -- line i, site j, rep r:
#   y_ijr = mu_j + rep_jr + QTL_ij + polyGxE_ij + epsilon_ijr
#
# rep_jr  : fixed block/rep effect within site (absorbed into fixed formula)
# epsilon : independent N(0, sigma2.e[j]) per observation

n.reps   <- 3              # complete blocks per site
sigma2.e <- c(0.30, 0.35)  # residual variances (site 1 and 2)

cat("0g. Replication: n.reps =", n.reps, "per site\n")
cat("    Residual variance Site1:", sigma2.e[1],
    " Site2:", sigma2.e[2], "\n\n")

# ---- 0h. Compute phenotypic values -----------------------------------------
#
# Genetic signal (QTL + polygenic G x E) is the same across reps within a
# site for each line.  Only the within-plot residual epsilon varies by rep.
#
# Y.rep[i, j, r] = mu[j] + rep_effect[j,r] + QTL[i,j] + polyGxE[i,j] + e[i,j,r]
#
# We absorb rep_effect into the fixed model (fixed = yield ~ Site + Site:Rep)
# and just add iid errors per observation.

mu       <- c(10, 12)     # overall site means
rep.effs <- list(
    Site1 = c( 0.3, -0.2,  0.0) * sqrt(sigma2.e[1]),   # 3 block effects, Site 1
    Site2 = c(-0.1,  0.2, -0.1) * sqrt(sigma2.e[2])    # 3 block effects, Site 2
)

# Compute genetic contribution once per line x site
genet <- matrix(0, n.lines, n.sites, dimnames = list(ids, c("Site1","Site2")))
for (j in 1:n.sites) {
    qtl.contrib <- rowSums(sweep(qtl.X, 2, qtl.effects[,j], "*"))
    genet[,j]   <- qtl.contrib + U.poly[,j]
}

# Build 3D array: lines x sites x reps
Y.rep <- array(
    0, dim = c(n.lines, n.sites, n.reps),
    dimnames = list(ids, c("Site1","Site2"), paste0("R", 1:n.reps))
)
for (j in 1:n.sites) {
    site.nm <- c("Site1","Site2")[j]
    for (r in 1:n.reps) {
        e.r         <- rnorm(n.lines, 0, sqrt(sigma2.e[j]))
        Y.rep[,j,r] <- mu[j] + rep.effs[[site.nm]][r] + genet[,j] + e.r
    }
}

cat("0h. Simulated phenotypes (sample statistics, averaged over reps):\n")
cat("    Site1: mean =", round(mean(Y.rep[,"Site1",]), 2),
    " sd =", round(sd(apply(Y.rep[,"Site1",], 1, mean)), 3), "(line means)\n")
cat("    Site2: mean =", round(mean(Y.rep[,"Site2",]), 2),
    " sd =", round(sd(apply(Y.rep[,"Site2",], 1, mean)), 3), "(line means)\n")
cat("    Between-site correlation of line means:",
    round(cor(apply(Y.rep[,"Site1",], 1, mean),
              apply(Y.rep[,"Site2",], 1, mean)), 3), "\n\n")

# ---- 0i. Build the phenotypic data frames ----------------------------------
#
# Long format: one row per line x site x rep combination.
# Column layout:
#   Line  -- DH line factor
#   Site  -- factor with levels "Site1" / "Site2"  (this is 'Trait' in qtlAim)
#   Rep   -- replicate/block factor, nested within Site
#   yield -- numeric response
#
# Total rows: n.lines x n.sites x n.reps = 100 x 2 x 3 = 600
#
# IMPORTANT: data must be sorted by Site for the dsum(~units | Site)
# residual structure.  We sort Site-first, then Line, then Rep.

pheno.long <- data.frame(
    Line  = factor(rep(rep(ids, times = n.reps), times = n.sites)),
    Site  = factor(rep(c("Site1","Site2"), each = n.lines * n.reps)),
    Rep   = factor(rep(rep(paste0("R", 1:n.reps), each = n.lines), times = n.sites)),
    yield = c(as.vector(Y.rep[,"Site1", ]),
              as.vector(Y.rep[,"Site2", ])),
    stringsAsFactors = FALSE
)
# Within each Site block the order is: R1 lines 1..100, R2 lines 1..100, ...
# Sort by Site, then by Line (to keep lines together within site is not required
# but tidy); the critical requirement is Site-contiguous for dsum.
pheno.long <- pheno.long[order(pheno.long$Site, pheno.long$Rep, pheno.long$Line), ]
rownames(pheno.long) <- NULL

# Per-site subsets for the univariate analyses
pheno.s1 <- pheno.long[pheno.long$Site == "Site1", c("Line", "Rep", "yield")]
pheno.s2 <- pheno.long[pheno.long$Site == "Site2", c("Line", "Rep", "yield")]

cat("0i. Long-format phenotypic data (with replication):\n")
cat("    Total rows:", nrow(pheno.long),
    " (", n.lines, "lines x", n.sites, "sites x", n.reps, "reps)\n")
cat("    Sites:", levels(pheno.long$Site), "\n")
cat("    Reps:", levels(pheno.long$Rep), "\n\n")

# =============================================================================
# SECTION 1 — Univariate qtlAim per site
# =============================================================================
cat("================================================================\n")
cat("  SECTION 1: Univariate qtlAim — one analysis per site\n")
cat("================================================================\n\n")

# ---- 1a. Site 1 ------------------------------------------------------------
#
# With replication, the fixed model includes Rep to absorb block effects.
# This allows ASReml to estimate Vg and Ve separately.
cat("1a. Site 1 — base model: yield ~ Rep, random = ~ Line\n")
base.s1 <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ Line,
    data      = pheno.s1,
    trace     = FALSE
)
cat("    Converged:", base.s1$converge, "\n")
cat("    Vg:", round(base.s1$vparameters["Line"], 3),
    " Ve:", round(base.s1$sigma2, 3), "\n")
cat("    h2 (site 1):", round(base.s1$vparameters["Line"] /
    (base.s1$vparameters["Line"] + 1), 3), "\n\n")

cat("    qtlAim() — Site 1\n")
qtl.s1 <- qtlAim(
    baseModel        = base.s1,
    genObj           = genObj,
    merge.by         = "Line",
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE
)
cat("    QTL found:", length(qtl.s1$QTL$qtl), "\n")

summ.s1 <- summary(qtl.s1, genObj = genObj)
if (!is.null(summ.s1)) {
    cat("\n    Summary — Site 1:\n")
    print(summ.s1)
} else {
    cat("    No significant QTL at TypeI = 0.05 for Site 1.\n")
}

# ---- 1b. Site 2 ------------------------------------------------------------
cat("\n1b. Site 2 — base model: yield ~ Rep, random = ~ Line\n")
base.s2 <- asreml(
    fixed     = yield ~ Rep,
    random    = ~ Line,
    data      = pheno.s2,
    trace     = FALSE
)
cat("    Converged:", base.s2$converge, "\n")
cat("    Vg:", round(base.s2$vparameters["Line"], 3),
    " Ve:", round(base.s2$sigma2, 3), "\n")
cat("    h2 (site 2):", round(base.s2$vparameters["Line"] /
    (base.s2$vparameters["Line"] + 1), 3), "\n\n")

cat("    qtlAim() — Site 2\n")
qtl.s2 <- qtlAim(
    baseModel        = base.s2,
    genObj           = genObj,
    merge.by         = "Line",
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE
)
cat("    QTL found:", length(qtl.s2$QTL$qtl), "\n")

summ.s2 <- summary(qtl.s2, genObj = genObj)
if (!is.null(summ.s2)) {
    cat("\n    Summary — Site 2:\n")
    print(summ.s2)
} else {
    cat("    No significant QTL at TypeI = 0.05 for Site 2.\n")
}

# ---- 1c. Univariate outlier plots ------------------------------------------
cat("\n1c. Outlier plots — univariate Site 1 and Site 2\n")

if (length(qtl.s1$QTL$qtl) > 0) {
    p.s1 <- plot(qtl.s1, genObj = genObj, type = "outlier", sig.col = "firebrick") +
        ggtitle("Site 1 — Univariate QTL Outlier Statistics")
    save_plot(p.s1, "UV01_Site1_outlier",
              height = 3 * length(qtl.s1$QTL$diag$oint))
}

if (length(qtl.s2$QTL$qtl) > 0) {
    p.s2 <- plot(qtl.s2, genObj = genObj, type = "outlier", sig.col = "steelblue") +
        ggtitle("Site 2 — Univariate QTL Outlier Statistics")
    save_plot(p.s2, "UV02_Site2_outlier",
              height = 3 * length(qtl.s2$QTL$diag$oint))
}

# ---- 1d. aimTable — combine both univariate summaries ----------------------
if (length(qtl.s1$QTL$qtl) > 0 && length(qtl.s2$QTL$qtl) > 0) {
    cat("\n1d. aimTable() — stacked summary of both univariate analyses\n")
    uv.table <- aimTable(qtl.s1, qtl.s2,
                         genObj = genObj,
                         labels = c("Site1", "Site2"))
    print(uv.table)
}

# =============================================================================
# SECTION 2 — Multivariate qtlAim with Trait = "Site"
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  SECTION 2: Multivariate qtlAim — Trait = 'Site'\n")
cat("================================================================\n\n")

# ---- 2a. Multivariate base model -------------------------------------------
#
# The base model must have a Line term that can absorb the polygenic G x E.
# We use diag(Site):Line which fits separate genetic variances per site
# with the same line-level effects (i.e., separate but not uncorrelated).
# This is equivalent to the 'diag' structure that .buildGenomeModel()
# will upgrade to when it adds :diag(Site) to the vm()/mbf() term.
#
# The residual uses dsum(~units | Site) for separate residual variances
# by site, consistent with our simulation (sigma2.e[1] != sigma2.e[2]).
#
# NOTE: In ASReml-R the data for a multivariate model must be sorted by
# the grouping factor (Site) to ensure correct blocking.

# pheno.long is already sorted by Site (done in 0i), so use it directly.
# Rename to pheno.mv for clarity and ensure factors.
pheno.mv       <- pheno.long
pheno.mv$Line  <- factor(pheno.mv$Line)
pheno.mv$Site  <- factor(pheno.mv$Site)
pheno.mv$Rep   <- factor(pheno.mv$Rep)

cat("2a. Multivariate base model:\n")
cat("    fixed    : yield ~ Site + Site:Rep  (site means + block effects)\n")
cat("    random   : ~ diag(Site):Line         (residual polygenic G x E)\n")
cat("    residual : ~ dsum(~units | Site)     (separate Ve per site)\n")
cat("\n")
cat("    Site:Rep captures block effects that differ between sites.\n")
cat("    After qtlAim() the random formula becomes:\n")
cat("      ~ diag(Site):vm(Line, covObj)  -- composite additive G x E\n")
cat("       + diag(Site):Line             -- residual polygenic G x E\n\n")

base.mv <- asreml(
    fixed    = yield ~ Site + Site:Rep,
    random   = ~ diag(Site):Line,
    residual = ~ dsum(~ units | Site),
    data     = pheno.mv,
    trace    = FALSE
)
cat("    Converged:", base.mv$converge, "\n")

# Extract and report genetic variance components per site
vp <- base.mv$vparameters
vp.site <- vp[grep("diag.*Site.*Line|Line.*Site.*diag", names(vp))]
cat("    Genetic variances (diag(Site):Line):\n")
for (i in seq_along(vp.site))
    cat("      ", names(vp.site)[i], ":",
        round(vp.site[i] * base.mv$sigma2, 4), "\n")
cat("    Residual sigma2:", round(base.mv$sigma2, 4), "\n\n")

# ---- 2b. qtlAim() — multivariate -------------------------------------------
cat("2b. qtlAim() — multivariate with Trait = 'Site'\n\n")

cat("    The engine will:\n")
cat("     - Fit vm(Line, covObj):diag(Site) as the genome-wide term\n")
cat("     - Use the multivariate quadratic-form outlier statistic:\n")
cat("       oint_i = t(atilde_i) %*% Ginv %*% atilde_i / vqtilde_i\n")
cat("     - Test significance via pchisq.mixture(stat, ntrait=2)\n")
cat("     - Add X.chr.idx + Site:X.chr.idx as fixed effects per iteration\n")
cat("     - Run waldTest() at the end to prune non-significant interactions\n\n")

cat("    The engine builds the genome-wide random formula as:\n")
cat("      diag(Site):vm(Line, covObj)   -- composite additive G x E term\n")
cat("     + diag(Site):Line              -- residual polygenic G x E term\n")
cat("    The base model's 'diag(Site):Line' term is stripped by .fixLines()\n")
cat("    (it contains merge.by = 'Line') and then re-added explicitly by\n")
cat("    .buildGenomeModel() as the residual polygenic term alongside the\n")
cat("    composite genomic term. This is the correct bivariate MET structure.\n\n")

qtl.mv <- qtlAim(
    baseModel        = base.mv,
    genObj           = genObj,
    merge.by         = "Line",
    Trait            = "Site",      # <-- the new argument
    TypeI            = 0.05,
    exclusion.window = 15,
    trace            = FALSE
)

cat("\n    Class:", paste(class(qtl.mv), collapse = ", "), "\n")
cat("    QTL found:", length(qtl.mv$QTL$qtl), "\n")
cat("    Trait:", qtl.mv$QTL$Trait, "\n")

if (!is.null(qtl.mv$QTL$is.interaction)) {
    cat("\n    Wald test — interaction vs main effect:\n")
    for (i in seq_along(qtl.mv$QTL$qtl)) {
        cat("      QTL", i, "(", qtl.mv$QTL$qtl[i], "): ",
            ifelse(qtl.mv$QTL$is.interaction[i], "INTERACTION", "MAIN EFFECT"),
            "\n")
    }
}

# ---- 2c. summary.qtlAim — multivariate -------------------------------------
cat("\n2c. summary() — multivariate (Trait-aware output)\n")
summ.mv <- summary(qtl.mv, genObj = genObj)
if (!is.null(summ.mv)) {
    cat("\n")
    print(summ.mv)
    cat("\n    Note: 'Trait' column shows 'MAIN' for main-effect QTL,\n")
    cat("          or the Site level name for interaction QTL.\n")
} else {
    cat("    No significant QTL detected.\n")
}

# ---- 2d. Multivariate outlier plot -----------------------------------------
if (length(qtl.mv$QTL$qtl) > 0) {
    cat("\n2d. Outlier plot — multivariate analysis\n")
    cat("    Note: the outlier statistic is now the multivariate quad-form,\n")
    cat("    NOT the univariate qtilde^2/vqtilde.\n")
    p.mv.oint <- plot(qtl.mv, genObj = genObj, type = "outlier",
                      sig.col = "darkorange") +
        ggtitle("Multivariate QTL Outlier Statistics (Site1 + Site2 joint)")
    save_plot(p.mv.oint, "MV01_multivariate_outlier",
              height = 3 * length(qtl.mv$QTL$diag$oint))
}

# ---- 2e. Heatmap comparison -------------------------------------------------
if (length(qtl.mv$QTL$qtl) > 0) {
    cat("\n2e. Heatmap — multivariate iterations\n")
    p.mv.heat <- plot(qtl.mv, genObj = genObj, type = "heatmap",
                      sig.col = "darkorange") +
        ggtitle("Multivariate QTL Heatmap (outlier statistic across iterations)")
    save_plot(p.mv.heat, "MV02_multivariate_heatmap",
              width = 16,
              height = 3 + length(qtl.mv$QTL$diag$oint) * 0.6)
}

# ---- 2f. n.fa demo — upgrade to corh (bivariate correlation structure) -----
cat("\n2f. Upgrading to corh(Site) variance structure (n.fa=1 for 2 traits)\n")
cat("    This estimates the genetic correlation between sites directly\n")
cat("    as part of the genome-wide composite term.\n\n")

qtl.mv.corh <- tryCatch({
    qtlAim(
        baseModel        = base.mv,
        genObj           = genObj,
        merge.by         = "Line",
        Trait            = "Site",
        n.fa             = 1L,       # upgrades diag -> corh for ntrait=2
        TypeI            = 0.05,
        exclusion.window = 15,
        trace            = FALSE
    )
}, error = function(e) {
    cat("    Note: corh model failed (may need better starting values):",
        conditionMessage(e), "\n")
    NULL
})

if (!is.null(qtl.mv.corh)) {
    cat("    corh model converged:", qtl.mv.corh$converge, "\n")
    cat("    QTL found:", length(qtl.mv.corh$QTL$qtl), "\n")

    # Extract the estimated genetic correlation from the final model
    vp.corh <- qtl.mv.corh$vparameters
    corh.idx <- grep("corh.*Site|Site.*corh", names(vp.corh))
    if (length(corh.idx)) {
        cat("    Estimated genetic correlation (corh):",
            round(vp.corh[corh.idx[length(corh.idx)]], 3), "\n")
        cat("    (True simulated r_g:", round(r.g, 3), ")\n")
    }
}

# =============================================================================
# SECTION 3 — Comparison: univariate vs multivariate
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  SECTION 3: Comparison — Univariate vs Multivariate\n")
cat("================================================================\n\n")

# ---- 3a. Which QTL did each approach find? ---------------------------------
cat("3a. QTL detection summary\n\n")

cat("  True QTL planted:\n")
cat("  QTL                Effect type           Site1   Site2\n")
cat("  -----------------  --------------------  ------  ------\n")
cat("  Chr.Chr1.7         Main effect           ",
    formatC(qtl.effects[1,1], 3, format="f"), " ",
    formatC(qtl.effects[1,2], 3, format="f"), "\n")
cat("  Chr.Chr2.10        Interaction (absent2) ",
    formatC(qtl.effects[2,1], 3, format="f"), " ",
    formatC(qtl.effects[2,2], 3, format="f"), "\n")
cat("  Chr.Chr3.4         Crossover GxE         ",
    formatC(qtl.effects[3,1], 3, format="f"), " ",
    formatC(qtl.effects[3,2], 3, format="f"), "\n\n")

cat("  Univariate Site1 detected:",
    if (length(qtl.s1$QTL$qtl) > 0) paste(qtl.s1$QTL$qtl, collapse=", ") else "none", "\n")
cat("  Univariate Site2 detected:",
    if (length(qtl.s2$QTL$qtl) > 0) paste(qtl.s2$QTL$qtl, collapse=", ") else "none", "\n")
cat("  Multivariate detected:    ",
    if (length(qtl.mv$QTL$qtl) > 0) paste(qtl.mv$QTL$qtl, collapse=", ") else "none", "\n\n")

# ---- 3b. Comparison plot: site-specific effects from multivariate summary --
if (!is.null(summ.mv) && nrow(summ.mv) > 0) {
    cat("3b. Effect size comparison plot\n")

    # Build a comparison data frame from the three summaries
    make_effect_df <- function(summ, site, type) {
        if (is.null(summ) || nrow(summ) == 0) return(NULL)
        nm <- if (type == "mv") as.character(summ$Trait) else site
        data.frame(
            QTL    = as.character(summ$Chromosome),
            Effect = as.numeric(summ$Size),
            Site   = nm,
            Type   = type,
            stringsAsFactors = FALSE
        )
    }

    df.s1 <- make_effect_df(summ.s1, "Site1", "univariate")
    df.s2 <- make_effect_df(summ.s2, "Site2", "univariate")
    df.mv <- make_effect_df(summ.mv, NA, "multivariate")

    df.comp <- do.call(rbind, list(df.s1, df.s2, df.mv))

    if (!is.null(df.comp) && nrow(df.comp) > 0) {
        p.comp <- ggplot(df.comp,
                         aes(x = QTL, y = Effect, fill = Site, alpha = Type)) +
            geom_col(position = position_dodge(width = 0.8), width = 0.7) +
            geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
            scale_fill_manual(values = c(
                Site1 = "steelblue", Site2 = "firebrick",
                MAIN  = "darkorange", Site1.1 = "steelblue", Site2.1 = "firebrick"
            )) +
            scale_alpha_manual(values = c(univariate = 0.6, multivariate = 1.0)) +
            labs(title = "QTL Effect Sizes: Univariate vs Multivariate",
                 subtitle = "Bars: univariate (transparent) and multivariate (solid)\nTrue planted effects: Chr1=main, Chr2=int(absent Site2), Chr3=crossover",
                 x = "Chromosome / QTL", y = "Estimated Effect",
                 fill = "Site", alpha = "Analysis") +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 30, hjust = 1))
        save_plot(p.comp, "CP01_effect_comparison", width = 10, height = 6)
    }
}

# ---- 3c. Power and precision commentary ------------------------------------
cat("3c. Why does the multivariate approach help?\n\n")
cat("  1. POWER for main-effect QTL:\n")
cat("     The joint model has n.lines x n.sites x n.reps =",
    n.lines * n.sites * n.reps, "observations.\n")
cat("     Univariate per-site uses only", n.lines * n.reps, "observations.\n")
cat("     For a main QTL with equal effects, the joint model is more powerful.\n\n")
cat("  2. POWER for interaction QTL:\n")
cat("     QTL 2 (Chr2) has a large effect in Site1 but is near-zero in Site2.\n")
cat("     The univariate Site1 analysis should detect it.\n")
cat("     The multivariate analysis detects it AND correctly classifies it as\n")
cat("     an interaction QTL via the post-loop waldTest().\n\n")
cat("  3. CROSSOVER GxE (QTL 3, Chr3):\n")
cat("     In Site1 the effect is", qtl.effects[3,1],
    "(negative) and in Site2 it is", qtl.effects[3,2], "(positive).\n")
cat("     The univariate analyses may detect it in one or both sites,\n")
cat("     but ONLY the multivariate approach correctly identifies the\n")
cat("     sign reversal as a significant Site x QTL interaction.\n\n")
cat("  4. GENETIC CORRELATION ESTIMATION:\n")
cat("     The n.fa=1 (corh) model estimates the between-site genetic\n")
cat("     correlation from the genome-wide composite term itself,\n")
cat("     exploiting the marker information across all loci simultaneously.\n\n")

# ---- 3d. LRT statistics from each approach ---------------------------------
cat("3d. LRT statistics per iteration\n\n")

print_lrt <- function(qtl.obj, label) {
    cat(" ", label, "\n")
    if (!is.null(qtl.obj$QTL$diag$lik.mat)) {
        lm <- qtl.obj$QTL$diag$lik.mat
        for (i in seq_len(nrow(lm)))
            cat(sprintf("    Iter %d: LRT stat = %.3f  p-value = %.4f  %s\n",
                        i, lm[i, "Statistic"], lm[i, "Pvalue"],
                        if (i <= length(qtl.obj$QTL$qtl)) paste("-> QTL:", qtl.obj$QTL$qtl[i]) else "(not significant)"))
    } else {
        cat("    No LRT info available.\n")
    }
    cat("\n")
}

print_lrt(qtl.s1, "Univariate Site 1:")
print_lrt(qtl.s2, "Univariate Site 2:")
print_lrt(qtl.mv, "Multivariate (joint):")

# =============================================================================
# SECTION 4 — Multivariate gwasAim demo (brief)
# =============================================================================
cat("================================================================\n")
cat("  SECTION 4: Multivariate gwasAim — brief demo\n")
cat("================================================================\n\n")

cat("4a. Building a simulated GWAS panel for the same 100 DH lines\n")
cat("    (In practice you'd use a diversity panel, but the same DH\n")
cat("    lines let us reuse the simulated phenotypes.)\n\n")

# Build a panel from the marker data already in the cross object
cat("    Extracting marker genotype matrix from DH cross...\n")
marker.mat.raw <- do.call(cbind,
    lapply(cross.sim$geno, function(g) g$data)
)
# Recode qtl cross codes 1/2 -> 0/2
marker.mat.012 <- (marker.mat.raw - 1) * 2   # 1->0, 2->2 (acts like 0/2 for DH)
rownames(marker.mat.012) <- ids

map.panel <- do.call(rbind, lapply(names(map.sim), function(chr) {
    data.frame(marker = names(map.sim[[chr]]),
               chr    = chr,
               pos    = as.numeric(map.sim[[chr]]),
               stringsAsFactors = FALSE)
}))

panel.gwas <- primePanel(
    geno     = marker.mat.012,
    map      = map.panel,
    id       = "id",
    map.id   = "marker",
    map.chr  = "chr",
    map.pos  = "pos",
    encoding = "012",
    maf      = 0.10,
    impute   = "none"
)
panel.gwas$pheno$id <- factor(ids)
cat("    Panel:", nrow(marker.mat.012), "lines,", ncol(marker.mat.012), "markers\n")
cat("    After MAF filter:", sum(sapply(panel.gwas$geno, function(g) ncol(g$imputed.data))),
    "markers\n\n")

# ---- 4b. Multivariate gwasAim base model ------------------------------------
cat("4b. Multivariate gwasAim base model: yield ~ Site + Site:Rep, ~ diag(Site):id\n")

# Rename Line -> id to match the panel's merge.by column
pheno.gwas.mv <- pheno.mv
names(pheno.gwas.mv)[names(pheno.gwas.mv) == "Line"] <- "id"
pheno.gwas.mv$id <- factor(pheno.gwas.mv$id)

base.gwas.mv <- asreml(
    fixed    = yield ~ Site + Site:Rep,
    random   = ~ diag(Site):id,
    residual = ~ dsum(~ units | Site),
    data     = pheno.gwas.mv,
    trace    = FALSE
)
cat("    Converged:", base.gwas.mv$converge, "\n\n")

cat("4c. gwasAim() — multivariate with Trait = 'Site'\n\n")

gwas.mv <- tryCatch(
    gwasAim(
        baseModel        = base.gwas.mv,
        genObj           = panel.gwas,
        merge.by         = "id",
        Trait            = "Site",
        TypeI            = 0.05,
        exclusion.window = 15,
        trace            = FALSE
    ),
    error = function(e) {
        cat("    gwasAim multivariate failed:", conditionMessage(e), "\n")
        NULL
    }
)

if (!is.null(gwas.mv)) {
    cat("    Markers found:", length(gwas.mv$QTL$qtl), "\n")
    if (!is.null(gwas.mv$QTL$is.interaction))
        cat("    Interaction markers:",
            sum(gwas.mv$QTL$is.interaction, na.rm = TRUE), "\n")

    summ.gwas.mv <- summary(gwas.mv, genObj = panel.gwas)
    if (!is.null(summ.gwas.mv)) {
        cat("\n    Multivariate GWAS summary:\n")
        print(summ.gwas.mv)
    }

    if (length(gwas.mv$QTL$qtl) > 0) {
        p.gwas.mv <- plot(gwas.mv, genObj = panel.gwas, type = "manhattan",
                          sig.col = "darkorange") +
            ggtitle("Multivariate GWAS Manhattan Plot")
        save_plot(p.gwas.mv, "GW01_multivariate_manhattan", height = 4)
    }
}

# =============================================================================
# Final summary
# =============================================================================
plots_made <- list.files(out_dir, pattern = "\\.png$")
cat("\n")
cat("================================================================\n")
cat("  Demo complete —", length(plots_made), "plots saved to", out_dir, "\n")
cat("================================================================\n\n")
cat("  Key takeaways:\n")
cat("  - Trait = NULL (default): identical to existing univariate analysis\n")
cat("  - Trait = 'Site':  engine upgrades genome term to vm():diag(Site)\n")
cat("                     uses multivariate outlier stat + pchisq.mixture LRT\n")
cat("                     post-loop waldTest prunes non-significant interactions\n")
cat("  - n.fa = 1L (2 traits): upgrades diag -> corh (estimates rg)\n")
cat("  - Both qtlAim() and gwasAim() support Trait\n")
cat("  - summary() output gains a leading 'Trait' column\n\n")

for (f in sort(plots_made))
    cat("  ", f, "\n")

# =============================================================================
# NOTES AND POTENTIAL ISSUES
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Notes on running this demo\n")
cat("================================================================\n\n")

cat("SIMULATION DESIGN:\n")
cat("  The polygenic G x E background is sampled as U ~ MN(0, K, Ga)\n")
cat("  using the SAME genomic relationship matrix K that qtlAim() will\n")
cat("  compute internally via .constructCM(). This makes the model\n")
cat("  correctly specified (data and model share the same K structure).\n\n")

cat("REPLICATION:\n")
cat("  Each site has n.reps = 3 complete blocks (RCBD).\n")
cat("  Without replication ASReml cannot separate Vg from Ve within a site\n")
cat("  and the model will be singular or fail to converge.\n")
cat("  The Rep factor is included in fixed effects as Site + Site:Rep,\n")
cat("  capturing block effects that differ between sites.\n\n")

cat("ASReml BASE MODEL:\n")
cat("  fixed    = yield ~ Site + Site:Rep  (site means + block effects)\n")
cat("  random   = ~ diag(Site):Line  -- residual polygenic G x E\n")
cat("    separate genetic variances per site, no between-site correlation\n")
cat("  residual = ~ dsum(~units | Site)  -- separate Ve per site\n")
cat("  Data MUST be sorted by Site before fitting the dsum residual (done above).\n")
cat("  After calling qtlAim(Trait='Site'), the random formula becomes:\n")
cat("    ~ diag(Site):vm(Line, covObj) + diag(Site):Line\n")
cat("  The first term is the composite additive G x E genome-wide term.\n")
cat("  The second is the residual polygenic G x E after QTL are extracted.\n\n")

cat("PREDICT() CLASSIFY TERM:\n")
cat("  In .qtlSelect(), the multivariate predict() call uses:\n")
cat("    classify = 'Line:Site'  (the two factor names)\n")
cat("    only     = 'diag(Site):vm(Line, covObj)'  (the full term label)\n")
cat("  ASReml-R returns predicted values for all Line x Site combinations.\n\n")

cat("MARKERS > LINES (vm path):\n")
cat("  The vm() path is used when intervals > lines.\n")
cat("  This demo uses 100 lines and 116 intervals (4 x 29).\n")
cat("  If you increase lines or reduce markers, the mbf() path will be used\n")
cat("  instead -- both paths give identical QTL selection results.\n\n")

cat("N.FA ARGUMENT:\n")
cat("  n.fa = 0 (default): uses diag(Site) -- separate, uncorrelated variances\n")
cat("  n.fa = 1, ntrait=2: upgrades to corh(Site) -- estimates rg\n")
cat("  n.fa = 1, ntrait>=3: upgrades to fa(Site,1) -- 1 factor\n")
cat("  The n.fa upgrade requires good starting values; if the model does not\n")
cat("  converge, fall back to n.fa=0 (diag) which is usually sufficient.\n\n")

cat("WALDT POST-LOOP PRUNING:\n")
cat("  After the selection loop, waldTest() tests each Trait:X.chr.idx\n")
cat("  interaction jointly. If p > TypeI the interaction is dropped and the\n")
cat("  QTL is retained as main-effect only. This is printed to console.\n")
cat("  Result stored in qtl.mv$QTL$is.interaction (TRUE = interaction QTL).\n\n")
