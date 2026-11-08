# =============================================================================
# scratch/smoke_gwasAim_mv.R
#
# End-to-end smoke test of the multi-environment GWAS pipeline on all 6 sites
# of phenoGWASsim.  Walks the base model up from diag -> fa1 -> fa2, then runs
# gwasAim() with Trait = "Site".
#
# Aim: verify base FA(2) converges, gwasAim runs to completion, and the 5
# planted QTL are recovered with correct MAIN vs INT classification.
#
# Data (from data-raw/simulate_GWASsim.R):
#   300 lines x 10 chr x 150 markers = 1500 markers, 6 sites x 2 reps = 3600 obs
#
# Planted QTL (effects on the biomAid yield scale, site_mean ~4500):
#   QTL 1  Chr2   45 cM   MAIN         +170 across all 6 sites
#   QTL 2  Chr5   90 cM   MAIN         -180 across all 6 sites
#   QTL 3  Chr7   30 cM   MAIN (weak)  +130 across all 6 sites
#   QTL 4  Chr3  110 cM   GxE fade     +250 (Env1) linearly fading to 0 (Env6)
#   QTL 5  Chr9   60 cM   GxE cross    -180 (Env1-3), +180 (Env4-6)
#
# Timing note: FA(2) + vm() genome-wide model over 3600 obs.  Expect each
# ASReml fit to take ~30-90 s; full script ~5-10 minutes on a reasonable
# machine.
# =============================================================================

suppressPackageStartupMessages({
    devtools::load_all()
    library(asreml)
})
asreml.options(trace       = TRUE,
                workspace   = "4gb",
                pworkspace  = "4gb")

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

# ---------------------------------------------------------------------------
# 1.  Panel prep (same as UV vignette)
# ---------------------------------------------------------------------------
cat("\n--- 1.  Panel prep ---\n")
chk   <- checkPanel(genoGWASraw, mapGWASsim, id = "id")
flt   <- filterPanel(chk)
panel <- primePanel(flt, impute = "knn", knn.k = 5)
cat(sprintf("Prepared panel: %d lines x %d markers\n",
            nrow(panel$pheno),
            sum(sapply(panel$geno, function(el) ncol(el$imputed.data)))))

# ---------------------------------------------------------------------------
# 2.  Phenotype prep: all sites, factors set up
# ---------------------------------------------------------------------------
pheno <- phenoGWASsim
pheno$id   <- factor(pheno$id)
# biomAid pads site names to minimum width -- for 6 sites -> Env1..Env6
n.sites    <- nlevels(pheno$Site)
site.width <- nchar(as.character(n.sites))
env.names  <- sprintf(paste0("Env%0", site.width, "d"), seq_len(n.sites))
pheno$Site <- factor(pheno$Site, levels = env.names)
pheno$Rep  <- factor(pheno$Rep)
# Row and Column left as integers -- we're not modelling spatial in this test
cat(sprintf("Phenotype  : %d rows across %d sites (%d lines x %d reps per site)\n",
            nrow(pheno), nlevels(pheno$Site),
            nlevels(pheno$id), nlevels(pheno$Rep)))

# ---------------------------------------------------------------------------
# 3.  Base model progression: diag -> fa(Site,1) -> fa(Site,2)
# ---------------------------------------------------------------------------
timing <- list()

fit_and_time <- function(label, expr) {
    t0 <- Sys.time()
    res <- eval(expr, envir = parent.frame())
    t1 <- Sys.time()
    timing[[label]] <<- as.numeric(difftime(t1, t0, units = "secs"))
    cat(sprintf("  %s: converged=%s  logL=%.3f  (%.1f s)\n",
                label, res$converge, res$loglik, timing[[label]]))
    res
}

pheno$yield <- pheno$yield/1000

cat("\n--- 2.  Base model: diag(Site):id ---\n")
base_diag <- fit_and_time("diag", quote({
    m <- asreml(
        fixed     = yield ~ Site + Site:Rep,
        random    = ~ diag(Site):id,
        residual  = ~ dsum(~units | Site),
        data      = pheno,
        na.action = na.method(x = "include"),
        maxit     = 100
    )
    m <- update(m)
    m
}))

cat("\n--- 3.  Base model: fa(Site,1):id ---\n")
base_fa1 <- fit_and_time("fa1", quote({
    m <- asreml(
        fixed     = yield ~ Site + Site:Rep,
        random    = ~ fa(Site, 1):id,
        residual  = ~ dsum(~units | Site),
        data      = pheno,
        na.action = na.method(x = "include"),
        maxit     = 100
    )
    m <- update(m)
    m
}))

cat("\n--- 4.  Base model: fa(Site,2):id ---\n")
base_fa2 <- fit_and_time("fa2", quote({
    m <- asreml(
        fixed     = yield ~ Site + Site:Rep,
        random    = ~ fa(Site, 2):id,
        residual  = ~ dsum(~units | Site),
        data      = pheno,
        na.action = na.method(x = "include"),
        maxit     = 100
    )
    m <- update(m)
    m
}))

# Model comparison table
cat("\n--- 5.  Model comparison ---\n")
# Parameter counts for reference (n.sites = number of environments):
#   diag(Site):id             ->  n.sites          (variance per site)
#   fa(Site, 1):id            ->  n.sites psi + n.sites loadings   = 2*n.sites
#   fa(Site, 2):id            ->  n.sites psi + 2*n.sites loadings = 3*n.sites
#   dsum(~units | Site) resid ->  n.sites          (variance per site, all models)
np <- c(diag = as.integer(   n.sites + n.sites),
        fa1  = as.integer(2L*n.sites + n.sites),
        fa2  = as.integer(3L*n.sites + n.sites))
comp <- data.frame(
    Model      = c("diag(Site):id", "fa(Site, 1):id", "fa(Site, 2):id"),
    Parameters = np,
    LogLik     = round(c(base_diag$loglik, base_fa1$loglik, base_fa2$loglik), 2),
    AIC        = round(c(-2*base_diag$loglik + 2*np["diag"],
                          -2*base_fa1$loglik  + 2*np["fa1"],
                          -2*base_fa2$loglik  + 2*np["fa2"]),   2),
    stringsAsFactors = FALSE
)
print(comp, row.names = FALSE)

# ---------------------------------------------------------------------------
# 4.  gwasAim MV run
# ---------------------------------------------------------------------------
cat("\n--- 6.  gwasAim MV (Trait = 'Site', str = NULL mirrors base) ---\n")
t0 <- Sys.time()
gwas_mv <- gwasAim(
    baseModel        = base_fa2,
    genObj           = panel,
    merge.by         = "id",
    Trait            = "Site",
    str              = NULL,            # mirror the FA(2) structure
    exclusion.window = 20,
    TypeI            = 0.05,
    trace            = TRUE,
    na.action        = na.method(x = "include"),
    maxit            = 100
)
t1 <- Sys.time()
timing[["gwasAim"]] <- as.numeric(difftime(t1, t0, units = "secs"))
cat(sprintf("\ngwasAim completed in %.1f s\n", timing[["gwasAim"]]))

# ---------------------------------------------------------------------------
# 5.  Results
# ---------------------------------------------------------------------------
cat("\n--- 7.  gwasAim results ---\n")
cat("class            :", paste(class(gwas_mv), collapse = ", "), "\n")
cat("engine path      :", gwas_mv$QTL$path, "\n")
cat("Trait            :", gwas_mv$QTL$Trait, "\n")
cat("n markers panel  :", gwas_mv$QTL$n.markers, "\n")
cat("QTL detected     :", length(gwas_mv$QTL$qtl), "\n")
cat("keys             :\n"); print(gwas_mv$QTL$qtl)
cat("is.interaction   :\n"); print(gwas_mv$QTL$is.interaction)

cat("\n--- 8.  print / summary ---\n")
print(gwas_mv, genObj = panel)

cat("\nSummary table:\n")
summ_mv <- summary(gwas_mv, genObj = panel)
print(summ_mv)

# ---------------------------------------------------------------------------
# 6.  Ground-truth recovery check
# ---------------------------------------------------------------------------
# The simulator picks planted-QTL markers dynamically (nearest by cM position
# on the target chromosome with MAF > 0.15), so marker names cannot be
# hardcoded.  We look up the nearest marker in the map for each planted
# position, and count a QTL as "recovered" if a detected marker lies within
# `window_cM` on the correct chromosome.
cat("\n--- 9.  Planted-QTL recovery check ---\n")

window_cM <- 5   # cM tolerance around planted position

# Planted-QTL truth table (chr / pos / effect pattern)
truth <- data.frame(
    chr    = c("Chr2", "Chr5", "Chr7", "Chr3", "Chr9"),
    pos    = c(   45,     90,     30,    110,     60),
    type   = c("MAIN  (+170 all sites)",
               "MAIN  (-180 all sites)",
               "MAIN  (+130 weak)",
               "GxE   (+250 fade to 0)",
               "GxE   (-/+180 crossover)"),
    expect = c("MAIN", "MAIN", "MAIN", "INT", "INT"),
    stringsAsFactors = FALSE
)

# Nearest map marker on the target chromosome (best guess at planted marker)
truth$nearest_marker <- mapply(function(ch, ps) {
    sub <- mapGWASsim[mapGWASsim$chr == ch, , drop = FALSE]
    sub$marker[which.min(abs(sub$pos - ps))]
}, truth$chr, truth$pos, USE.NAMES = FALSE)

# Assemble detected-QTL table: marker, chr, cM, MAIN/INT tag
if (length(gwas_mv$QTL$qtl) > 0L) {
    keys   <- gwas_mv$QTL$qtl
    is_int <- gwas_mv$QTL$is.interaction
    det_chr <- vapply(keys, function(k) strsplit(k, "\\.")[[1L]][2L],
                      character(1L))
    det_idx <- vapply(keys, function(k) as.integer(strsplit(k, "\\.")[[1L]][3L]),
                      integer(1L))
    det_mk  <- mapply(function(ch, ix) names(panel$geno[[ch]]$map)[ix],
                       det_chr, det_idx, USE.NAMES = FALSE)
    det_pos <- mapply(function(ch, ix) unname(panel$geno[[ch]]$map[ix]),
                       det_chr, det_idx, USE.NAMES = FALSE)
    det <- data.frame(
        marker = det_mk,
        chr    = det_chr,
        pos    = round(det_pos, 2),
        tag    = ifelse(is_int, "INT", "MAIN"),
        stringsAsFactors = FALSE
    )
} else {
    det <- data.frame(marker = character(0), chr = character(0),
                      pos = numeric(0), tag = character(0),
                      stringsAsFactors = FALSE)
}

cat("\nDetected QTL:\n"); print(det, row.names = FALSE)

# For each planted QTL, find the closest detected marker on the same
# chromosome; recovered if within window_cM.
truth$hit_marker <- NA_character_
truth$hit_pos    <- NA_real_
truth$hit_tag    <- NA_character_
truth$dist_cM    <- NA_real_
for (i in seq_len(nrow(truth))) {
    same_chr <- det[det$chr == truth$chr[i], , drop = FALSE]
    if (nrow(same_chr) == 0L) next
    d  <- abs(same_chr$pos - truth$pos[i])
    ix <- which.min(d)
    if (d[ix] <= window_cM) {
        truth$hit_marker[i] <- same_chr$marker[ix]
        truth$hit_pos[i]    <- round(same_chr$pos[ix], 2)
        truth$hit_tag[i]    <- same_chr$tag[ix]
        truth$dist_cM[i]    <- round(d[ix], 2)
    }
}
truth$recovered <- !is.na(truth$hit_marker)
truth$tag_ok    <- !is.na(truth$hit_tag) & truth$hit_tag == truth$expect

cat(sprintf("\nRecovery (window = %.1f cM):\n", window_cM))
print(truth[, c("chr","pos","type","expect","nearest_marker",
                "hit_marker","hit_pos","dist_cM","hit_tag",
                "recovered","tag_ok")],
      row.names = FALSE)

cat(sprintf("\nSummary: %d / %d planted QTL recovered; %d / %d with correct MAIN/INT tag.\n",
            sum(truth$recovered), nrow(truth),
            sum(truth$tag_ok),    nrow(truth)))

# ---------------------------------------------------------------------------
# 7.  Timing report
# ---------------------------------------------------------------------------
cat("\n--- 10.  Timings (seconds) ---\n")
print(round(unlist(timing), 1))
cat(sprintf("Total: %.1f s (%.1f min)\n",
            sum(unlist(timing)), sum(unlist(timing))/60))
