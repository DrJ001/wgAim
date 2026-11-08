# =============================================================================
# demo_diag1.R  — Demo 1 diagnostics
# Runs single-trial qtlAim + gwasAim and produces:
#   - aimTrace() LRT plots
#   - plot.qtlAim() / plot.gwasAim() outlier plots
# =============================================================================
devtools::load_all()
library(asreml)
library(qtl)
library(MASS)
library(ggplot2)

set.seed(101)

out_dir <- file.path("scratch", "plots_diag")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp <- function(p, name, w = 14, h = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggsave(path, plot = p, width = min(w, 48), height = min(h, 48),
           dpi = 150, bg = "white", limitsize = FALSE)
    message("Saved: ", path)
}

# =============================================================================
# Section 0 — Simulate (same seed as demo_mv_large.R)
# =============================================================================
n.lines <- 200; n.chr <- 5; n.mar <- 200; mar.step <- 0.5
map.sim <- sim.map(len = rep((n.mar - 1) * mar.step, n.chr),
                   n.mar = rep(n.mar, n.chr),
                   eq.spac = TRUE, include.x = FALSE, sex.sp = FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

cross.sim <- sim.cross(map.sim, n.ind = n.lines, type = "bc", error.prob = 0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)

genObj <- primeCross(cross.sim, type = "interval", id = "id")

# Relationship matrix K (same as engine builds internally)
M  <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K  <- tcrossprod(M) / ncol(M)
K  <- K / mean(diag(K))

# primePanel from imputed.data (pm1 coding)
marker.mat.pm1 <- do.call("cbind", lapply(genObj$geno, function(g) g$imputed.data))
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

# Planted QTL effects
qtl.info <- list(list(chr="Chr1", int=50L, eff= 0.80),
                 list(chr="Chr2", int=80L, eff=-0.65),
                 list(chr="Chr3", int=30L, eff= 0.55),
                 list(chr="Chr4", int=120L,eff=-0.50))
n.qtl <- length(qtl.info)
qtl.X <- matrix(0, nrow = n.lines, ncol = n.qtl, dimnames = list(ids, NULL))
for (i in seq_len(n.qtl))
    qtl.X[, i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]
qtl.main <- rowSums(sweep(qtl.X, 2, sapply(qtl.info, "[[", "eff"), "*"))

# Polygenic background kept small so the LRT drops after extracting
# the 4 true QTL.  With sigma2.g=0.60 and sigma2.e=0.40, the infinitesimal
# background carries ~0.60 units of variance across all 995 intervals and
# the composite LRT keeps passing well beyond 4 iterations (false positives).
# Using sigma2.g=0.05, sigma2.e=0.80 gives h² ≈ 0.67 (QTL-driven) and the
# polygenic signal is negligible relative to the 4 planted QTL effects.
sigma2.g.uni <- 0.05
sigma2.e.uni <- 0.80
u.uni <- MASS::mvrnorm(1, rep(0, n.lines), sigma2.g.uni * K)
names(u.uni) <- ids

# Phenotypic data — single trial, 2 reps
pheno.uni <- do.call(rbind, lapply(1:2, function(r)
    data.frame(id    = factor(ids),
               Rep   = factor(r),
               yield = 10 + c(0.25, -0.25)[r] + qtl.main + u.uni +
                       rnorm(n.lines, 0, sqrt(sigma2.e.uni)))))
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep), ]

cat("True QTL planted:\n")
for (i in seq_len(n.qtl))
    cat("  ", qtl.info[[i]]$chr, "interval", qtl.info[[i]]$int,
        " effect =", qtl.info[[i]]$eff, "\n")
cat("\n")

# =============================================================================
# Demo 1a — qtlAim
# =============================================================================
cat("Fitting base model (qtlAim)...\n")
base.uni <- asreml(fixed = yield ~ Rep, random = ~ id, data = pheno.uni,
                   trace = FALSE, na.action = asreml::na.method(x = "include"))

cat("Running qtlAim()...\n")
qtl.1a <- qtlAim(base.uni, genObj = genObj, merge.by = "id",
                 TypeI = 0.05, exclusion.window = 20,
                 trace = "scratch/trace_diag1a.txt",
                 na.action = asreml::na.method(x = "include"))

cat("\n--- qtlAim results ---\n")
cat("QTL detected:", length(qtl.1a$QTL$qtl), "\n")
cat("QTL locations:", paste(qtl.1a$QTL$qtl, collapse = ", "), "\n\n")
print(summary(qtl.1a, genObj = genObj))

# aimTrace LRT plot
cat("\nGenerating aimTrace lrt plot...\n")
p.lrt.qtl <- aimTrace(qtl.1a, plot = "lrt")
sp(p.lrt.qtl, "D1a_qtlAim_lrt",
   w = 8, h = 4 + 0.3 * length(qtl.1a$QTL$diag$lik))

# aimTrace stability plot
cat("Generating aimTrace stability plot...\n")
p.stab.qtl <- aimTrace(qtl.1a, plot = "stability")
sp(p.stab.qtl, "D1a_qtlAim_stability",
   w = 10, h = 3 + length(qtl.1a$QTL$qtl))

# outlier plot
cat("Generating outlier plot...\n")
p.oint.qtl <- plot(qtl.1a, genObj = genObj, type = "outlier", sig.col = "firebrick")
sp(p.oint.qtl, "D1a_qtlAim_outlier",
   w = 14, h = 3 * length(qtl.1a$QTL$diag$oint))

# blups plot
cat("Generating blups plot...\n")
p.blup.qtl <- plot(qtl.1a, genObj = genObj, type = "blups", sig.col = "firebrick")
sp(p.blup.qtl, "D1a_qtlAim_blups",
   w = 14, h = 3 * length(qtl.1a$QTL$diag$blups))

# =============================================================================
# Demo 1b — gwasAim
# =============================================================================
cat("\nFitting base model (gwasAim)...\n")
base.gwas <- asreml(fixed = yield ~ Rep, random = ~ id, data = pheno.uni,
                    trace = FALSE, na.action = asreml::na.method(x = "include"))

cat("Running gwasAim()...\n")
gwas.1b <- gwasAim(base.gwas, genObj = panel, merge.by = "id",
                   TypeI = 0.05, exclusion.window = 20,
                   trace = "scratch/trace_diag1b.txt",
                   na.action = asreml::na.method(x = "include"))

cat("\n--- gwasAim results ---\n")
cat("Markers detected:", length(gwas.1b$QTL$qtl), "\n")
cat("Marker locations:", paste(gwas.1b$QTL$qtl, collapse = ", "), "\n\n")
print(summary(gwas.1b, genObj = panel))

# aimTrace LRT plot
cat("\nGenerating aimTrace lrt plot...\n")
p.lrt.gwas <- aimTrace(gwas.1b, plot = "lrt")
sp(p.lrt.gwas, "D1b_gwasAim_lrt",
   w = 8, h = 4 + 0.3 * length(gwas.1b$QTL$diag$lik))

# aimTrace stability plot
cat("Generating aimTrace stability plot...\n")
p.stab.gwas <- aimTrace(gwas.1b, plot = "stability")
sp(p.stab.gwas, "D1b_gwasAim_stability",
   w = 10, h = 3 + length(gwas.1b$QTL$qtl))

# manhattan plot
cat("Generating manhattan plot...\n")
p.mh <- plot(gwas.1b, genObj = panel, type = "manhattan", sig.col = "firebrick")
sp(p.mh, "D1b_gwasAim_manhattan",
   w = 14, h = 3 * length(gwas.1b$QTL$diag$oint))

# outlier plot
cat("Generating outlier plot...\n")
p.oint.gwas <- plot(gwas.1b, genObj = panel, type = "outlier", sig.col = "firebrick")
sp(p.oint.gwas, "D1b_gwasAim_outlier",
   w = 14, h = 3 * length(gwas.1b$QTL$diag$oint))

# =============================================================================
# LRT statistic summary — print the table aimTrace would show
# =============================================================================
cat("\n=== LRT table — qtlAim ===\n")
aimTrace(qtl.1a, lik.out = TRUE, plot = FALSE)

cat("\n=== LRT table — gwasAim ===\n")
aimTrace(gwas.1b, lik.out = TRUE, plot = FALSE)

plots_made <- list.files(out_dir, pattern = "\\.png$")
cat("\n", length(plots_made), "plots saved to", out_dir, "\n")
for (f in sort(plots_made)) cat("  ", f, "\n")
