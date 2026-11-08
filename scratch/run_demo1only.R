# =============================================================================
# run_demo1only.R  —  Single-trial qtlAim and gwasAim verification
# =============================================================================
devtools::load_all()
library(asreml)
library(qtl)
library(MASS)
library(ggplot2)

set.seed(101)

out_dir <- file.path("scratch", "plots_demo1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sp <- function(p, name, w = 14, h = 6) {
    path <- file.path(out_dir, paste0(name, ".png"))
    ggplot2::ggsave(path, plot = p,
                    width  = min(w, 48), height = min(h, 48),
                    dpi = 150, bg = "white", limitsize = FALSE)
    message("  Saved: ", path)
}

# =============================================================================
# Section 0 — Simulate (same as run_demo12only.R)
# =============================================================================
n.lines <- 200; n.chr <- 5; n.mar <- 80; mar.step <- 1.0
map.sim <- sim.map(len = rep((n.mar - 1) * mar.step, n.chr),
                   n.mar = rep(n.mar, n.chr),
                   eq.spac = TRUE, include.x = FALSE, sex.sp = FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))

cross.sim <- sim.cross(map.sim, n.ind = n.lines, type = "bc", error.prob = 0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)

genObj <- primeCross(cross.sim, type = "interval", id = "id")
M  <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K  <- tcrossprod(M) / ncol(M);  K <- K / mean(diag(K))

# Panel
marker.mat.pm1 <- do.call("cbind", lapply(genObj$geno, function(g) g$imputed.data))
rownames(marker.mat.pm1) <- ids
map.panel <- do.call(rbind, lapply(names(genObj$geno), function(chr)
    data.frame(marker = colnames(genObj$geno[[chr]]$imputed.data),
               chr    = chr,
               pos    = as.numeric(genObj$geno[[chr]]$map),
               stringsAsFactors = FALSE)))
panel <- primePanel(marker.mat.pm1, map.panel,
                    id = "id", map.id = "marker", map.chr = "chr", map.pos = "pos",
                    encoding = "pm1", impute = "none")
panel$pheno$id <- factor(ids)

# QTL
qtl.info <- list(list(chr = "Chr1", int = 20L, eff =  0.80),
                 list(chr = "Chr2", int = 55L, eff = -0.65),
                 list(chr = "Chr3", int = 30L, eff =  0.55),
                 list(chr = "Chr4", int = 60L, eff = -0.50))
n.qtl <- length(qtl.info)
qtl.X <- matrix(0, nrow = n.lines, ncol = n.qtl, dimnames = list(ids, NULL))
for (i in seq_len(n.qtl))
    qtl.X[, i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]
qtl.main <- rowSums(sweep(qtl.X, 2, sapply(qtl.info, "[[", "eff"), "*"))

u.uni <- MASS::mvrnorm(1, rep(0, n.lines), 0.05 * K)
names(u.uni) <- ids

pheno.uni <- do.call(rbind, lapply(1:2, function(r)
    data.frame(id    = factor(ids),
               Rep   = factor(r),
               yield = 10 + c(0.25, -0.25)[r] + qtl.main + u.uni +
                       rnorm(n.lines, 0, sqrt(0.80)))))
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep), ]

cat("True QTL planted:\n")
for (i in seq_len(n.qtl))
    cat("  Chr", qtl.info[[i]]$chr, "int", qtl.info[[i]]$int,
        " effect =", qtl.info[[i]]$eff, "\n")

# =============================================================================
# Demo 1a — qtlAim
# =============================================================================
cat("\n------ Demo 1a: qtlAim() ------\n")
base.uni <- asreml(fixed = yield ~ Rep, random = ~ id, data = pheno.uni,
                   trace = FALSE, na.action = asreml::na.method(x = "include"))
cat("base converged:", base.uni$converge,
    " Vg:", round(base.uni$vparameters["id"] * base.uni$sigma2, 3),
    " Ve:", round(base.uni$sigma2, 3), "\n\n")

qtl.1a <- qtlAim(base.uni, genObj = genObj, merge.by = "id",
                 TypeI = 0.05, exclusion.window = 20,
                 trace = "scratch/trace_demo1a.txt",
                 na.action = asreml::na.method(x = "include"))

cat("\nQTL detected:", length(qtl.1a$QTL$qtl), "\n")
cat("Locations:", paste(qtl.1a$QTL$qtl, collapse = ", "), "\n\n")

cat("summary.qtlAim():\n")
print(summary(qtl.1a, genObj = genObj))

cat("\naimTrace() LRT table:\n")
aimTrace(qtl.1a, lik.out = TRUE, plot = FALSE)

p.lrt  <- aimTrace(qtl.1a, plot = "lrt")
p.stab <- aimTrace(qtl.1a, plot = "stability")
p.oint <- plot(qtl.1a, genObj = genObj, type = "outlier",  sig.col = "firebrick")
p.blup <- plot(qtl.1a, genObj = genObj, type = "blups",    sig.col = "firebrick")
p.heat <- plot(qtl.1a, genObj = genObj, type = "heatmap",  sig.col = "firebrick")
p.eff  <- plot(qtl.1a, genObj = genObj, type = "effects")

sp(p.lrt,  "D1a_lrt",      w = 9,  h = 5)
sp(p.stab, "D1a_stability", w = 12, h = 4 + length(qtl.1a$QTL$qtl))
sp(p.oint, "D1a_outlier",   w = 14, h = 3 * length(qtl.1a$QTL$diag$oint))
sp(p.blup, "D1a_blups",     w = 14, h = 3 * length(qtl.1a$QTL$diag$blups))
sp(p.heat, "D1a_heatmap",   w = 14, h = 3 + length(qtl.1a$QTL$diag$oint) * 0.6)
sp(p.eff,  "D1a_effects",   w = 10, h = 4)

# =============================================================================
# Demo 1b — gwasAim
# =============================================================================
cat("\n------ Demo 1b: gwasAim() ------\n")
base.gwas <- asreml(fixed = yield ~ Rep, random = ~ id, data = pheno.uni,
                    trace = FALSE, na.action = asreml::na.method(x = "include"))

gwas.1b <- gwasAim(base.gwas, genObj = panel, merge.by = "id",
                   TypeI = 0.05, exclusion.window = 20,
                   trace = "scratch/trace_demo1b.txt",
                   na.action = asreml::na.method(x = "include"))

cat("\nMarkers detected:", length(gwas.1b$QTL$qtl), "\n")
cat("Locations:", paste(gwas.1b$QTL$qtl, collapse = ", "), "\n\n")

cat("summary.gwasAim():\n")
print(summary(gwas.1b, genObj = panel))

cat("\naimTrace() LRT table:\n")
aimTrace(gwas.1b, lik.out = TRUE, plot = "both", print.out = FALSE)

p.glrt  <- aimTrace(gwas.1b, plot = "lrt")
p.gstab <- aimTrace(gwas.1b, plot = "stability")
p.mh    <- plot(gwas.1b, genObj = panel, type = "manhattan", sig.col = "firebrick")
p.goint <- plot(gwas.1b, genObj = panel, type = "outlier",   sig.col = "firebrick")
p.gheat <- plot(gwas.1b, genObj = panel, type = "heatmap",   sig.col = "firebrick")
p.geff  <- plot(gwas.1b, genObj = panel, type = "effects")

sp(p.glrt,  "D1b_lrt",       w = 9,  h = 5)
sp(p.gstab, "D1b_stability",  w = 12, h = 4 + length(gwas.1b$QTL$qtl))
sp(p.mh,    "D1b_manhattan",  w = 14, h = 3 * length(gwas.1b$QTL$diag$oint))
sp(p.goint, "D1b_outlier",    w = 14, h = 3 * length(gwas.1b$QTL$diag$oint))
sp(p.gheat, "D1b_heatmap",    w = 14, h = 3 + length(gwas.1b$QTL$diag$oint) * 0.6)
sp(p.geff,  "D1b_effects",    w = 10, h = 4)

# =============================================================================
# Summary
# =============================================================================
cat("\n================================================================\n")
cat("  True QTL:\n")
for (i in seq_len(n.qtl))
    cat("    Chr", qtl.info[[i]]$chr, "int", qtl.info[[i]]$int,
        " effect =", sprintf("%+.2f", qtl.info[[i]]$eff), "\n")

cat("\n  qtlAim  detected:", length(qtl.1a$QTL$qtl),
    "—", paste(qtl.1a$QTL$qtl, collapse = ", "), "\n")
cat("  gwasAim detected:", length(gwas.1b$QTL$qtl),
    "—", paste(gwas.1b$QTL$qtl, collapse = ", "), "\n")

plots_made <- list.files(out_dir, pattern = "\\.png$")
cat("\n ", length(plots_made), "plots saved to", out_dir, "\n")
for (f in sort(plots_made)) cat("   ", f, "\n")
cat("================================================================\n")
