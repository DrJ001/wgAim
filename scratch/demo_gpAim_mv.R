# =============================================================================
# demo_gpAim_mv.R
# Multivariate genomic prediction demo for gpAim().
#
# Three GP analyses demonstrating per-BLUP accuracy (Mrode 2014) and
# generalised H² (Cullis 2006), added to the GEBV data frame in all paths:
#
#   Section 1: Univariate G-BLUP       -- Trial 1 only
#   Section 2: Bivariate G-BLUP        -- Trials 1 & 2, corh(Trial):id
#   Section 3: Eight-trial G-BLUP      -- all 8 trials, fa(Trial,2):id
#   Section 4: Accuracy & H² summary   -- compare across analyses
#   Section 5: GEBV rank comparison    -- univariate vs 8-trial MV
#
# Variance structure on the genomic vm term is determined by the base model:
#   bare id          -> bare vm()         (univariate)
#   corh(Trial):id   -> corgh(Trial):vm() (bivariate, mirror)
#   fa(Trial,2):id   -> fa(Trial,2):vm()  (8-trial, mirror)
#
# Outputs: scratch/plots_demo_gp/  (created automatically)
# =============================================================================

devtools::load_all()
library(asreml)
library(qtl)

outdir <- "scratch/plots_demo_gp"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 8, h = 6) {
    ggplot2::ggsave(file.path(outdir, paste0(name, ".png")),
                    plot = p, width = w, height = h, dpi = 150)
    cat("  saved:", file.path(outdir, paste0(name, ".png")), "\n")
    invisible(p)
}

# =============================================================================
# Section 0 — Simulate 8-trial multi-environment data
#
#   120 DH lines, 5 chromosomes x 80 markers @ 1 cM (79 intervals each)
#   8 trials, 2 reps per trial  ->  1920 observations
#   4 planted QTL
# =============================================================================
cat("Section 0: Simulating 8-trial data...\n")
set.seed(202)

n.lines  <- 120
n.chr    <- 5
n.mar    <- 80
n.trials <- 8
n.reps   <- 2

map.sim <- sim.map(len = rep((n.mar - 1) * 1.0, n.chr),
                   n.mar = n.mar, include.x = FALSE, eq.spacing = TRUE)
names(map.sim) <- paste0("C", seq_len(n.chr))   # "C1" .. "C5"
cross.sim <- sim.cross(map.sim, n.ind = n.lines, type = "bc",
                       error.prob = 0, missing.prob = 0)
ids <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
cross.sim$pheno$id <- ids

genObj <- primeCross(cross.sim, type = "interval", id = "id",
                     impute = "Martinez")
cat(sprintf("  genObj: %d lines, %d chromosomes\n", n.lines, n.chr))

M <- do.call(cbind, lapply(genObj$geno, function(g) g$interval.data))
K <- tcrossprod(M) / mean(diag(tcrossprod(M)))

qtl.pos  <- list(c("C1", 20), c("C2", 55), c("C3", 30), c("C4", 60))
qtl.vecs <- lapply(qtl.pos, function(qp) {
    chr <- qp[[1]]; pos <- as.numeric(qp[[2]])
    im  <- genObj$geno[[chr]]$inferred.map
    as.numeric(genObj$geno[[chr]]$interval.data[, which.min(abs(im - pos))])
})
qtl.X <- matrix(unlist(qtl.vecs), nrow = n.lines,
                dimnames = list(ids, paste0("Q", seq_len(4))))
colnames(qtl.X) <- paste0("Q", seq_len(4))

Ga.true <- rWishart(1, df = n.trials + 2,
                    Sigma = diag(0.005, n.trials))[,, 1]
Ga.true <- Ga.true / mean(diag(Ga.true)) * 0.04

qtl.eff         <- matrix(0, 4, n.trials)
qtl.eff[1, ]    <- 0.25
qtl.eff[2, ]    <- -0.20
qtl.eff[3, 1:4] <- 0.18
qtl.eff[4, ]    <- seq(0.20, -0.20, length.out = n.trials)

ch    <- chol(Ga.true)
Znorm <- matrix(rnorm(n.lines * n.trials), n.lines, n.trials)
U     <- Znorm %*% t(ch)
colnames(U) <- paste0("Trial", seq_len(n.trials))

sigma2.e <- rep(0.70, n.trials)
mu.site  <- runif(n.trials, 9.0, 13.0)

pheno <- do.call(rbind, lapply(seq_len(n.trials), function(t) {
    do.call(rbind, lapply(seq_len(n.reps), function(r) {
        data.frame(
            id    = factor(ids, levels = ids),
            Trial = factor(paste0("Trial", t),
                           levels = paste0("Trial", seq_len(n.trials))),
            Rep   = factor(r),
            yield = mu.site[t] + rnorm(1, 0, 0.1) +
                    U[, t] + as.vector(qtl.X %*% qtl.eff[, t]) +
                    rnorm(n.lines, 0, sqrt(sigma2.e[t])),
            stringsAsFactors = FALSE
        )
    }))
}))
pheno <- pheno[order(pheno$Trial, pheno$id, pheno$Rep), ]
cat(sprintf("  pheno: %d obs (%d lines x %d trials x %d reps)\n",
            nrow(pheno), n.lines, n.trials, n.reps))

# =============================================================================
# Section 1 — Univariate G-BLUP (Trial 1)
# =============================================================================
cat("\n--- Section 1: Univariate G-BLUP (Trial 1) ---\n")

pheno1   <- droplevels(pheno[pheno$Trial == "Trial1", ])
base.uni <- asreml(yield ~ Rep, random = ~ id,
                   residual = ~ units, data = pheno1, trace = FALSE)

gp.uni <- gpAim(base.uni, genObj = genObj, merge.by = "id",
                trace = "scratch/trace_gp_univariate.txt")
print(gp.uni)

cat("\nTop 10 lines by GEBV (Trial 1):\n")
print(head(summary(gp.uni), 10))

cat(sprintf("\nAccuracy range : %.4f -- %.4f  (mean %.4f)\n",
            min(gp.uni$GP$gebv$Accuracy),
            max(gp.uni$GP$gebv$Accuracy),
            mean(gp.uni$GP$gebv$Accuracy)))
cat(sprintf("gen.H2         : %.4f\n", gp.uni$GP$gen.H2))
cat(sprintf("naive h2       : %.4f\n", gp.uni$GP$heritability))

save_plot(plot(gp.uni, type = "caterpillar"),              "gp1_caterpillar")
save_plot(plot(gp.uni, type = "density", prop.select = .1),"gp1_density")
save_plot(plot(gp.uni, type = "heatmap"),                  "gp1_heatmap")
save_plot(plot(gp.uni, type = "blups", genObj = genObj),   "gp1_blups")

# =============================================================================
# Section 2 — Bivariate G-BLUP (Trials 1 & 2)
#
# corh(Trial):id in the base model mirrors to corgh(Trial):vm() on the
# genomic term.  Per-trial Accuracy and gen.H2 in the GEBV data frame.
# =============================================================================
cat("\n--- Section 2: Bivariate G-BLUP (Trials 1 & 2) ---\n")

pheno2  <- droplevels(pheno[pheno$Trial %in% c("Trial1", "Trial2"), ])
base2   <- asreml(yield ~ Trial + Trial:Rep,
                  random   = ~ corh(Trial):id,
                  residual = ~ dsum(~units | Trial),
                  data = pheno2, maxit = 40, trace = FALSE)

gp2 <- gpAim(base2, genObj = genObj, merge.by = "id",
             Trait = "Trial",
             trace = "scratch/trace_gp_bivariate.txt")
print(gp2)

cat("\nTop 8 lines per trial:\n")
summ2 <- summary(gp2)
print(do.call(rbind, lapply(levels(summ2$Trial), function(t)
    head(summ2[summ2$Trial == t, ], 8))))

cat("\nAccuracy per trial:\n")
acc2 <- tapply(gp2$GP$gebv$Accuracy, gp2$GP$gebv$Trial, function(a)
    c(min = round(min(a), 4), mean = round(mean(a), 4), max = round(max(a), 4)))
print(do.call(rbind, acc2))

cat("\nGeneralised H2 per trial:\n")
print(round(gp2$GP$gen.H2, 4))

cat("\nGenetic correlation (Gcor):\n")
print(round(gp2$GP$Gcor, 4))

save_plot(plot(gp2, type = "caterpillar"), "gp2_caterpillar", w = 10, h = 5)
save_plot(plot(gp2, type = "density", prop.select = .1), "gp2_density", w = 10, h = 5)
save_plot(plot(gp2, type = "heatmap"), "gp2_heatmap")

# =============================================================================
# Section 3 — Eight-trial G-BLUP with fa(Trial, 2)
#
# fa(Trial,2):id in the base model mirrors to fa(Trial,2):vm().  With 8
# trials and 2 factors the structure is parsimonious (16 loadings + 8
# specifics vs 36 unstructured parameters).
# =============================================================================
cat("\n--- Section 3: Eight-trial G-BLUP (fa(Trial,2)) ---\n")

base8 <- asreml(yield ~ Trial + Trial:Rep,
                random   = ~ fa(Trial, 2):id,
                residual = ~ dsum(~units | Trial),
                data = pheno, maxit = 80, trace = FALSE)

gp8 <- gpAim(base8, genObj = genObj, merge.by = "id",
             Trait = "Trial",
             trace = "scratch/trace_gp_8trial.txt")
print(gp8)

cat("\nTop 5 lines per trial:\n")
summ8 <- summary(gp8)
print(do.call(rbind, lapply(levels(summ8$Trial), function(t)
    head(summ8[summ8$Trial == t, ], 5))))

cat("\nAccuracy per trial:\n")
acc8 <- tapply(gp8$GP$gebv$Accuracy, gp8$GP$gebv$Trial, function(a)
    c(min = round(min(a), 4), mean = round(mean(a), 4), max = round(max(a), 4)))
print(do.call(rbind, acc8))

cat("\nGeneralised H2 per trial:\n")
print(round(gp8$GP$gen.H2, 4))

save_plot(plot(gp8, type = "caterpillar"), "gp8_caterpillar", w = 14, h = 8)
save_plot(plot(gp8, type = "density", prop.select = .1), "gp8_density", w = 14, h = 8)
save_plot(plot(gp8, type = "heatmap"), "gp8_heatmap")

# =============================================================================
# Section 4 — Accuracy & H² summary across all three analyses
# =============================================================================
cat("\n--- Section 4: Accuracy and H2 summary ---\n")

# Univariate
acc.summary <- data.frame(
    Analysis  = c("Univariate (Trial 1)",
                  rep("Bivariate", 2),
                  paste0("8-trial FA(2): ", paste0("Trial", 1:8))),
    Trial     = c("Trial1",
                  paste0("Trial", 1:2),
                  paste0("Trial", 1:8)),
    mean.Acc  = c(
        mean(gp.uni$GP$gebv$Accuracy),
        tapply(gp2$GP$gebv$Accuracy, gp2$GP$gebv$Trial, mean)[paste0("Trial", 1:2)],
        tapply(gp8$GP$gebv$Accuracy, gp8$GP$gebv$Trial, mean)[paste0("Trial", 1:8)]
    ),
    gen.H2    = c(
        gp.uni$GP$gen.H2,
        gp2$GP$gen.H2[paste0("Trial", 1:2)],
        gp8$GP$gen.H2[paste0("Trial", 1:8)]
    ),
    naive.h2  = c(
        gp.uni$GP$heritability,
        gp2$GP$heritability[paste0("Trial", 1:2)],
        gp8$GP$heritability[paste0("Trial", 1:8)]
    ),
    stringsAsFactors = FALSE
)
acc.summary[, 3:5] <- round(acc.summary[, 3:5], 4)
cat("\n")
print(acc.summary, row.names = FALSE)

# Plot: mean accuracy per trial for 8-trial model
acc8.df <- data.frame(
    Trial    = names(tapply(gp8$GP$gebv$Accuracy, gp8$GP$gebv$Trial, mean)),
    mean.Acc = as.numeric(tapply(gp8$GP$gebv$Accuracy, gp8$GP$gebv$Trial, mean)),
    gen.H2   = gp8$GP$gen.H2[paste0("Trial", 1:8)],
    stringsAsFactors = FALSE
)
acc8.df$Trial <- factor(acc8.df$Trial, levels = paste0("Trial", 1:8))

p.acc <- ggplot2::ggplot(acc8.df,
    ggplot2::aes(x = .data$Trial, group = 1)) +
    ggplot2::geom_line(ggplot2::aes(y = .data$mean.Acc),
                       colour = "steelblue", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = .data$mean.Acc),
                        colour = "steelblue", size = 3) +
    ggplot2::geom_line(ggplot2::aes(y = .data$gen.H2),
                       colour = "firebrick", linewidth = 1, linetype = "dashed") +
    ggplot2::geom_point(ggplot2::aes(y = .data$gen.H2),
                        colour = "firebrick", size = 3, shape = 17) +
    ggplot2::annotate("text", x = 7.5, y = max(acc8.df$mean.Acc) * 1.01,
                      label = "mean Accuracy", colour = "steelblue",
                      hjust = 1, size = 3.5) +
    ggplot2::annotate("text", x = 7.5, y = max(acc8.df$gen.H2) * 0.97,
                      label = "gen. H\u00b2", colour = "firebrick",
                      hjust = 1, size = 3.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Value",
                  title = "Mean BLUP Accuracy and Generalised H\u00b2 per Trial",
                  subtitle = "8-trial FA(2) G-BLUP model") +
    theme_scatter()
save_plot(p.acc, "acc_H2_per_trial", w = 9, h = 5)

# =============================================================================
# Section 5 — GEBV rank comparison: univariate vs 8-trial MV (Trial 1)
# =============================================================================
cat("\n--- Section 5: GEBV comparison (univariate vs 8-trial, Trial 1) ---\n")

gebv.uni <- gp.uni$GP$gebv[, c("id", "GEBV", "Accuracy")]
names(gebv.uni)[2:3] <- c("GEBV_uni", "Acc_uni")

gebv.mv8 <- gp8$GP$gebv[gp8$GP$gebv$Trial == "Trial1", c("id", "GEBV", "Accuracy")]
names(gebv.mv8)[2:3] <- c("GEBV_mv8", "Acc_mv8")

comp <- merge(gebv.uni, gebv.mv8, by = "id")
r.pearson  <- round(cor(comp$GEBV_uni, comp$GEBV_mv8), 4)
r.spearman <- round(cor(comp$GEBV_uni, comp$GEBV_mv8, method = "spearman"), 4)

cat(sprintf("  Pearson  r (GEBV uni vs MV-8, Trial 1): %.4f\n", r.pearson))
cat(sprintf("  Spearman r (GEBV uni vs MV-8, Trial 1): %.4f\n", r.spearman))
cat(sprintf("  Mean accuracy uni : %.4f\n", mean(comp$Acc_uni)))
cat(sprintf("  Mean accuracy MV-8: %.4f\n", mean(comp$Acc_mv8)))

p.comp <- ggplot2::ggplot(comp,
    ggplot2::aes(x = .data$GEBV_uni, y = .data$GEBV_mv8,
                 colour = .data$Acc_mv8)) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_colour_gradient(low = "grey70", high = "steelblue",
                                   name = "MV-8\nAccuracy") +
    ggplot2::geom_abline(intercept = 0, slope = 1,
                         colour = "firebrick", linetype = "dashed",
                         linewidth = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = TRUE,
                         colour = "grey30", fill = "grey85",
                         linewidth = 0.5) +
    ggplot2::labs(
        x       = "Univariate GEBV (Trial 1)",
        y       = "8-trial FA(2) GEBV (Trial 1)",
        title   = "G-BLUP comparison: univariate vs 8-trial MV",
        caption = sprintf("Pearson r = %.4f   Spearman r = %.4f",
                          r.pearson, r.spearman)) +
    ggplot2::theme_bw()
save_plot(p.comp, "comparison_uni_vs_mv8", w = 7, h = 6)

cat("\n=== Demo complete. Plots saved to:", outdir, "===\n")
