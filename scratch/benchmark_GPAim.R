# =============================================================================
# benchmark_GPAim.R
#
# Tests GPAim.asreml() on simulated polygenic data.
# A polygenic trait (many small effects) is simulated so true GEBVs are known.
# Checks:
#   1. GPAim runs end-to-end
#   2. GEBVs are extracted correctly
#   3. Heritability estimate is sensible
#   4. Rank correlation between true and predicted GEBVs is reasonable
#   5. print(), summary(), plot() work
#
# Run from the wgaim project root:
#   Rscript scratch/benchmark_GPAim.R
# =============================================================================

library(asreml)
devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim", quiet = TRUE)

set.seed(101)

cat("\n============================================================\n")
cat(" wgAim GPAim Test: Simulated Polygenic Data\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. Simulate polygenic GWAS data
#    200 lines, 300 markers across 5 chromosomes
#    Polygenic trait: 30 markers each contribute a small effect
#    True h2 ~ 0.5
# -----------------------------------------------------------------------------
cat("Simulating polygenic data...\n")
n.lines   <- 200
n.chr     <- 5
n.per.chr <- 60
n.markers <- n.chr * n.per.chr

# Random 0/1/2 genotypes
geno.mat <- matrix(rbinom(n.lines * n.markers, 2, 0.5),
                   nrow = n.lines, ncol = n.markers)
line.ids     <- paste0("Line_", sprintf("%03d", 1:n.lines))
marker.ids   <- paste0("M", sprintf("%04d", 1:n.markers))
rownames(geno.mat) <- line.ids
colnames(geno.mat) <- marker.ids

# Genetic map
map.df <- data.frame(
    marker = marker.ids,
    chr    = rep(paste0("Chr", 1:n.chr), each = n.per.chr),
    pos    = rep(seq(1, n.per.chr, by = 1), times = n.chr),
    stringsAsFactors = FALSE
)

# Polygenic effects: 30 markers with small effects
n.causal    <- 30
causal.idx  <- sample(n.markers, n.causal)
effects     <- rnorm(n.causal, 0, 0.15)

# ±1 coding for true genetic values
geno.pm1    <- geno.mat - 1
true.gebv   <- as.numeric(geno.pm1[, causal.idx] %*% effects)
var.g.true  <- var(true.gebv)
var.e       <- var.g.true          # set Ve = Vg so h2 ~ 0.5
residual    <- rnorm(n.lines, 0, sqrt(var.e))

# 2 reps so base model is estimable
pheno.df <- do.call(rbind, lapply(1:2, function(r)
    data.frame(line  = factor(line.ids),
               rep   = factor(r),
               yield = 10 + true.gebv + rnorm(n.lines, 0, sqrt(var.e)),
               stringsAsFactors = FALSE)))

cat(sprintf("  Lines     : %d (x2 reps)\n", n.lines))
cat(sprintf("  Markers   : %d across %d chromosomes\n", n.markers, n.chr))
cat(sprintf("  Causal    : %d markers\n", n.causal))
cat(sprintf("  True Vg   : %.4f\n", var.g.true))
cat(sprintf("  True Ve   : %.4f\n", var.e))
cat(sprintf("  True h2   : %.4f\n\n", var.g.true / (var.g.true + var.e)))

# -----------------------------------------------------------------------------
# 2. Build panel object
# -----------------------------------------------------------------------------
cat("Building panel object...\n")
panel <- makePanel(geno = geno.mat, map = map.df,
                   encoding = "012", maf = 0.05, impute = "none")
panel$pheno$line <- factor(line.ids)

# -----------------------------------------------------------------------------
# 3. Fit base model
# -----------------------------------------------------------------------------
cat("Fitting base model...\n")
base.mod <- asreml(yield ~ rep,
                   random   = ~ line,
                   data     = pheno.df)
cat(sprintf("  Base logLik: %.3f\n\n", base.mod$loglik))

# -----------------------------------------------------------------------------
# 4. Run GPAim
# -----------------------------------------------------------------------------
cat("Running GPAim()...\n")
t.gp <- system.time({
    gp.fit <- GPAim(base.mod,
                    genoObj  = panel,
                    merge.by = "line",
                    gen.type = "marker",
                    trace    = "scratch/trace_GPAim.txt",
                    na.action = na.method(x = "include"))
})
cat(sprintf("\n  Elapsed: %.2f sec\n", t.gp["elapsed"]))

# -----------------------------------------------------------------------------
# 5. Results
# -----------------------------------------------------------------------------
cat("\n============================================================\n")
cat(" Results\n")
cat("============================================================\n")

cat("\nClass:", paste(class(gp.fit), collapse = ", "), "\n")

# print method
cat("\n--- print(gp.fit) ---\n")
print(gp.fit)

# Rank correlation with true GEBVs
pred.gebv  <- gp.fit$GP$gebv
true.df    <- data.frame(line = line.ids, true.GEBV = true.gebv)
combined   <- merge(pred.gebv, true.df, by.x = "line", by.y = "line")
rank.cor   <- cor(combined$GEBV, combined$true.GEBV, method = "spearman")
pearson.cor <- cor(combined$GEBV, combined$true.GEBV)

cat(sprintf("\nPrediction accuracy:\n"))
cat(sprintf("  Pearson correlation  (predicted vs true GEBV): %.4f\n", pearson.cor))
cat(sprintf("  Spearman correlation (predicted vs true GEBV): %.4f\n", rank.cor))

# Heritability check
cat(sprintf("\nHeritability:\n"))
cat(sprintf("  True h2     : %.4f\n", var.g.true / (var.g.true + var.e)))
cat(sprintf("  Estimated h2: %.4f\n", gp.fit$GP$heritability))

# summary method
cat("\n--- summary(gp.fit) [top 10 lines] ---\n")
st <- summary(gp.fit)
print(head(st, 10))

# All three plot types
for (ptype in c("caterpillar", "heatmap", "density")) {
    fname <- sprintf("scratch/plot_GPAim_%s.png", ptype)
    cat(sprintf("\n--- %s plot -> %s ---\n", ptype, fname))
    tryCatch({
        p <- plot(gp.fit, type = ptype)
        ggplot2::ggsave(fname, p, width = 8, height = 5, dpi = 150)
        cat("  Saved.\n")
    }, error = function(e) cat("  Error:", conditionMessage(e), "\n"))
}

cat("\n============================================================\n")
cat(" Test 2: mbf path (lines > markers)\n")
cat("============================================================\n\n")

# Simulate with more lines than markers to trigger mbf path
set.seed(202)
n.lines2  <- 200
n.markers2 <- 80    # lines > markers -> mbf path

geno2  <- matrix(rbinom(n.lines2 * n.markers2, 2, 0.5),
                 nrow = n.lines2, ncol = n.markers2)
ids2   <- paste0("Line_", sprintf("%03d", 1:n.lines2))
mids2  <- paste0("M", sprintf("%04d", 1:n.markers2))
rownames(geno2) <- ids2; colnames(geno2) <- mids2

map2 <- data.frame(
    marker = mids2,
    chr    = rep(paste0("Chr", 1:4), each = n.markers2 / 4),
    pos    = rep(seq(1, n.markers2 / 4), times = 4),
    stringsAsFactors = FALSE
)
geno2.pm1  <- geno2 - 1
causal2    <- sample(n.markers2, 10)
eff2       <- rnorm(10, 0, 0.2)
true.g2    <- as.numeric(geno2.pm1[, causal2] %*% eff2)
pheno2     <- do.call(rbind, lapply(1:2, function(r)
    data.frame(line  = factor(ids2), rep = factor(r),
               yield = 10 + true.g2 + rnorm(n.lines2, 0, sd(true.g2)),
               stringsAsFactors = FALSE)))

panel2          <- makePanel(geno2, map2, encoding = "012", maf = 0.05, impute = "none")
panel2$pheno$line <- factor(ids2)

base2 <- asreml(yield ~ rep, random = ~ line, data = pheno2)

cat(sprintf("  Lines: %d  Markers: %d  -> expects mbf path\n\n", n.lines2, n.markers2))
gp2 <- GPAim(base2, genoObj = panel2, merge.by = "line",
             trace = "scratch/trace_GPAim2.txt",
             na.action = na.method(x = "include"))

cat("\n--- print(gp2) ---\n")
print(gp2)

pred.pearson2 <- cor(gp2$GP$gebv$GEBV,
                     true.g2[match(gp2$GP$gebv$line, ids2)])
cat(sprintf("\nPearson accuracy (mbf path): %.4f\n", pred.pearson2))

cat("\n============================================================\n")
cat(" Done.\n")
cat("============================================================\n\n")
