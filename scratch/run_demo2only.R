# =============================================================================
# run_demo2only.R
# Demo 2 only: two field trials, corh(Trial):id base random structure.
#
#  Demo 2a -- qtlAim()  two trials, Trait = "Trial"
#  Demo 2b -- gwasAim() two trials, Trait = "Trial"
#
# Population / genetic map identical to run_demo12only.R:
#   set.seed(101), 200 DH lines, 5 chr x 80 markers @ 1 cM.
# =============================================================================

devtools::load_all()
library(asreml)
library(qtl)
library(MASS)

set.seed(101)

try_run <- function(label, expr) {
    cat("\n  [", label, "]\n", sep = "")
    tryCatch(expr, error = function(e) {
        cat("  ERROR in", label, ":", conditionMessage(e), "\n"); NULL
    })
}

# =============================================================================
# SECTION 0 — population (must be identical to run_demo12only.R)
# =============================================================================
cat("\n------ SECTION 0: population ------\n")
n.lines <- 200; n.chr <- 5; n.mar <- 80; mar.step <- 1.0
map.sim <- sim.map(len=rep((n.mar-1)*mar.step, n.chr), n.mar=rep(n.mar, n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)

genObj <- primeCross(cross.sim, type="interval", id="id")
cat("  intervals:", sum(sapply(genObj$geno, function(g) ncol(g$interval.data))), "\n")

M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

marker.mat.pm1 <- do.call(cbind, lapply(genObj$geno, function(g) g$imputed.data))
rownames(marker.mat.pm1) <- ids
map.panel <- do.call(rbind, lapply(names(genObj$geno), function(chr)
    data.frame(marker=colnames(genObj$geno[[chr]]$imputed.data),
               chr=chr, pos=as.numeric(genObj$geno[[chr]]$map),
               stringsAsFactors=FALSE)))
panel <- primePanel(marker.mat.pm1, map.panel, id="id", map.id="marker",
                    map.chr="chr", map.pos="pos", encoding="pm1",
                    maf=0.05, impute="none")
panel$pheno$id <- factor(ids)
cat("  panel markers:", sum(sapply(panel$geno, function(g) ncol(g$imputed.data))), "\n")

qtl.info <- list(list(chr="Chr1",int=20L,eff=0.80), list(chr="Chr2",int=55L,eff=-0.65),
                 list(chr="Chr3",int=30L,eff=0.55),  list(chr="Chr4",int=60L,eff=-0.50))
n.qtl <- length(qtl.info)
qtl.X <- matrix(0, nrow=n.lines, ncol=n.qtl,
                dimnames=list(ids, sapply(qtl.info, function(q)
                    paste("Chr",q$chr,q$int,sep="."))))
for (i in seq_len(n.qtl))
    qtl.X[,i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]
cat("  QTL planted: Chr1/20, Chr2/55, Chr3/30, Chr4/60\n\n")

# =============================================================================
# DEMO 2 — two trials, corh
# =============================================================================
cat("================================================================\n")
cat("  DEMO 2: Two field trials — corh(Trial):id\n")
cat("  Trait = 'Trial'  |  200 lines x 2 trials x 2 reps\n")
cat("================================================================\n\n")

n.t2 <- 2L
sigma2.g2 <- c(0.02, 0.02); r.g2 <- 0.75
Ga2 <- matrix(c(sigma2.g2[1],
                r.g2 * sqrt(prod(sigma2.g2)),
                r.g2 * sqrt(prod(sigma2.g2)),
                sigma2.g2[2]), 2, 2)

# QTL effects: QTL 1+2 main, QTL 3+4 interaction / crossover
qtl.eff.2 <- matrix(c(0.80, 0.80, -0.65, -0.65, 0.55, 0.15, -0.50, 0.45),
                    nrow=n.qtl, ncol=n.t2, byrow=TRUE)
cat("Two-trial QTL effects:\n")
cat("                    Trial1   Trial2   Type\n")
types <- c("Main","Main","Interaction","Crossover GxE")
for (i in seq_len(n.qtl))
    cat(sprintf("  QTL %d (%-12s) %6.3f   %6.3f   %s\n",
                i, colnames(qtl.X)[i],
                qtl.eff.2[i,1], qtl.eff.2[i,2], types[i]))
cat("\n")

# Simulate polygenic background and phenotypes
Sigma2 <- kronecker(Ga2, K)
U2 <- matrix(MASS::mvrnorm(1, rep(0, n.lines*n.t2), Sigma2),
             nrow=n.lines, ncol=n.t2,
             dimnames=list(ids, paste0("Trial", 1:n.t2)))

pheno.2 <- do.call(rbind, lapply(seq_len(n.t2), function(j) {
    qc <- rowSums(sweep(qtl.X, 2, qtl.eff.2[,j], "*"))
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id=factor(ids), Trial=factor(paste0("Trial",j)),
                   Rep=factor(paste0("R",r)),
                   yield=c(10,11.5)[j] + qc + U2[,j] +
                         rnorm(n.lines, 0, sqrt(0.70)))))
}))
pheno.2 <- pheno.2[order(pheno.2$Trial, pheno.2$Rep, pheno.2$id),]
cat("Phenotypic data:", nrow(pheno.2), "rows\n")
cat("yield range:", round(range(pheno.2$yield), 2), "\n\n")

# ---- Base model ------------------------------------------------------------
cat("Fitting base model  corh(Trial):id ...\n")
cat("  fixed    = yield ~ Trial + Trial:Rep\n")
cat("  random   = ~ corh(Trial):id\n")
cat("  residual = ~ dsum(~units | Trial)\n\n")

base.2 <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                 residual=~dsum(~units|Trial), data=pheno.2, trace=FALSE,
                 na.action=asreml::na.method(x="include"), maxit=40)
if (!base.2$converge) { cat("  updating...\n"); base.2 <- update(base.2) }
cat("Base model converged:", base.2$converge, "\n")
vp2 <- base.2$vparameters * base.2$sigma2
cat("Variance components:\n")
for (nm in names(vp2)[grep("id", names(vp2))])
    cat(sprintf("  %-40s %.4f\n", nm, vp2[nm]))
cat(sprintf("  %-40s %.4f\n", "sigma2 (Ve)", base.2$sigma2), "\n")

# ---- Demo 2a: qtlAim -------------------------------------------------------
cat("================================================================\n")
cat("  Demo 2a: qtlAim() — two trials, Trait = 'Trial'\n")
cat("  Engine: diag(Trial):vm()  ->  corgh(Trial):vm()\n")
cat("================================================================\n\n")

qtl.2a <- try_run("Demo 2a — qtlAim corh", {
    fit <- qtlAim(base.2, genObj=genObj, merge.by="id", Trait="Trial",
                  TypeI=0.05, exclusion.window=20,
                  trace="scratch/trace_demo2a.txt",
                  na.action=asreml::na.method(x="include"))
    cat("\n  QTL detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        for (i in seq_along(fit$QTL$qtl))
            cat(sprintf("    %d. %-22s  %s\n", i, fit$QTL$qtl[i],
                        ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN")))
        cat("\n  summary.qtlAim():\n")
        print(summary(fit, genObj=genObj))
    }
    fit
})

# ---- Demo 2b: gwasAim ------------------------------------------------------
cat("\n")
cat("================================================================\n")
cat("  Demo 2b: gwasAim() — two trials, Trait = 'Trial'\n")
cat("  Engine: diag(Trial):vm()  ->  corgh(Trial):vm()\n")
cat("================================================================\n\n")

base.2.gwas <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                      residual=~dsum(~units|Trial), data=pheno.2, trace=FALSE,
                      na.action=asreml::na.method(x="include"), maxit=40)
if (!base.2.gwas$converge) { cat("  updating...\n"); base.2.gwas <- update(base.2.gwas) }
cat("gwasAim base model converged:", base.2.gwas$converge, "\n\n")

gwas.2b <- try_run("Demo 2b — gwasAim corh", {
    fit <- gwasAim(base.2.gwas, genObj=panel, merge.by="id", Trait="Trial",
                   TypeI=0.05, exclusion.window=20,
                   trace="scratch/trace_demo2b.txt",
                   na.action=asreml::na.method(x="include"))
    cat("\n  Markers detected:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        for (i in seq_along(fit$QTL$qtl))
            cat(sprintf("    %d. %-22s  %s\n", i, fit$QTL$qtl[i],
                        ifelse(fit$QTL$is.interaction[i], "INTERACTION", "MAIN")))
        cat("\n  summary.gwasAim():\n")
        print(summary(fit, genObj=panel))
    }
    fit
})

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("  Demo 2 — summary\n")
cat("================================================================\n\n")
cat("Planted QTL:\n")
cat("  QTL 1 — Chr1 int 20 — MAIN       (+0.80 both trials)\n")
cat("  QTL 2 — Chr2 int 55 — MAIN       (-0.65 both trials)\n")
cat("  QTL 3 — Chr3 int 30 — INTERACTION (+0.55 Trial1, +0.15 Trial2)\n")
cat("  QTL 4 — Chr4 int 60 — CROSSOVER  (-0.50 Trial1, +0.45 Trial2)\n\n")

report <- function(label, obj) {
    n <- if (!is.null(obj) && length(obj$QTL$qtl)) length(obj$QTL$qtl) else 0L
    cat(sprintf("  %-36s %d QTL\n", label, n))
    if (n > 0)
        for (i in seq_len(n))
            cat(sprintf("      %d. %-22s %s\n", i, obj$QTL$qtl[i],
                        ifelse(obj$QTL$is.interaction[i], "INTERACTION", "MAIN")))
}
report("Demo 2a  qtlAim  corh:", qtl.2a)
report("Demo 2b  gwasAim corh:", gwas.2b)
cat("\nTrace files: scratch/trace_demo2a.txt  scratch/trace_demo2b.txt\n\n")
