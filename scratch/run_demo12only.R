# Run only Section 0 + Demos 1 and 2 to verify they work before Demo 3
# Paste the demo script contents up to (but not including) DEMO 3.

devtools::load_all()
library(asreml)
library(qtl)
library(MASS)
library(ggplot2)

set.seed(101)

# ---- try_run ---------------------------------------------------------------
# Convenience wrapper: catch hard R errors only; let ASReml warnings pass
# through unmuffled. ASReml emits many expected warnings during fitting
# (singularities at late iterations, log-lik decreases) that should not
# abort the run — muffling them was incorrectly aborting qtlAim/gwasAim.
try_run <- function(label, expr) {
    cat("\n  [", label, "]\n", sep = "")
    tryCatch(
        expr,
        error = function(e) {
            cat("  ERROR in", label, ":", conditionMessage(e), "\n")
            NULL
        }
    )
}

# =============================================================================
# SECTION 0
# =============================================================================
n.lines  <- 200; n.chr <- 5; n.mar <- 80; mar.step <- 1.0
cat("0a. sim.cross bc, 200 lines, 5 chr x 80 markers @ 1 cM\n")
cat("    395 intervals > 200 lines => vm() path\n")
map.sim <- sim.map(len=rep((n.mar-1)*mar.step, n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)
cat("    markers:", sum(nmar(cross.sim)), "\n")

cat("0b. primeCross\n")
genObj <- primeCross(cross.sim, type="interval", id="id")
cat("    intervals:", sum(sapply(genObj$geno, function(g) ncol(g$interval.data))), "\n")

cat("0c. K matrix\n")
M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

cat("0d. primePanel (pm1 from imputed.data)\n")
marker.mat.pm1 <- do.call(cbind, lapply(genObj$geno, function(g) g$imputed.data))
rownames(marker.mat.pm1) <- ids
map.panel <- do.call(rbind, lapply(names(genObj$geno), function(chr) {
    data.frame(marker=colnames(genObj$geno[[chr]]$imputed.data),
               chr=chr, pos=as.numeric(genObj$geno[[chr]]$map),
               stringsAsFactors=FALSE)
}))
panel <- primePanel(marker.mat.pm1, map.panel, id="id", map.id="marker",
                    map.chr="chr", map.pos="pos", encoding="pm1",
                    maf=0.05, impute="none")
panel$pheno$id <- factor(ids)
cat("    markers after MAF:", sum(sapply(panel$geno, function(g) ncol(g$imputed.data))), "\n")

cat("0e. QTL setup\n")
# QTL interval indices must be within 1..(n.mar-1)=79.
# Positions chosen to be well-separated and on different chromosomes.
qtl.info <- list(list(chr="Chr1",int=20L,eff=0.80), list(chr="Chr2",int=55L,eff=-0.65),
                 list(chr="Chr3",int=30L,eff=0.55),  list(chr="Chr4",int=60L,eff=-0.50))
n.qtl <- length(qtl.info)
qtl.X <- matrix(0, nrow=n.lines, ncol=n.qtl,
                dimnames=list(ids, sapply(qtl.info, function(q)
                    paste("Chr",q$chr,q$int,sep="."))))
for (i in seq_len(n.qtl))
    qtl.X[,i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[, qtl.info[[i]]$int]
qtl.main <- rowSums(sweep(qtl.X, 2, sapply(qtl.info,"[[","eff"), "*"))

cat("0f. Polygenic effects (single trial)\n")
# Keep polygenic background small so the LRT drops after the 4 true QTL.
# The 4 QTL explain ~1.6 units of variance; with sigma2.g=0.05 the residual
# polygenic signal is negligible and the algorithm stops cleanly.
# sigma2.e=0.80 gives h2 ~ 1.6/(1.6+0.05+0.80) ~ 0.65 (QTL-driven).
sigma2.g.uni <- 0.05
sigma2.e.uni <- 0.80
u.uni <- MASS::mvrnorm(1, rep(0,n.lines), sigma2.g.uni * K)
names(u.uni) <- ids

# =============================================================================
# DEMO 1 — single trial
# =============================================================================
cat("\n---- DEMO 1: single trial ----\n")
pheno.uni <- do.call(rbind, lapply(1:2, function(r)
    data.frame(id=factor(ids), Rep=factor(r),
               yield=10 + c(0.25,-0.25)[r] + qtl.main + u.uni +
                     rnorm(n.lines,0,sqrt(sigma2.e.uni)))))
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep),]

base.uni <- asreml(fixed=yield~Rep, random=~id, data=pheno.uni, trace=FALSE,
                   na.action=asreml::na.method(x="include"))
cat("base.uni converged:", base.uni$converge,
    " Vg:", round(base.uni$vparameters["id"]*base.uni$sigma2,3),
    " Ve:", round(base.uni$sigma2,3), "\n")

qtl.1a <- try_run("Demo 1a qtlAim", {
    fit <- qtlAim(base.uni, genObj=genObj, merge.by="id",
                  TypeI=0.05, exclusion.window=20,
                  trace="scratch/trace_demo1a.txt",
                  na.action=asreml::na.method(x="include"))
    cat("  QTL:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        print(summary(fit, genObj=genObj))
    }
    fit
})

base.gwas.uni <- asreml(fixed=yield~Rep, random=~id, data=pheno.uni, trace=FALSE,
                        na.action=asreml::na.method(x="include"))
gwas.1b <- try_run("Demo 1b gwasAim", {
    fit <- gwasAim(base.gwas.uni, genObj=panel, merge.by="id",
                   TypeI=0.05, exclusion.window=20,
                   trace="scratch/trace_demo1b.txt",
                   na.action=asreml::na.method(x="include"))
    cat("  Markers:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) print(summary(fit, genObj=panel))
    fit
})

# =============================================================================
# DEMO 2 — two trials, corh
# =============================================================================
cat("\n---- DEMO 2: two trials, corh ----\n")
# Polygenic variances large enough for ASReml to estimate the corh structure
# reliably, but small enough that the LRT drops after the true QTL.
# sigma2.g2 = c(0.10, 0.10) with r.g2 = 0.75 gives a plausible MET background.
# QTL effects explain ~1.6 (Trial1) and ~1.3 (Trial2) units; polygenic adds
# 0.10 per trial. sigma2.e = 0.70 per trial.
# h2 (trial1) ~ 1.6 / (1.6+0.10+0.70) ~ 0.67; h2 (trial2) ~ 0.62.
n.t2 <- 2L
sigma2.g2 <- c(0.02, 0.02); r.g2 <- 0.75
Ga2 <- matrix(c(sigma2.g2[1], r.g2*sqrt(prod(sigma2.g2)),
                r.g2*sqrt(prod(sigma2.g2)), sigma2.g2[2]), 2, 2)
qtl.eff.2 <- matrix(c(0.80,0.80, -0.65,-0.65, 0.55,0.15, -0.50,0.45),
                    nrow=n.qtl, ncol=n.t2, byrow=TRUE)
Sigma2 <- kronecker(Ga2, K)
U2 <- matrix(MASS::mvrnorm(1, rep(0,n.lines*n.t2), Sigma2),
             nrow=n.lines, ncol=n.t2, dimnames=list(ids, paste0("Trial",1:n.t2)))
pheno.2 <- do.call(rbind, lapply(seq_len(n.t2), function(j) {
    qc <- rowSums(sweep(qtl.X,2,qtl.eff.2[,j],"*"))
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id=factor(ids), Trial=factor(paste0("Trial",j)),
                   Rep=factor(paste0("R",r)),
                   yield=c(10,11.5)[j]+qc+U2[,j]+rnorm(n.lines,0,sqrt(0.70))))
    )
}))
pheno.2 <- pheno.2[order(pheno.2$Trial, pheno.2$Rep, pheno.2$id),]

base.2 <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                 residual=~dsum(~units|Trial), data=pheno.2, trace=FALSE,
                 na.action=asreml::na.method(x="include"), maxit=40)
if (!base.2$converge) base.2 <- update(base.2)
cat("base.2 converged:", base.2$converge, "\n")

qtl.2a <- try_run("Demo 2a qtlAim corh", {
    fit <- qtlAim(base.2, genObj=genObj, merge.by="id", Trait="Trial",
                  TypeI=0.05, exclusion.window=20,
                  trace="scratch/trace_demo2a.txt",
                  na.action=asreml::na.method(x="include"))
    cat("  QTL:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        for (i in seq_along(fit$QTL$qtl))
            cat("   ", fit$QTL$qtl[i], ":",
                ifelse(fit$QTL$is.interaction[i], "INTERACTION","MAIN"), "\n")
        print(summary(fit, genObj=genObj))
    }
    fit
})

base.2.gwas <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                      residual=~dsum(~units|Trial), data=pheno.2, trace=FALSE,
                      na.action=asreml::na.method(x="include"), maxit=40)
if (!base.2.gwas$converge) base.2.gwas <- update(base.2.gwas)
gwas.2b <- try_run("Demo 2b gwasAim corh", {
    fit <- gwasAim(base.2.gwas, genObj=panel, merge.by="id", Trait="Trial",
                   TypeI=0.05, exclusion.window=20,
                   trace="scratch/trace_demo2b.txt",
                   na.action=asreml::na.method(x="include"))
    cat("  Markers:", length(fit$QTL$qtl), "\n")
    if (length(fit$QTL$qtl)) {
        for (i in seq_along(fit$QTL$qtl))
            cat("   ", fit$QTL$qtl[i], ":",
                ifelse(fit$QTL$is.interaction[i], "INTERACTION","MAIN"), "\n")
        print(summary(fit, genObj=panel))
    }
    fit
})

cat("\n==== Demos 1 and 2 complete ====\n")
cat("  1a qtlAim  single:", if (!is.null(qtl.1a)  && length(qtl.1a$QTL$qtl))  length(qtl.1a$QTL$qtl)  else 0, "\n")
cat("  1b gwasAim single:", if (!is.null(gwas.1b) && length(gwas.1b$QTL$qtl)) length(gwas.1b$QTL$qtl) else 0, "\n")
cat("  2a qtlAim  corh:  ", if (!is.null(qtl.2a)  && length(qtl.2a$QTL$qtl))  length(qtl.2a$QTL$qtl)  else 0, "\n")
cat("  2b gwasAim corh:  ", if (!is.null(gwas.2b) && length(gwas.2b$QTL$qtl)) length(gwas.2b$QTL$qtl) else 0, "\n")
