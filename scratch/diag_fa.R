# diag_fa.R — investigate fa(Trial,2) vparameter structure and diag-downgrade
devtools::load_all()
library(asreml); library(qtl); library(MASS)
set.seed(101)

# Replicate Section 0 of run_demo3only.R (150 lines, same seed)
n.lines <- 150; n.chr <- 5; n.mar <- 80; mar.step <- 1.0
map.sim <- sim.map(len=rep((n.mar-1)*mar.step,n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)
genObj <- primeCross(cross.sim, type="interval", id="id")
M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

# Minimal 6-trial phenotype (noise only — no planted QTL, just structure)
n.t <- 6L
set.seed(42)
pheno.6 <- do.call(rbind, lapply(seq_len(n.t), function(j)
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id    = factor(ids),
                   Trial = factor(paste0("Trial", j), levels=paste0("Trial", 1:n.t)),
                   Rep   = factor(paste0("R", r)),
                   yield = 10 + rnorm(n.lines))))))
pheno.6 <- pheno.6[order(pheno.6$Trial, pheno.6$Rep, pheno.6$id),]
phenoData <- pheno.6  # needed for ASReml formula resolution

# -----------------------------------------------------------
# 1. Fit base model and inspect vparameter names
# -----------------------------------------------------------
cat("=== 1. Base model fa(Trial,2):id ===\n\n")
base.6 <- asreml(fixed=yield~Trial+Trial:Rep, random=~fa(Trial,2):id,
                 residual=~dsum(~units|Trial), data=pheno.6, trace=FALSE,
                 na.action=asreml::na.method(x="include"), maxit=60)
if (!base.6$converge) base.6 <- update(base.6)
cat("Converged:", base.6$converge, "  sigma2:", round(base.6$sigma2, 6), "\n\n")

vp.base <- base.6$vparameters
cat("All vparameter names (base model):\n")
for (nm in names(vp.base))
    cat(sprintf("  %-55s = %.6f  con=%s\n", nm, vp.base[nm],
                base.6$vparameters.con[nm]))

# -----------------------------------------------------------
# 2. Build genome model: diag(Trial):vm then fa(Trial,2):vm
# -----------------------------------------------------------
cat("\n=== 2. Genome model: diag -> fa(2) step-up ===\n\n")
cov.env <- wgAim:::.constructCM(M)
covObj  <- cov.env$relm
assign("covObj",    covObj,    envir=globalenv())
assign("phenoData", phenoData, envir=globalenv())

Trait <- "Trial"
resid.term <- "fa(Trial, 2):id"

# Step 1: diag
ran.diag <- as.formula(paste0("~ diag(Trial):vm(id, covObj) + ", resid.term))
qtlM.diag <- update(base.6, random.=ran.diag, data=quote(phenoData),
                    na.action=asreml::na.method(x="include"), trace=FALSE)
cat("diag vm model converged:", qtlM.diag$converge, "\n\n")

# Step 2: fa(Trial,1):vm
ran.fa1 <- as.formula(paste0("~ fa(Trial,1):vm(id, covObj) + ", resid.term))
qtlM.fa1 <- update(qtlM.diag, random.=ran.fa1,
                   na.action=asreml::na.method(x="include"), trace=FALSE)
cat("fa(1) vm model converged:", qtlM.fa1$converge, "\n\n")

# Step 3: fa(Trial,2):vm
ran.fa2 <- as.formula(paste0("~ fa(Trial,2):vm(id, covObj) + ", resid.term))
qtlM.fa2 <- update(qtlM.fa1, random.=ran.fa2,
                   na.action=asreml::na.method(x="include"), trace=FALSE)
cat("fa(2) vm model converged:", qtlM.fa2$converge, "\n")
cat("sigma2:", round(qtlM.fa2$sigma2, 6), "\n\n")

vp.fa2 <- qtlM.fa2$vparameters
cat("All vparameter names (fa(2) vm model):\n")
for (nm in names(vp.fa2))
    cat(sprintf("  %-65s = %.6f  con=%s\n", nm, vp.fa2[nm],
                qtlM.fa2$vparameters.con[nm]))

# -----------------------------------------------------------
# 3. Dissect the vm-term vparameters — sterms grep
# -----------------------------------------------------------
cat("\n=== 3. sterms grep on fa(2) vm model ===\n\n")
sterms <- "vm.*covObj|mbf.*ints"
vm.nm <- names(vp.fa2)[grep(sterms, names(vp.fa2))]
cat("vm-term parameters (", length(vm.nm), "total):\n", sep="")
for (nm in vm.nm)
    cat(sprintf("  %-65s = %.6f  con=%s\n", nm, vp.fa2[nm],
                qtlM.fa2$vparameters.con[nm]))

# Check: which are !var (Psi), which are !fa1, !fa2
is.var <- grepl("!var$",  vm.nm)
is.fa1 <- grepl("!fa1$",  vm.nm)
is.fa2 <- grepl("!fa2$",  vm.nm)
is.other <- !(is.var | is.fa1 | is.fa2)
cat("\n  !var  (Psi, specific variances):", sum(is.var),  "params\n")
cat("  !fa1  (Lambda col 1)           :", sum(is.fa1), "params\n")
cat("  !fa2  (Lambda col 2)           :", sum(is.fa2), "params\n")
cat("  other                          :", sum(is.other), "params\n")
if (any(is.other)) {
    cat("  Other param names:\n")
    for (nm in vm.nm[is.other]) cat("   ", nm, "\n")
}

# -----------------------------------------------------------
# 4. Ga construction — current code vs corrected
# -----------------------------------------------------------
cat("\n=== 4. Ga construction check ===\n\n")
ntrait <- n.t
sigma2 <- qtlM.fa2$sigma2
vp <- vp.fa2

# Current code (multiplies ALL vm params by sigma2)
vpars.all <- vp[grep(sterms, names(vp))] * sigma2
n.fa.val  <- 2L
psi.cur   <- vpars.all[1:ntrait]
lam.cur   <- matrix(vpars.all[(ntrait+1L):((n.fa.val+1L)*ntrait)], ncol=n.fa.val)
Ga.cur    <- lam.cur %*% t(lam.cur) + diag(psi.cur)
cat("Current Ga (all params * sigma2):\n")
cat("  diagonal:", round(diag(Ga.cur), 4), "\n")
cat("  off-diag range:", round(range(Ga.cur[lower.tri(Ga.cur)]), 4), "\n")
cat("  is PD:", all(eigen(Ga.cur)$values > 0), "\n\n")

# Corrected code (only !var params * sigma2; !fa loadings in raw scale)
vpars.raw <- vp[grep(sterms, names(vp))]
vp.nms    <- names(vp)[grep(sterms, names(vp))]
is.loading <- grepl("!fa[0-9]+$", vp.nms)
psi.cor   <- pmax(vpars.raw[!is.loading] * sigma2, 0)
lam.raw   <- vpars.raw[is.loading]
cat("Corrected: n psi params =", sum(!is.loading), "  n loading params =", sum(is.loading), "\n")
lam.cor   <- matrix(lam.raw, ncol=n.fa.val)
Ga.cor    <- lam.cor %*% t(lam.cor) + diag(as.numeric(psi.cor))
cat("Corrected Ga (psi * sigma2, loadings raw):\n")
cat("  diagonal:", round(diag(Ga.cor), 4), "\n")
cat("  off-diag range:", round(range(Ga.cor[lower.tri(Ga.cor)]), 4), "\n")
cat("  is PD:", all(eigen(Ga.cor)$values > 0), "\n\n")

cat("sigma2 =", round(sigma2, 6), "  (if ~1, both approaches give similar Ga)\n\n")

# -----------------------------------------------------------
# 5. Diag-downgrade mechanics for fa term
# -----------------------------------------------------------
cat("=== 5. Diag-downgrade term extraction ===\n\n")
lrt.all.rt <- attr(terms(qtlM.fa2$call$random), "term.labels")
cat("term.labels from fa(2) model:\n")
for (i in seq_along(lrt.all.rt)) cat(sprintf("  [%d] %s\n", i, lrt.all.rt[i]))

lrt.vm.idx <- grep(sterms, lrt.all.rt)
cat("\ngrep(sterms) index:", lrt.vm.idx, "\n")

if (length(lrt.vm.idx) > 0) {
    lrt.cur.term <- lrt.all.rt[lrt.vm.idx[1L]]
    cat("lrt.cur.term       :", lrt.cur.term, "\n")
    lrt.diag.pfx <- paste0("diag(", Trait, "):")
    cat("starts with diag() :", startsWith(lrt.cur.term, lrt.diag.pfx), "\n")
    lrt.vm.core    <- sub("^[^:]+:", "", lrt.cur.term)
    cat("lrt.vm.core        :", lrt.vm.core, "\n")
    lrt.diag.gterm <- paste0(lrt.diag.pfx, lrt.vm.core)
    cat("lrt.diag.gterm     :", lrt.diag.gterm, "\n\n")

    # Actually fit the downgraded model and compare logliks
    cat("Fitting diag-downgraded LRT model...\n")
    lrt.all.rt[lrt.vm.idx] <- lrt.diag.gterm
    lrt.diag.ran <- as.formula(paste("~", paste(lrt.all.rt, collapse=" + ")))
    qtlForLRT <- update(qtlM.fa2, random.=lrt.diag.ran,
                        na.action=asreml::na.method(x="include"), trace=FALSE)
    cat("LRT model converged:", qtlForLRT$converge, "\n")
    stat <- 2 * (qtlForLRT$loglik - base.6$loglik)
    pval <- 1 - pchisq.mixture(stat, ntrait=ntrait)
    cat(sprintf("LRT stat: %.4f  p-value (mixture, ntrait=%d): %.4f\n", stat, ntrait, pval))
    cat(sprintf("fa(2) vm loglik:  %.4f\n", qtlM.fa2$loglik))
    cat(sprintf("diag  vm loglik:  %.4f\n", qtlForLRT$loglik))
    cat(sprintf("base  loglik:     %.4f\n", base.6$loglik))
    cat(sprintf("stat using fa(2): %.4f  (would use this without downgrade)\n",
                2*(qtlM.fa2$loglik - base.6$loglik)))
} else {
    cat("ERROR: grep found nothing — diag-downgrade would silently fail!\n")
}
