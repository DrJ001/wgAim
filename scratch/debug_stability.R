# Debug: print raw coef.list to check ordering and signs
devtools::load_all(); library(asreml); library(qtl); library(MASS)
set.seed(101)

n.lines <- 200; n.chr <- 5; n.mar <- 200; mar.step <- 0.5
map.sim <- sim.map(len=rep((n.mar-1)*mar.step,n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno[["id"]] <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno[["id"]])
genObj <- primeCross(cross.sim, type="interval", id="id")
M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

qtl.eff <- c(0.80, -0.65, 0.55, -0.50)
qtl.chr <- c("Chr1","Chr2","Chr3","Chr4")
qtl.int <- c(50L, 80L, 30L, 120L)
qtl.X <- matrix(0, nrow=n.lines, ncol=4, dimnames=list(ids,NULL))
for (i in 1:4)
    qtl.X[,i] <- genObj$geno[[qtl.chr[i]]]$interval.data[,qtl.int[i]]
qtl.main <- as.numeric(qtl.X %*% qtl.eff)
u.uni <- as.numeric(MASS::mvrnorm(1, rep(0,n.lines), 0.05*K))
names(u.uni) <- ids

rn <- function(r) 10 + c(0.25,-0.25)[r] + qtl.main + u.uni + rnorm(n.lines,0,sqrt(0.80))
pheno.uni <- do.call(rbind, lapply(1:2, function(r)
    data.frame(id=factor(ids), Rep=factor(r), yield=rn(r))))
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep),]

base.uni <- asreml(fixed=yield~Rep, random=~id, data=pheno.uni, trace=FALSE,
                   na.action=asreml::na.method(x="include"))
fit <- qtlAim(base.uni, genObj=genObj, merge.by="id", TypeI=0.05,
              exclusion.window=20, trace=FALSE,
              na.action=asreml::na.method(x="include"))

cat("QTL found:", fit$QTL$qtl, "\n\n")
cat("coef.list (newest-first via rev()):\n")
for (j in seq_along(fit$QTL$diag$coef.list)) {
    cl_j <- fit$QTL$diag$coef.list[[j]]
    cat("  iter", j, ":", paste(sprintf("%s=%.4f", names(cl_j), cl_j), collapse=", "), "\n")
}

cat("\n.build_stability_df extraction (k=QTL rank, j=iter, idx=j-k+1):\n")
cl     <- fit$QTL$diag$coef.list
n_qtl  <- length(fit$QTL$qtl)
sigma2 <- fit$sigma2
for (k in seq_len(n_qtl)) {
    true_eff <- if (k <= 4) qtl.eff[k] else NA
    cat("  QTL", k, "(", fit$QTL$qtl[k], ") true eff =", true_eff, "\n")
    for (j in k:n_qtl) {
        idx <- j - k + 1L
        eff <- cl[[j]][idx]
        cat("    iter", j, "idx", idx, ": eff =", round(eff,4),
            " name =", names(cl[[j]])[idx], "\n")
    }
}
