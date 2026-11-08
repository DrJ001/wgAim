# Debug: check what ASReml names the corgh(Trial):vm() term internally
devtools::load_all()
library(asreml); library(qtl); library(MASS)
set.seed(101)

n.lines <- 200; n.chr <- 5; n.mar <- 200; mar.step <- 0.5
map.sim <- sim.map(len=rep((n.mar-1)*mar.step,n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)

genObj <- primeCross(cross.sim, type="interval", id="id")
M  <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K  <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

# Build a tiny 2-trial pheno
qtl.info <- list(list(chr="Chr1",int=50L,eff=0.80), list(chr="Chr2",int=80L,eff=-0.65))
n.qtl <- 2L
qtl.X <- matrix(0, nrow=n.lines, ncol=n.qtl, dimnames=list(ids, c("q1","q2")))
for (i in 1:n.qtl)
    qtl.X[,i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[,qtl.info[[i]]$int]
qtl.main <- rowSums(sweep(qtl.X,2,c(0.80,-0.65),"*"))

n.t2 <- 2L; Ga2 <- matrix(c(0.55,0.35,0.35,0.50),2,2)
U2 <- matrix(MASS::mvrnorm(1,rep(0,n.lines*n.t2),kronecker(Ga2,K)),
             nrow=n.lines, ncol=n.t2, dimnames=list(ids,paste0("Trial",1:2)))
pheno.2 <- do.call(rbind, lapply(1:n.t2, function(j)
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id=factor(ids), Trial=factor(paste0("Trial",j)),
                   Rep=factor(paste0("R",r)),
                   yield=c(10,11.5)[j]+qtl.main+U2[,j]+rnorm(n.lines,0,0.5))
    ))
))
pheno.2 <- pheno.2[order(pheno.2$Trial,pheno.2$Rep,pheno.2$id),]

# Fit base corh model
base.2 <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                 residual=~dsum(~units|Trial), data=pheno.2, trace=FALSE)
cat("base.2 converged:", base.2$converge, "\n")

# --- Now manually build the genome model exactly as .buildGenomeModel does ---
genoData <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
nint <- lapply(genObj$geno, function(el) 1:ncol(el$interval.data))
lint <- unlist(lapply(nint, length))
mnams <- paste("Chr", rep(names(genObj$geno), times=lint), unlist(nint), sep=".")
dimnames(genoData) <- list(ids, mnams)

rm.env <- wgAim:::.constructCM(genoData)
covObj <- rm.env$relm
attr(genObj, "env") <- rm.env

# Build the vm model manually (what .buildGenomeModel does with Trait="Trial")
ran.1 <- as.formula("~ diag(Trial):vm(id, covObj) + corh(Trial):id")
cat("\nFitting diag(Trial):vm(id,covObj) + corh(Trial):id ...\n")
m1 <- update(base.2, random.=ran.1)
cat("m1 converged:", m1$converge, "\n")

# Check vparameters names
cat("\nm1 vparameters names:\n")
print(names(m1$vparameters))

cat("\nm1 random terms (from call):\n")
rterms <- attr(terms(m1$call$random), "term.labels")
print(rterms)

cat("\nm1 G.param names:\n")
print(names(m1$G.param))

# Now upgrade to corgh
ran.2 <- as.formula("~ corgh(Trial):vm(id, covObj) + corh(Trial):id")
cat("\nFitting corgh(Trial):vm(id,covObj) + corh(Trial):id ...\n")
m2 <- update(m1, random.=ran.2)
cat("m2 converged:", m2$converge, "\n")

cat("\nm2 vparameters names:\n")
print(names(m2$vparameters))
cat("\nm2 random terms:\n")
print(attr(terms(m2$call$random), "term.labels"))

# Now try predict with the vmterm.raw as stored
vmterm.raw <- attr(terms(m2$call$random), "term.labels")[1]
cat("\nvmterm.raw =", vmterm.raw, "\n")
cat("Attempting predict(m2, classify='id:Trial', only=vmterm.raw)...\n")
pv <- tryCatch(
    predict(m2, classify="id:Trial", only=vmterm.raw, vcov=TRUE, data=pheno.2),
    error=function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(pv)) cat("predict() succeeded, nrow pvals:", nrow(pv$pvals), "\n")
