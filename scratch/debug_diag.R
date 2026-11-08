devtools::load_all(); library(asreml); library(qtl); library(MASS)
set.seed(101)
n.lines <- 200; n.chr <- 5; n.mar <- 80; mar.step <- 1.0
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

qtl.info <- list(list(chr="Chr1",int=20L,eff=0.80), list(chr="Chr2",int=55L,eff=-0.65),
                 list(chr="Chr3",int=30L,eff=0.55),  list(chr="Chr4",int=60L,eff=-0.50))
n.qtl <- 4L
qtl.X <- matrix(0, nrow=n.lines, ncol=n.qtl, dimnames=list(ids,NULL))
for (i in 1:n.qtl)
    qtl.X[,i] <- genObj$geno[[qtl.info[[i]]$chr]]$interval.data[,qtl.info[[i]]$int]
qtl.main <- as.numeric(qtl.X %*% sapply(qtl.info,"[[","eff"))
u.uni <- MASS::mvrnorm(1, rep(0,n.lines), 0.05*K); names(u.uni) <- ids

Ga2 <- matrix(c(0.02,0.75*sqrt(0.02*0.02),0.75*sqrt(0.02*0.02),0.02),2,2)
qtl.eff.2 <- matrix(c(0.80,0.80,-0.65,-0.65,0.55,0.15,-0.50,0.45),nrow=4,byrow=TRUE)
U2 <- matrix(MASS::mvrnorm(1,rep(0,n.lines*2),kronecker(Ga2,K)),
             nrow=n.lines, ncol=2, dimnames=list(ids,c("Trial1","Trial2")))
phenoMV <- do.call(rbind, lapply(1:2, function(j) {
    qc <- rowSums(sweep(qtl.X,2,qtl.eff.2[,j],"*"))
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id=factor(ids),Trial=factor(paste0("Trial",j)),
                   Rep=factor(paste0("R",r)),
                   yield=c(10,11.5)[j]+qc+U2[,j]+rnorm(n.lines,0,sqrt(0.70)))))
}))
phenoMV <- phenoMV[order(phenoMV$Trial, phenoMV$Rep, phenoMV$id),]

base.2 <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                 residual=~dsum(~units|Trial), data=phenoMV, trace=FALSE,
                 na.action=asreml::na.method(x="include"), maxit=40)
cat("base converged:", base.2$converge, "\n")

# Manually call buildGenomeModel to get the initial qtlModel
# then test buildDiagVmModel directly
genoData <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
mnams <- paste("Chr", rep(names(genObj$geno), times=sapply(genObj$geno, function(g) ncol(g$interval.data))),
               unlist(lapply(genObj$geno, function(g) seq_len(ncol(g$interval.data)))), sep=".")
dimnames(genoData) <- list(ids, mnams)

plines <- phenoMV[["id"]]
glines <- genObj$pheno[["id"]]
genoData <- genoData[rownames(genoData) %in% as.character(plines),,drop=FALSE]

# Build initial vm model manually
cov.env <- wgAim:::.constructCM(genoData)
covObj  <- cov.env$relm
gterm   <- "diag(Trial):vm(id, covObj)"
resid.term <- "corh(Trial):id"
ran.form <- as.formula(paste("~", gterm, "+", resid.term))
qtlModel <- base.2
qtlModel$call$data <- quote(phenoMV)
qtlModel <- update(qtlModel, random. = ran.form,
                   na.action=asreml::na.method(x="include"))
cat("diag vm model converged:", qtlModel$converge, "\n")
cat("vmterms[1]:", gterm, "\n")

# Upgrade to corgh
ran.corgh <- as.formula(paste("~ corgh(Trial):vm(id, covObj) +", resid.term))
qtlCorgh <- update(qtlModel, random. = ran.corgh,
                   na.action=asreml::na.method(x="include"))
cat("corgh vm model converged:", qtlCorgh$converge, "\n")
vmterms <- c("corgh(Trial):vm(id, covObj)", resid.term)

# Now test buildDiagVmModel
cat("\nTesting .buildDiagVmModel()...\n")
cat("vmterms[1]:", vmterms[1], "\n")
cat("startsWith diag:", startsWith(vmterms[1], "diag(Trial):"), "\n")
diag.back <- tryCatch(
    wgAim:::.buildDiagVmModel(qtlCorgh, vmterms, "Trial",
                               na.action=asreml::na.method(x="include")),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(diag.back))
    cat("Back to diag converged:", diag.back$converge,
        " loglik:", diag.back$loglik, "\n")
