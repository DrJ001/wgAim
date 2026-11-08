devtools::load_all(); library(asreml); library(qtl); library(MASS)
set.seed(101)
n.lines <- 50; n.chr <- 2; n.mar <- 20; mar.step <- 2.0
map.sim <- sim.map(len=rep((n.mar-1)*mar.step,n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr",1:n.chr)
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("L",sprintf("%02d",1:n.lines))
ids <- as.character(cross.sim$pheno$id)
genObj <- primeCross(cross.sim, type="interval", id="id")
M <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
rownames(M) <- ids
K <- tcrossprod(M)/ncol(M); K <- K/mean(diag(K))

qtl.X <- genObj$geno[["Chr1"]]$interval.data[,10L,drop=FALSE]
qtl.eff2 <- c(0.8,0.8)
U2 <- matrix(MASS::mvrnorm(1,rep(0,n.lines*2),kronecker(diag(0.02,2),K)),
             nrow=n.lines, ncol=2, dimnames=list(ids,c("Trial1","Trial2")))
pheno2 <- do.call(rbind, lapply(1:2, function(j)
    do.call(rbind, lapply(1:2, function(r)
        data.frame(id=factor(ids), Trial=factor(paste0("Trial",j)),
                   Rep=factor(paste0("R",r)),
                   yield=c(10,11)[j]+qtl.eff2[j]*qtl.X[,1]+U2[,j]+rnorm(n.lines,0,.8))))))
pheno2 <- pheno2[order(pheno2$Trial,pheno2$Rep,pheno2$id),]

base2 <- asreml(fixed=yield~Trial+Trial:Rep, random=~corh(Trial):id,
                residual=~dsum(~units|Trial), data=pheno2, trace=FALSE,
                na.action=asreml::na.method(x="include"))
cat("base converged:", base2$converge, "\n")

covObj <- wgAim:::.constructCM(M)$relm

# Simulate what happens in the loop: qtlModel has corgh vm term + QTL fixed effects
ran1 <- as.formula("~ diag(Trial):vm(id, covObj) + corh(Trial):id")
m1 <- update(base2, random.=ran1, na.action=asreml::na.method(x="include"))
cat("diag vm converged:", m1$converge, "\n")
ran2 <- as.formula("~ corgh(Trial):vm(id, covObj) + corh(Trial):id")
m2 <- update(m1, random.=ran2, na.action=asreml::na.method(x="include"))
cat("corgh vm converged:", m2$converge, "\n")

vmterms <- c("corgh(Trial):vm(id, covObj)", "corh(Trial):id")
cat("Testing buildDiagVmModel from within a nested function...\n")

# Simulate being called from inside qtlAim loop
test_fn <- function(qtlModel, vmterms, Trait, ...) {
    cat("Inside test_fn, phenoData in parent.frame:", exists("pheno2", envir=parent.frame()), "\n")
    # This simulates what .buildDiagVmModel does
    diag.prefix <- paste0("diag(", Trait, "):")
    vm.core     <- sub("^[^:]+:", "", vmterms[1L])
    diag.gterm  <- paste0(diag.prefix, vm.core)
    all.rterms  <- attr(terms(qtlModel$call$random), "term.labels")
    vm.idx      <- grep("vm.*covObj|mbf.*ints", all.rterms)
    all.rterms[vm.idx] <- diag.gterm
    diag.ran    <- as.formula(paste("~", paste(all.rterms, collapse=" + ")))
    cat("Formula:", deparse(diag.ran), "\n")
    cat("qtlModel$call$data:", deparse(qtlModel$call$data), "\n")
    qtlModel$call$data <- quote(pheno2)
    tryCatch(update(qtlModel, random.=diag.ran, ...),
             error=function(e) cat("ERROR:", conditionMessage(e), "\n"))
}
res <- test_fn(m2, vmterms, "Trial", na.action=asreml::na.method(x="include"))
if (!is.null(res) && !inherits(res,"character"))
    cat("Result converged:", res$converge, " loglik:", res$loglik, "\n")

# Now wrap in another level of indirection (like the engine does)
outer_fn <- function(qtlModel, vmterms, Trait, ...) {
    phenoData <- pheno2  # simulate qtlAim.asreml's local phenoData
    cat("Calling with phenoData in local scope...\n")
    qtlModel$call$data <- quote(phenoData)
    test_fn(qtlModel, vmterms, Trait, ...)
}
res2 <- outer_fn(m2, vmterms, "Trial", na.action=asreml::na.method(x="include"))
if (!is.null(res2) && !inherits(res2,"character"))
    cat("Result2 converged:", res2$converge, " loglik:", res2$loglik, "\n")
