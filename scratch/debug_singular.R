# Debug: confirm singularity is fixed by na.action and check identical lines
devtools::load_all(); library(asreml); library(qtl); library(MASS)
set.seed(101)
n.lines <- 200; n.chr <- 5; n.mar <- 200; mar.step <- 0.5
map.sim <- sim.map(len=rep((n.mar-1)*mar.step,n.chr), n.mar=rep(n.mar,n.chr),
                   eq.spac=TRUE, include.x=FALSE, sex.sp=FALSE)
names(map.sim) <- paste0("Chr", seq_len(n.chr))
cross.sim <- sim.cross(map.sim, n.ind=n.lines, type="bc", error.prob=0.005)
cross.sim$pheno$id <- paste0("DH", sprintf("%03d", seq_len(n.lines)))
ids <- as.character(cross.sim$pheno$id)
genObj <- primeCross(cross.sim, type="interval", id="id")

# Check for identical lines in imputed.data (marker matrix)
imp.mat <- do.call("cbind", lapply(genObj$geno, function(el) el$imputed.data))
n.unique <- nrow(unique(imp.mat))
cat("Total lines:", n.lines, " | Unique lines (imputed.data):", n.unique, "\n")
dup.rows <- which(duplicated(imp.mat))
if (length(dup.rows)) cat("Duplicate line indices:", dup.rows, "\n") else cat("No duplicate lines\n")

# Check interval data too
int.mat <- do.call("cbind", lapply(genObj$geno, function(el) el$interval.data))
n.unique.int <- nrow(unique(round(int.mat, 8)))
cat("Unique lines (interval.data):", n.unique.int, "\n")

# Confirm: without na.action -> singularity; with na.action -> OK
pheno.uni <- do.call(rbind, lapply(1:2, function(r)
    data.frame(id=factor(ids), Rep=factor(r),
               yield=10 + rnorm(n.lines,0,1),
               stringsAsFactors=FALSE)))
pheno.uni <- pheno.uni[order(pheno.uni$id, pheno.uni$Rep),]
base.uni <- asreml(fixed=yield~Rep, random=~id, data=pheno.uni, trace=FALSE)

cat("\n--- Without na.action ---\n")
tryCatch(
    qtlAim(base.uni, genObj=genObj, merge.by="id", TypeI=0.05,
           breakout=1, trace=FALSE),
    error=function(e) cat("ERROR:", conditionMessage(e), "\n")
)

cat("\n--- With na.action ---\n")
tryCatch(
    qtlAim(base.uni, genObj=genObj, merge.by="id", TypeI=0.05,
           breakout=1, trace=FALSE,
           na.action=asreml::na.method(x="include")),
    error=function(e) cat("ERROR:", conditionMessage(e), "\n")
)
cat("Done\n")
