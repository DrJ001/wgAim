devtools::load_all(quiet = TRUE)
library(asreml)
data(genoRxK, package = "wgAim")
data(phenoRxK, package = "wgAim")
rxk_int <- primeCross(genoRxK, type = "interval", id = "Genotype")
pheno <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)
base_yld <- asreml::asreml(
    fixed = yld ~ Range, random = ~ Genotype + Row, data = pheno,
    na.action = asreml::na.method(x = "include"), trace = FALSE)
qtl_yld <- qtlAim(baseModel = base_yld, genObj = rxk_int,
                   merge.by = "Genotype", TypeI = 0.05, trace = FALSE)

cat("call$data:", deparse(qtl_yld$call$data), "\n")
fe <- environment(qtl_yld$call$fixed)
cat("formula env objects:", paste(ls(fe), collapse = ", "), "\n")
cat("exists 'phenoData' in formula env:", exists("phenoData", envir = fe, inherits = FALSE), "\n")
cat("exists 'yld.data'  in .GlobalEnv:", exists("yld.data",  envir = .GlobalEnv, inherits = FALSE), "\n")
cat("exists 'phenoData' in .GlobalEnv:", exists("phenoData", envir = .GlobalEnv, inherits = FALSE), "\n")

# What does the model frame look like?
cat("\nModel frame class:", class(qtl_yld$mf), "\n")
cat("mf dim:", dim(qtl_yld$mf), "\n")
cat("mf head names:", paste(head(names(qtl_yld$mf), 6), collapse = ", "), "\n")
