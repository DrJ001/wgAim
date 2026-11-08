devtools::load_all(quiet = TRUE)
library(asreml)
data(genoRxK, package = "wgAim"); data(phenoRxK, package = "wgAim")
rxk_int <- primeCross(genoRxK, type = "interval", id = "Genotype")
pheno <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row <- factor(pheno$Row)
pheno$Range <- factor(pheno$Range)
base_yld <- asreml::asreml(fixed = yld ~ Range, random = ~ Genotype + Row,
    data = pheno, na.action = asreml::na.method(x = "include"), trace = FALSE)
qtl_yld <- qtlAim(baseModel = base_yld, genObj = rxk_int, merge.by = "Genotype",
                   TypeI = 0.05, trace = FALSE)

state <- qtl_yld$QTL$diag$state
cat("Length of state:", length(state), "\n")
cat("First 6 state names:", paste(head(names(state), 6), collapse = ", "), "\n\n")

# How many markers / intervals per chr?
chr1a <- rxk_int$geno[["1A"]]
cat("Chr 1A: markers =", length(chr1a$map),
    "  interval.data cols =", ncol(chr1a$interval.data),
    "  inferred.map length =", length(chr1a$inferred.map), "\n")

# So state is named "Chr.<chr>.<idx>" — interval indices
# The genoData we need for fineMap's covObj rebuild is interval.data, not imputed.data
cat("\nTotal interval cols across all chrs:",
    sum(sapply(rxk_int$geno, function(el) ncol(el$interval.data))), "\n")
cat("Total marker cols across all chrs:",
    sum(sapply(rxk_int$geno, function(el) ncol(el$imputed.data))), "\n")
cat("State length (should match interval cols):", length(state), "\n")
