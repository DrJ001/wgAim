suppressPackageStartupMessages(devtools::load_all())

data(genoGWASraw,  package = "wgAim")
data(mapGWASsim,   package = "wgAim")
data(phenoGWASsim, package = "wgAim")

cat("\n=========================== checkPanel ==============================\n")
chk <- checkPanel(genoGWASraw, mapGWASsim, id = "id",
                   map.id = "marker", map.chr = "chr", map.pos = "pos",
                   encoding = "012")
print(chk)

cat("\n=========================== filterPanel =============================\n")
flt <- filterPanel(chk)   # default steps
print(flt)

cat("\n=========================== primePanel (knn) ========================\n")
panel <- primePanel(flt, impute = "knn", knn.k = 5)
cat("\npanel classes:", paste(class(panel), collapse = ", "), "\n")
cat("panel chromosomes:", length(panel$geno), "\n")
n.mk <- sum(sapply(panel$geno, function(el) ncol(el$imputed.data)))
cat("panel markers   :", n.mk, "\n")
cat("panel lines     :", nrow(panel$pheno), "\n")
cat("any NA remaining:", anyNA(do.call("cbind",
       lapply(panel$geno, function(el) el$imputed.data))), "\n")

cat("\n============================ Phenotype ==============================\n")
cat("Structure:\n")
str(phenoGWASsim)
cat("\nSites x Reps table:\n")
print(table(phenoGWASsim$Site, phenoGWASsim$Rep))
cat("\nPer-site yield summary:\n")
ysumm <- t(sapply(split(phenoGWASsim$yield, phenoGWASsim$Site),
                   function(y) c(mean = mean(y), sd = sd(y),
                                  min = min(y), max = max(y))))
print(round(ysumm))

cat("\n=================== QTL markers survived pipeline ==================\n")
qtl.expected <- c("SNP_C02_0095", "SNP_C05_0183", "SNP_C07_0060",
                   "SNP_C03_0214", "SNP_C09_0123")
panel.mk <- unlist(lapply(panel$geno,
                            function(el) colnames(el$imputed.data)),
                     use.names = FALSE)
for (m in qtl.expected)
    cat(sprintf("  %-16s in panel: %s\n", m, m %in% panel.mk))
