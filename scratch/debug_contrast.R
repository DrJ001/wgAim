devtools::load_all("C:/jtaylor/Research/Computing/R Packages/Claude/wgaim", quiet=TRUE)
library(asreml)
asreml.options(trace=FALSE)
data(genoRxK, package="wgAim")
data(phenoRxK, package="wgAim")
pheno <- phenoRxK
pheno$Genotype <- factor(pheno$Genotype)
pheno$Row      <- factor(pheno$Row)
pheno$Range    <- factor(pheno$Range)
rxk <- primeCross(genoRxK, type="interval", id="Genotype")
base_yld <- asreml(fixed=yld~Range, random=~Genotype+Row, data=pheno,
                   na.action=na.method(x="include"), trace=FALSE)
base_yld <- update(base_yld, trace=FALSE)
qtl_yld  <- qtlAim(base_yld, genObj=rxk, merge.by="Genotype",
                   exclusion.window=50, TypeI=0.05, trace=FALSE)

# phenoData is assigned to the caller env by qtlAim() — check it exists
cat("phenoData exists:", exists("phenoData"), "\n")
cat("yld.data exists:", exists("yld.data"), "\n")
data <- as.data.frame(phenoData)
cat("phenoData class:", paste(class(phenoData), collapse=", "), "\n")
cat("phenoData nrow:", nrow(data), "\n")
cat("phenoData names (first 15):", paste(head(names(data), 15), collapse=", "), "\n")

effects <- qtl_yld$QTL$effects
cat("\neffects names:", paste(names(effects), collapse=", "), "\n")

qtl_cols <- gsub("Chr\\.", "X.", names(effects))
qtl_cols_present <- qtl_cols[qtl_cols %in% names(data)]
cat("qtl_cols present in phenoData:", paste(qtl_cols_present, collapse=", "), "\n")

id_col   <- qtl_yld$QTL$diag$genetic.term
cat("id_col:", id_col, "\n")
line_ids <- as.character(data[[id_col]])
cat("line_ids length:", length(line_ids), "\n")
cat("NAs in line_ids:", sum(is.na(line_ids)), "\n")

qtl_mat <- as.matrix(data[, qtl_cols_present, drop=FALSE])
cat("qtl_mat dim:", nrow(qtl_mat), "x", ncol(qtl_mat), "\n")
cat("NAs in qtl_mat:", sum(is.na(qtl_mat)), "\n")

eff_vec     <- effects[match(qtl_cols_present, gsub("Chr\\.", "X.", names(effects)))]
qtl_contrib <- as.vector(qtl_mat %*% eff_vec)
cat("qtl_contrib length:", length(qtl_contrib), "\n")
cat("NAs in qtl_contrib:", sum(is.na(qtl_contrib)), "\n")

cat("\nAttempting aggregate...\n")
tryCatch({
    contrib_df <- stats::aggregate(
        data.frame(qtl_contrib=qtl_contrib),
        by=list(line=line_ids),
        FUN=mean)
    cat("aggregate OK, nrow:", nrow(contrib_df), "\n")
}, error=function(e) cat("aggregate ERROR:", conditionMessage(e), "\n"))

cat("\nAttempting full contrast plot...\n")
tryCatch({
    p <- plot(qtl_yld, genObj=rxk, type="contrast", data=phenoData)
    cat("plot OK, class:", paste(class(p), collapse=", "), "\n")
}, error=function(e) cat("plot ERROR:", conditionMessage(e), "\n"))
