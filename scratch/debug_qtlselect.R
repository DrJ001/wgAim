# Instrument .qtlSelect() inside gwasAim() to inspect Ga, predict output, oint
args   <- commandArgs(trailingOnly = TRUE)
struct <- if (length(args) >= 1L) args[1L] else "corgh"
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })

src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
src <- sub('^struct <- .*$', sprintf('struct <- "%s"', struct), src)
stop.at <- grep("^gw <- gwasAim", src) - 1L
invisible(capture.output(eval(parse(text = src[1:stop.at]))))
planted <- colnames(M)[qtl.cols]
# marker index within chromosome (state names are Chr.<chr>.<idx within chr>)
planted.key <- vapply(qtl.cols, function(j) {
    ch <- rep(names(panel$geno), vapply(panel$geno, function(g) ncol(g$imputed.data), 1L))[j]
    paste0("Chr.", ch, ".", j - min(which(rep(names(panel$geno),
        vapply(panel$geno, function(g) ncol(g$imputed.data), 1L)) == ch)) + 1L)
}, "")

ns <- asNamespace("wgAim")
orig <- get(".qtlSelect", ns)
dbg <- function(asm, phenoData, intervalObj, gen.type, selection,
                exclusion.window, state, verboseLev, merge.by = NULL,
                Trait = NULL, ntrait = 1L) {
    cat("\n================ DEBUG .qtlSelect ================\n")
    rterms.all <- attr(terms.formula(asm$call$random), "term.labels")
    cat("random term labels:\n"); print(rterms.all)
    vp <- asm$vparameters
    sel <- grep("vm.*covObj|mbf.*ints", names(vp))
    cat("vm vparameters:\n"); print(vp[sel])
    cat("G.param names:\n"); print(names(asm$G.param))
    vmterm.raw <- rterms.all[grep("vm.*covObj", rterms.all)]
    cat("vmterm.raw:", vmterm.raw, "\n")
    pv <- tryCatch(predict(asm, classify = paste(Trait, merge.by, sep = ":"),
                           only = vmterm.raw, vcov = TRUE, data = phenoData, maxit = 1),
                   error = function(e) { cat("predict ERROR:", conditionMessage(e), "\n"); NULL })
    if (!is.null(pv)) {
        cat("pv names:", names(pv), "\n")
        cat("pvals dim:", dim(pv$pvals), " vcov is NULL:", is.null(pv$vcov), "\n")
        print(head(pv$pvals)); print(summary(pv$pvals$predicted.value))
    }
    res <- orig(asm, phenoData, intervalObj, gen.type, selection, exclusion.window,
                state, verboseLev, merge.by, Trait, ntrait)
    cat("oint summary:\n"); print(summary(res$oint[as.logical(state)]))
    cat("selected:", res$qtl, "\n")
    oi <- res$oint; rk <- rank(-oi, ties.method = "min")
    pl <- planted.key
    cat("planted:", paste0(planted, " -> ", pl, " oint=", round(oi[pl], 2), " rank=", rk[pl]), sep = "\n  ")
    top <- head(order(-oi), 5)
    cat("top 5:", paste0(names(oi)[top], " (", round(oi[top], 2), ")"), "\n")
    print(round(Ga_dbg <- wgAim:::.extractVmGa(asm, Trait, levels(phenoData[[Trait]]), 1), 4))
    cat("==================================================\n")
    res
}
environment(dbg) <- ns
unlockBinding(".qtlSelect", ns); assign(".qtlSelect", dbg, envir = ns)

# Also print Ga via trace on MASS::ginv


gw <- gwasAim(baseModel = base, genObj = panel, merge.by = "id",
              Trait = "Trial", str = str.arg, TypeI = 0.05,
              trace = FALSE, breakout = 3, maxit = 50)
