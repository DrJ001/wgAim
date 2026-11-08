# =============================================================================
# scratch/probe_gpAim_upgrade.R
# Verify gpAim() now upgrades the MV genomic term to the mirrored / requested
# structure (previously it always stayed diag(Trait):vm).
# Run: Rscript scratch/probe_gpAim_upgrade.R
# =============================================================================
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L

cases <- list(c("diag", "NULL", "3"), c("corgh", "NULL", "3"), c("us", "NULL", "3"),
              c("fa(Trial,1)", "NULL", "3"), c("diag", "us", "3"),
              c("us", "NULL", "4"), c("fa(Trial,1)", "NULL", "4"))
for (cs in cases) {
    commandArgs <- function(...) cs
    invisible(capture.output(eval(parse(text = src[1:stop.at]))))
    rm(commandArgs)
    rf <- as.formula(if (grepl("\\(", struct)) sprintf("~ %s:id", struct)
                     else sprintf("~ %s(Trial):id", struct))
    base <- suppressWarnings(asreml(yield ~ Trial + Trial:Rep, random = rf,
                   residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                   data = pheno, workspace = "1gb", trace = FALSE))
    base <- suppressWarnings(update(base))   # embeds the formula in base$call
    cat(sprintf("\n===== base=%s  str=%s  trials=%s =====\n", struct,
                if (is.null(str.arg)) "NULL" else str.arg, n.t))
    gp <- tryCatch(suppressWarnings(suppressMessages({
              out <- NULL
              invisible(capture.output(out <- gpAim(baseModel = base, genObj = panel,
                  merge.by = "id", Trait = "Trial", str = str.arg, trace = FALSE, maxit = 50)))
              out })),
          error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
    if (is.null(gp)) next
    mod <- gp; res <- gp$GP
    cat("random formula:", paste(deparse(mod$call$random, width.cutoff = 500L), collapse = ""), "\n")
    vpn <- sub("^[^!]*", "", names(mod$vparameters)[grep("vm.*covObj", names(mod$vparameters))])
    cat(sprintf("vm params (%d): %s\n", length(vpn), paste(head(vpn, 6), collapse = " ")))
    cat("path:", res$path, " var.genetic:", round(res$var.genetic, 4), "\n")
    if (!is.null(res$Gcor)) { cat("Gcor (genomic):\n"); print(round(res$Gcor, 3)) }
}
