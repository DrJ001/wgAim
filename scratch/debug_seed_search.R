# Search sim seeds / Ga_L settings for which diag(Trial):vm(Gsave) + corgh(Trial):Gsave converges
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
src <- readLines("C:/jtaylor/Claude/wgAim/scratch/repro_corgh_3trial.R")
stop.at <- grep("^# Base model \\+ gwasAim", src) - 2L

for (cfg in list(list(2026, "c(0.06,0.05,0.07)", 0.5), list(11, "c(0.06,0.05,0.07)", 0.5),
                 list(42, "c(0.08,0.06,0.09)", 0.3), list(7, "c(0.05,0.05,0.05)", 0.6),
                 list(99, "c(0.10,0.08,0.12)", 0.4), list(123, "c(0.06,0.05,0.07)", 0.5))) {
    s <- src[1:stop.at]
    s <- sub("set.seed\\(2026\\)", sprintf("set.seed(%d)", cfg[[1]]), s)
    s <- sub("^Ga_L <- .*$", sprintf("Ga_L <- mkcov(%s, %s)", cfg[[2]], cfg[[3]]), s)
    invisible(capture.output(suppressMessages(eval(parse(text = s)))))
    covObj <- K; dimnames(covObj) <- list(geno.ids, geno.ids)
    ph <- pheno; whg <- levels(ph$id) %in% geno.ids
    ph$Gomit <- ph$Gsave <- ph$id
    levels(ph$Gsave)[!whg] <- NA; levels(ph$Gomit)[whg] <- "GEN"
    b0 <- tryCatch(suppressWarnings(asreml(yield ~ Trial:Gomit + Trial + Trial:Rep,
                    random = ~ corgh(Trial):Gsave, residual = ~ dsum(~units | Trial),
                    na.action = na.method(x = "include"), data = ph, workspace = "1gb", trace = FALSE)),
                   error = function(e) NULL)
    m <- tryCatch(suppressWarnings(update(b0, random. = ~ diag(Trial):vm(Gsave, covObj) + corgh(Trial):Gsave, maxit = 30)),
                  error = function(e) NULL)
    cat(sprintf("seed=%-5d GaL=%s r=%.1f : %s\n", cfg[[1]], cfg[[2]], cfg[[3]],
                if (is.null(m)) "ABORT" else sprintf("converge=%s LRT=%.1f", m$converge, 2*(m$loglik-b0$loglik))))
}
