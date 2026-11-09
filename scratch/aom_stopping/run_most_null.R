# =============================================================================
# scratch/aom_stopping/run_most_null.R
#
# STANDALONE MOST, step 1, under the global null (NO planted QTL).  No LRT is
# used to decide anything; the legacy LRT is computed on the same fits only
# for comparison.
#
# Protocol (per replicate):
#   1. simulate phenotypes: K-structured background (Ga_K x kscale) + IID
#      line effects (Ga_L x lscale) + residual; kscale = 0 and/or lscale = 0
#      allowed (no genetic signal at all when both are 0)
#   2. base model: UV  yield ~ Rep, random ~ id
#                  MV  us(Trial):id (falls back to diag(Trial):id if the us
#                      base fit fails -- user's base-model choice, not MOST)
#   3. engine phases 1-3 via wgAim internals (.validateModel, .buildGenoData,
#      .fixLines, .buildGenomeModel(str = "diag")): genome-wide term fitted as
#      diag(Trait):vm / vm; continued (maxit 100, twice) if not converged
#   4. trials whose genomic variance is at the boundary (con "B" or ~0) are
#      dropped; if ALL are bounded -> stop ("no-signal")
#   5. MOST on the remaining trials: max outlier statistic vs the per-trait
#      ("diag") studentised simulated null (B draws)
#   Legacy LRT for comparison: 2(L1 - L0) on the same diag vm / base fits,
#   pchisq.mixture(., ntrait) as in the package.
#
# Run : Rscript scratch/aom_stopping/run_most_null.R <pop> <ntrait> <kscale> <lscale> <nrep> <seed0> [B] [out.sub]
# Out : scratch/aom_stopping/results/<out.sub>/<pop>_t<nt>_k<k>_l<l>.csv
# =============================================================================
args   <- commandArgs(trailingOnly = TRUE)
pop    <- args[1L]; n.t <- as.integer(args[2L])
kscale <- as.numeric(args[3L]); lscale <- as.numeric(args[4L])
nrep   <- as.integer(args[5L]); seed0 <- as.integer(args[6L])
B      <- if (length(args) >= 7L) as.integer(args[7L]) else 1000L
out.sub <- if (length(args) >= 8L) args[8L] else "most_null"
stopifnot(pop %in% c("gwas", "dh"), n.t %in% c(1L, 3L))
out.dir <- file.path("scratch/aom_stopping/results", out.sub)
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
out.csv <- file.path(out.dir, sprintf("%s_t%d_k%s_l%s.csv", pop, n.t, kscale, lscale))

.libPaths(c("scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })
source("scratch/aom_stopping/aom_null_diag.R")
ns <- asNamespace("wgAim")
for (f in c(".validateModel", ".buildGenoData", ".fixLines", ".buildGenomeModel",
            "pchisq.mixture")) assign(f, get(f, envir = ns))

# ---- genotypes ----------------------------------------------------------------
if (pop == "gwas") {
    data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
    flt <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                            map.id = "marker", map.chr = "chr", map.pos = "pos"))
    genObj <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                               map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
    M <- do.call(cbind, lapply(genObj$geno, function(el) el$imputed.data))
    n.extra <- 9L; gen.type <- "marker"
} else {
    data(genoMVsim, package = "wgAim")
    genObj <- primeCross(genoMVsim, type = "interval", id = "id")
    gm0 <- primeCross(genoMVsim, type = "marker", id = "id")
    M <- do.call(cbind, lapply(gm0$geno, function(el) el$imputed.data))
    n.extra <- 0L; gen.type <- "interval"
}
geno.ids <- rownames(M); n.g <- nrow(M); p <- ncol(M)
Mc <- sweep(M, 2L, colMeans(M))
K  <- tcrossprod(Mc) / p; K <- K / mean(diag(K)) + diag(1e-4, n.g)
cK <- t(chol(K))
trials <- paste0("T", seq_len(n.t))
Trait  <- if (n.t > 1L) "Trial" else NULL
mkcov <- function(v, r) { v <- rep_len(v, n.t); S <- matrix(r, n.t, n.t); diag(S) <- 1
                          D <- diag(sqrt(v), n.t); D %*% S %*% D }
Ga_K <- mkcov(c(0.12, 0.08, 0.15), 0.6); Ga_L <- mkcov(0.05, 0.6)
mvn <- function(nr, S, sc, left = NULL) {
    if (sc == 0) return(matrix(0, nr, n.t))
    Z <- matrix(rnorm(nr * n.t), nr)
    if (!is.null(left)) Z <- left %*% Z
    Z %*% chol(sc * S)
}

sim_pheno <- function(seed) {
    set.seed(seed)
    U <- mvn(n.g, Ga_K, kscale, cK)                     # NO QTL
    L <- mvn(n.g, Ga_L, lscale)
    G <- U + L; rownames(G) <- geno.ids
    if (n.extra > 0L) {
        Ge <- if (kscale + lscale == 0) matrix(0, n.extra, n.t)
              else mvn(n.extra, kscale * Ga_K + lscale * Ga_L, 1)
        rownames(Ge) <- sprintf("NoGeno%02d", seq_len(n.extra)); G <- rbind(G, Ge)
    }
    ids <- rownames(G)
    ph <- expand.grid(id = ids, Rep = c("1", "2"), Trial = trials, stringsAsFactors = FALSE)
    ph$yield <- setNames(c(4.5, 4.8, 4.2)[seq_len(n.t)], trials)[ph$Trial] +
        G[cbind(match(ph$id, ids), match(ph$Trial, trials))] +
        ifelse(ph$Rep == "2", 0.1, 0) + rnorm(nrow(ph), sd = sqrt(0.15))
    ph$id <- factor(ph$id); ph$Rep <- factor(ph$Rep); ph$Trial <- factor(ph$Trial)
    ph
}

fit_base <- function(ph) {
    assign("pheno", ph, envir = globalenv())
    if (n.t == 1L) {
        b <- suppressWarnings(eval(quote(asreml(yield ~ Rep, random = ~ id,
                 na.action = na.method(x = "include"), data = pheno, trace = FALSE)),
                 envir = globalenv()))
        return(list(base = b, str = "uv"))
    }
    for (s in c("us", "diag")) {
        assign("rf.", as.formula(sprintf("~ %s(Trial):id", s)), envir = globalenv())
        b <- tryCatch(suppressWarnings(eval(quote(asreml(yield ~ Trial + Trial:Rep, random = rf.,
                 residual = ~ dsum(~units | Trial), na.action = na.method(x = "include"),
                 data = pheno, workspace = "1gb", maxit = 50, trace = FALSE)),
                 envir = globalenv())), error = function(e) NULL)
        if (!is.null(b)) {
            assign("base", b, envir = globalenv())
            b <- tryCatch(suppressWarnings(eval(quote(update(base, trace = FALSE)),
                          envir = globalenv())), error = function(e) b)
            return(list(base = b, str = s))
        }
    }
    stop("base model failed for us and diag")
}

most_step1 <- function(base) {
    vd <- .validateModel(base, "id", "fixed", "interval", -1)
    ph <- vd$phenoData
    gd <- .buildGenoData(genObj, gen.type, genObj$pheno[, "id"], ph[, "id"])
    fl <- invisible(.fixLines(vd$baseModel, ph, gd$genoData, "id", ph[, "id"], TRUE,
                              Trait = Trait))
    gm <- suppressWarnings(.buildGenomeModel(fl$baseModel, gd$genoData, fl$phenoData,
               fl$merge.by, genObj, FALSE, fl$rterms, globalenv(), Trait = Trait,
               str = "diag"))
    qm <- gm$qtlModel
    assign("phenoData", fl$phenoData, envir = globalenv())
    assign("qm.", qm, envir = globalenv())
    n.cont <- 0L
    while (!isTRUE(qm$converge) && n.cont < 2L) {
        qm <- suppressWarnings(eval(quote(update(qm., maxit = 100, trace = FALSE)),
                                    envir = globalenv()))
        assign("qm.", qm, envir = globalenv()); n.cont <- n.cont + 1L
    }
    # legacy LRT on the same fits (comparison only)
    stat  <- 2 * (qm$loglik - fl$baseModel$loglik)
    p.lrt <- if (n.t == 1L) (1 - pchisq(stat, 1)) / 2 else 1 - pchisq.mixture(stat, n.t)
    # bounded trials
    s2 <- qm$sigma2; if (qm$vparameters.con[length(qm$vparameters.con)] == 4) s2 <- 1
    sel <- grep("vm.*covObj", names(qm$vparameters))
    val <- as.numeric(qm$vparameters[sel]) * s2
    con <- as.character(qm$vparameters.con[sel])
    lev <- if (n.t == 1L) 1L else match(sub(".*!Trial_", "", names(qm$vparameters)[sel]), trials)
    bnd <- con %in% c("B", "S", "F") | val <= 1e-6 * max(abs(val), 1e-12) | val <= 0
    keep <- sort(lev[!bnd])
    res <- list(conv = isTRUE(qm$converge), n.cont = n.cont, stat = stat, p.lrt = p.lrt,
                n.bounded = sum(bnd), d = NULL)
    if (!length(keep)) return(res)
    res$d <- aomNullDiag(qm, fl$phenoData, gm$intervalObj, gd$state, fl$merge.by,
                         Trait, n.t, B = B,
                         keep.traits = if (length(keep) < n.t) keep else NULL)$summary
    res
}

cat(sprintf("STANDALONE MOST NULL: pop=%s ntrait=%d kscale=%s lscale=%s nrep=%d seeds %d..%d B=%d\n",
            pop, n.t, kscale, lscale, nrep, seed0, seed0 + nrep - 1L, B))
for (r in seq_len(nrep)) {
    seed <- seed0 + r - 1L; t0 <- Sys.time()
    row <- tryCatch({
        fb <- fit_base(sim_pheno(seed))
        s1 <- most_step1(fb$base)
        d  <- s1$d
        data.frame(pop = pop, ntrait = n.t, kscale = kscale, lscale = lscale, seed = seed,
                   base.str = fb$str, vm.conv = s1$conv, n.cont = s1$n.cont,
                   n.bounded = s1$n.bounded,
                   status = if (is.null(d)) "no-signal" else "ok",
                   lrt.stat = s1$stat, p.lrt = s1$p.lrt,
                   max.oint = if (is.null(d)) NA else d$max.oint,
                   p.most = if (is.null(d)) 1 else d$p.fwe.studD,
                   p.orig = if (is.null(d)) 1 else d$p.fwe,
                   p.studM = if (is.null(d)) 1 else d$p.fwe.studM,
                   p.bonf = if (is.null(d)) 1 else d$p.bonf,
                   thr.most = if (is.null(d)) NA else d$fwe.thr.studD,
                   diag.pin.dev = if (is.null(d)) NA else d$diag.pin.dev,
                   Meff = if (is.null(d)) NA else d$Meff,
                   null.mean = if (is.null(d)) NA else d$null.mean.oint)
    }, error = function(e)
        data.frame(pop = pop, ntrait = n.t, kscale = kscale, lscale = lscale, seed = seed,
                   base.str = NA, vm.conv = NA, n.cont = NA, n.bounded = NA,
                   status = paste("error:", substr(conditionMessage(e), 1, 80)),
                   lrt.stat = NA, p.lrt = NA, max.oint = NA, p.most = NA, p.orig = NA,
                   p.studM = NA, p.bonf = NA, thr.most = NA, diag.pin.dev = NA,
                   Meff = NA, null.mean = NA))
    row$secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    write.table(row, out.csv, sep = ",", row.names = FALSE,
                col.names = !file.exists(out.csv), append = file.exists(out.csv))
    cat(sprintf("rep %3d seed %d  %-10s base %-4s conv %s bnd %s  LRT p %.4f  MOST p %.4f  orig %.4f  pin.dev %.4f  [%.0fs]\n",
                r, seed, substr(row$status, 1, 10), row$base.str, row$vm.conv,
                as.character(row$n.bounded), row$p.lrt, row$p.most, row$p.orig,
                row$diag.pin.dev, row$secs))
}
cat("DONE\n")
