# =============================================================================
# scratch/check_corgh_order.R
#
# Colleague's question: ASReml lists corgh correlations row-wise
# (2,1),(3,1),(3,2),(4,1),...; the OLD wgAim loop filled Ga column-wise
# (2,1),(3,1),(4,1),(3,2),...  -- identical for t <= 3, different for t >= 4.
#
# This script, for t = 4 and 5 trials with clearly DISTINCT true correlations:
#   1. fits corgh(Trial):vm(id, covObj) + diag(Trial):id
#   2. prints the .cor parameter names in ASReml's order
#   3. rebuilds Ga (a) with the OLD positional column-wise fill,
#                  (b) with the NEW name-based .extractVmGa()
#   4. compares both against the us(Trial):vm fit of the same data
#      (us estimates covariances directly, so it is an independent reference)
# =============================================================================
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages({ library(asreml); library(wgAim) })

data(genoGWASraw, package = "wgAim"); data(mapGWASsim, package = "wgAim")
flt   <- suppressMessages(filterPanel(genoGWASraw, mapGWASsim, id = "id",
                     map.id = "marker", map.chr = "chr", map.pos = "pos"))
panel <- suppressMessages(primePanel(flt, id = "id", map.id = "marker",
                    map.chr = "chr", map.pos = "pos", impute = "knn", knn.k = 5L))
M <- do.call(cbind, lapply(panel$geno, function(el) el$imputed.data))
ids <- rownames(M); n.g <- nrow(M)
Mc <- sweep(M, 2L, colMeans(M)); K <- tcrossprod(Mc) / ncol(M)
K <- K / mean(diag(K)) + diag(1e-4, n.g)
covObj <- K; dimnames(covObj) <- list(ids, ids)

old_fill <- function(model, ntrait) {
    # verbatim logic of the pre-fix corgh branch in .qtlSelect()
    vp <- model$vparameters; sel <- grep("vm.*covObj", names(vp))
    nms <- names(vp)[sel]; raw <- vp[sel]
    is.cor <- grepl("\\.cor$", nms)
    v <- pmax(raw[!is.cor], 0); r <- raw[is.cor]; s <- sqrt(v)
    Ga <- diag(v); idx <- 0L
    for (col in seq_len(ntrait - 1L)) for (row in (col + 1L):ntrait) {
        idx <- idx + 1L
        Ga[row, col] <- Ga[col, row] <- r[idx] * s[row] * s[col]
    }
    unname(Ga)
}
cor_of <- function(G) { d <- 1 / sqrt(diag(G)); G * tcrossprod(d) }

runs <- list(c(4L, 1L), c(4L, 2L), c(4L, 3L), c(5L, 105L))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) runs <- lapply(strsplit(args, ":"), as.integer)
for (rs in runs) {
    nt <- rs[1L]
    cat(sprintf("\n######################## %d trials (seed %d) ########################\n", nt, rs[2L]))
    set.seed(rs[2L])
    tl <- paste0("T", seq_len(nt))
    # Distinct true correlations: AR(1)-like in trial index, so (4,1) is
    # much weaker than (3,2) -- misplacement is easy to see.
    Rt <- 0.9^abs(outer(seq_len(nt), seq_len(nt), "-"))^1.5
    Rt <- 0.95 ^ (abs(outer(seq_len(nt), seq_len(nt), "-")) * 3)
    vt <- rep_len(c(0.30, 0.20, 0.40, 0.25, 0.35), nt)
    Gt <- Rt * tcrossprod(sqrt(vt))
    U  <- t(chol(K)) %*% matrix(rnorm(n.g * nt), n.g) %*% chol(Gt)
    L  <- matrix(rnorm(n.g * nt), n.g) * sqrt(0.05)
    ph <- expand.grid(id = ids, Rep = c("1", "2"), Trial = tl, stringsAsFactors = FALSE)
    ph$yield <- 4.5 + (U + L)[cbind(match(ph$id, ids), match(ph$Trial, tl))] +
        rnorm(nrow(ph), sd = sqrt(0.12))
    ph$id <- factor(ph$id); ph$Rep <- factor(ph$Rep); ph$Trial <- factor(ph$Trial)

    fit <- function(s) tryCatch(suppressWarnings(asreml(yield ~ Trial + Trial:Rep,
        random = as.formula(sprintf("~ %s(Trial):vm(id, covObj) + diag(Trial):id", s)),
        residual = ~ dsum(~units | Trial), data = ph, workspace = "1gb",
        trace = FALSE, maxit = 50)), error = function(e) { cat(s, "ERROR:", conditionMessage(e), "\n"); NULL })
    m.cg <- fit("corgh"); m.us <- fit("us")
    if (is.null(m.cg) || is.null(m.us)) next
    cat("converged: corgh =", m.cg$converge, "  us =", m.us$converge, "\n")

    vc <- summary(m.cg)$varcomp
    cat("\nASReml order of corgh .cor parameters (summary(model)$varcomp):\n")
    cn <- rownames(vc)[grepl("vm.*\\.cor$", rownames(vc))]
    pairs <- sub(".*!Trial!(T[0-9]+):!Trial!(T[0-9]+)\\.cor$", "(\\1,\\2)", cn)
    cat(" ", paste(gsub("T", "", pairs), collapse = " "), "\n")
    old.order <- unlist(lapply(seq_len(nt - 1L), function(c) sprintf("(%d,%d)", (c + 1L):nt, c)))
    cat("Order assumed by the OLD loop (column-wise):\n ", paste(old.order, collapse = " "), "\n")

    G.old <- old_fill(m.cg, nt)
    G.new <- unname(wgAim:::.extractVmGa(m.cg, "Trial", tl, 1))
    G.us  <- unname(wgAim:::.extractVmGa(m.us, "Trial", tl, 1))

    lt <- lower.tri(G.new)
    tab <- data.frame(pair   = sprintf("(%d,%d)", row(G.new)[lt], col(G.new)[lt]),
                      true   = round(Rt[lt], 3),
                      us.fit = round(cor_of(G.us)[lt], 3),
                      NEW    = round(cor_of(G.new)[lt], 3),
                      OLD    = round(cor_of(G.old)[lt], 3))
    tab$OLD.wrong <- ifelse(abs(tab$OLD - tab$NEW) > 1e-8, "<-- misplaced", "")
    cat("\nGenetic correlations of the vm term, lower triangle:\n")
    print(tab, row.names = FALSE)
    cat(sprintf("\nmax |Ga.new - Ga.us| = %.4f   max |Ga.old - Ga.us| = %.4f\n",
                max(abs(G.new - G.us)), max(abs(G.old - G.us))))
    cat("old Ga positive definite:", all(eigen(G.old, only.values = TRUE)$values > 0),
        "  new Ga positive definite:", all(eigen(G.new, only.values = TRUE)$values > 0), "\n")
}
