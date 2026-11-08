# =============================================================================
# research/aom_stopping/aom_null_diag.R
#
# STANDALONE DIAGNOSTIC -- not part of the package.  Nothing in R/ calls this.
#
# Null distribution of the AOM outlier statistic, as a candidate stopping rule.
# See research/aom_stopping/aom_stopping_rule.qmd for the mathematics.
#
# Under the AOM null for marker j (its effect is just another draw from the
# polygenic distribution), treating the variance parameters as known, the
# marker BLUPs are jointly Gaussian:
#
#     g~  ~ N(0, Vg),  Vg = kron(Ga, K) - PEV(g~)       (line level, trait-major)
#     q~_k = T %*% g~_k                                 (marker level, per trait)
#
#   vm path : T = cov.env$trans,  K = cov.env$relm   (package .constructCM)
#   mbf path: T = M^+ (Moore-Penrose), K = M M'       (M = current marker matrix)
#
# Statistics (as in .qtlSelect):
#   marker      : O_j = q~_j' Ginv q~_j / tr(Ginv V_j)
#                 (UV: q~_j^2 / v_j ~ chi^2_1; MV: sum lambda chi^2_1, mean 1)
#   chromosome  : O_c = sum_{j in c} q~_j' Ginv q~_j / sum_{j in c} tr(Ginv V_j)
#
# Returned (current iteration):
#   * observed statistics (to check against the package)
#   * per-marker p-value of the max (UV chi^2_1, MV Imhof), Bonferroni
#   * family-wise (max-statistic) null by simulation of g~* ~ N(0, Vg)
#   * the same for the chromosome statistic
#   * the standardised alternative S_j = q~_j' V_j^- q~_j ~ chi^2_r
#   * timings
# =============================================================================

# Imhof (1961): P(sum lambda_k chi^2_1 > x), lambda > 0
.imhof_upper <- function(x, lambda) {
    lambda <- lambda[lambda > 1e-12 * max(lambda)]
    if (length(lambda) == 1L) return(pchisq(x / lambda, 1, lower.tail = FALSE))
    f <- function(u) {
        if (u < 1e-12) return(0.5 * (sum(lambda) - x))   # limit u -> 0
        th  <- 0.5 * sum(atan(lambda * u)) - 0.5 * x * u
        rho <- prod((1 + lambda^2 * u^2)^0.25)
        sin(th) / (u * rho)
    }
    fv  <- Vectorize(f)
    # Oscillatory integrand: integrate half-period chunks of sin(x u / 2),
    # truncated where the envelope 1/(u prod(1+l^2u^2)^(1/4)) is negligible.
    env <- function(u) 1 / (u * prod((1 + lambda^2 * u^2)^0.25))
    U   <- 1; while (env(U) * U > 1e-10 && U < 1e8) U <- U * 2
    w   <- 2 * pi / max(x, 1e-3)
    br  <- seq(0, min(U, 4000 * w), by = w)
    val <- 0
    for (i in seq_len(length(br) - 1L))
        val <- val + tryCatch(integrate(fv, br[i], br[i + 1L], rel.tol = 1e-8,
                                        abs.tol = 1e-15, subdivisions = 200L)$value,
                              error = function(e) NA_real_)
    if (is.na(val)) return(NA_real_)
    min(max(0.5 + val / pi, 0), 1)
}

# Square-root factor L (N x r) with L L' = A (A symmetric PSD, possibly
# singular / with tiny negative eigenvalues from PEV round-off).
#   "chol" : pivoted Cholesky, truncated at the numerical rank (~10x faster
#            than eigen with reference BLAS; N = 4800: 14 s vs 138 s)
#   "eigen": eigen-decomposition, negative/negligible eigenvalues dropped
.sqrt_factor <- function(A, method = c("chol", "eigen")) {
    method <- match.arg(method)
    if (method == "chol") {
        R <- tryCatch(suppressWarnings(chol(A, pivot = TRUE,
                                            tol = 1e-10 * max(diag(A)))),
                      error = function(e) NULL)
        if (!is.null(R)) {
            r <- attr(R, "rank")
            return(t(R[seq_len(r), order(attr(R, "pivot")), drop = FALSE]))
        }
    }
    ev   <- eigen(A, symmetric = TRUE)
    keep <- ev$values > 1e-10 * max(ev$values)
    ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), sum(keep))
}

# -----------------------------------------------------------------------------
# .aom_inputs: Ga, Ginv, T, K, g~ (n x t), PEV(g~) for the current model,
# for the vm or mbf path, UV or MV.  Mirrors .qtlSelect() exactly.
# -----------------------------------------------------------------------------
.aom_inputs <- function(asm, phenoData, intervalObj, merge.by, Trait, ntrait) {
    ns <- asNamespace("wgAim")
    sigma2 <- asm$sigma2
    if (asm$vparameters.con[length(asm$vparameters.con)] == 4) sigma2 <- 1
    vp <- asm$vparameters
    trait.levels <- if (ntrait > 1L) levels(phenoData[[Trait]]) else "UV"

    Ga <- if (ntrait == 1L)
        matrix(vp[grep("vm.*covObj|mbf.*ints", names(vp))] * sigma2, 1L, 1L)
    else
        get(".extractVmGa", envir = ns)(asm, Trait, trait.levels, sigma2)

    cov.env <- attr(intervalObj, "env")
    if (!is.null(cov.env)) {
        # ---------------- vm path ----------------
        path  <- "vm"
        T     <- cov.env$trans; K <- cov.env$relm
        vmkey <- names(asm$G.param)[grep("vm.*covObj", names(asm$G.param))]
        if (ntrait == 1L) {
            pv <- predict(asm, classify = vmkey, only = vmkey, vcov = TRUE,
                          data = phenoData, maxit = 1)
            g  <- matrix(pv$pvals[, "predicted.value"], ncol = 1L)  # package order
            C  <- as.matrix(pv$vcov)
        } else {
            pv  <- predict(asm, classify = paste(Trait, merge.by, sep = ":"),
                           only = vmkey, vcov = TRUE, data = phenoData, maxit = 1)
            n.l <- nrow(pv$pvals) %/% ntrait
            ord <- order(pv$pvals[[Trait]], pv$pvals[[merge.by]])
            g   <- matrix(pv$pvals[ord, "predicted.value"], nrow = n.l, ncol = ntrait)
            C   <- as.matrix(pv$vcov)[ord, ord]
        }
    } else {
        # ---------------- mbf path ----------------
        path <- "mbf"
        M <- attr(intervalObj, "mbf.geno")
        if (is.null(M)) stop("aomNullDiag: no vm env and no mbf.geno on intervalObj.")
        key    <- asm$call$mbf$ints$key[1L]
        mbfkey <- names(asm$G.param)[grep("mbf.*ints", names(asm$G.param))]
        cls    <- if (ntrait == 1L) key else paste(Trait, key, sep = ":")
        pv <- predict(asm, classify = cls, only = mbfkey, vcov = TRUE,
                      data = phenoData, maxit = 1)
        if (is.null(pv$vcov)) stop("aomNullDiag: predict() returned no vcov (mbf).")
        lines.M <- rownames(M)
        pl   <- as.character(pv$pvals[[key]])
        keep <- which(pl %in% lines.M)                     # drops ".NoGeno"
        ord  <- if (ntrait == 1L) keep[match(lines.M, pl[keep])]
                else keep[order(match(as.character(pv$pvals[[Trait]][keep]), trait.levels),
                                match(pl[keep], lines.M))]
        if (length(ord) != length(lines.M) * ntrait || anyNA(ord))
            stop("aomNullDiag: mbf line predictions do not match the marker matrix.")
        g <- matrix(pv$pvals$predicted.value[ord], nrow = length(lines.M), ncol = ntrait)
        C <- as.matrix(pv$vcov)[ord, ord]
        T <- MASS::ginv(M); dimnames(T) <- list(colnames(M), rownames(M))
        K <- tcrossprod(M)
    }
    g[is.na(g)] <- 0; C[is.na(C)] <- 0
    if (nrow(g) != ncol(T)) stop("aomNullDiag: line dimension mismatch.")
    Ginv <- if (ntrait == 1L) matrix(1 / Ga[1L, 1L], 1L, 1L) else MASS::ginv(Ga)
    list(path = path, Ga = Ga, Ginv = Ginv, T = T, K = K, g = g, C = C,
         trait.levels = trait.levels)
}

# -----------------------------------------------------------------------------
# aomNullDiag
# -----------------------------------------------------------------------------
aomNullDiag <- function(asm, phenoData, intervalObj, state, merge.by,
                        Trait = NULL, ntrait = 1L, selection = "interval",
                        B = 2000L, alpha = 0.05, batch = 250L, seed = 1L,
                        factor.method = "chol") {
    tm  <- c(inputs = 0, blocks = 0, factor = 0, simulate = 0, pvals = 0)
    t0  <- proc.time()[["elapsed"]]
    inp <- .aom_inputs(asm, phenoData, intervalObj, merge.by, Trait, ntrait)
    tm["inputs"] <- proc.time()[["elapsed"]] - t0
    T <- inp$T; Ginv <- as.matrix(inp$Ginv); t <- ntrait
    m <- nrow(T); n <- ncol(T)
    Vg  <- kronecker(inp$Ga, inp$K) - inp$C
    blk <- function(k) ((k - 1L) * n + 1L):(k * n)

    # ---- per-marker BLUP variance blocks V_j (t x t) -----------------------
    t0 <- proc.time()[["elapsed"]]
    V <- array(0, c(m, t, t))
    for (k in seq_len(t)) for (l in k:t) {
        v <- rowSums((T %*% Vg[blk(k), blk(l)]) * T)
        V[, k, l] <- V[, l, k] <- v
    }
    den <- if (t == 1L) V[, 1, 1] * Ginv[1, 1] else
        apply(V, 1L, function(Vj) sum(Ginv * Vj))
    den <- pmax(den, 0)                       # package floors vqtilde at 0
    tm["blocks"] <- proc.time()[["elapsed"]] - t0

    # ---- observed statistics ------------------------------------------------
    Q    <- T %*% inp$g                                       # m x t
    num  <- rowSums((Q %*% Ginv) * Q)
    oint <- ifelse(den > 0, num / den, 0)
    names(oint) <- rownames(T)
    jmax <- which.max(oint)
    chr  <- sapply(strsplit(rownames(T), "\\."), "[", 2L)
    ochr <- tapply(num, chr, sum) / tapply(den, chr, sum)
    ochr[!is.finite(ochr)] <- 0
    cmax <- names(ochr)[which.max(ochr)]

    # standardised alternative for the selected marker
    Vj   <- matrix(V[jmax, , ], t, t)
    maha <- sum(Q[jmax, ] * (MASS::ginv(Vj) %*% Q[jmax, ]))
    ev   <- eigen(Vj, symmetric = TRUE, only.values = TRUE)$values
    rnk  <- sum(ev > 1e-8 * max(ev))

    # ---- per-marker null p-value (UV chi^2_1, MV Imhof) ---------------------
    t0 <- proc.time()[["elapsed"]]
    Gh <- with(eigen(Ginv, symmetric = TRUE),
               vectors %*% diag(sqrt(pmax(values, 0)), t) %*% t(vectors))
    lam_j <- function(j) {
        l <- eigen(Gh %*% matrix(V[j, , ], t, t) %*% Gh, symmetric = TRUE,
                   only.values = TRUE)$values
        pmax(l, 0) / den[j]
    }
    pmark <- function(j, x) {
        if (den[j] <= 0) return(1)
        if (t == 1L) return(pchisq(x, 1, lower.tail = FALSE))
        .imhof_upper(x, lam_j(j))
    }
    p.max <- pmark(jmax, oint[jmax])
    tm["pvals"] <- proc.time()[["elapsed"]] - t0

    # ---- family-wise null by simulation -------------------------------------
    t0 <- proc.time()[["elapsed"]]
    L  <- .sqrt_factor(Vg, factor.method)
    tm["factor"] <- proc.time()[["elapsed"]] - t0
    t0 <- proc.time()[["elapsed"]]
    set.seed(seed)
    maxnull <- maxchr <- numeric(0); nullsum <- 0
    chr.f <- factor(chr)
    den.c <- tapply(den, chr.f, sum)
    while (length(maxnull) < B) {
        bb <- min(batch, B - length(maxnull))
        Z  <- L %*% matrix(rnorm(ncol(L) * bb), ncol(L), bb)
        Qs <- lapply(seq_len(t), function(k) T %*% Z[blk(k), , drop = FALSE])
        numS <- matrix(0, m, bb)
        for (k in seq_len(t)) for (l in seq_len(t))
            numS <- numS + Ginv[k, l] * Qs[[k]] * Qs[[l]]
        st <- numS / ifelse(den > 0, den, Inf)
        maxnull <- c(maxnull, apply(st, 2L, max))
        nullsum <- nullsum + sum(st[den > 0, ])
        ch <- rowsum(numS, chr.f) / as.vector(den.c)
        ch[!is.finite(ch)] <- 0
        maxchr <- c(maxchr, apply(ch, 2L, max))
    }
    tm["simulate"] <- proc.time()[["elapsed"]] - t0

    thr   <- unname(quantile(maxnull, 1 - alpha))
    p.fwe <- (1 + sum(maxnull >= oint[jmax])) / (1 + B)
    thr.c <- unname(quantile(maxchr, 1 - alpha))
    p.fwe.c <- (1 + sum(maxchr >= ochr[cmax])) / (1 + B)
    p.star <- pmark(jmax, thr)
    Meff   <- log(1 - alpha) / log(1 - p.star)

    list(summary = data.frame(
             path = inp$path, ntrait = t, n.markers = m, n.lines = n,
             marker = names(oint)[jmax], max.oint = oint[jmax],
             p.marker = p.max, p.bonf = min(1, p.max * m),
             fwe.thr = thr, p.fwe = p.fwe, Meff = Meff,
             maha.max = maha, maha.df = rnk,
             p.maha.bonf = min(1, m * pchisq(maha, rnk, lower.tail = FALSE)),
             null.mean.oint = nullsum / (sum(den > 0) * B),
             chr.max = cmax, max.ochr = unname(ochr[cmax]),
             fwe.thr.chr = thr.c, p.fwe.chr = p.fwe.c,
             stringsAsFactors = FALSE),
         oint = oint, ochr = ochr, maxnull = maxnull, maxchr = maxchr,
         Q = Q, den = den, Ga = inp$Ga,
         timing = tm)
}
