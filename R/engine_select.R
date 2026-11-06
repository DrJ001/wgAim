# =============================================================================
# engine_select.R
# Internal functions for interval/marker selection and significance testing.
# =============================================================================

#' @keywords internal
# Trait    : character name of the trait factor column (NULL = univariate)
# ntrait   : integer, length(levels(phenoData[[Trait]])); 1 when Trait = NULL
.qtlSelect <- function(asm, phenoData, intervalObj, gen.type, selection,
                       exclusion.window, state, verboseLev, merge.by = NULL,
                       Trait = NULL, ntrait = 1L) {
    sigma2 <- asm$sigma2
    if (asm$vparameters.con[length(asm$vparameters.con)] == 4)
        sigma2 <- 1

    # ------------------------------------------------------------------
    # Compute Ga (ntrait x ntrait genetic covariance matrix) and Ginv.
    # For ntrait == 1: Ga = avar (scalar wrapped in 1x1 matrix).
    # For ntrait > 1: Ga extracted from vparameters based on variance
    #   structure (diag / corh / fa(k)).
    # ------------------------------------------------------------------
    sterms <- "vm.*covObj|mbf.*ints"
    vp     <- asm$vparameters

    # Trait levels used for Ga extraction and column names of scaled.mat
    trait.levels <- if (!is.null(Trait)) levels(phenoData[[Trait]]) else NULL

    if (ntrait == 1L) {
        avar <- vp[grep(sterms, names(vp))] * sigma2
        Ga   <- matrix(avar, 1L, 1L)
    } else {
        # Structure-agnostic extraction keyed on vparameter names
        # (diag / corh / corgh / us / fa(k)).  When fa specific variances hit
        # the boundary (zero), Ga is rank-deficient; MASS::ginv() then gives
        # the correct Moore-Penrose generalised inverse Ga^-.
        Ga   <- .extractVmGa(asm, Trait, trait.levels, sigma2)
        avar <- NA_real_   # not used directly when ntrait > 1
    }
    Ginv <- if (ntrait == 1L) matrix(1 / Ga[1L, 1L], 1L, 1L) else MASS::ginv(Ga)

    # scaled.mat: nmarkers x ntrait signed z-scores; populated inside vm/mbf branches
    scaled.mat <- NULL

    # ------------------------------------------------------------------
    # Predict conditional BLUPs (atilde) and compute posterior variance
    # ------------------------------------------------------------------
    cat(" Predict step for outlier statistics \n")
    cat("=====================================\n")

    if (!is.null(cov.env <- attr(intervalObj, "env"))) {
        # vm path
        # vmterm.only: the G.param key of the vm term, used for predict(only=)
        # in both UV and MV paths.  ASReml strips non-fa variance-model
        # wrappers from term names, e.g.
        #   corgh(Trait):vm(id, covObj)  -> "Trait:vm(id, covObj)"
        #   diag/corh/us(Trait):vm(...)  -> "Trait:vm(id, covObj)"
        #   fa(Trait, k):vm(id, covObj)  -> "fa(Trait, k):vm(id, covObj)"
        # Passing the raw formula string for a non-fa structure does NOT match
        # any model term: predict() then silently ignores `only` (returning
        # full fitted values incl. fixed effects) or, on some ASReml builds,
        # returns no vcov -- both of which break the outlier statistics.
        vmterm.only <- names(asm$G.param)[grep("vm.*covObj", names(asm$G.param))]

        if (ntrait == 1L) {
            pv      <- predict(asm, classify = vmterm.only,
                               only = vmterm.only, vcov = TRUE, data = phenoData, maxit = 1)
            atilde  <- pv$pvals[, "predicted.value"]
            qtilde  <- as.vector(cov.env$trans %*% atilde)
            vatilde <- avar * cov.env$relm - as.matrix(pv$vcov)
            qhalf   <- cov.env$trans %*% vatilde
            vqtilde <- colSums(t(qhalf) * t(cov.env$trans))
        } else {
            classify.term <- paste(Trait, merge.by, sep = ":")
            pv   <- predict(asm, classify = classify.term,
                            only = vmterm.only, vcov = TRUE, data = phenoData, maxit = 1)
            if (is.null(pv$vcov))
                stop("predict() returned no variance matrix for the genomic term '",
                     vmterm.only, "'. Cannot compute outlier statistics.")
            # Sort Trait-major: all lines for Trait 1 first, then Trait 2, ...
            # This matches kronecker(Ga, relm) which is block-structured as
            # Ga[i,j] * relm in Trait-major order.
            # byrow=FALSE then correctly fills atilde: column k = all lines for trait k.
            nlines <- nrow(pv$pvals) %/% ntrait
            ord    <- order(pv$pvals[[Trait]], pv$pvals[[merge.by]])
            atilde <- matrix(pv$pvals[ord, "predicted.value"],
                             nrow = nlines, ncol = ntrait, byrow = FALSE)
            atilde[is.na(atilde)] <- 0
            pev    <- as.matrix(pv$vcov)[ord, ord]
            pev[is.na(pev)] <- 0
            # Back-transform: qtilde.mat is nmarkers x ntrait
            qtilde.mat <- cov.env$trans %*% atilde
            vatilde    <- kronecker(Ga, cov.env$relm) - pev

            # Use the Rcpp implementation for speed and correctness.
            vqtilde    <- compute_vqtilde(cov.env$trans, Ginv,
                                          as.matrix(vatilde), ntrait)
            # Scalar outlier stat: t(qtilde_i) %*% Ginv %*% qtilde_i
            qtilde <- apply(qtilde.mat, 1L, function(q) sum(q * (Ginv %*% q)))
            # Signed scaled BLUPs (nmarkers x ntrait): per-trial z-scores.
            # For trial j, posterior variance = diag(trans %*% vatilde_jj %*% t(trans))
            # where vatilde_jj is the j-th nlines x nlines diagonal block of vatilde.
            scaled.mat <- matrix(0, nrow = nrow(qtilde.mat), ncol = ntrait,
                                 dimnames = list(rownames(qtilde.mat), trait.levels))
            for (j in seq_len(ntrait)) {
                rows_j  <- ((j - 1L) * nlines + 1L):(j * nlines)
                v_j     <- rowSums((cov.env$trans %*% vatilde[rows_j, rows_j]) *
                                    cov.env$trans)
                v_j     <- pmax(v_j, 0)
                scaled.mat[, j] <- ifelse(v_j > 0, qtilde.mat[, j] / sqrt(v_j), 0)
            }
        }
    } else {
        # mbf path
        if (ntrait == 1L) {
            avar    <- vp[grep("mbf.*ints", names(vp))] * sigma2
            mbf     <- grep("mbf", rownames(asm$coefficients$random))
            qtilde  <- asm$coefficients$random[mbf, 1L]
            pevar   <- sigma2 * asm$vcoeff$random[mbf]
            vqtilde <- avar - pevar
        } else {
            # mbf multivariate.
            #
            # The denominator needs, per marker, the FULL ntrait x ntrait
            # posterior block Ga - PEV(q_i).  ASReml's vcoeff only gives the
            # PEV diagonal; combining it with the full Ga ignores the large
            # cross-trait PEV covariances, drives every vqtilde <= 0 and zeroes
            # all outlier statistics (forward selection then just takes the
            # first marker left on Chr1).
            #
            # Exact route, mirroring the vm path: predict the line-level values
            # u = M q of the mbf term per trait with a full vcov, then map to
            # marker space with trans = M^+ (Moore-Penrose).  q lies in the row
            # space of M, so q = M^+ u exactly for any lines/markers ratio
            # (verified: max |M^+ u - q_BLUP| ~ 1e-5).  The prior of u is
            # Ga (x) MM' with Ga the per-marker covariance.
            M <- attr(intervalObj, "mbf.geno")
            if (is.null(M))
                stop("Internal error: marker matrix for the mbf path not found.")
            mbf.key  <- asm$call$mbf$ints$key[1L]
            mbf.only <- names(asm$G.param)[grep("mbf.*ints", names(asm$G.param))]
            pv <- predict(asm, classify = paste(Trait, mbf.key, sep = ":"),
                          only = mbf.only, vcov = TRUE, data = phenoData, maxit = 1)
            if (is.null(pv$vcov))
                stop("predict() returned no variance matrix for the genomic term '",
                     mbf.only, "'. Cannot compute outlier statistics.")
            lines.M <- rownames(M)
            pl      <- as.character(pv$pvals[[mbf.key]])
            keep    <- which(pl %in% lines.M)       # drops ".NoGeno"
            ord     <- keep[order(match(as.character(pv$pvals[[Trait]][keep]), trait.levels),
                                  match(pl[keep], lines.M))]
            nlines  <- length(lines.M)
            if (length(ord) != nlines * ntrait)
                stop("Internal error: mbf line predictions do not match the marker matrix.")
            atilde  <- matrix(pv$pvals$predicted.value[ord], nrow = nlines, ncol = ntrait)
            atilde[is.na(atilde)] <- 0
            pev     <- as.matrix(pv$vcov)[ord, ord]
            pev[is.na(pev)] <- 0
            trans   <- MASS::ginv(M)                 # nmarkers x nlines
            dimnames(trans) <- list(colnames(M), rownames(M))
            relm    <- tcrossprod(M)                 # nlines x nlines
            qtilde.mat <- trans %*% atilde
            vatilde    <- kronecker(Ga, relm) - pev
            vqtilde    <- compute_vqtilde(trans, Ginv, as.matrix(vatilde), ntrait)
            qtilde     <- apply(qtilde.mat, 1L, function(q) sum(q * (Ginv %*% q)))
            scaled.mat <- matrix(0, nrow = nrow(qtilde.mat), ncol = ntrait,
                                 dimnames = list(NULL, trait.levels))
            for (j in seq_len(ntrait)) {
                rows_j <- ((j - 1L) * nlines + 1L):(j * nlines)
                v_j    <- pmax(rowSums((trans %*% vatilde[rows_j, rows_j]) * trans), 0)
                scaled.mat[, j] <- ifelse(v_j > 0, qtilde.mat[, j] / sqrt(v_j), 0)
            }
        }
    }

    # ------------------------------------------------------------------
    # Outlier statistics (scalar per interval/marker regardless of ntrait)
    # ------------------------------------------------------------------
    # vqtilde is the posterior variance diagonal: prior - PEV.  At late
    # iterations, after most markers are excluded, the remaining covObj can
    # be nearly rank-deficient and PEV can numerically exceed the prior,
    # yielding negative vqtilde values.  Floor at zero so that oint is
    # always non-negative (a negative outlier stat is meaningless).
    gnams <- names(state)[as.logical(state)]
    vqtilde <- pmax(vqtilde, 0)
    names(qtilde) <- names(vqtilde) <- gnams
    # Outlier statistic: two forms depending on ntrait.
    #
    # ntrait == 1: qtilde is a signed scalar BLUP; squaring gives a
    #   chi-squared-scale Wald statistic: oint = qtilde^2 / vqtilde.
    #
    # ntrait > 1: qtilde = t(q_i) %*% Ginv %*% q_i is already a quadratic
    #   form (>= 0, chi-squared-scale); vqtilde = tr(Ginv Ti vatilde Ti') is
    #   its variance.  No further squaring: oint = qtilde / vqtilde.
    oint <- if (ntrait == 1L)
        ifelse(vqtilde > 0, qtilde^2 / vqtilde, 0)
    else
        ifelse(vqtilde > 0, qtilde  / vqtilde, 0)
    names(oint) <- gnams
    ochr <- NULL
    if (selection == "chromosome") {
        chr.names <- names(intervalObj$geno)
        nochr     <- length(chr.names)
        allc      <- sapply(strsplit(gnams, "\\."), "[", 2)
        ochr      <- numeric(nochr)
        for (c in seq_len(nochr)) {
            whc     <- allc %in% chr.names[c]
            if (ntrait == 1L) {
                nums <- qtilde[whc]^2
            } else {
                nums <- qtilde[whc]   # already quadratic forms
            }
            dens    <- vqtilde[whc]
            ochr[c] <- ifelse(!is.na(sum(nums) / sum(dens)), sum(nums) / sum(dens), 0)
        }
        names(ochr) <- chr.names
        mchr <- chr.names[ochr == max(ochr)]
        cint <- allc %in% mchr
        chri <- oint[cint]
        mint <- (1:length(chri))[chri == max(chri)]
        qtl  <- names(chri)[mint]
        if (verboseLev > 0) {
            cgen <- ifelse(gen.type == "marker", "Marker", "Interval")
            cat("\n Selection of chromosome using the AOM statistic\n")
            cat("=============================================== \n")
            for (i in seq_len(nochr))
                cat(" Chromosome ", chr.names[i], "Outlier Statistic ", ochr[i], "\n")
            cat("============================================= \n\n")
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for (i in seq_along(chri))
                cat(cgen, names(chri)[i], "Outlier Statistic ", chri[i], "\n")
            cat("=============================================== \n\n")
        }
    } else {
        qtl <- names(oint)[oint == max(oint)]
        qsp <- unlist(strsplit(qtl, split = "\\."))
        mint <- as.numeric(qsp[3])
        mchr <- qsp[2]
        if (verboseLev > 0) {
            cgen <- ifelse(gen.type == "marker", "Marker", "Interval")
            cat(cgen, "outlier statistics \n")
            cat("=============================================== \n")
            for (i in seq_along(oint))
                cat(cgen, names(oint)[i], "Outlier Statistic ", oint[i], "\n")
            cat("=============================================== \n\n")
        }
    }

    # Take first if tied; fill out full state vectors; apply exclusion window
    qtl  <- qtl[1]
    tint <- state
    tint[as.logical(state)] <- oint
    # blups: signed z-scores for BLUP/contrast plots.
    #   ntrait == 1: named vector (length = all markers); zero where vqtilde == 0.
    #   ntrait >  1: named matrix (all markers x ntrait); zero rows for excluded markers.
    if (ntrait == 1L) {
        blups <- state
        blups[as.logical(state)] <- ifelse(
            vqtilde > 0,
            qtilde / sqrt(vqtilde),
            0
        )
    } else {
        rownames(scaled.mat) <- gnams
        blups <- matrix(0, nrow = length(state), ncol = ntrait,
                        dimnames = list(names(state), trait.levels))
        blups[as.logical(state), ] <- scaled.mat
    }
    oint <- tint
    schr  <- sapply(strsplit(names(state), "\\."), "[", 2)
    wnams <- names(state)[schr %in% mchr]
    inums <- as.numeric(sapply(strsplit(wnams, "\\."), "[", 3))
    if (gen.type == "marker")
        dists <- intervalObj$geno[[mchr]]$map
    else
        dists <- intervalObj$geno[[mchr]]$inferred.map
    dists <- dists[inums]
    exc   <- wnams[abs(dists - dists[mint]) <= exclusion.window]
    state[exc] <- 0
    list(state = state, qtl = qtl, ochr = ochr, oint = oint, blups = blups)
}

#' @keywords internal
# Boundary LRT comparing qtlModel (with composite vm/mbf term) against
# baseModel (same fixed QTL effects and residual structure, no vm/mbf term).
#
# ntrait == 1 : one-sided test of a scalar variance parameter at its lower
#   boundary.  Null: 1/2*chi^2(0) + 1/2*chi^2(1).
#
# ntrait > 1  : boundary test of the ntrait x ntrait genetic covariance matrix
#   Ga of the vm term being zero.  Null: pchisq.mixture(stat, ntrait).
#
#   Justification: under H0 (Ga = 0), only the ntrait diagonal variance
#   parameters are at their boundary (>= 0).  Correlation and factor-loading
#   parameters are non-identifiable at Ga = 0 and do not contribute degrees
#   of freedom to the null distribution (Self & Liang 1987, non-regular case).
#   Consequently pchisq.mixture(stat, ntrait) is the exact null distribution
#   for diag(Trait):vm() and an asymptotically valid approximation for
#   corgh(Trait):vm() and fa(Trait,k):vm().
#
#   For ntrait > 1 the caller is responsible for passing a qtlModel whose
#   vm/mbf term uses diag(Trait): (built via .buildDiagVmModel()).  This
#   removes correlation/loading parameters that are non-identifiable under H0
#   asymptotically but weakly identifiable in finite samples, which would
#   otherwise inflate the LRT statistic and cause anti-conservative tests.
#   baseModel (no vm term, same residual structure and QTL fixed effects)
#   is the correct null and requires no modification.
.lrtTest <- function(qtlModel, baseModel, TypeI, ntrait = 1L) {
    baseLogL <- baseModel$loglik
    stat     <- 2 * (qtlModel$loglik - baseLogL)
    if (ntrait == 1L) {
        pvalue <- (1 - pchisq(stat, 1)) / 2
        pass   <- stat >= qchisq(1 - 2 * TypeI, 1)
    } else {
        pvalue <- 1 - pchisq.mixture(stat, ntrait = ntrait)
        pass   <- pvalue <= TypeI
    }
    list(baseLogL = baseLogL, stat = stat, pvalue = pvalue, pass = pass)
}

# =============================================================================
# .mbfEffectRows
#
# Row indices into model$coefficients$random of the per-level marker effects
# of a multivariate mbf('ints') term, as an nmarkers x ntrait integer matrix
# (column j = trait.levels[j], rows in the order of `markers`).
#
# ASReml-R 4.2 names these "<prefix>_<level>:mbf(\"ints\")_<marker>", where
# <prefix> is the Trait name for diag/corh/corgh/us and "fa(Trait, k)" for
# fa(Trait, k).  fa terms also carry "fa(Trait, k)_Comp<k>" factor-score
# rows, which match no level and are excluded.  Rows are matched by level and
# marker name, so no ordering is assumed.
# =============================================================================
#' @keywords internal
.mbfEffectRows <- function(rn, Trait, trait.levels, markers) {
    rows   <- grep("mbf", rn)
    nms    <- rn[rows]
    pre    <- sub(":mbf.*$", "", nms)
    tr.pfx <- paste0(Trait, "_")
    lev    <- ifelse(grepl("^fa\\(", pre),
                     sub("^fa\\([^)]*\\)_", "", pre),
                     ifelse(startsWith(pre, tr.pfx),
                            substring(pre, nchar(tr.pfx) + 1L), NA_character_))
    mark   <- sub("^.*mbf\\([^)]*\\)_", "", nms)
    idx <- vapply(trait.levels, function(l) {
        s <- which(lev == l)
        rows[s][match(markers, mark[s])]
    }, integer(length(markers)))
    idx <- matrix(idx, nrow = length(markers), dimnames = list(markers, trait.levels))
    if (anyNA(idx))
        stop("Could not match the mbf marker effects of the fitted model to the ",
             "genotype matrix for ", Trait, " level(s): ",
             paste(trait.levels[colSums(is.na(idx)) > 0L], collapse = ", "), ".")
    idx
}

# =============================================================================
# .extractVmGa
#
# Build the ntrait x ntrait genetic covariance matrix Ga of the multivariate
# genomic term (vm(...covObj) or mbf('ints')) from a fitted ASReml model.
#
# The structure is identified from the vparameter NAMES rather than the
# formula text, because ASReml strips non-fa wrappers from term names
# (corgh(Trait):vm(id, covObj) has G.param key "Trait:vm(id, covObj)").
# Each parameter is placed by exact name, so no ordering is assumed.
# ASReml 4.2 naming (suffix after the term key), for Trait = "Trial":
#   diag  : !Trial_T1                         variances
#   corh  : !Trial!cor + !Trial_T1 ...         one correlation + variances
#   corgh : !Trial!T2:!Trial!T1.cor + !Trial_T1 ...
#   us    : !Trial_T1:T1, !Trial_T2:T1, ...    lower triangle incl. diagonal
#   fa(k) : !T1!var (Psi), !T1!fa1 ... !T1!fak (Lambda)
#
# Variances / covariances are ratios to sigma2 and are scaled by sigma2;
# correlations and fa loadings are used as returned (consistent with the
# previous implementation).  Negative variances are floored at zero.
#
# Returns a numeric matrix with dimnames = list(trait.levels, trait.levels).
# =============================================================================
#' @keywords internal
# `term` : regular expression selecting the vparameters of the term to read.
#   Default = the genomic vm/mbf term.  For the residual genetic term use
#   paste0(":", gterm, "!") -- matches "Trial:Gsave!..." / "fa(Trial, 1):Gsave!..."
#   but not "Trial:vm(Gsave, covObj)!...".
.extractVmGa <- function(model, Trait, trait.levels, sigma2 = 1,
                         term = "vm.*covObj|mbf.*ints") {
    ntrait <- length(trait.levels)
    vp     <- model$vparameters
    sel    <- grep(term, names(vp))
    if (!length(sel))
        stop("No genomic (vm/mbf) variance parameters found in the model.")
    vals <- as.numeric(vp[sel])
    # Strip the term key (everything before the first "!") to get the suffix
    sfx  <- sub("^[^!]*", "", names(vp)[sel])
    pget <- function(s) vals[match(s, sfx)]
    lev  <- trait.levels

    fa.var <- paste0("!", lev, "!var")
    us.dg  <- paste0("!", Trait, "_", lev, ":", lev)
    v.nm   <- paste0("!", Trait, "_", lev)

    if (all(fa.var %in% sfx)) {
        # ---- fa(Trait, k) ------------------------------------------------
        k   <- max(as.integer(sub(".*!fa([0-9]+)$", "\\1",
                                  sfx[grepl("!fa[0-9]+$", sfx)])))
        psi <- pmax(pget(fa.var) * sigma2, 0)
        Lam <- vapply(seq_len(k), function(j) pget(paste0("!", lev, "!fa", j)),
                      numeric(ntrait))
        Lam <- matrix(Lam, nrow = ntrait)
        Lam[is.na(Lam)] <- 0
        Ga  <- Lam %*% t(Lam) + diag(psi, ntrait)
    } else if (all(us.dg %in% sfx)) {
        # ---- us(Trait) ---------------------------------------------------
        Ga <- matrix(0, ntrait, ntrait)
        for (i in seq_len(ntrait)) for (j in seq_len(i)) {
            Ga[i, j] <- Ga[j, i] <-
                pget(paste0("!", Trait, "_", lev[i], ":", lev[j])) * sigma2
        }
        diag(Ga) <- pmax(diag(Ga), 0)
    } else if (all(v.nm %in% sfx)) {
        # ---- diag / corh / corgh -----------------------------------------
        var.pars <- pmax(pget(v.nm) * sigma2, 0)
        sds      <- sqrt(var.pars)
        R        <- diag(ntrait)
        corh.nm  <- paste0("!", Trait, "!cor")
        if (corh.nm %in% sfx) {
            R[] <- pget(corh.nm)
            diag(R) <- 1
        } else {
            for (i in seq_len(ntrait)) for (j in seq_len(i - 1L)) {
                r <- pget(paste0("!", Trait, "!", lev[i], ":!", Trait, "!",
                                lev[j], ".cor"))
                if (!is.na(r)) R[i, j] <- R[j, i] <- r
            }
        }
        Ga <- R * tcrossprod(sds)
    } else {
        stop("Unrecognised variance structure for the genomic term. ",
             "Supported: diag, corh, corgh, us, fa(k). Parameters found: ",
             paste(names(vp)[sel], collapse = ", "))
    }
    dimnames(Ga) <- list(trait.levels, trait.levels)
    Ga
}

# =============================================================================
# .mvPercVarComponents
#
# Per-trial variance components for the multivariate Perc.Var denominator in
# summary.qtlAim() / summary.gwasAim():
#   var.mark : diag(Ga) of the genomic term / rel.scale   (length ntrait)
#   var.res  : diag(G)  of the residual genetic term      (length ntrait)
#
# Previously the summaries averaged ALL vparameters of these terms, which
# includes correlations (corh/corgh), covariances (us) and loadings (fa), so
# Perc.Var depended on the parameterisation of the same model.  Returns NULL
# if a structure cannot be identified (caller keeps its fallback).
# =============================================================================
#' @keywords internal
.mvPercVarComponents <- function(object, gterm, scale, sigma2) {
    Trait <- object$QTL$Trait
    lev   <- object$QTL$trait.levels
    if (is.null(Trait) || is.null(lev)) return(NULL)
    mb.esc <- gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", gterm)
    tryCatch(list(
        var.mark = diag(.extractVmGa(object, Trait, lev, sigma2)) / scale,
        var.res  = diag(.extractVmGa(object, Trait, lev, sigma2,
                                     term = paste0(":", mb.esc, "!")))
    ), error = function(e) NULL)
}
