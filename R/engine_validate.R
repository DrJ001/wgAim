# =============================================================================
# engine_validate.R
# Internal shared validation for wgAim analyses.
# =============================================================================

#' @keywords internal
# .resolveCallFormulae
#
# When a base model is fitted with a formula held in a variable, e.g.
#   rf <- ~ us(Trial):id;  asreml(y ~ Trial, random = rf, ...)
# the stored call contains the SYMBOL `rf`, not the formula.  Downstream code
# parses model$call$random / $fixed / $residual directly (.rhs_terms(),
# terms.formula(), update(fixed. = . ~ ...)), which fails with "object of
# type 'symbol' is not subsettable" -- and the symbol may not even resolve
# from inside wgAim.  ASReml keeps the evaluated formulae in model$formulae,
# so substitute those into the call for any component that is not a literal
# `~` formula.  Verified with ASReml-R 4.2: update() on the substituted model
# reproduces the original fit exactly.
.resolveCallFormulae <- function(model) {
    for (nm in c("fixed", "random", "residual", "sparse")) {
        el <- model$call[[nm]]
        if (is.null(el) || inherits(el, "formula")) next
        if (is.call(el) && identical(el[[1L]], as.name("~"))) next
        f <- model$formulae[[nm]]
        if (is.null(f))
            stop("Cannot resolve the '", nm, "' formula of the base model. ",
                 "Refit it with the formula written directly in the asreml() call.")
        # Strip the terms attributes: keep a plain formula, as when written
        # literally in the call.
        env <- environment(f)
        if (is.null(env)) env <- globalenv()
        f <- stats::as.formula(paste(deparse(f, width.cutoff = 500L), collapse = " "),
                               env = env)
        model$call[[nm]] <- f
    }
    model
}

#' @keywords internal
.validateModel <- function(baseModel, merge.by, method, selection, breakout) {
    baseModel <- .resolveCallFormulae(baseModel)
    if (!baseModel$converge) {
        cat("Warning: Base model has not converged. Updating base model\n")
        baseModel <- update(baseModel)
        if (!baseModel$converge)
            stop("Base model not converged: Check base model before proceeding.")
    }
    asremlEnv <- lapply(baseModel$formulae, function(el) attr(el, ".Environment"))
    phenoData <- eval(baseModel$call$data)
    if (is.null(merge.by))
        stop("merge.by: name of the column linking phenotypic and genotypic data is required.")
    if (!(method %in% c("fixed", "random")))
        stop("method must be either \"fixed\" or \"random\".")
    if (!(selection %in% c("interval", "chromosome")))
        stop("selection must be either \"interval\" or \"chromosome\".")
    if (!is.numeric(breakout) | breakout < -1 | breakout == 0)
        stop("breakout must be -1 or a positive integer.")
    list(baseModel = baseModel, asremlEnv = asremlEnv, phenoData = phenoData)
}
