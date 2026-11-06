# =============================================================================
# test-engine-helpers.R
# Pure helpers (no ASReml needed):
#   .resolveCallFormulae()  -- formulae passed to asreml() as variables
#   .mbfKey() / .mbfCovObj() -- mbf key + covariate table (NA-key workaround)
#   .mbfEffectRows()        -- MV mbf coefficient matching (incl. fa naming)
# =============================================================================

.mock_call_model <- function(call, formulae) {
    m <- list(call = call, formulae = formulae)
    class(m) <- c("asreml", "list")
    m
}

test_that(".resolveCallFormulae substitutes formulae for symbols in the call", {
    rf <- ~ us(Trial):id
    fx <- yield ~ Trial
    cl <- quote(asreml(fixed = fx, random = rf, residual = ~ units, data = dat))
    m  <- .mock_call_model(cl, list(fixed = fx, random = rf, residual = ~ units))
    out <- .resolveCallFormulae(m)
    expect_true(inherits(out$call$random, "formula"))
    expect_true(inherits(out$call$fixed, "formula"))
    expect_equal(.rhs_terms(out$call$random), "us(Trial):id")
    expect_equal(deparse(out$call$fixed), "yield ~ Trial")
    # literal residual formula left untouched
    expect_identical(out$call$residual, cl$residual)
})

test_that(".resolveCallFormulae leaves literal formulae unchanged", {
    cl <- quote(asreml(fixed = yield ~ Trial, random = ~ diag(Trial):id, data = dat))
    m  <- .mock_call_model(cl, list())
    out <- .resolveCallFormulae(m)
    expect_identical(out$call, cl)
})

test_that(".resolveCallFormulae errors clearly when a formula cannot be resolved", {
    cl <- quote(asreml(fixed = yield ~ 1, random = rf, data = dat))
    m  <- .mock_call_model(cl, list())
    expect_error(.resolveCallFormulae(m), "Cannot resolve the 'random' formula")
})

test_that(".mbfKey uses Gmbf only after fixLines (merge.by == 'Gsave')", {
    expect_identical(.mbfKey("Gsave"), "Gmbf")
    expect_identical(.mbfKey("id"), "id")
})

test_that(".mbfCovObj adds an all-zero .NoGeno row for the Gmbf key", {
    G <- matrix(c(1, -1, 0, 1, 1, -1), nrow = 3,
                dimnames = list(c("L1", "L2", "L3"), c("Chr.1.1", "Chr.1.2")))
    co <- .mbfCovObj(G, "Gmbf")
    expect_equal(names(co), c("Gmbf", "Chr.1.1", "Chr.1.2"))
    expect_equal(co$Gmbf, c("L1", "L2", "L3", ".NoGeno"))
    expect_equal(unname(unlist(co[4L, -1L])), c(0, 0))
    expect_equal(unname(as.matrix(co[1:3, -1L])), unname(G))
    co2 <- .mbfCovObj(G, "id")
    expect_equal(nrow(co2), 3L)
    expect_equal(names(co2)[1L], "id")
})

test_that(".mbfEffectRows matches diag/us naming by level and marker", {
    mk <- c("Chr.1.1", "Chr.1.2", "Chr.1.3")
    rn <- c("Trial_T1:Gsave_L1",
            paste0("Trial_T1:mbf(\"ints\")_", mk),
            paste0("Trial_T2:mbf(\"ints\")_", rev(mk)))   # scrambled within T2
    idx <- .mbfEffectRows(rn, "Trial", c("T1", "T2"), mk)
    expect_equal(dim(idx), c(3L, 2L))
    expect_equal(rn[idx[, "T1"]], paste0("Trial_T1:mbf(\"ints\")_", mk))
    expect_equal(rn[idx[, "T2"]], paste0("Trial_T2:mbf(\"ints\")_", mk))
})

test_that(".mbfEffectRows handles fa naming and excludes Comp rows", {
    mk <- c("Chr.1.1", "Chr.1.2")
    rn <- c(paste0("fa(Trial, 1)_Comp1:mbf(\"ints\")_", mk),
            paste0("fa(Trial, 1)_T1:mbf(\"ints\")_", mk),
            paste0("fa(Trial, 1)_T2:mbf(\"ints\")_", mk))
    idx <- .mbfEffectRows(rn, "Trial", c("T1", "T2"), mk)
    expect_equal(as.vector(idx), 3:6)
})

test_that(".mbfEffectRows handles a Trait name containing '_' and errors on gaps", {
    mk <- c("m1", "m2")
    rn <- paste0("My_Env_", rep(c("A", "B"), each = 2), ":mbf(\"ints\")_", mk)
    idx <- .mbfEffectRows(rn, "My_Env", c("A", "B"), mk)
    expect_equal(as.vector(idx), 1:4)
    expect_error(.mbfEffectRows(rn[-4], "My_Env", c("A", "B"), mk), "level\\(s\\): B")
})
