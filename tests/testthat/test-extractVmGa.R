# =============================================================================
# test-extractVmGa.R
# .extractVmGa(): structure-agnostic Ga extraction for the MV genomic term.
# Mock vparameter names are copied verbatim from ASReml-R 4.2 fits of
#   <str>(Trial):vm(id, covObj) + diag(Trial):id
# =============================================================================

.mock_vm <- function(vals, sigma2 = 1) {
    m <- list(vparameters = c(vals, "Trial_T1!R" = 0.15), sigma2 = sigma2)
    class(m) <- c("asreml", "list")
    m
}
key <- "Trial:vm(id, covObj)"
tl3 <- c("T1", "T2", "T3")

test_that(".extractVmGa: diag", {
    m  <- .mock_vm(setNames(c(0.2, 0.1, 0.3), paste0(key, "!Trial_", tl3)))
    Ga <- .extractVmGa(m, "Trial", tl3)
    expect_equal(unname(Ga), diag(c(0.2, 0.1, 0.3)))
    expect_equal(dimnames(Ga), list(tl3, tl3))
})

test_that(".extractVmGa: diag variances scaled by sigma2 and floored at 0", {
    m  <- .mock_vm(setNames(c(0.2, -0.01, 0.3), paste0(key, "!Trial_", tl3)))
    Ga <- .extractVmGa(m, "Trial", tl3, sigma2 = 2)
    expect_equal(unname(diag(Ga)), c(0.4, 0, 0.6))
})

test_that(".extractVmGa: corh uses the single '!Trial!cor' correlation", {
    v  <- c(0.26, 0.13, 0.28)
    m  <- .mock_vm(setNames(c(0.745, v),
                            c(paste0(key, "!Trial!cor"), paste0(key, "!Trial_", tl3))))
    Ga <- .extractVmGa(m, "Trial", tl3)
    expect_equal(unname(diag(Ga)), v)
    expect_equal(Ga[2, 1], 0.745 * sqrt(v[2] * v[1]))
    expect_equal(Ga[3, 2], 0.745 * sqrt(v[3] * v[2]))
    expect_false(anyNA(Ga))
})

test_that(".extractVmGa: corgh places correlations by name (3 trials)", {
    v  <- c(0.27, 0.13, 0.28)
    cn <- paste0(key, c("!Trial!T2:!Trial!T1.cor", "!Trial!T3:!Trial!T1.cor",
                        "!Trial!T3:!Trial!T2.cor"))
    m  <- .mock_vm(setNames(c(0.80, 0.64, 0.88, v), c(cn, paste0(key, "!Trial_", tl3))))
    Ga <- .extractVmGa(m, "Trial", tl3)
    expect_equal(Ga["T2", "T1"], 0.80 * sqrt(v[2] * v[1]))
    expect_equal(Ga["T3", "T1"], 0.64 * sqrt(v[3] * v[1]))
    expect_equal(Ga["T3", "T2"], 0.88 * sqrt(v[3] * v[2]))
    expect_equal(Ga, t(Ga))
})

test_that(".extractVmGa: corgh with 4 trials is independent of parameter order", {
    tl4 <- paste0("T", 1:4)
    v   <- c(0.1, 0.2, 0.3, 0.4)
    R   <- matrix(c(1, .1, .2, .3,  .1, 1, .4, .5,  .2, .4, 1, .6,  .3, .5, .6, 1), 4)
    pairs <- list(c(2,1), c(3,1), c(3,2), c(4,1), c(4,2), c(4,3))
    cn  <- vapply(pairs, function(p) sprintf("%s!Trial!T%d:!Trial!T%d.cor", key, p[1], p[2]), "")
    cv  <- vapply(pairs, function(p) R[p[1], p[2]], 0)
    o   <- rev(seq_along(cn))        # scramble order
    m   <- .mock_vm(setNames(c(cv[o], v), c(cn[o], paste0(key, "!Trial_", tl4))))
    Ga  <- .extractVmGa(m, "Trial", tl4)
    expect_equal(unname(Ga), R * tcrossprod(sqrt(v)))
})

test_that(".extractVmGa: us reads full lower triangle (incl. covariances)", {
    nm <- paste0(key, c("!Trial_T1:T1", "!Trial_T2:T1", "!Trial_T2:T2",
                        "!Trial_T3:T1", "!Trial_T3:T2", "!Trial_T3:T3"))
    m  <- .mock_vm(setNames(c(0.27, 0.15, 0.13, 0.18, 0.17, 0.28), nm))
    Ga <- .extractVmGa(m, "Trial", tl3)
    expect_equal(unname(Ga),
                 matrix(c(.27, .15, .18,  .15, .13, .17,  .18, .17, .28), 3))
})

test_that(".extractVmGa: fa(k) builds Lam Lam' + Psi (incl. fixed zero loading)", {
    fk  <- "fa(Trial, 2):vm(id, covObj)"
    psi <- c(0.075, 0, 0.046)
    L1  <- c(0.437, 0.339, 0.403); L2 <- c(0, 0.113, 0.271)
    nm  <- c(paste0(fk, "!", tl3, "!var"), paste0(fk, "!", tl3, "!fa1"),
             paste0(fk, "!", tl3, "!fa2"))
    m   <- .mock_vm(setNames(c(psi, L1, L2), nm))
    Ga  <- .extractVmGa(m, "Trial", tl3)
    Lam <- cbind(L1, L2)
    expect_equal(unname(Ga), Lam %*% t(Lam) + diag(psi))
})

test_that(".extractVmGa: errors on unrecognised structure / missing term", {
    m <- .mock_vm(setNames(0.1, paste0(key, "!weird")))
    expect_error(.extractVmGa(m, "Trial", tl3), "Unrecognised variance structure")
    m2 <- .mock_vm(c("Trial:id!Trial_T1" = 0.1))
    expect_error(.extractVmGa(m2, "Trial", tl3), "No genomic")
})
