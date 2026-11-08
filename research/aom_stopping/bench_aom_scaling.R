# Scaling of the AOM max-statistic simulation (no ASReml).
#   factor   : square-root factor of Vg (N = n * t): eigen vs pivoted Cholesky
#   simulate : T (m x n) %*% Z (n x B) for each of t traits, B = 2000
cat("BLAS:", if (!is.null(extSoftVersion()["BLAS"])) extSoftVersion()["BLAS"] else "?", "\n")
cat("La_library:", La_library(), "\n\n")
set.seed(1)
mkpsd <- function(N, r = N) { X <- matrix(rnorm(N * r), N, r); tcrossprod(X) / r }
tf <- function(expr) system.time(expr)[["elapsed"]]

cat("Factorisation of Vg (N x N):\n")
res <- data.frame()
for (N in c(600, 1600, 3200, 4800)) {
    A  <- mkpsd(N)
    te <- tf(eigen(A, symmetric = TRUE))
    tc <- tf(suppressWarnings(chol(A, pivot = TRUE)))
    res <- rbind(res, data.frame(N = N, eigen.s = te, chol.pivot.s = tc,
                                 mem.MB = round(8 * N^2 / 2^20)))
    rm(A); gc(FALSE)
}
print(res, row.names = FALSE)

cat("\nSimulation multiply per trait, B = 2000 (batches of 250):\n")
res2 <- data.frame()
for (cfg in list(c(m = 1500, n = 200), c(m = 1500, n = 600), c(m = 6000, n = 400),
                 c(m = 12000, n = 400), c(m = 12000, n = 800))) {
    m <- cfg[["m"]]; n <- cfg[["n"]]
    T <- matrix(rnorm(m * n), m, n)
    Z <- matrix(rnorm(n * 250), n, 250)
    t1 <- tf(for (i in 1:8) T %*% Z)          # 8 x 250 = 2000 draws
    res2 <- rbind(res2, data.frame(m = m, n = n, secs.per.trait = t1))
}
print(res2, row.names = FALSE)
