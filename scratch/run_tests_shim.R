# Minimal testthat shim (testthat not installed in this R library).
# Usage: Rscript scratch/run_tests_shim.R tests/testthat/test-xxx.R [...]
.libPaths(c("C:/jtaylor/Claude/wgAim/scratch/rlib", .libPaths()))
suppressPackageStartupMessages(library(wgAim))
env <- new.env(parent = asNamespace("wgAim"))
n.pass <- 0L; n.fail <- 0L; n.skip <- 0L
env$test_that <- function(desc, code) {
    res <- tryCatch({ suppressWarnings(capture.output(eval(substitute(code), new.env(parent = env)))); "pass" },
                    error = function(e) {
                        if (identical(conditionMessage(e), "SKIP")) return("skip")
                        cat("FAIL:", desc, "\n   ", conditionMessage(e), "\n"); "fail" })
    if (res == "pass") n.pass <<- n.pass + 1L
    else if (res == "skip") n.skip <<- n.skip + 1L
    else n.fail <<- n.fail + 1L
}
chk <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
env$expect_equal <- function(a, b, ..., info = NULL, label = NULL) { r <- all.equal(a, b, ...); chk(isTRUE(r), paste(r, collapse = "; ")) }
env$expect_true  <- function(x, info = NULL, label = NULL) chk(isTRUE(x), "not TRUE")
env$expect_false <- function(x, info = NULL, label = NULL) chk(identical(x, FALSE), "not FALSE")
env$expect_error <- function(expr, pattern = NULL, ...) {
    e <- tryCatch({ expr; NULL }, error = function(e) e)
    chk(!is.null(e), "no error thrown")
    if (!is.null(pattern)) chk(grepl(pattern, conditionMessage(e)), paste("message mismatch:", conditionMessage(e)))
}
env$expect_s3_class <- function(x, cls) chk(inherits(x, cls), "class mismatch")
env$skip_if_not_installed <- function(pkg, ...) if (!requireNamespace(pkg, quietly = TRUE)) stop("SKIP")
env$skip_if_not <- function(cond, ...) if (!isTRUE(cond)) stop("SKIP")
env$skip_if <- function(cond, ...) if (isTRUE(cond)) stop("SKIP")
env$skip <- function(...) stop("SKIP")
env$skip_on_cran <- function() invisible()
env$expect_gte <- function(a, b) chk(all(a >= b), "not >=")
env$expect_lte <- function(a, b) chk(all(a <= b), "not <=")
env$expect_gt  <- function(a, b) chk(all(a > b), "not >")
env$expect_lt  <- function(a, b) chk(all(a < b), "not <")
env$expect_length <- function(x, n) chk(length(x) == n, sprintf("length %d != %d", length(x), n))
env$expect_null <- function(x) chk(is.null(x), "not NULL")
env$expect_identical <- function(a, b) chk(identical(a, b), "not identical")
env$expect_type <- function(x, t) chk(typeof(x) == t, "type mismatch")
env$expect_is <- function(x, cls) chk(inherits(x, cls), "class mismatch")
env$expect_named <- function(x, nms = NULL, ...) chk(!is.null(names(x)) && (is.null(nms) || setequal(names(x), nms)), "names mismatch")
env$expect_match <- function(x, p, ...) chk(all(grepl(p, x, ...)), "no match")
env$expect_no_error <- function(expr, ...) invisible(expr)
env$expect_silent <- function(expr) invisible(expr)
env$expect_invisible <- function(expr) invisible(expr)
env$expect_visible <- function(expr) invisible(expr)
env$expect_setequal <- function(a, b) chk(setequal(a, b), "sets differ")
env$expect_warning <- function(expr, pattern = NULL, ...) {
    w <- NULL; withCallingHandlers(expr, warning = function(c) { w <<- c; invokeRestart("muffleWarning") })
    chk(!is.null(w), "no warning")
}
env$expect_message <- function(expr, pattern = NULL, ...) {
    m <- NULL; withCallingHandlers(expr, message = function(c) { m <<- c; invokeRestart("muffleMessage") })
    chk(!is.null(m), "no message")
}
env$expect_output <- function(expr, pattern = NULL, ...) {
    out <- capture.output(expr); if (!is.null(pattern)) chk(any(grepl(pattern, out, ...)), "output mismatch")
}
env$local_mocked_bindings <- function(...) stop("SKIP")
env$with_mocked_bindings <- function(...) stop("SKIP")
for (f in commandArgs(trailingOnly = TRUE)) {
    cat("==>", f, "\n")
    sys.source(f, envir = env)
}
cat(sprintf("\nPASS: %d  FAIL: %d  SKIP: %d\n", n.pass, n.fail, n.skip))
