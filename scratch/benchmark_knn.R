# =============================================================================
# benchmark_knn.R
# Compare original vs vectorised KNN imputation across chromosome sizes.
#
# Two implementations compared:
#
#   ORIGINAL:
#     - Distance: dist() on mean-substituted matrix
#     - Impute:   loop over every missing (line, marker) pair; calls order()
#                 inside the loop for every pair
#
#   VECTORISED:
#     - Distance: tcrossprod() via BLAS (row outer-product formula);
#                 significantly faster than dist() for large p
#     - Impute:   pre-sort the n x n neighbour matrix once; per-marker step
#                 uses vectorised logical indexing then vapply over missing
#                 lines only
#
# Panel sizes:  n_markers = 5000, 10000, 15000 per chromosome
# Line counts:  n_lines   = 100, 200, 500
# Missing rate: 5%,  k = 5
# =============================================================================

# ---- Original implementation ------------------------------------------------
knn_impute_original <- function(geno_chr, k = 5L) {

    na_cols   <- which(colSums(is.na(geno_chr)) > 0L)
    col_means <- colMeans(geno_chr, na.rm = TRUE)

    # Step A: dist() on mean-substituted matrix
    geno_sub <- geno_chr
    for (j in na_cols)
        geno_sub[is.na(geno_sub[, j]), j] <- col_means[j]

    t_dist <- system.time(d_mat <- as.matrix(dist(geno_sub)))["elapsed"]

    # Step B: loop — order() called for every missing (line, marker) pair
    result <- geno_chr
    t_imp  <- system.time({
        for (j in na_cols) {
            miss_lines <- which(is.na(geno_chr[, j]))
            if (!length(miss_lines)) next
            nonmiss <- which(!is.na(geno_chr[, j]))
            for (i in miss_lines) {
                d_i  <- d_mat[i, nonmiss]
                nn   <- nonmiss[order(d_i)[seq_len(min(k, length(nonmiss)))]]
                result[i, j] <- mean(geno_chr[nn, j])
            }
        }
    })["elapsed"]

    list(t_dist = t_dist, t_imp = t_imp, t_total = t_dist + t_imp)
}

# ---- Vectorised implementation ----------------------------------------------
knn_impute_vectorised <- function(geno_chr, k = 5L) {

    na_cols   <- which(colSums(is.na(geno_chr)) > 0L)
    col_means <- colMeans(geno_chr, na.rm = TRUE)
    n         <- nrow(geno_chr)

    # Step A: distance via tcrossprod() (BLAS inner-product formula)
    #   dist^2(i,j) = ||x_i||^2 + ||x_j||^2 - 2 * <x_i, x_j>
    geno_sub <- geno_chr
    for (j in na_cols)
        geno_sub[is.na(geno_sub[, j]), j] <- col_means[j]

    t_dist <- system.time({
        row_ss   <- rowSums(geno_sub^2)                     # n-vector
        cp       <- tcrossprod(geno_sub)                    # n x n via BLAS
        dist_sq  <- outer(row_ss, row_ss, "+") - 2 * cp
        d_mat    <- sqrt(pmax(dist_sq, 0))
        diag(d_mat) <- 0
    })["elapsed"]

    # Step B: pre-sort neighbour matrix once; impute with vectorised indexing
    result <- geno_chr
    t_imp  <- system.time({
        # nn_rank[r, i] = index of r-th nearest neighbour of line i
        nn_rank <- apply(d_mat, 1L, order)          # n x n; done ONCE

        for (j in na_cols) {
            miss_j    <- which( is.na(geno_chr[, j]))
            nonmiss_j <- which(!is.na(geno_chr[, j]))
            if (!length(miss_j) || !length(nonmiss_j)) next

            # Logical flag: which neighbour slots are non-missing at marker j
            nonmiss_flag          <- logical(n)
            nonmiss_flag[nonmiss_j] <- TRUE

            # Sub-matrix of pre-sorted neighbours for missing lines only
            nn_sub    <- nn_rank[, miss_j, drop = FALSE]    # n x n_miss
            valid_mat <- matrix(nonmiss_flag[nn_sub],
                                nrow = nrow(nn_sub))         # n x n_miss (logical)

            # Impute each missing line from its k nearest valid neighbours
            result[miss_j, j] <- vapply(seq_along(miss_j), function(idx) {
                vr <- which(valid_mat[, idx])
                if (!length(vr)) return(col_means[j])
                nn_idx <- nn_sub[vr[seq_len(min(k, length(vr)))], idx]
                mean(geno_chr[nn_idx, j])
            }, numeric(1L))
        }
    })["elapsed"]

    list(t_dist = t_dist, t_imp = t_imp, t_total = t_dist + t_imp)
}

# =============================================================================
# Benchmark grid
# =============================================================================
n_markers_vec <- c(5000L, 10000L, 15000L)
n_lines_vec   <- c(100L, 200L, 500L)
miss_rate     <- 0.05
k             <- 5L

cat("\n")
cat(strrep("=", 78), "\n")
cat("  KNN imputation benchmark — original vs vectorised\n")
cat(sprintf("  Missing rate: %.0f%%   k = %d\n", miss_rate * 100, k))
cat(strrep("=", 78), "\n")
cat(sprintf("  %-6s  %-8s  %10s  %10s  %12s  %10s\n",
            "lines", "markers", "orig (s)", "vect (s)", "speedup (x)", "ratio dist"))
cat(strrep("-", 78), "\n")

all_res <- list()
for (nl in n_lines_vec) {
    for (nm in n_markers_vec) {
        set.seed(42)
        geno_chr <- matrix(
            sample(c(-1L, 0L, 1L), nl * nm, replace = TRUE,
                   prob = c(0.35, 0.30, 0.35)),
            nrow = nl,
            dimnames = list(paste0("L", seq_len(nl)),
                            paste0("M", seq_len(nm)))
        )
        geno_chr[sample(length(geno_chr),
                         round(miss_rate * length(geno_chr)))] <- NA

        t_orig <- knn_impute_original(geno_chr, k = k)
        t_vect <- knn_impute_vectorised(geno_chr, k = k)

        speedup <- t_orig$t_total / max(t_vect$t_total, 1e-6)

        cat(sprintf("  %-6d  %-8d  %10.2f  %10.2f  %12.1fx  %9.1f%%\n",
                    nl, nm,
                    t_orig$t_total, t_vect$t_total,
                    speedup,
                    100 * t_vect$t_dist / max(t_vect$t_total, 1e-6)))

        all_res[[length(all_res) + 1L]] <- list(
            nl = nl, nm = nm,
            orig = t_orig, vect = t_vect, speedup = speedup)
    }
    cat(strrep("-", 78), "\n")
}

# ---- Per-step breakdown ------------------------------------------------------
cat("\n")
cat(strrep("=", 78), "\n")
cat("  Vectorised: time breakdown by step\n")
cat(strrep("=", 78), "\n")
cat(sprintf("  %-6s  %-8s  %14s  %14s  %12s\n",
            "lines", "markers", "dist (s)", "impute (s)", "total (s)"))
cat(strrep("-", 78), "\n")
for (r in all_res) {
    cat(sprintf("  %-6d  %-8d  %14.2f  %14.2f  %12.2f\n",
                r$nl, r$nm,
                r$vect$t_dist, r$vect$t_imp, r$vect$t_total))
}
cat(strrep("-", 78), "\n")

# ---- Whole-genome estimate (10 chromosomes) ----------------------------------
cat("\n")
cat(strrep("=", 78), "\n")
cat("  Estimated whole-genome time — vectorised (10 equal chromosomes)\n")
cat(strrep("=", 78), "\n")
cat(sprintf("  %-6s  %-14s  %12s  %12s\n",
            "lines", "total markers", "time (s)", "time (min)"))
cat(strrep("-", 78), "\n")
for (r in all_res) {
    total_s <- r$vect$t_total * 10
    cat(sprintf("  %-6d  %-14d  %12.1f  %12.2f\n",
                r$nl, r$nm * 10L, total_s, total_s / 60))
}
cat(strrep("=", 78), "\n\n")
