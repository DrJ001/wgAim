# =============================================================================
# accuracy.R  —  Mrode (2014) model-based BLUP accuracy for ASReml-R V4
# =============================================================================
# Computes per-group mean accuracy:  r_j = mean[ sqrt(1 - PEV_ij / G_jj) ]
# Optionally also returns Cullis (2006) H² per group.
#
# Supports:  fa(Group, k):id(Var),  diag(Group):id(Var),
#            corgh(Group):id(Var),  corh(Group):id(Var),
#            us(Group):id(Var),     and id(Var) alone.
#
# predict.asreml() only= string used per structure:
#   fa(Group, k):id(Var)  ->  "fa(Group, k):Var"   (space after comma required)
#   diag / corgh / corh / us  ->  "Group:Var"
#   id(Var) alone         ->  "Var"
#
# Usage:
#   accuracy(model)                              # both metrics (default)
#   accuracy(model, metric = "accuracy")         # Mrode accuracy only
#   accuracy(model, metric = "gen.H2")           # Cullis (2006) generalised H2 only
#   accuracy(model, by_variety = TRUE)           # per-variety rows
# =============================================================================


#' Model-Based BLUP Accuracy for ASReml-R V4 Multi-Environment Trials
#'
#' Computes per-group mean Mrode (2014) accuracy and/or Cullis (2006) H\eqn{^2}
#' from an ASReml-R V4 mixed model. Supports \code{fa()}, \code{diag()},
#' \code{corgh()}, \code{corh()}, \code{us()}, and single-environment
#' \code{id()} random structures for the variety term.
#'
#' @param model      A fitted \code{asreml} model object.
#' @param term       Character. The full random-effect interaction term for the
#'   variety effect, written in the same format as the \code{random =} formula
#'   argument, e.g. \code{"fa(Site, 2):id(Variety)"},
#'   \code{"corgh(Site):vm(Variety, giv1)"},
#'   or \code{"vm(Variety, giv1)"} for a single-environment model.
#'   When \code{NULL} (default) the term is auto-detected from the model's
#'   random formula.
#' @param metric     Character vector: \code{"accuracy"}, \code{"gen.H2"},
#'   or both (default). Controls which columns appear in the output.
#' @param pworkspace Character. Passed to \code{predict.asreml()}.
#'   Default \code{"2gb"}.
#' @param by_variety Logical. If \code{TRUE}, return one row per variety x
#'   group instead of group-level summaries. Default \code{FALSE}.
#'
#' @return A data frame with columns \code{group}, \code{n_vars}, \code{G_jj},
#'   and one or both of \code{mean_acc}/\code{sd_acc} (Mrode accuracy) and
#'   \code{gen.H2} (Cullis H\eqn{^2}). When \code{by_variety = TRUE} the
#'   \code{variety}, \code{pev}, and \code{present} columns replace
#'   \code{n_vars}/\code{sd_acc}. The \code{present} column is \code{TRUE}
#'   for variety x group combinations observed in the data and \code{FALSE}
#'   for unobserved variety x group combinations (e.g. unobserved variety
#'   x site pairs in factor-analytic models).
#'
#' @examples
#' \dontrun{
#' acc    <- accuracy(model_fa)
#' acc_bv <- accuracy(model_fa, by_variety = TRUE)
#' accuracy(model_fa, metric = "gen.H2")
#' # Override auto-detection with the full term string
#' accuracy(model_fa, term = "fa(Site, 2):id(Variety)")
#' accuracy(model_vm, term = "corgh(Site):vm(Variety, giv1)")
#' }
#'
#' @export
accuracy <- function(model,
                     term       = NULL,
                     metric     = c("accuracy", "gen.H2"),
                     pworkspace = "2gb",
                     by_variety = FALSE) {

  metric <- match.arg(metric, several.ok = TRUE)
  want_H2  <- "gen.H2"   %in% metric
  want_acc <- "accuracy" %in% metric

  p <- if (!is.null(term)) .parse_acc_term(term) else .parse_met_formula(model)

  dat <- .get_model_data(model)
  gl  <- if (!is.null(p$group_var))
           levels(droplevels(dat[[p$group_var]])) else NULL

  # When there is no group structure (single-environment), treat the whole
  # model as one group labelled "all" so it runs through the same loop.
  single_env <- is.null(gl)
  if (single_env) gl <- "all"

  G <- .extract_G_diag(model, if (single_env) NULL else gl, p$group_var)
  if (is.null(G)) stop("[accuracy] Cannot extract G_jj from model.")
  if (single_env) G <- setNames(G, "all")

  # predict(only=) returns pure BLUP SEs — avoids fixed-effect inflation
  # that causes PEV > G_jj and spurious zero accuracies
  need_sed <- want_H2
  out <- predict(model, classify = p$classify, only = p$only_term,
                 sed = need_sed, pworkspace = pworkspace)
  pv  <- out$pvals
  names(pv) <- tolower(names(pv))
  pv$pev    <- pv$std.error^2
  pv$.row   <- seq_len(nrow(pv))   # preserve original row order for SED lookup

  # Mark which variety x group combinations are observed in the data.
  # pv$present = TRUE  : combination observed in the training data
  # pv$present = FALSE : combination not observed in the data (unobserved)
  # pv$finite          : predict() returned a usable (finite, positive) SE
  gc <- if (!single_env) tolower(p$group_var) else NULL
  vc <- tolower(p$by_var)

  gd <- if (!single_env) names(dat)[match(tolower(p$group_var), tolower(names(dat)))] else NULL
  vd <- names(dat)[match(tolower(p$by_var), tolower(names(dat)))]

  pv$finite <- is.finite(pv$std.error) & pv$std.error > 0

  if (!single_env && !is.null(gd) && !is.null(vd)) {
    obs        <- paste(dat[[gd]], dat[[vd]], sep = "\001")
    pv$present <- paste(pv[[gc]], pv[[vc]], sep = "\001") %in% obs
  } else {
    # Single-environment: all finite rows are "present"
    pv$present <- pv$finite
  }

  # ---- Unified loop (works for single- and multi-environment) ---------------
  rows <- lapply(gl, function(grp) {

    # Group-level: count only observed varieties (keeps n_vars meaningful).
    # by_variety: return all finite-SE rows; observed status carried in $present.
    if (by_variety) {
      idx <- if (single_env) which(pv$finite)
             else            which(pv[[gc]] == grp & pv$finite)
    } else {
      idx <- if (single_env) which(pv$present & pv$finite)
             else            which(pv[[gc]] == grp & pv$present & pv$finite)
    }
    gjj <- unname(G[grp])   # unname() prevents row-name warnings in data.frame

    if (length(idx) < 2L || is.na(gjj) || gjj <= 0) {
      if (by_variety) {
        # Return zero-row frame with the correct columns
        row <- data.frame(group = character(0), variety = character(0),
                          G_jj = numeric(0),   pev     = numeric(0),
                          present = logical(0))
        if (want_acc) row$accuracy  <- numeric(0)
        if (want_H2)  row$gen.H2   <- numeric(0)
      } else {
        row <- data.frame(group = grp, n_vars = 0L, G_jj = gjj)
        if (want_acc)  { row$mean_acc <- NA_real_; row$sd_acc <- NA_real_ }
        if (want_H2)     row$gen.H2  <- NA_real_
      }
      return(row)
    }

    # Mrode accuracy
    a <- sqrt(pmax(0, 1 - pv$pev[idx] / gjj))

    # Cullis (2006) generalised H²: 1 - (mean_SED)² / (2 * G_jj)
    h2 <- if (want_H2) {
      if (!is.null(out$sed) && max(idx) <= nrow(out$sed)) {
        ut <- out$sed[idx, idx]
        ut <- ut[upper.tri(ut)]
        ut <- ut[is.finite(ut) & ut > 0]
        if (length(ut)) max(0, min(1, 1 - mean(ut)^2 / (2 * gjj))) else NA_real_
      } else NA_real_
    } else NULL

    if (by_variety) {
      out_bv <- data.frame(group    = grp,
                           variety  = as.character(pv[[vc]][idx]),
                           G_jj     = gjj,
                           pev      = pv$pev[idx],
                           present  = pv$present[idx],
                           stringsAsFactors = FALSE)
      if (want_acc) out_bv$accuracy <- a
      if (want_H2)  out_bv$gen.H2  <- h2   # group-level statistic (Cullis H²)
      return(out_bv)
    }

    row <- data.frame(group  = grp,
                      n_vars = length(idx),
                      G_jj   = gjj)
    if (want_acc) {
      row$mean_acc <- mean(a, na.rm = TRUE)
      row$sd_acc   <- sd(a,   na.rm = TRUE)
    }
    if (want_H2) row$gen.H2 <- h2
    row
  })

  do.call(rbind, c(rows, make.row.names = FALSE))
}

# =============================================================================
# HELPERS
# =============================================================================

.get_model_data <- function(model) {
  dat <- tryCatch(eval(model$call$data), error = function(e) NULL)
  if (is.data.frame(dat)) return(dat)
  dat <- tryCatch(eval(model$call$data, envir = parent.frame(2L)),
                  error = function(e) NULL)
  if (is.data.frame(dat)) return(dat)
  stop("[accuracy] Cannot retrieve data from model$call$data.")
}

.extract_G_diag <- function(model, group_levels, group_var = NULL) {

  vc  <- summary(model)$varcomp
  nms <- rownames(vc)
  cmp <- vc$component

  if (is.null(group_levels)) {
    ok <- which(!grepl("units|!R$|R!|residual", nms, ignore.case = TRUE))
    if (length(ok) >= 1L) return(setNames(max(0, cmp[ok[1L]]), "all"))
    return(NULL)
  }

  # ---- FA: specific variances (Psi) + factor loadings (Lambda) --------------
  # Psi row names end with "!GroupLevel!var" (FA notation)
  vi_fa <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("!", g, "!var")))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_fa))) {
    Psi <- pmax(0, cmp[vi_fa])

    Lc <- list()
    for (f in seq_len(20L)) {
      fi <- sapply(group_levels, function(g) {
        i <- which(endsWith(nms, paste0("!", g, "!fa", f)))
        if (length(i) == 1L) i else NA_integer_
      })
      if (any(is.na(fi))) break
      Lc[[f]] <- cmp[fi]
    }

    G <- if (length(Lc)) rowSums(do.call(cbind, Lc)^2) + Psi else Psi
    return(setNames(G, group_levels))
  }

  # ---- diag / corgh / corh: specific variances -------------------------------
  # Row names end with "!Env_GroupLevel"  (e.g. "Env:Variety!Env_E1")
  vi_diag <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("_", g)))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_diag))) {
    G <- pmax(0, cmp[vi_diag])
    return(setNames(G, group_levels))
  }

  # ---- us: unstructured — diagonal entries are "!Env_Gj:Gj" ----------------
  # Off-diagonal entries are "!Env_Gj:Gi" (j > i); we want only j == i.
  vi_us <- sapply(group_levels, function(g) {
    i <- which(endsWith(nms, paste0("_", g, ":", g)))
    if (length(i) == 1L) i else NA_integer_
  })

  if (!all(is.na(vi_us))) {
    G <- pmax(0, cmp[vi_us])
    return(setNames(G, group_levels))
  }

  NULL
}

# Split string on sep at parenthesis depth 0
.split_top <- function(s, sep = "+") {
  res <- character(0L); depth <- 0L; start <- 1L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if      (ch == "(")           depth <- depth + 1L
    else if (ch == ")")           depth <- depth - 1L
    else if (ch == sep && !depth) { res <- c(res, substr(s, start, i-1L)); start <- i+1L }
  }
  c(res, substr(s, start, nchar(s)))
}

# Position of first ":" at depth 0
.top_colon <- function(s) {
  depth <- 0L
  for (i in seq_len(nchar(s))) {
    ch <- substr(s, i, i)
    if      (ch == "(") depth <- depth + 1L
    else if (ch == ")") depth <- depth - 1L
    else if (ch == ":" && !depth) return(i)
  }
  NA_integer_
}

# Content inside outermost parens, e.g. "id(Var)" -> "Var"
.inside <- function(s) {
  op <- regexpr("\\(", s)
  if (op < 0L) return(s)
  substr(s, op + 1L, nchar(s) - 1L)
}

# Parse a full random-effect interaction term string into the components needed
# by accuracy(): group_var, by_var, n_fa, only_term, classify.
#
# Accepts the same term formats as randomRegress() plus single-environment
# terms (no colon):
#   Multi-env:   "fa(Site, 2):id(Variety)"  "corgh(Site):vm(Variety, giv1)"
#   Single-env:  "id(Variety)"  "vm(Variety, giv1)"  "Variety"
#
# Rules for classify vs only=:
#   classify  — always bare variable names: "Site:Variety" or "Variety"
#   only      — bare group name (or fa spec) : rhs
#               id()  stripped; vm()/ide() kept; bare name used as-is
#' @noRd
.parse_acc_term <- function(term_str) {

  # Preserve original spacing before stripping whitespace.
  # ASReml-R's predict(only=) must match term labels exactly — including spaces
  # inside wrappers such as vm(Variety, giv1).  We use the stripped string only
  # for parsing (field names, structure detection), and the original for only=.
  cp_orig <- .top_colon(term_str)

  # ---- Single-environment (no top-level colon) ----------------------------
  if (is.na(cp_orig)) {
    rhs_orig <- term_str                         # preserve original spacing
    term_str <- gsub("\\s+", "", term_str)
    rhs      <- term_str
    if (grepl("^id\\(", rhs)) {
      by_var   <- .inside(rhs)
      only_rhs <- by_var
    } else if (grepl("^[A-Za-z][A-Za-z0-9_]*\\(", rhs)) {
      by_var   <- trimws(strsplit(.inside(rhs), ",")[[1L]][1L])
      only_rhs <- rhs_orig                       # use original spacing for only=
    } else {
      by_var   <- rhs
      only_rhs <- rhs
    }
    return(list(group_var  = NULL,
                by_var     = by_var,
                n_fa       = NULL,
                only_term  = only_rhs,
                classify   = by_var))
  }

  # ---- Multi-environment --------------------------------------------------
  rgt_orig <- substr(term_str, cp_orig + 1L, nchar(term_str))  # rhs original spacing
  term_str  <- gsub("\\s+", "", term_str)
  cp        <- .top_colon(term_str)

  lft  <- substr(term_str, 1L, cp - 1L)
  rgt  <- substr(term_str, cp + 1L, nchar(term_str))
  type <- sub("\\(.*", "", lft)

  if (!type %in% c("fa", "us", "corgh", "corh", "diag"))
    stop("Unsupported variance structure '", type, "' in term. ",
         "Supported: fa, us, corgh, corh, diag.")

  grp <- trimws(strsplit(.inside(lft), ",")[[1L]][1L])
  nfa <- if (type == "fa")
    suppressWarnings(as.integer(trimws(strsplit(.inside(lft), ",")[[1L]][2L])))
  else NULL

  # Right-hand side: id() stripped; vm()/ide() kept (with original spacing); bare as-is
  if (grepl("^id\\(", rgt)) {
    by_var   <- .inside(rgt)
    only_rhs <- by_var
  } else if (grepl("^[A-Za-z][A-Za-z0-9_]*\\(", rgt)) {
    by_var   <- trimws(strsplit(.inside(rgt), ",")[[1L]][1L])
    only_rhs <- rgt_orig                         # use original spacing for only=
  } else {
    by_var   <- rgt
    only_rhs <- rgt
  }

  # Build only= string (fa requires "fa(Group, k):rhs"; space after comma required)
  only <- if (type == "fa" && !is.null(nfa))
    sprintf("fa(%s, %d):%s", grp, nfa, only_rhs)
  else
    paste0(grp, ":", only_rhs)

  list(group_var  = grp,
       by_var     = by_var,
       n_fa       = nfa,
       only_term  = only,
       classify   = paste0(grp, ":", by_var))
}

# Extract the relevant term string from the model's random formula and
# delegate to .parse_acc_term() for the actual parsing.
#' @noRd
.parse_met_formula <- function(model) {

  raw <- tryCatch(paste(deparse(model$call$random), collapse = ""),
                  error = function(e)
                    stop("[accuracy] Cannot read model$call$random."))

  # Strip only the leading "~" and any surrounding whitespace, then split on
  # top-level "+" operators.  Do NOT collapse all internal whitespace — we need
  # to preserve spaces inside wrappers (e.g. "vm(Variety, giv1)") so that
  # the only= string passed to predict.asreml() matches the stored term label.
  rnd   <- trimws(sub("^[[:space:]]*~[[:space:]]*", "", raw))
  terms <- .split_top(rnd, "+")
  terms <- trimws(terms)   # trim whitespace around each top-level term

  # Prefer MET terms (supported multi-env variance structures); fall back to
  # the first term for single-environment models.
  met      <- grep("^(fa|diag|corh|corgh|us)[[:space:]]*[(]", terms, value = TRUE)
  term_str <- if (length(met)) met[1L] else terms[1L]

  .parse_acc_term(term_str)
}
