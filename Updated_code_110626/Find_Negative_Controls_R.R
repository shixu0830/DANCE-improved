# DANCE negative control triplet discovery (R-only implementation)
# For marginal DANCE: baseline = character(0), conditional = FALSE.
# For conditional DANCE: baseline = X_names, conditional = TRUE.

fisher_pval_vec <- function(r, n) {
  m <- length(r)
  if (n <= 3) return(rep(1, m))
  s <- sqrt(n - 3)
  r <- pmax(-0.999999, pmin(0.999999, r))
  z <- 0.5 * log((1 + r) / (1 - r)) * s
  2 * pnorm(-abs(z))
}

vt_pval <- function(x, y, z, w, Co, N) {
  D <- Co[x, z] * Co[y, w] - Co[x, w] * Co[y, z]
  
  det_xy <- Co[x, x] * Co[y, y] - Co[x, y]^2
  det_zw <- Co[z, z] * Co[w, w] - Co[z, w]^2
  
  # det of full 4x4 submatrix
  idx4 <- c(x, y, z, w)
  det4 <- det(Co[idx4, idx4])
  
  # Wishart variance estimate
  v <- ((N + 1) / (N - 1) * det_xy * det_zw - det4) / (N - 2)
  if (!is.finite(v)) return(0)     
  if (v <= 0) v <- abs(v)          
  if (v < 1e-12) v <- 1e-12        
  
  zstat <- D / sqrt(v)
  2 * pnorm(-abs(zstat))
}

triplet_valid_minp <- function(a, b, cc, t_idx, o_idx, Co, N, vtalpha) {
  p1 <- vt_pval(a,  b,  cc, t_idx, Co, N); if (p1 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  p2 <- vt_pval(a,  cc, b,  t_idx, Co, N); if (p2 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  p3 <- vt_pval(b,  cc, a,  t_idx, Co, N); if (p3 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  p4 <- vt_pval(a,  b,  cc, o_idx, Co, N); if (p4 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  p5 <- vt_pval(a,  cc, b,  o_idx, Co, N); if (p5 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  p6 <- vt_pval(b,  cc, a,  o_idx, Co, N); if (p6 <= vtalpha) return(list(valid=FALSE, minp=NA_real_))
  
  list(valid = TRUE, minp = min(p1, p2, p3, p4, p5, p6))
}

# Tests VanTetrad({NC_i, NC_j}, {T, O}) — this is NOT one of the 6 DNCT tetrads.
# It catches candidates that violate the NC model structure entirely.
precompute_pair_ok <- function(cidx, t_idx, o_idx, Co, N, vtalpha) {
  m <- length(cidx)
  pair_ok <- matrix(0L, m, m)
  if (m < 2) return(pair_ok)
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      p <- vt_pval(cidx[i], cidx[j], t_idx, o_idx, Co, N)
      v <- as.integer(p > vtalpha)   
      pair_ok[i, j] <- v
      pair_ok[j, i] <- v
    }
  }
  pair_ok
}

# pair_ok is indexed by POSITION within cidx (1..m).
# Returns $idx as a K x 3 matrix of 1-BASED Co column indices.
screen_triplets_with_p <- function(cidx, t_idx, o_idx, Co, N, vtalpha,
                                   pair_ok, use_pair) {
  m <- length(cidx)
  if (m < 3) return(list(idx = matrix(nrow = 0, ncol = 3), pmin = numeric(0)))
  
  res_idx  <- vector("list", 0)
  res_pmin <- numeric(0)
  
  for (i in 1:(m - 2)) {
    for (j in (i + 1):(m - 1)) {
      if (use_pair && pair_ok[i, j] == 0L) next
      for (k in (j + 1):m) {
        if (use_pair && (pair_ok[i, k] == 0L || pair_ok[j, k] == 0L)) next
        
        res <- triplet_valid_minp(cidx[i], cidx[j], cidx[k],
                                  t_idx, o_idx, Co, N, vtalpha)
        if (res$valid) {
          res_idx[[length(res_idx) + 1L]] <- c(cidx[i], cidx[j], cidx[k])
          res_pmin[length(res_pmin) + 1L] <- res$minp
        }
      }
    }
  }
  
  if (length(res_idx) == 0L) {
    return(list(idx = matrix(nrow = 0, ncol = 3), pmin = numeric(0)))
  }
  list(idx = do.call(rbind, res_idx), pmin = res_pmin)
}


# Helper: consistent empty output
empty_triplet_matrix <- function() {
  matrix(character(0), nrow = 3L, ncol = 0L)
}

# Helper: check that a matrix can be used for tetrad screening
is_bad_matrix <- function(M) {
  is.null(M) || !is.matrix(M) || nrow(M) != ncol(M) ||
    nrow(M) < 2L || any(!is.finite(M))
}

# Helper: marginal or conditional correlation via Schur complement
conditional_cor_schur <- function(df, V, B = NULL) {
  B <- if (is.null(B)) character(0) else B
  D <- df[, c(V, B), drop = FALSE]
  cc <- stats::complete.cases(D)
  D <- D[cc, , drop = FALSE]
  n <- nrow(D)
  
  pV <- length(V)
  pB <- length(B)
  N_eff <- n - pB
  if (!is.finite(N_eff) || N_eff < 4L) {
    return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                reason = "too_few_complete_rows"))
  }
  
  # A covariance matrix is defined even if a candidate has zero variance.
  S <- tryCatch(stats::cov(D), error = function(e) NULL)
  if (is_bad_matrix(S)) {
    return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                reason = "bad_covariance_matrix"))
  }
  
  if (pB == 0L) {
    S_VgB <- S
  } else {
    SVV <- S[seq_len(pV), seq_len(pV), drop = FALSE]
    SVB <- S[seq_len(pV), pV + seq_len(pB), drop = FALSE]
    SBB <- S[pV + seq_len(pB), pV + seq_len(pB), drop = FALSE]
    
    if (any(diag(SBB) <= 0)) {
      return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                  reason = "zero_variance_baseline"))
    }
    
    SBBinv <- tryCatch(
      solve(SBB),
      error = function(e) tryCatch(MASS::ginv(SBB), error = function(e2) NULL)
    )
    if (is.null(SBBinv) || any(!is.finite(SBBinv))) {
      return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                  reason = "bad_baseline_inverse"))
    }
    
    S_VgB <- tryCatch(SVV - SVB %*% SBBinv %*% t(SVB),
                      error = function(e) NULL)
    if (is_bad_matrix(S_VgB)) {
      return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                  reason = "bad_schur_complement"))
    }
  }
  
  # Check the covariance matrix used for correlations (after Schur if conditional).
  # Only exactly zero/nonpositive or nonfinite variances are rejected.
  variances <- diag(S_VgB)
  bad_vars <- V[which(!is.finite(variances) | variances <= 0)]
  if (length(bad_vars)) {
    return(list(Co = NULL, N_eff = N_eff, bad_vars = bad_vars,
                reason = "nonpositive_variance"))
  }
  
  sdv <- sqrt(variances)
  Dn <- diag(1 / sdv, nrow = length(sdv))
  Co <- Dn %*% S_VgB %*% Dn
  if (is_bad_matrix(Co)) {
    return(list(Co = NULL, N_eff = N_eff, bad_vars = character(0),
                reason = "bad_correlation_matrix"))
  }
  dimnames(Co) <- list(V, V)
  list(Co = Co, N_eff = N_eff, bad_vars = character(0), reason = NULL)
}

FindNegativeControls <- function(data, tname, oname,
                                 baseline = character(0), # Vector of baseline covariate names
                                 vtalpha  = 0.05,         # Pairwise screening / tetrad p-value threshold
                                 independence = TRUE,   # Whether to pre-filter candidates using independence screen
                                 exclude_rare = TRUE,   # Exclude rare binary NCs
                                 alpha_corr   = 0.2,     # P-value threshold for independence screen
                                 use_pair = FALSE,       # Whether to use pairwise pre-screening
                                 FindBest = FALSE,       # Whether to select only the best triplet
                                 conditional = TRUE,    # Whether to condition on baseline variables
                                 return_diagnostics = FALSE) { # Attach screening diagnostics to output
  
  stopifnot(is.data.frame(data))
  if (!(tname %in% names(data))) stop("treatment name not in data")
  if (!(oname %in% names(data))) stop("outcome name not in data")
  if (tname == oname) stop("treatment and outcome must differ")
  
  # Keep the original 3-row triplet output; optionally attach diagnostics.
  excluded_rare <- character(0)
  excluded_zero_variance <- character(0)
  finish <- function(out, status = "success") {
    if (isTRUE(return_diagnostics)) {
      attr(out, "status") <- status
      attr(out, "excluded_rare") <- excluded_rare
      attr(out, "excluded_zero_variance") <- excluded_zero_variance
    }
    out
  }
  fail0 <- function(message) {
    warning(message)
    finish(empty_triplet_matrix(), "calculation_failed")
  }
  no_triplets <- function() finish(empty_triplet_matrix(), "no_triplets")
  
  # Step 1: Data preparation
  data_work <- data
  Bnames <- intersect(baseline, names(data_work))
  
  if (!conditional && length(baseline) > 0L) {
    stop("For completely marginal DANCE, use baseline = character(0).")
  }
  if (conditional && length(Bnames) == 0L) {
    stop("Conditional DANCE requires baseline X variables present in the data.")
  }
  
  df_num <- data_work[, c(tname, oname,
                          setdiff(names(data_work), c(tname, oname))), drop = FALSE]
  num_idx <- vapply(df_num, is.numeric, logical(1))
  df_num <- df_num[, num_idx, drop = FALSE]
  if (!(tname %in% names(df_num)) || !(oname %in% names(df_num))) {
    return(fail0("Treatment or outcome is nonnumeric."))
  }
  
  if (conditional && any(!vapply(data_work[, Bnames, drop = FALSE],
                                 is.numeric, logical(1)))) {
    return(fail0("Conditional baseline variables must be numeric."))
  }
  
  # For plasmode simulations, exclude simulated unmeasured confounders (U variables)
  # from the negative control candidate pool.
  # For applications to real-world data, review or remove this exclusion rule,
  # as variables starting with "U" may be legitimate NC candidates.
  Uvars <- grep("^U", colnames(df_num), value = TRUE)
  cand_all <- setdiff(colnames(df_num), c(tname, oname, Bnames, Uvars))
  if (length(cand_all) < 3L) return(no_triplets())
  
  # Step 2: Optional 1% prevalence screen 
  if (isTRUE(exclude_rare) && length(cand_all) > 0L) {
    is_ok <- vapply(cand_all, function(v) {
      x <- df_num[[v]]
      x <- x[!is.na(x)]
      ux <- sort(unique(x))
      if (length(ux) != 2L || !all(ux %in% c(0, 1))) return(TRUE)
      p1 <- mean(x == 1)
      p0 <- 1 - p1
      (p1 >= 0.01) && (p0 >= 0.01)
    }, logical(1))
    excluded_rare <- cand_all[!is_ok]
    cand_all <- cand_all[is_ok]
    if (length(cand_all) < 3L) return(no_triplets())
  }
  
  # Step 3: Compute correlations; remove any zero-variance candidate and retry.
  # Complete cases are recalculated inside conditional_cor_schur() each time.
  B_for_cc <- if (conditional) Bnames else NULL
  repeat {
    V <- c(tname, oname, cand_all)
    cc <- conditional_cor_schur(data_work, V = V, B = B_for_cc)
    
    if (length(cc$bad_vars)) {
      if (any(cc$bad_vars %in% c(tname, oname))) {
        return(fail0("Treatment or outcome has zero/nonpositive variance."))
      }
      bad_cands <- intersect(cc$bad_vars, cand_all)
      if (!length(bad_cands)) {
        return(fail0("Zero/nonpositive variance could not be attributed to candidates."))
      }
      excluded_zero_variance <- unique(c(excluded_zero_variance, bad_cands))
      cand_all <- setdiff(cand_all, bad_cands)
      if (length(cand_all) < 3L) return(no_triplets())
      next
    }
    
    if (is.null(cc$Co)) {
      if (identical(cc$reason, "too_few_complete_rows")) {
        return(fail0("Too few effective degrees of freedom after conditioning (N_eff < 4)."))
      }
      return(fail0(sprintf("Correlation construction failed: %s", cc$reason)))
    }
    break
  }
  
  Co_all <- cc$Co
  N_eff <- cc$N_eff
  
  # Step 4: Optional independence screening (uses same matrix as tetrad testing)
  cand0 <- cand_all
  if (independence) {
    rT <- as.numeric(Co_all[tname, cand0])
    rO <- as.numeric(Co_all[oname, cand0])
    pT <- fisher_pval_vec(rT, N_eff)
    pO <- fisher_pval_vec(rO, N_eff)
    keep <- is.finite(pT) & is.finite(pO) &
      (pT <= alpha_corr) & (pO <= alpha_corr)
    cand0 <- cand0[keep]
    if (length(cand0) < 3L) return(no_triplets())
  }
  
  # Step 5: Correlation matrix for selected candidates
  nm <- c(tname, oname, cand0)
  Co <- Co_all[nm, nm, drop = FALSE]
  storage.mode(Co) <- "double"
  if (is_bad_matrix(Co)) return(fail0("Final correlation matrix is invalid."))
  if ((ncol(Co) - 2L) < 3L) return(no_triplets())
  
  # Step 6: Prepare indices for R triplet screening
  cidx <- 3L:ncol(Co) # 1-based indices
  t_idx <- 1L
  o_idx <- 2L
  pair_ok <- if (use_pair) {
    precompute_pair_ok(cidx, t_idx, o_idx, Co, N_eff, vtalpha)
  } else {
    matrix(1L, nrow = length(cidx), ncol = length(cidx))
  }
  
  # Step 7: Triplet screen
  res <- tryCatch(
    screen_triplets_with_p(cidx, t_idx, o_idx, Co, N_eff, vtalpha,
                                pair_ok, use_pair),
    error = function(e) {
      warning(sprintf("screen_triplets_with_p failed: %s", e$message))
      NULL
    }
  )
  if (is.null(res)) return(finish(empty_triplet_matrix(), "calculation_failed"))
  
  trips <- res$idx
  pmins <- res$pmin
  if (!nrow(trips)) return(no_triplets())
  if (length(pmins) != nrow(trips) || any(!is.finite(pmins))) {
    return(fail0("Invalid R output: triplet counts or pmin values."))
  }
  
  # Step 8: Select best triplet or rank by minimum tetrad p-value
  if (isTRUE(FindBest)) {
    best <- which.max(pmins)
    trips <- matrix(trips[best, ], nrow = 1L)
  } else {
    ord <- order(pmins, decreasing = TRUE)
    trips <- trips[ord, , drop = FALSE]
  }
  
  # Step 9: Map R indices back to candidate names
  out <- matrix(nm[trips], ncol = 3L)
  finish(t(out))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
