Rcpp::sourceCpp("~/Dance_Plasmode/Simulation_code_1109/tetrad_helpers.cpp")

# --- helper: conditional correlation via Schur complement ---
conditional_cor_schur <- function(df, V, B) {
  if (length(B) == 0L) {
    D  <- df[, V, drop = FALSE]
    cc <- stats::complete.cases(D)
    D  <- D[cc, , drop = FALSE]
    n  <- nrow(D)
    Co <- stats::cor(D)                         # standard correlation (no baseline)
    dimnames(Co) <- list(V, V)
    return(list(Co = Co, N_eff = n))
  }
  D  <- df[, c(V, B), drop = FALSE]
  cc <- stats::complete.cases(D)
  D  <- D[cc, , drop = FALSE]
  n  <- nrow(D)
  
  pV <- length(V); pB <- length(B)
  S  <- stats::cov(D)                          
  
  SVV <- S[seq_len(pV),             seq_len(pV), drop = FALSE]
  SVB <- S[seq_len(pV),             pV + seq_len(pB), drop = FALSE]
  SBB <- S[pV + seq_len(pB),        pV + seq_len(pB), drop = FALSE]
  
  SBBinv <- tryCatch(solve(SBB), error = function(e) MASS::ginv(SBB))
  S_VgB  <- SVV - SVB %*% SBBinv %*% t(SVB)     # Schur complement Σ_{V|B} (up to scale)
  
  sdv <- sqrt(pmax(diag(S_VgB), 0))
  Dn  <- diag(ifelse(sdv > 0, 1/sdv, 0), nrow = length(sdv))
  Co  <- Dn %*% S_VgB %*% Dn                    # Corr_{V|B}
  dimnames(Co) <- list(V, V)
  
  N_eff <- n - pB                          # effective df for Wishart part
  list(Co = Co, N_eff = N_eff)
}

FindNegativeControls <- function(data, tname, oname,
                                 baseline = baseline,     # Vector of baseline covariate names
                                 vtalpha  = 0.05,        # Pairwise screening / tetrad p-value threshold
                                 independence = FALSE,   # Whether to pre-filter candidates using independence screen
                                 exclude_rare = FALSE,   # Exclude rare binary NCs
                                 alpha_corr   = 0.2,     # P-value threshold for independence screen
                                 use_pair = FALSE,       # Whether to use pairwise pre-screening
                                 FindBest = FALSE){       # Whether to select only the best triplet
  stopifnot(is.data.frame(data))
  if (!(tname %in% names(data))) stop("treatment name not in data")
  if (!(oname %in% names(data))) stop("outcome name not in data")
  if (tname == oname) stop("treatment and outcome must differ")
  
  # Step 1: Data preparation
  data_work <- data
  Bnames <- intersect(baseline, names(data_work))  # Keep only baseline columns present in data
  df_num <- data_work[, c(tname, oname, setdiff(names(data_work), c(tname, oname))), drop = FALSE] 
  num_idx <- vapply(df_num, is.numeric, logical(1)) # Identify numeric candidate variables
  df_num <- df_num[, num_idx, drop = FALSE]
  
  if (!(tname %in% names(df_num)) || !(oname %in% names(df_num))) {
    return(matrix(nrow = 3, ncol = 0))
  }
  
  cand_all <- setdiff(colnames(df_num), c(tname, oname, Bnames))
  if (length(cand_all) < 3L) return(matrix(nrow = 3, ncol = 0))
  
  # Step 2: Exclude rare binary NCs
  if (isTRUE(exclude_rare) && length(cand_all) > 0L) {
    is_ok <- vapply(cand_all, function(v) {
      x <- df_num[[v]]
      x <- x[!is.na(x)]
      ux <- sort(unique(x))
      
      # Only treat as binary if coded 0/1 (numeric)
      if (length(ux) != 2L || !all(ux %in% c(0, 1))) {
        return(TRUE)  # keep non-binary or non-0/1 variables
      }
      
      p1 <- mean(x == 1)
      p0 <- 1 - p1
      # keep if both levels are at least rare_threshold
      (p1 >= 0.01) && (p0 >= 0.01)
    }, logical(1))
    
    cand_all <- cand_all[is_ok]
    if (length(cand_all) < 3L) {
      # not enough candidates left to form a triplet
      return(matrix(nrow = 3, ncol = 0))
    }
  }
  
  
  # Step 3: Compute (conditional) correlation matrix
  B_for_cc <- if (conditional && length(Bnames) > 0) Bnames else NULL # Pass baseline to conditional correlation only if conditional = TRUE
  V <- c(tname, oname, cand_all)
  cc <- conditional_cor_schur(data_work, V = V, B = B_for_cc)
  Co_all <- cc$Co
  N_eff  <- cc$N_eff
  
  
  # Step 4: Optional independence screen step
  cand0 <- cand_all
  if (independence) {
    # --- (a) Marginal screen ---
    rT <- as.numeric(Co_all[tname, cand0])
    rO <- as.numeric(Co_all[oname, cand0])
    pT <- rcpp_fisher_pval_vec(rT, N_eff)
    pO <- rcpp_fisher_pval_vec(rO, N_eff)
    
    # Require association with both T and O
    keep <- (pT <= alpha_corr) & (pO <= alpha_corr)
    
    # --- (b) Additional conditional-on-baseline screen ---
    if (length(Bnames) > 0L && any(keep)) {
      V_cond <- c(tname, oname, cand0)
      cc_cond <- conditional_cor_schur(data_work, V = V_cond, B = Bnames)
      Co_cond <- cc_cond$Co
      N_eff_cond <- cc_cond$N_eff
      
      if (is.finite(N_eff_cond) && N_eff_cond > 3L) {
        rT_c <- as.numeric(Co_cond[tname, cand0])
        rO_c <- as.numeric(Co_cond[oname, cand0])
        pT_c <- rcpp_fisher_pval_vec(rT_c, N_eff_cond)
        pO_c <- rcpp_fisher_pval_vec(rO_c, N_eff_cond)
        
        # Require association with both T and O after adjusting for baselines
        keep <- keep & (pT_c <= alpha_corr) & (pO_c <= alpha_corr)
      }
    }
    cand0 <- cand0[keep]
    if (length(cand0) < 3L) return(matrix(nrow = 3, ncol = 0))
  }
  
  
  # Step 5: Final correlation matrix among treatment, outcome, and filtered candidates
  nm <- c(tname, oname, cand0)
  Co <- Co_all[nm, nm, drop = FALSE]
  storage.mode(Co) <- "double"
  
  if ((ncol(Co) - 2L) < 3L) return(matrix(nrow = 3, ncol = 0))
  
  # Step 6: Prepare indices for C++ triplet screening after optional pairwise screening
  cidx_cpp <- as.integer(2L:(ncol(Co) - 1L))  # 0-based indices
  t_idx <- 0L
  o_idx <- 1L
  
  pair_ok <- if (use_pair) {
    rcpp_precompute_pair_ok(cidx_cpp, t_idx, o_idx, Co, N_eff, vtalpha)
  } else {
    matrix(1L, nrow = length(cidx_cpp), ncol = length(cidx_cpp))
  }
  
  # Step 7: Triplet screen
  res <- tryCatch(
    rcpp_screen_triplets_with_p(cidx_cpp, t_idx, o_idx, Co, N_eff, vtalpha, pair_ok, use_pair),
    error = function(e) {
      warning(sprintf("rcpp_screen_triplets_with_p failed: %s", e$message))
      return(list(idx = matrix(nrow = 0, ncol = 3), pmin = numeric(0)))
    }
  )
  
  trips <- res$idx
  pmins <- res$pmin
  if (!nrow(trips)) return(matrix(nrow = 3, ncol = 0))
  
  # Step 8: Select triplets based on FindBest or order by p-values
  if (isTRUE(FindBest)) {
    best <- which.max(pmins)
    trips <- matrix(trips[best, ], nrow = 1)
  } else {
    ord <- order(pmins, decreasing = TRUE)
    trips <- trips[ord, , drop = FALSE]
  }
  
  # Step 9: Map indices back to variable names
  out <- matrix(nm[trips + 1L], ncol = 3)
  t(out)
}


`%||%` <- function(x, y) if (is.null(x)) y else x