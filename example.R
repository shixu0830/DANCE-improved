## ============================================================
## Example ROC simulation for DANCE-improved
## ============================================================

## ----------------- Load core implementation -----------------

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Package 'Rcpp' is required. Please install it first.")
}
if (!requireNamespace("MESS", quietly = TRUE)) {
  stop("Package 'MESS' is required for AUC. Please install it.")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("Package 'RColorBrewer' is required for plotting. Please install it.")
}

Rcpp::sourceCpp("tetrad_helpers.cpp")   # C++ helpers
source("find_negative_controls.R")        # your main NC search

## Use conditional correlation on baseline X by default
conditional <- TRUE

## ---------------- Data generating process (continuous A) ----------------
## A, Y, Z, X all continuous.
## Z-types created:
##   - disconnected      : U, X -> Z   (valid NCs given X)
##   - independent       : X -> Z
##   - instrumental_A    : U, X, A -> Z
##   - instrumental_Y    : U, X -> Z,  Z -> Y
##   - confounder        : U, X, A -> Z, Z -> Y
##   - chain             : pair (Z_chain_ia, Z_chain_ib):
##                         U, X -> Z1; (Z1, U, X) -> Z2

generate_example_data <- function(
    n             = 1000,
    p_x           = 5,
    n_disc        = 20,
    n_indep       = 20,
    n_instr_A     = 20,
    n_instr_Y     = 20,
    n_confounder  = 20,
    n_chain_pairs = 10,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  ## Baseline covariates X
  X <- matrix(rnorm(n * p_x), nrow = n, ncol = p_x)
  colnames(X) <- paste0("X", seq_len(p_x))
  
  ## Latent confounder
  U <- rnorm(n)
  
  ## Continuous treatment A | X, U
  alpha_x <- runif(p_x, -0.5, 0.5)
  alpha_u <- 1.0
  A <- 0.5 + as.numeric(X %*% alpha_x) + alpha_u * U + rnorm(n)
  
  ## Outcome Y | A, X, U (before Z effects)
  beta_a <- 1.0
  beta_x <- runif(p_x, -0.5, 0.5)
  beta_u <- 1.0
  Y <- 0.5 + beta_a * A + as.numeric(X %*% beta_x) + beta_u * U + rnorm(n)
  
  ## Container for Z's + types
  z_list <- list()
  z_type <- character(0)
  add_var <- function(name, vec, lbl) {
    z_list[[name]] <<- vec
    z_type[name]   <<- lbl
  }
  
  ## 1) Valid NCs: disconnected (U, X -> Z)
  if (n_disc > 0) {
    for (j in seq_len(n_disc)) {
      gx <- runif(p_x, -0.5, 0.5)
      gu <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Zj <- as.numeric(X %*% gx + gu * U + rnorm(n))
      add_var(sprintf("Z_disc_%d", j), Zj, "disconnected")
    }
  }
  
  ## 2) Independent-type: X -> Z (no U, no A)
  if (n_indep > 0) {
    for (j in seq_len(n_indep)) {
      gx <- runif(p_x, -0.5, 0.5)
      Zj <- as.numeric(X %*% gx + rnorm(n))
      add_var(sprintf("Z_indep_%d", j), Zj, "independent")
    }
  }
  
  ## 3) Instrumental for A: U, X, A -> Z (no direct Y effect)
  if (n_instr_A > 0) {
    for (j in seq_len(n_instr_A)) {
      gx <- runif(p_x, -0.5, 0.5)
      gu <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      ga <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Zj <- as.numeric(X %*% gx + gu * U + ga * A + rnorm(n))
      add_var(sprintf("Z_instrA_%d", j), Zj, "instrumental_A")
    }
  }
  
  ## 4) Instrumental for Y: U, X -> Z; Z -> Y
  if (n_instr_Y > 0) {
    for (j in seq_len(n_instr_Y)) {
      gx <- runif(p_x, -0.5, 0.5)
      gu <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Zj <- as.numeric(X %*% gx + gu * U + rnorm(n))
      gy <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Y  <- Y + gy * Zj
      add_var(sprintf("Z_instrY_%d", j), Zj, "instrumental_Y")
    }
  }
  
  ## 5) Confounders: U, X, A -> Z; Z -> Y
  if (n_confounder > 0) {
    for (j in seq_len(n_confounder)) {
      gx <- runif(p_x, -0.5, 0.5)
      gu <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      ga <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Zj <- as.numeric(X %*% gx + gu * U + ga * A + rnorm(n))
      gy <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Y  <- Y + gy * Zj
      add_var(sprintf("Z_conf_%d", j), Zj, "confounder")
    }
  }
  
  ## 6) Chain pairs: U, X -> Z1; (Z1, U, X) -> Z2
  if (n_chain_pairs > 0) {
    for (j in seq_len(n_chain_pairs)) {
      gx1 <- runif(p_x, -0.5, 0.5)
      gu1 <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Z1  <- as.numeric(X %*% gx1 + gu1 * U + rnorm(n))
      
      gx2 <- runif(p_x, -0.5, 0.5)
      gu2 <- runif(1, 0.5, 1.5) * sample(c(-1, 1), 1)
      Z2  <- as.numeric(0.8 * Z1 + X %*% gx2 + gu2 * U + rnorm(n))
      
      add_var(sprintf("Z_chain_%da", j), Z1, "chain")
      add_var(sprintf("Z_chain_%db", j), Z2, "chain")
    }
  }
  
  ## Assemble observed data: U is latent (not included)
  df <- data.frame(
    A = A,
    Y = as.numeric(Y),
    as.data.frame(z_list),
    X,
    check.names = FALSE
  )
  
  Z_types <- data.frame(
    variable = names(z_type),
    type     = unname(z_type),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  list(
    data    = df,
    Z_types = Z_types,
    X_names = colnames(X)
  )
}

## ---------------- Helper utilities for ROC -------------------

# Map Z names to integer indices
z_indexer <- function(z_names) {
  idx <- seq_along(z_names)
  names(idx) <- z_names
  idx
}

# Encode unordered triplet (a, b, c) as a single integer key
encode_triplets <- function(a, b, c, N) {
  i1 <- pmin(a, pmin(b, c))
  i3 <- pmax(a, pmax(b, c))
  i2 <- a + b + c - i1 - i3
  i1 * (N^2) + i2 * N + i3
}

# All "true" valid triplets: choose 3 from disconnected Z's
make_key_valid_idx <- function(disc_idx, N) {
  if (length(disc_idx) < 3L) return(integer(0))
  C <- utils::combn(disc_idx, 3L)
  encode_triplets(C[1, ], C[2, ], C[3, ], N)
}

# Compare discovered triplets (by names) to answer keys (encoded)
test_check_int <- function(result_names, answer_keys, z_names, z_idx) {
  to_idx <- function(M, z_idx) {
    if (is.null(M) || !length(M)) return(matrix(integer(0), 3L, 0L))
    if (is.null(dim(M))) M <- matrix(M, nrow = 3L)
    out <- z_idx[M]
    dim(out) <- dim(M)
    unname(out)
  }
  R <- to_idx(result_names, z_idx)
  if (ncol(R) == 0L) return(list(tp = 0L, fp = 0L))
  
  N       <- length(z_names)
  res_key <- encode_triplets(R[1, ], R[2, ], R[3, ], N)
  idx     <- match(res_key, answer_keys, nomatch = 0L)
  tp      <- sum(idx > 0L)
  fp      <- length(res_key) - tp
  list(tp = tp, fp = fp)
}

# Run FindNegativeControls across vtalpha grid
discover_nc <- function(df, tname, oname, multiple,
                        baseline,
                        independence = TRUE,
                        alpha_corr  = 0.2,
                        use_pair    = TRUE) {
  n <- nrow(df)
  lapply(multiple, function(m) {
    vtalpha <- m / n
    FindNegativeControls(
      data         = df,
      tname        = tname,
      oname        = oname,
      baseline     = baseline,
      vtalpha      = vtalpha,
      independence = independence,
      exclude_rare = FALSE,
      alpha_corr   = alpha_corr,
      use_pair     = use_pair,
      FindBest     = FALSE
    )
  })
}

## -------------------- ROC for one sample size -----------------

sim_roc <- function(size = 1000, ntry = 200,
                    tname = "A", oname = "Y",
                    independence = TRUE,
                    alpha_corr = 0.2,
                    use_pair = TRUE) {
  
  # NEW: vtalpha grid from very strict to extremely loose
  vt_grid <- seq(0.001, 1, length.out = 12)
  L <- length(vt_grid)
  
  tp <- fp <- fn <- tn <- na <- numeric(L)
  trials_used <- 0L
  
  for (j in seq_len(ntry)) {
    set.seed(j)
    out <- generate_example_data(n = size, seed = j)
    df      <- out$data
    Z_types <- out$Z_types
    X_names <- out$X_names
    
    z_names <- Z_types$variable
    N       <- length(z_names)
    if (N < 3L) next
    
    z_idx    <- z_indexer(z_names)
    disc_nm  <- Z_types$variable[Z_types$type == "disconnected"]
    disc_idx <- unname(z_idx[disc_nm])
    answer   <- make_key_valid_idx(disc_idx, N)
    dim_ans  <- length(answer)
    if (dim_ans == 0L) next
    
    tot <- choose(N, 3L)
    
    varError <- FALSE
    res_list <- tryCatch(
      lapply(vt_grid, function(vtalpha) {
        FindNegativeControls(
          data         = df,
          tname        = tname,
          oname        = oname,
          baseline     = X_names,
          vtalpha      = vtalpha,   # use vtalpha directly
          independence = independence,
          exclude_rare = FALSE,
          alpha_corr   = alpha_corr,
          use_pair     = use_pair,
          FindBest     = FALSE
        )
      }),
      warning = function(w) { varError <<- TRUE; NULL },
      error   = function(e) { varError <<- TRUE; NULL }
    )
    if (varError || is.null(res_list)) next
    trials_used <- trials_used + 1L
    
    for (i in seq_along(vt_grid)) {
      res_mat <- res_list[[i]]
      
      if (is.null(res_mat) ||
          !length(res_mat) ||
          (is.matrix(res_mat) && ncol(res_mat) == 0L)) {
        fn[i] <- fn[i] + dim_ans
        tn[i] <- tn[i] + (tot - dim_ans)
        na[i] <- na[i] + 1
        next
      }
      
      chk <- test_check_int(res_mat, answer,
                            z_names = z_names, z_idx = z_idx)
      tp[i] <- tp[i] + chk$tp
      fp[i] <- fp[i] + chk$fp
      fn[i] <- fn[i] + dim_ans - chk$tp
      tn[i] <- tn[i] + (tot - dim_ans) - chk$fp
    }
  }
  
  tpr <- tp / pmax(tp + fn, 1)
  fpr <- fp / pmax(fp + tn, 1)
  
  list(tpr = tpr,
       fpr = fpr,
       vt_grid = vt_grid,
       trials_used = trials_used)
}

## ----------------- Multiple sample sizes + one ROC plot -------

s   <- c(100, 500, 1000, 2000, 10000)
auc <- numeric(length(s))
roc <- vector("list", length(s))

for (i in seq_along(s)) {
  cat("Running sim_roc for n =", s[i], "...\n")
  roc[[i]] <- sim_roc(size = s[i], ntry = 20)
  
  # AUC prefers FPR in ascending order
  ord <- order(roc[[i]]$fpr)
  auc[i] <- round(
    MESS::auc(roc[[i]]$fpr[ord], roc[[i]]$tpr[ord],
              rule = 2, from = 0, to = 1),
    3
  )
}

library(RColorBrewer)
col <- brewer.pal(6, "Set1")

par(mar = c(3, 3, 1, 1))
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "",
     type = "n", font.lab = 1,
     cex.lab = 0.7, cex.axis = 0.7)
title(ylab = "True positive rate", cex.lab = 0.7, mgp = c(2, 2, 0))
title(xlab = "False positive rate", cex.lab = 0.7, mgp = c(2, 2, 0))
lines(c(0, 1), c(0, 1), lwd = 2)

for (i in seq_along(s)) {
  ord <- order(roc[[i]]$fpr)
  lines(c(0, roc[[i]]$fpr[ord], 1),
        c(0, roc[[i]]$tpr[ord], 1),
        lty = i, col = col[i], lwd = 2)
  
  # highlight mid-grid thresholds if they exist
  if (length(roc[[i]]$fpr) >= 7) {
    points(roc[[i]]$fpr[6], roc[[i]]$tpr[6],
           col = col[i], lwd = 1, pch = 19)
    points(roc[[i]]$fpr[7], roc[[i]]$tpr[7],
           col = col[i], lwd = 1, pch = 1)
  }
}

legend("bottomright", inset = 0.05,
       legend = paste("n =", s, ", AUC:", sprintf("%.3f", auc)),
       lty = 1:length(s), lwd = 3, col = col,
       bty = "n", cex = 0.6, pt.cex = 0.5)

