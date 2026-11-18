# Load ONE DANCE implementation below; both are called from R.
# Update the file paths to match your environment.
# R-only: requires MASS; no C++ compiler is needed.
source("Find_Negative_Controls_R.R")

# Rcpp: requires MASS, Rcpp, RcppParallel, and a C++11-compatible compiler.
# Also update the Rcpp::sourceCpp() path inside that R file to tetrad_helpers.cpp.
# To use Rcpp, comment out the R-only source() above and uncomment the line below.
# source("Find_Negative_Controls_Rcpp.R")
# Do not source both: the second would replace FindNegativeControls().

# ---- 5. DATA GENERATION IN DANCE PAPER ----

sim.data <- function(cont, graph, n, c, d) {
  if (graph == 1 & cont == FALSE) {
    U  <- rbinom(n, 1, 0.5)
    A  <- rbinom(n, 1, plogis(-1 + c[2]*U))
    Y  <- rbinom(n, 1, plogis(-1 + c[3]*U + c[1]*A + c[4]*U*A))
    Z1 <- rbinom(n, 1, plogis(-1 + c[5]*U))
    Z2 <- rbinom(n, 1, plogis(-1 + c[6]*U + d[1]*Z1 + d[2]*U*Z1))
    Z3 <- rbinom(n, 1, plogis(-1 + c[7]*U))
    Z4 <- rbinom(n, 1, plogis(-1 + c[8]*U))
    data.frame(A, Y, Z1, Z2, Z3, Z4)
  } else if (graph == 1 & cont == TRUE) {
    U  <- rnorm(n, 0, 2)
    A  <- c[2]*U + rnorm(n)
    Y  <- c[1]*A + c[3]*U + rnorm(n)
    Z1 <- c[4]*U + rnorm(n)
    Z2 <- c[5]*U + d[1]*Z1 + rnorm(n)
    Z3 <- c[6]*U + rnorm(n)
    Z4 <- c[7]*U + rnorm(n)
    data.frame(A, Y, Z1, Z2, Z3, Z4)
  } else if (graph == 2 & cont == FALSE) {
    U  <- rbinom(n, 1, 0.5)
    A  <- rbinom(n, 1, plogis(-1 + c[2]*U))
    Y  <- rbinom(n, 1, plogis(-1 + c[3]*U + c[1]*A + c[4]*U*A))
    Z1 <- rbinom(n, 1, plogis(-1 + c[5]*U))
    Z2 <- rbinom(n, 1, plogis(-1 + c[6]*U + d[1]*Z1 + d[2]*U*Z1))
    Z3 <- rbinom(n, 1, plogis(-1 + c[7]*U))
    Z4 <- rbinom(n, 1, plogis(-1 + c[8]*U + d[3]*Z3 + d[4]*U*Z3))
    Z5 <- rbinom(n, 1, plogis(-1 + c[9]*U + d[5]*Z3 + d[6]*Z4 + d[7]*U*Z3 + d[8]*U*Z4 + d[9]*Z3*Z4 + d[10]*U*Z3*Z4))
    Z6 <- rbinom(n, 1, plogis(-1 + c[10]*U))
    Z7 <- rbinom(n, 1, plogis(-1 + c[11]*U + d[11]*Z6 + d[12]*U*Z6))
    data.frame(A, Y, Z1, Z2, Z3, Z4, Z5, Z6, Z7)
  } else {
    # graph == 2, cont == TRUE, with baseline X
    X  <- rbinom(n, 1, 0.5)
    U  <- X + rnorm(n, 0, 2)
    A  <- -1 - 2*X + c[2]*U + rnorm(n)
    Y  <- 3*X + c[1]*A + c[3]*U + rnorm(n)
    Z1 <- -X + c[4]*U + rnorm(n)
    Z2 <- -2*X + c[5]*U + d[1]*Z1 + rnorm(n)
    Z3 <- X + c[6]*U + rnorm(n)
    Z4 <- 2*X + c[7]*U + d[2]*Z3 + rnorm(n)
    Z5 <- -1.5*X + c[8]*U + d[3]*Z3 + d[4]*Z4 + rnorm(n)
    Z6 <- 1.5*X + c[9]*U + rnorm(n)
    Z7 <- c[10]*U + d[5]*Z6 + rnorm(n)
    data.frame(A, Y, X, Z1, Z2, Z3, Z4, Z5, Z6, Z7)
  }
}

answer.key <- function(graph) {
  if (graph == 1) {
    matrix(c("Z1","Z3","Z4","Z2","Z3","Z4"), 3, 2)
  } else {
    matrix(c("Z1","Z3","Z6", "Z1","Z3","Z7", "Z1","Z4","Z6", "Z1","Z4","Z7",
             "Z1","Z5","Z6", "Z1","Z5","Z7", "Z2","Z3","Z6", "Z2","Z3","Z7",
             "Z2","Z4","Z6", "Z2","Z4","Z7", "Z2","Z5","Z6", "Z2","Z5","Z7"), 3, 12)
  }
}

test.check <- function(result, answer) {
  dim1 <- dim(result)[2]
  dim2 <- dim(answer)[2]
  tp <- 0; fp <- 0
  count <- rep(0, dim2)
  if (dim1 == 0) return(list(tp = tp, fp = fp, count = count))
  for (i in 1:dim1) {
    for (j in 1:dim2) {
      if (length(setdiff(answer[, j], result[, i])) == 0) {
        tp <- tp + 1
        count[j] <- count[j] + 1
      }
    }
  }
  fp <- dim1 - tp
  list(tp = tp, fp = fp, count = count)
}

discover.nc <- function(df, tname, oname, multiple, baseline) {
  res <- vector("list", length(multiple))
  n <- nrow(df)
  for (i in seq_along(multiple)) {
    vtalpha <- multiple[i] / n
    res[[i]] <- FindNegativeControls(
      data = df, tname = tname, oname = oname,
      vtalpha = vtalpha,
      baseline = baseline,
      independence = TRUE, alpha_corr = 0.2,
      use_pair = TRUE, FindBest = FALSE
    )
  }
  res
}

sim.roc <- function(cont, size, graph, c, d, ntry) {
  multiple <- 2^((1:12) - 6)
  na <- rep(0, length(multiple))
  answer <- answer.key(graph)
  count <- rep(list(matrix(0, 3, dim(answer)[2])), length(multiple))
  tp <- fp <- fn <- tn <- rep(0, length(multiple))
  tot <- ifelse(graph == 1, 4, 35)
  noError <- 0
  
  for (j in 1:ntry) {
    set.seed(j)
    df <- sim.data(cont = cont, graph = graph, n = size, c = c, d = d)
    varError <- FALSE
    tryCatch({
      res_list <- discover.nc(df, "A", "Y", multiple, "X")
    }, warning = function(w) { varError <<- TRUE },
    error   = function(e) { varError <<- TRUE })
    
    if (!varError) {
      noError <- noError + 1
      for (i in seq_along(multiple)) {
        if (dim(res_list[[i]])[2] == 0) {
          fn[i] <- fn[i] + dim(answer)[2]
          tn[i] <- tn[i] + tot - dim(answer)[2]
          na[i] <- na[i] + 1
        } else {
          test <- test.check(res_list[[i]], answer)
          tp[i] <- tp[i] + test$tp
          fp[i] <- fp[i] + test$fp
          fn[i] <- fn[i] + dim(answer)[2] - test$tp
          tn[i] <- tn[i] + tot - dim(answer)[2] - test$fp
          count[[i]] <- count[[i]] + test$count
        }
      }
    }
  }
  
  tpr <- tp / (tp + fn)
  fpr <- fp / (fp + tn)
  na  <- round(na / noError, 2)
  
  list(tpr = tpr, fpr = fpr, na = na,
       summary = round(rbind(tpr, fpr, na), 2))
}

c.wk <- c(0.6262472, 0.3817754, 0.6795356, 0.3837477, 0.3565017, 0.4012941,
          0.3450560, 0.5648041, 0.3724304, 0.5551002, 0.4074619, 0.3746760,
          0.4442358, 0.4431560, 0.3439058, 0.5940274, 0.6612201, 0.3517061,
          0.3082211, 0.5345576)
d.wk <- c(1.413934, 1.245254, 1.270487, 1.479389, 1.652445, 1.836359,
          1.466640, 1.542527, 1.048164, 1.589371, 1.960653, 1.841541,
          1.371016, 1.488179, 1.361607, 1.816975, 1.833493, 1.337110,
          1.171592, 1.558444)
c.str <- c(0.8800869, 0.8755867, 0.8124299, 0.6633984, 0.6223731, 0.7401199,
           0.9992366, 0.7725701, 0.6719096, 0.8698893, 0.7598463, 0.6975442,
           0.7107355, 0.7997793, 0.9001751, 0.6447044, 0.6683586, 0.8351065,
           0.6506573, 0.6790947)
d.str <- c(3.749785, 3.370727, 2.306108, 2.881296, 2.347439, 3.622491,
           2.573287, 2.680017, 2.166247, 3.917405, 2.155028, 3.238999,
           3.962036, 3.788603, 3.248937, 2.461236, 2.460915, 2.980951,
           2.005269, 2.548388)
c.bin <- c(1.380936, 1.488693, 1.650510, 1.315629, 1.465796, 1.928048,
           1.597286, 1.364279, 1.049533, 1.485732, 1.064882, 1.005135,
           1.401728, 1.382683, 1.389307, 1.499075, 1.243281, 1.732134,
           1.546704, 1.167175, 1.975921, 1.160478, 1.270171, 1.082588,
           1.908118, 1.869197, 1.626949, 1.442009, 1.158431, 1.539326)
d.bin <- c(1.226158, 1.969087, 1.901274, 1.847480, 1.946074, 1.583777,
           1.314048, 1.967713, 1.434478, 1.527845, 1.155267, 1.267736,
           1.044821, 1.084514, 1.404189, 1.067947, 1.868043, 1.500284,
           1.093246, 1.784776)


# ---- 8. TRAPEZOIDAL AUC (replaces MESS::auc dependency) ----

trap.auc <- function(fpr, tpr) {
  # Prepend (0,0) and append (1,1) for a complete ROC curve
  x <- c(0, fpr, 1)
  y <- c(0, tpr, 1)
  # Sort by x (then by y for ties)
  ord <- order(x, y)
  x <- x[ord]; y <- y[ord]
  # Remove duplicates in x, keeping max y
  ux <- unique(x)
  uy <- numeric(length(ux))
  for (i in seq_along(ux)) uy[i] <- max(y[x == ux[i]])
  # Trapezoidal rule
  sum(diff(ux) * (head(uy, -1) + tail(uy, -1)) / 2)
}


# =========================================================
# PART A: ROC CURVE GENERATION
# =========================================================

run_roc_simulation <- function(cont = TRUE, graph = 2,
                               sizes = c(10, 30, 100, 300, 1000),
                               ntry = 200, cc = c.wk, dd = d.wk,
                               outfile = "roc_curve.png") {
  ns <- length(sizes)
  auc_vals <- rep(0, ns)
  roc <- vector("list", ns)
  
  cat("Running ROC simulation...\n")
  cat(sprintf("  cont=%s, graph=%d, ntry=%d\n", cont, graph, ntry))
  t0 <- proc.time()
  
  for (i in 1:ns) {
    cat(sprintf("  N = %d ...\n", sizes[i]))
    roc[[i]] <- sim.roc(cont = cont, size = sizes[i], graph = graph,
                        c = cc, d = dd, ntry = ntry)
    auc_vals[i] <- round(trap.auc(roc[[i]]$fpr, roc[[i]]$tpr), 3)
  }
  
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("Done in %.1f seconds.\n\n", elapsed))
  
  # ---- Plot ----
  cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
  
  png(outfile, width = 4, height = 3, units = "in", res = 200)
  par(mar = c(3, 3, 1, 1))
  plot(0, 0, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "", type = "n",
       font.lab = 1, cex.lab = 0.7, cex.axis = 0.7)
  title(ylab = "True positive rate", cex.lab = 0.7, mgp = c(2, 2, 0))
  title(xlab = "False positive rate", cex.lab = 0.7, mgp = c(2, 2, 0))
  lines(c(0, 1), c(0, 1), lwd = 2)
  
  for (i in 1:ns) {
    # User's plotting style: line from (1,1) through data to (0,0)
    lines(c(1, roc[[i]]$fpr, 0), c(1, roc[[i]]$tpr, 0),
          lty = i, col = cols[i], lwd = 2)
    # Solid circle at gamma = 1/n (multiplier index 6)
    points(x = roc[[i]]$fpr[6], y = roc[[i]]$tpr[6],
           col = cols[i], lwd = 1, pch = 19)
    # Open circle at gamma = 2/n (multiplier index 7)
    points(x = roc[[i]]$fpr[7], y = roc[[i]]$tpr[7],
           col = cols[i], lwd = 1, pch = 1)
  }
  
  legend("bottomright", inset = 0.05,
         legend = paste0("AUC: ", sprintf("%.3f", auc_vals)),
         lty = 1:ns, lwd = 3, col = cols[1:ns],
         bty = "n", cex = 0.6, pt.cex = 0.5)
  dev.off()
  cat(sprintf("ROC plot saved to: %s\n", outfile))
  
  # ---- Print summary table ----
  cat("\nSummary:\n")
  for (i in 1:ns) {
    cat(sprintf("  N=%5d  AUC=%.3f  ", sizes[i], auc_vals[i]))
    cat(sprintf("TPR[gamma=1/n]=%.2f  FPR[gamma=1/n]=%.2f\n",
                roc[[i]]$tpr[6], roc[[i]]$fpr[6]))
  }
  
  invisible(list(roc = roc, auc = auc_vals, sizes = sizes))
}

# =========================================================
# PART C: ATE ESTIMATION DEMO
# =========================================================

run_estimation_demo <- function(n = 3000, cc = c.wk, dd = d.wk,
                                n_mc = 50, n_mc_obs = 1000) {
  
  true_ate <- cc[1]
  cat("============================================================\n")
  cat("ATE Estimation Demo\n")
  cat(sprintf("True ATE = %.4f\n", true_ate))
  cat("============================================================\n\n")
  
  # ---- Single dataset ----
  cat(sprintf("--- Single dataset (n=%d) ---\n", n))
  set.seed(123)
  # Generate graph-2 continuous data (no baseline X)
  df <- sim.data(cont = TRUE, graph = 1, n = n, c = cc, d = dd)
  set.seed(123)
  U  <- rnorm(n, 0, 2)
  A  <- cc[2]*U + rnorm(n); Y  <- cc[1]*A + cc[3]*U + rnorm(n)
  Z1 <- cc[4]*U + rnorm(n); Z2 <- cc[5]*U + dd[1]*Z1 + rnorm(n)
  Z3 <- cc[6]*U + rnorm(n); Z4 <- cc[7]*U + dd[2]*Z3 + rnorm(n)
  Z5 <- cc[8]*U + dd[3]*Z3 + dd[4]*Z4 + rnorm(n)
  Z6 <- cc[9]*U + rnorm(n); Z7 <- cc[10]*U + dd[5]*Z6 + rnorm(n)
  df <- data.frame(A=A, Y=Y, Z1=Z1, Z2=Z2, Z3=Z3,
                   Z4=Z4, Z5=Z5, Z6=Z6, Z7=Z7)
  
  naive <- coef(lm(Y ~ A, data = df))["A"]
  cat(sprintf("Naive regression ATE: %.4f  (bias=%.4f, %.1f%%)\n",
              naive, naive - true_ate, 100*(naive - true_ate)/true_ate))
  
  ate_est <- effect.est(df, "A", "Y")
  cat(sprintf("DANCE ATE:            %.4f  (bias=%.4f, %.1f%%)\n\n",
              ate_est, ate_est - true_ate, 100*(ate_est - true_ate)/true_ate))
  
  # ---- Monte Carlo ----
  cat(sprintf("--- Monte Carlo (%d reps, n=%d each) ---\n", n_mc, n_mc_obs))
  naive_ates <- dance_ates <- numeric(n_mc)
  dance_ok <- logical(n_mc)
  
  for (r in 1:n_mc) {
    set.seed(r)
    U  <- rnorm(n_mc_obs, 0, 2)
    A  <- cc[2]*U + rnorm(n_mc_obs); Y <- cc[1]*A + cc[3]*U + rnorm(n_mc_obs)
    Z1 <- cc[4]*U + rnorm(n_mc_obs)
    Z2 <- cc[5]*U + dd[1]*Z1 + rnorm(n_mc_obs)
    Z3 <- cc[6]*U + rnorm(n_mc_obs)
    Z4 <- cc[7]*U + dd[2]*Z3 + rnorm(n_mc_obs)
    Z5 <- cc[8]*U + dd[3]*Z3 + dd[4]*Z4 + rnorm(n_mc_obs)
    Z6 <- cc[9]*U + rnorm(n_mc_obs)
    Z7 <- cc[10]*U + dd[5]*Z6 + rnorm(n_mc_obs)
    dfi <- data.frame(A=A, Y=Y, Z1=Z1, Z2=Z2, Z3=Z3,
                      Z4=Z4, Z5=Z5, Z6=Z6, Z7=Z7)
    
    naive_ates[r] <- coef(lm(Y ~ A, data = dfi))["A"]
    ate_r <- tryCatch(
      suppressWarnings(effect.est(dfi, "A", "Y")),
      error = function(e) NA_real_
    )
    dance_ates[r] <- ate_r
    dance_ok[r]   <- !is.na(ate_r)
    if (r %% 10 == 0) cat(sprintf("  rep %d/%d\n", r, n_mc))
  }
  
  ok <- dance_ok
  cat(sprintf("\nResults (true ATE = %.3f):\n", true_ate))
  cat(sprintf("  Naive:  mean=%.4f  bias=%.4f (%.1f%%)  SD=%.4f\n",
              mean(naive_ates), mean(naive_ates) - true_ate,
              100*(mean(naive_ates) - true_ate)/true_ate, sd(naive_ates)))
  if (sum(ok) > 0) {
    cat(sprintf("  DANCE:  mean=%.4f  bias=%.4f (%.1f%%)  SD=%.4f  (%d/%d succeeded)\n",
                mean(dance_ates[ok]), mean(dance_ates[ok]) - true_ate,
                100*(mean(dance_ates[ok]) - true_ate)/true_ate,
                sd(dance_ates[ok]), sum(ok), n_mc))
  } else {
    cat("  DANCE:  all replications failed\n")
  }
}


# =========================================================
# EXECUTION
# =========================================================

# ---- A) ROC 
run_roc_simulation(cont = TRUE, graph = 2,
                   sizes = c(30, 50, 100, 500, 1000),
                   ntry = 200, cc = c.str, dd = d.str,
                   outfile = "roc_user_style.png")

# ---- B) ATE estimation

run_estimation_demo(n = 3000, cc = c.wk, dd = d.wk,
                    n_mc = 50, n_mc_obs = 1000)

