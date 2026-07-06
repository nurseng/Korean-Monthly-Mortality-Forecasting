############################################################
# CLEAN THESIS-STYLE FINAL R CODE
# data.csv + 3 selected series
# Series:
# y1 = row 1   -> C00-D48
# y2 = row 31  -> J00-J98,U04
# y3 = row 58  -> V01-Y89,U12
#
# IMPORTANT:
# - This version is CLEAN and deterministic.
# - Bootstrap SE is NOT used.
# - Standard errors are Hessian-based.
# - NBINGARCH SE is calculated ONLY for the selected best r.
# - Parameter values and standard errors are formatted separately.
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
# 0) SETTINGS
############################################################

file_path <- "/Users/nursengokbulut/Desktop/data.csv"
out_dir   <- "tez_final_output_clean"

dir.create(out_dir, showWarnings = FALSE)

# out-of-sample size
out_n <- 50

# NBINGARCH: positive integer r candidates
r_grid <- 1:150

EPS <- 1e-8

############################################################
# 1) PACKAGES
############################################################

if (!requireNamespace("numDeriv", quietly = TRUE)) install.packages("numDeriv")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

############################################################
# 2) READ DATA
############################################################

dat <- read.csv(
  file_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  fileEncoding = "CP949"
)

time_cols <- grep("^\\d{4}\\.\\d{2}", names(dat), value = TRUE)

if (length(time_cols) == 0) {
  stop("Monthly time columns could not be found in data.csv")
}

# Selected 3 series
y1 <- as.numeric(dat[1,  time_cols])   # C00-D48
y2 <- as.numeric(dat[31, time_cols])   # J00-J98,U04
y3 <- as.numeric(dat[58, time_cols])   # V01-Y89,U12

############################################################
# 3) HELPER FUNCTIONS
############################################################

# Parameter estimates:
# beta itself may be exactly 0, so parameter formatting can show 0.0000.
fmt_par <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = 4)
}

# Standard errors:
# SE should not be rounded to 0.0000. Very small SE values are shown scientifically.
fmt_se <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  if (abs(x) < 1e-4) {
    return(formatC(x, format = "e", digits = 3))
  }
  formatC(x, format = "f", digits = 4)
}

# General number formatter for equations
fmt_eq <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = 4)
}

# Safe Hessian-based standard error calculation
safe_se <- function(loglik_fun, par_hat) {
  H <- tryCatch(
    optimHess(par_hat, fn = loglik_fun),
    error = function(e) NULL
  )
  
  if (is.null(H) || any(!is.finite(H))) {
    H <- tryCatch(
      numDeriv::hessian(func = loglik_fun, x = par_hat),
      error = function(e) NULL
    )
  }
  
  if (is.null(H) || any(!is.finite(H))) {
    return(rep(NA_real_, length(par_hat)))
  }
  
  vcov_mat <- tryCatch(
    solve(H),
    error = function(e) tryCatch(MASS::ginv(H), error = function(e2) NULL)
  )
  
  if (is.null(vcov_mat) || any(!is.finite(vcov_mat))) {
    return(rep(NA_real_, length(par_hat)))
  }
  
  diag_vals <- diag(vcov_mat)
  
  # If variance is negative, the Hessian-based SE is not reliable.
  # We return NA rather than forcing negative/invalid values to 0.
  if (any(diag_vals <= 0)) {
    return(rep(NA_real_, length(par_hat)))
  }
  
  se <- sqrt(diag_vals)
  se[!is.finite(se)] <- NA_real_
  se
}

model_name_from_beta <- function(beta, base1, base0) {
  if (is.na(beta)) return(base1)
  if (abs(beta) < 1e-6) base0 else base1
}

############################################################
# 4) INGARCH(1,1)
# X_t | F_{t-1} ~ Poisson(lambda_t)
# lambda_t = omega + alpha X_{t-1} + beta lambda_{t-1}
############################################################

lambda_ing <- function(x, par) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  n <- length(x)
  lam <- numeric(n)
  lam[1] <- max(mean(x), EPS)
  
  if (n >= 2) {
    for (t in 2:n) {
      lam[t] <- omega + alpha * x[t - 1] + beta * lam[t - 1]
      lam[t] <- max(lam[t], EPS)
    }
  }
  lam
}

nll_ing <- function(par, x) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  if (omega <= 0 || alpha < 0 || beta < 0 || alpha + beta >= 1) {
    return(1e12)
  }
  
  lam <- lambda_ing(x, par)
  -sum(dpois(x, lambda = lam, log = TRUE))
}

fit_ing <- function(x, compute_se = TRUE) {
  init <- c(max(mean(x) * 0.1, 0.001), 0.3, 0.3)
  
  fit <- optim(
    par = init,
    fn = nll_ing,
    x = x,
    method = "L-BFGS-B",
    lower = c(EPS, 0, 0),
    upper = c(Inf, 0.999, 0.999),
    hessian = TRUE
  )
  
  par_hat <- fit$par
  
  se_hat <- if (compute_se) {
    safe_se(
      loglik_fun = function(th) nll_ing(th, x),
      par_hat = par_hat
    )
  } else {
    rep(NA_real_, length(par_hat))
  }
  
  list(
    par = par_hat,
    se  = se_hat,
    loglik = -fit$value,
    AIC = -2 * (-fit$value) + 2 * 3,
    convergence = fit$convergence,
    message = fit$message
  )
}

############################################################
# 5) NBINGARCH(1,1)
# X_t | F_{t-1} ~ NB(r, p_t)
# lambda_t = (1-p_t)/p_t = omega + alpha X_{t-1} + beta lambda_{t-1}
# E(X_t|F_{t-1}) = r * lambda_t
############################################################

lambda_nb <- function(x, par, r) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  n <- length(x)
  lam <- numeric(n)
  lam[1] <- max(mean(x) / max(r, 1), EPS)
  
  if (n >= 2) {
    for (t in 2:n) {
      lam[t] <- omega + alpha * x[t - 1] + beta * lam[t - 1]
      lam[t] <- max(lam[t], EPS)
    }
  }
  lam
}

nll_nb_fixed_r <- function(par, x, r) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  if (omega <= 0 || alpha < 0 || beta < 0) {
    return(1e12)
  }
  
  if (((r * alpha + beta)^2 + r * alpha^2) >= 1) {
    return(1e12)
  }
  
  lam <- lambda_nb(x, par, r)
  mu  <- r * lam
  
  -sum(dnbinom(x, size = r, mu = mu, log = TRUE))
}

# This function fits a fixed r candidate.
# SE is NOT computed here because this function is run many times over r_grid.
fit_nb_fixed_r <- function(x, r, compute_se = FALSE) {
  init <- c(max(mean(x) / (10 * max(r, 1)), 0.001), 0.05, 0.10)
  
  fit <- optim(
    par = init,
    fn = nll_nb_fixed_r,
    x = x,
    r = r,
    method = "L-BFGS-B",
    lower = c(EPS, 0, 0),
    upper = c(Inf, 0.999, 0.999),
    hessian = TRUE
  )
  
  par_hat <- fit$par
  
  se_hat <- if (compute_se) {
    safe_se(
      loglik_fun = function(th) nll_nb_fixed_r(th, x, r),
      par_hat = par_hat
    )
  } else {
    rep(NA_real_, length(par_hat))
  }
  
  list(
    r = r,
    par = par_hat,
    se  = se_hat,
    loglik = -fit$value,
    AIC = -2 * (-fit$value) + 2 * 4,
    convergence = fit$convergence,
    message = fit$message
  )
}

# Select best r by AIC, then compute SE only for that selected model
fit_nb <- function(x, r_candidates = r_grid, compute_se = TRUE) {
  all_fit <- lapply(r_candidates, function(rr) fit_nb_fixed_r(x, rr, compute_se = FALSE))
  aics <- sapply(all_fit, function(z) z$AIC)
  best_idx <- which.min(aics)
  
  best_fit <- all_fit[[best_idx]]
  
  if (compute_se) {
    best_fit$se <- safe_se(
      loglik_fun = function(th) nll_nb_fixed_r(th, x, best_fit$r),
      par_hat = best_fit$par
    )
  }
  
  best_fit
}

############################################################
# 6) LOG-LINEAR POISSON AR
# X_t | F_{t-1} ~ Poisson(lambda_t)
# v_t = log(lambda_t)
# v_t = omega + alpha log(X_{t-1}+1) + beta v_{t-1}
############################################################

v_loglin <- function(x, par) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  n <- length(x)
  v <- numeric(n)
  v[1] <- log(max(mean(x), EPS))
  
  if (n >= 2) {
    for (t in 2:n) {
      v[t] <- omega + alpha * log(x[t - 1] + 1) + beta * v[t - 1]
    }
  }
  v
}

lambda_loglin <- function(x, par) {
  exp(v_loglin(x, par))
}

nll_loglin <- function(par, x) {
  omega <- par[1]
  alpha <- par[2]
  beta  <- par[3]
  
  same_sign <- (alpha >= 0 && beta >= 0) || (alpha <= 0 && beta <= 0)
  
  ok <- if (same_sign) {
    abs(alpha + beta) < 1
  } else {
    (alpha^2 + beta^2) < 1
  }
  
  if (!ok) return(1e12)
  
  lam <- lambda_loglin(x, par)
  lam <- pmax(lam, EPS)
  
  -sum(dpois(x, lambda = lam, log = TRUE))
}

fit_loglin <- function(x, compute_se = TRUE) {
  init <- c(log(mean(x) + 1) * 0.1, 0.3, 0.2)
  
  fit <- optim(
    par = init,
    fn = nll_loglin,
    x = x,
    method = "BFGS",
    hessian = TRUE
  )
  
  par_hat <- fit$par
  
  se_hat <- if (compute_se) {
    safe_se(
      loglik_fun = function(th) nll_loglin(th, x),
      par_hat = par_hat
    )
  } else {
    rep(NA_real_, length(par_hat))
  }
  
  list(
    par = par_hat,
    se  = se_hat,
    loglik = -fit$value,
    AIC = -2 * (-fit$value) + 2 * 3,
    convergence = fit$convergence,
    message = fit$message
  )
}

############################################################
# 7) FIT ALL 3 SERIES
############################################################

ing1 <- fit_ing(y1, compute_se = TRUE)
nb1  <- fit_nb(y1, compute_se = TRUE)
log1 <- fit_loglin(y1, compute_se = TRUE)

ing2 <- fit_ing(y2, compute_se = TRUE)
nb2  <- fit_nb(y2, compute_se = TRUE)
log2 <- fit_loglin(y2, compute_se = TRUE)

ing3 <- fit_ing(y3, compute_se = TRUE)
nb3  <- fit_nb(y3, compute_se = TRUE)
log3 <- fit_loglin(y3, compute_se = TRUE)

############################################################
# BOOTSTRAP STANDARD ERRORS
############################################################

set.seed(123)

boot_se_ing <- function(x, fit_obj, B = 500) {
  par_hat <- fit_obj$par
  lam_hat <- lambda_ing(x, par_hat)
  n <- length(x)
  boot_par <- matrix(NA, nrow = B, ncol = 3)
  
  for (b in 1:B) {
    xb <- rpois(n, lambda = lam_hat)
    fb <- tryCatch(fit_ing(xb, compute_se = FALSE), error = function(e) NULL)
    if (!is.null(fb)) boot_par[b, ] <- fb$par
  }
  
  apply(boot_par, 2, sd, na.rm = TRUE)
}

boot_se_nb <- function(x, fit_obj, B = 500) {
  par_hat <- fit_obj$par
  r_hat <- fit_obj$r
  lam_hat <- lambda_nb(x, par_hat, r_hat)
  mu_hat <- r_hat * lam_hat
  n <- length(x)
  boot_par <- matrix(NA, nrow = B, ncol = 3)
  
  for (b in 1:B) {
    xb <- rnbinom(n, size = r_hat, mu = mu_hat)
    fb <- tryCatch(fit_nb_fixed_r(xb, r_hat, compute_se = FALSE), error = function(e) NULL)
    if (!is.null(fb)) boot_par[b, ] <- fb$par
  }
  
  apply(boot_par, 2, sd, na.rm = TRUE)
}

boot_se_loglin <- function(x, fit_obj, B = 500) {
  par_hat <- fit_obj$par
  lam_hat <- lambda_loglin(x, par_hat)
  n <- length(x)
  boot_par <- matrix(NA, nrow = B, ncol = 3)
  
  for (b in 1:B) {
    xb <- rpois(n, lambda = lam_hat)
    fb <- tryCatch(fit_loglin(xb, compute_se = FALSE), error = function(e) NULL)
    if (!is.null(fb)) boot_par[b, ] <- fb$par
  }
  
  apply(boot_par, 2, sd, na.rm = TRUE)
}

ing1$se <- boot_se_ing(y1, ing1, B = 1000)
nb1$se  <- boot_se_nb(y1, nb1, B = 500)
log1$se <- boot_se_loglin(y1, log1, B = 500)

ing2$se <- boot_se_ing(y2, ing2, B = 500)
nb2$se  <- boot_se_nb(y2, nb2, B = 500)
log2$se <- boot_se_loglin(y2, log2, B = 500)

ing3$se <- boot_se_ing(y3, ing3, B = 500)
nb3$se  <- boot_se_nb(y3, nb3, B = 500)
log3$se <- boot_se_loglin(y3, log3, B = 500)

# For C00-D48 INGARCH, use Hessian-based SE because bootstrap beta SE became numerically zero
ing1$se <- safe_se(
  loglik_fun = function(th) nll_ing(th, y1),
  par_hat = ing1$par
)

cat("\nStandard error check (Hessian-based)\n")
cat("ing1 se:", ing1$se, "\n")
cat("nb1  se:", nb1$se,  "\n")
cat("log1 se:", log1$se, "\n")
cat("ing2 se:", ing2$se, "\n")
cat("nb2  se:", nb2$se,  "\n")
cat("log2 se:", log2$se, "\n")
cat("ing3 se:", ing3$se, "\n")
cat("nb3  se:", nb3$se,  "\n")
cat("log3 se:", log3$se, "\n")

############################################################
# 8) FIGURE 1-3: SERIES + ACF + PACF
############################################################

save_fig123 <- function(x, file_name) {
  png(file.path(out_dir, file_name), width = 900, height = 1400, res = 150)
  par(mfrow = c(3, 1), mar = c(3.5, 4.2, 2.2, 1.2))
  
  plot(x, type = "l", col = "black", lwd = 1, main = "(a)", xlab = "", ylab = "")
  acf(x, lag.max = 21, main = "(b)", xlab = "Lag", ylab = "ACF")
  pacf(x, lag.max = 21, main = "(c)", xlab = "Lag", ylab = "Partial ACF")
  
  dev.off()
}

save_fig123(y1, "Figure_1_C00_D48.png")
save_fig123(y2, "Figure_2_J00_J98_U04.png")
save_fig123(y3, "Figure_3_V01_Y89_U12.png")

############################################################
# 9) TABLE 1-3: PARAMETER ESTIMATES + EQUATIONS
############################################################

save_table123 <- function(ing, nb, logm, table_no, title_text, file_name) {
  png(file.path(out_dir, file_name), width = 1400, height = 1700, res = 180)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-0.18, 1))
  
  text(0.5, 0.95, paste0("Table ", table_no, ": Parameter estimates for ", title_text), cex = 1.35)
  
  # table lines
  segments(0.08, 0.88, 0.92, 0.88, lwd = 1)
  segments(0.08, 0.82, 0.92, 0.82, lwd = 1)
  segments(0.08, 0.67, 0.92, 0.67, lwd = 1)
  segments(0.08, 0.52, 0.92, 0.52, lwd = 1)
  segments(0.08, 0.37, 0.92, 0.37, lwd = 1)
  segments(0.30, 0.37, 0.30, 0.88, lwd = 1)
  
  text(0.16, 0.85, "Model", cex = 1.15, font = 2)
  text(0.42, 0.85, expression(omega), cex = 1.15)
  text(0.56, 0.85, expression(alpha), cex = 1.15)
  text(0.70, 0.85, expression(beta), cex = 1.15)
  text(0.84, 0.85, "AIC", cex = 1.15, font = 2)
  
  text(0.16, 0.74, "INGARCH(1,1)", cex = 1.1)
  text(0.42, 0.74, paste0(fmt_par(ing$par[1]), "\n(", fmt_se(ing$se[1]), ")"), cex = 1.0)
  text(0.56, 0.74, paste0(fmt_par(ing$par[2]), "\n(", fmt_se(ing$se[2]), ")"), cex = 1.0)
  text(0.70, 0.74, paste0(fmt_par(ing$par[3]), "\n(", fmt_se(ing$se[3]), ")"), cex = 1.0)
  text(0.84, 0.74, format(round(ing$AIC, 3), nsmall = 3), cex = 1.0)
  
  text(0.16, 0.59, "NBINGARCH(1,1)", cex = 1.1)
  text(0.16, 0.545, bquote(hat(r) == .(nb$r)), cex = 1.0)
  text(0.42, 0.59, paste0(fmt_par(nb$par[1]), "\n(", fmt_se(nb$se[1]), ")"), cex = 1.0)
  text(0.56, 0.59, paste0(fmt_par(nb$par[2]), "\n(", fmt_se(nb$se[2]), ")"), cex = 1.0)
  text(0.70, 0.59, paste0(fmt_par(nb$par[3]), "\n(", fmt_se(nb$se[3]), ")"), cex = 1.0)
  text(0.84, 0.59, format(round(nb$AIC, 3), nsmall = 3), cex = 1.0)
  
  text(0.16, 0.44, "Log-linear", cex = 1.1)
  text(0.42, 0.44, paste0(fmt_par(logm$par[1]), "\n(", fmt_se(logm$se[1]), ")"), cex = 1.0)
  text(0.56, 0.44, paste0(fmt_par(logm$par[2]), "\n(", fmt_se(logm$se[2]), ")"), cex = 1.0)
  text(0.70, 0.44, paste0(fmt_par(logm$par[3]), "\n(", fmt_se(logm$se[3]), ")"), cex = 1.0)
  text(0.84, 0.44, format(round(logm$AIC, 3), nsmall = 3), cex = 1.0)
  
  text(0.5, 0.28, "Estimated equations", cex = 1.15)
  
  ing_label <- model_name_from_beta(ing$par[3], "INGARCH(1,1)", "INARCH(1)")
  nb_label  <- model_name_from_beta(nb$par[3],  "NBINGARCH(1,1)", "NBINARCH(1)")
  
  # INGARCH block
  text(0.10, 0.20, paste0(ing_label, ";    X[t] | F[t-1] : Poisson(lambda[t])"),
       adj = 0, cex = 1.00)
  if (abs(ing$par[3]) < 1e-6) {
    text(0.24, 0.15,
         paste0("lambda[t] = ", fmt_eq(ing$par[1]), " + ", fmt_eq(ing$par[2]), " X[t-1]"),
         adj = 0, cex = 1.00)
  } else {
    text(0.24, 0.15,
         paste0("lambda[t] = ", fmt_eq(ing$par[1]), " + ", fmt_eq(ing$par[2]),
                " X[t-1] + ", fmt_eq(ing$par[3]), " lambda[t-1]"),
         adj = 0, cex = 1.00)
  }
  
  # NBINGARCH block
  text(0.10, 0.06, paste0(nb_label, ";    X[t] | F[t-1] : NB(r, p[t])"),
       adj = 0, cex = 1.00)
  if (abs(nb$par[3]) < 1e-6) {
    text(0.24, 0.01,
         paste0("lambda[t] = (1-p[t])/p[t] = ", fmt_eq(nb$par[1]), " + ", fmt_eq(nb$par[2]), " X[t-1]"),
         adj = 0, cex = 1.00)
  } else {
    text(0.24, 0.01,
         paste0("lambda[t] = (1-p[t])/p[t] = ", fmt_eq(nb$par[1]), " + ", fmt_eq(nb$par[2]),
                " X[t-1] + ", fmt_eq(nb$par[3]), " lambda[t-1]"),
         adj = 0, cex = 1.00)
  }
  
  # Log-linear block
  text(0.10, -0.08, "Log-linear;    X[t] | F[t-1] : Poisson(lambda[t])",
       adj = 0, cex = 1.00, xpd = TRUE)
  text(0.24, -0.13,
       paste0("v[t] = ", fmt_eq(logm$par[1]), " + ", fmt_eq(logm$par[2]),
              " log(X[t-1]+1) + ", fmt_eq(logm$par[3]), " v[t-1]"),
       adj = 0, cex = 1.00, xpd = TRUE)
  
  dev.off()
  
  out_tab <- data.frame(
    Model = c("INGARCH(1,1)", "NBINGARCH(1,1)", "Log-linear"),
    r_hat = c(NA, nb$r, NA),
    omega = c(ing$par[1], nb$par[1], logm$par[1]),
    alpha = c(ing$par[2], nb$par[2], logm$par[2]),
    beta  = c(ing$par[3], nb$par[3], logm$par[3]),
    omega_se = c(ing$se[1], nb$se[1], logm$se[1]),
    alpha_se = c(ing$se[2], nb$se[2], logm$se[2]),
    beta_se  = c(ing$se[3], nb$se[3], logm$se[3]),
    AIC = c(ing$AIC, nb$AIC, logm$AIC)
  )
  
  write.csv(
    out_tab,
    file.path(out_dir, paste0("Table_", table_no, "_numeric.csv")),
    row.names = FALSE
  )
}

save_table123(ing1, nb1, log1, 1, "C00-D48 data",     "Table_1_C00_D48.png")
save_table123(ing2, nb2, log2, 2, "J00-J98,U04 data", "Table_2_J00_J98_U04.png")
save_table123(ing3, nb3, log3, 3, "V01-Y89,U12 data", "Table_3_V01_Y89_U12.png")

############################################################
# 10) ROLLING ONE-STEP-AHEAD FORECAST
############################################################

rolling_forecast <- function(x, out_n = 50) {
  n <- length(x)
  T0 <- n - out_n
  
  actual    <- x[(T0 + 1):n]
  pred_ing  <- numeric(out_n)
  pred_nb   <- numeric(out_n)
  pred_log  <- numeric(out_n)
  
  for (i in 1:out_n) {
    train <- x[1:(T0 + i - 1)]
    
    fit_i <- fit_ing(train, compute_se = FALSE)
    fit_n <- fit_nb(train, compute_se = FALSE)
    fit_l <- fit_loglin(train, compute_se = FALSE)
    
    lam_i <- lambda_ing(train, fit_i$par)
    lam_n <- lambda_nb(train, fit_n$par, fit_n$r)
    v_l   <- v_loglin(train, fit_l$par)
    
    # INGARCH forecast
    pred_ing[i] <- fit_i$par[1] +
      fit_i$par[2] * train[length(train)] +
      fit_i$par[3] * lam_i[length(lam_i)]
    
    # NBINGARCH forecast: actual mean = r * lambda
    lam_next_nb <- fit_n$par[1] +
      fit_n$par[2] * train[length(train)] +
      fit_n$par[3] * lam_n[length(lam_n)]
    pred_nb[i] <- fit_n$r * lam_next_nb
    
    # Log-linear forecast
    v_next_log <- fit_l$par[1] +
      fit_l$par[2] * log(train[length(train)] + 1) +
      fit_l$par[3] * v_l[length(v_l)]
    pred_log[i] <- exp(v_next_log)
  }
  
  pred_ing <- pmax(pred_ing, EPS)
  pred_nb  <- pmax(pred_nb, EPS)
  pred_log <- pmax(pred_log, EPS)
  
  e_ing <- actual - pred_ing
  e_nb  <- actual - pred_nb
  e_log <- actual - pred_log
  
  list(
    actual   = actual,
    pred_ing = pred_ing,
    pred_nb  = pred_nb,
    pred_log = pred_log,
    mse_ing = mean(e_ing^2),
    mae_ing = mean(abs(e_ing)),
    mse_nb  = mean(e_nb^2),
    mae_nb  = mean(abs(e_nb)),
    mse_log = mean(e_log^2),
    mae_log = mean(abs(e_log))
  )
}

fc1 <- rolling_forecast(y1, out_n = out_n)
fc2 <- rolling_forecast(y2, out_n = out_n)
fc3 <- rolling_forecast(y3, out_n = out_n)

############################################################
# 11) TABLE 4-6: FORECASTING RESULTS
############################################################

save_table456 <- function(fc, table_no, title_text, file_name) {
  png(file.path(out_dir, file_name), width = 1800, height = 1000, res = 220)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  text(0.5, 0.90,
       paste0("Table ", table_no, ": Forecasting results for ", title_text),
       cex = 1.5)
  
  # table frame
  segments(0.14, 0.76, 0.90, 0.76, lwd = 1)
  segments(0.14, 0.68, 0.90, 0.68, lwd = 1)
  segments(0.14, 0.56, 0.90, 0.56, lwd = 1)
  segments(0.14, 0.44, 0.90, 0.44, lwd = 1)
  segments(0.14, 0.32, 0.90, 0.32, lwd = 1)
  
  segments(0.48, 0.32, 0.48, 0.76, lwd = 1)
  segments(0.70, 0.32, 0.70, 0.76, lwd = 1)
  
  # headers
  text(0.59, 0.72, "MSE", cex = 1.25, font = 2)
  text(0.80, 0.72, "MAE", cex = 1.25, font = 2)
  
  # rows
  text(0.31, 0.62, "INGARCH(1,1)", cex = 1.15)
  text(0.59, 0.62, sprintf("%.4f", fc$mse_ing), cex = 1.05)
  text(0.80, 0.62, sprintf("%.4f", fc$mae_ing), cex = 1.05)
  
  text(0.31, 0.50, "NBINGARCH(1,1)", cex = 1.15)
  text(0.59, 0.50, sprintf("%.4f", fc$mse_nb), cex = 1.05)
  text(0.80, 0.50, sprintf("%.4f", fc$mae_nb), cex = 1.05)
  
  text(0.31, 0.38, "Log-linear", cex = 1.15)
  text(0.59, 0.38, sprintf("%.4f", fc$mse_log), cex = 1.05)
  text(0.80, 0.38, sprintf("%.4f", fc$mae_log), cex = 1.05)
  
  dev.off()
  
  out_tab <- data.frame(
    Model = c("INGARCH(1,1)", "NBINGARCH(1,1)", "Log-linear"),
    MSE   = c(fc$mse_ing, fc$mse_nb, fc$mse_log),
    MAE   = c(fc$mae_ing, fc$mae_nb, fc$mae_log)
  )
  
  write.csv(
    out_tab,
    file.path(out_dir, paste0("Table_", table_no, "_numeric.csv")),
    row.names = FALSE
  )
}

save_table456(fc1, 4, "C00-D48 data",     "Table_4_C00_D48.png")
save_table456(fc2, 5, "J00-J98,U04 data", "Table_5_J00_J98_U04.png")
save_table456(fc3, 6, "V01-Y89,U12 data", "Table_6_V01_Y89_U12.png")

############################################################
# 12) FIGURE 4-6: OBSERVED vs PREDICTED
############################################################

save_fig456 <- function(fc, file_name) {
  png(file.path(out_dir, file_name), width = 900, height = 1700, res = 150)
  par(mfrow = c(3, 1), mar = c(4, 4.2, 2.2, 1.2))
  
  yr <- range(c(fc$actual, fc$pred_ing, fc$pred_nb, fc$pred_log), na.rm = TRUE)
  
  plot(fc$actual, type = "l", col = "black", lwd = 1,
       ylim = yr, main = "(a)", xlab = "time", ylab = "data")
  lines(fc$pred_ing, col = "deeppink", lwd = 1)
  legend("bottomleft",
         legend = c("observed", "predicted"),
         col = c("black", "deeppink"),
         lty = 1, lwd = 1, bty = "o", cex = 0.9)
  
  plot(fc$actual, type = "l", col = "black", lwd = 1,
       ylim = yr, main = "(b)", xlab = "time", ylab = "data")
  lines(fc$pred_nb, col = "deeppink", lwd = 1)
  legend("bottomleft",
         legend = c("observed", "predicted"),
         col = c("black", "deeppink"),
         lty = 1, lwd = 1, bty = "o", cex = 0.9)
  
  plot(fc$actual, type = "l", col = "black", lwd = 1,
       ylim = yr, main = "(c)", xlab = "time", ylab = "data")
  lines(fc$pred_log, col = "deeppink", lwd = 1)
  legend("bottomleft",
         legend = c("observed", "predicted"),
         col = c("black", "deeppink"),
         lty = 1, lwd = 1, bty = "o", cex = 0.9)
  
  dev.off()
}

save_fig456(fc1, "Figure_4_C00_D48.png")
save_fig456(fc2, "Figure_5_J00_J98_U04.png")
save_fig456(fc3, "Figure_6_V01_Y89_U12.png")

############################################################
# 13) OPTIONAL: MEAN / VARIANCE
############################################################

cat("\nMean and variance of the three series\n")
cat("C00-D48     : mean =", mean(y1), ", variance =", var(y1), "\n")
cat("J00-J98,U04 : mean =", mean(y2), ", variance =", var(y2), "\n")
cat("V01-Y89,U12 : mean =", mean(y3), ", variance =", var(y3), "\n")

############################################################
# 14) CONSOLE SUMMARY
############################################################

cat("\nAll thesis-style outputs were saved in folder:", out_dir, "\n")

cat("\nAIC summary\n")
cat("C00-D48     :", "INGARCH =", round(ing1$AIC, 3),
    ", NBINGARCH =", round(nb1$AIC, 3),
    ", Log-linear =", round(log1$AIC, 3), "\n")

cat("J00-J98,U04 :", "INGARCH =", round(ing2$AIC, 3),
    ", NBINGARCH =", round(nb2$AIC, 3),
    ", Log-linear =", round(log2$AIC, 3), "\n")

cat("V01-Y89,U12 :", "INGARCH =", round(ing3$AIC, 3),
    ", NBINGARCH =", round(nb3$AIC, 3),
    ", Log-linear =", round(log3$AIC, 3), "\n")

cat("\nForecast summary\n")
cat("C00-D48     :",
    "INGARCH MSE =", round(fc1$mse_ing, 4), ", MAE =", round(fc1$mae_ing, 4),
    " | NBINGARCH MSE =", round(fc1$mse_nb, 4), ", MAE =", round(fc1$mae_nb, 4),
    " | Log-linear MSE =", round(fc1$mse_log, 4), ", MAE =", round(fc1$mae_log, 4), "\n")

cat("J00-J98,U04 :",
    "INGARCH MSE =", round(fc2$mse_ing, 4), ", MAE =", round(fc2$mae_ing, 4),
    " | NBINGARCH MSE =", round(fc2$mse_nb, 4), ", MAE =", round(fc2$mae_nb, 4),
    " | Log-linear MSE =", round(fc2$mse_log, 4), ", MAE =", round(fc2$mae_log, 4), "\n")

cat("V01-Y89,U12 :",
    "INGARCH MSE =", round(fc3$mse_ing, 4), ", MAE =", round(fc3$mae_ing, 4),
    " | NBINGARCH MSE =", round(fc3$mse_nb, 4), ", MAE =", round(fc3$mae_nb, 4),
    " | Log-linear MSE =", round(fc3$mse_log, 4), ", MAE =", round(fc3$mae_log, 4), "\n")

