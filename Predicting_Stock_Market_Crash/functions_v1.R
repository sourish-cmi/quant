gk_volatility <- function(open, high, low, close) {
  u <- log(high / open)
  d <- log(low / open)
  r <- log(close / open)
  sigma2 <- 0.5 * (u - d)^2 - (2 * log(2) - 1) * r^2
  return(sigma2)
}


compute_aic_bic <- function(par_hat, y, X, mu, neg_log_lik_fn) {
  # Compute negative log-likelihood (without prior)
  nll_value <- neg_log_lik_fn(par = par_hat, y = y, X = X, mu = mu)
  
  # Number of estimated parameters
  k <- length(par_hat)
  
  # Number of observations
  n <- length(y)
  
  # AIC and BIC
  aic <- 2 * k + 2 * nll_value
  bic <- log(n) * k + 2 * nll_value
  
  return(list(AIC = aic, BIC = bic))
}

# Function to simulate from GPD regression
simulate_gpd_regression <- function(n, beta, mu = 5, xi = 0.1, X = NULL, sigma = 1, rho = 0.5) {
  library(evd)       # for rgpd
  library(mvtnorm)   # for rmvnorm
  
  p <- length(beta)
  
  # Covariance matrix with AR(1) structure
  Sigma <- matrix(NA, p, p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i, j] <- sigma * rho^(abs(i - j))
    }
  }
  
  mu_X <- rep(0, p)  # mean of X (not GPD threshold!)
  
  # Generate design matrix if not provided
  if (is.null(X)) {
    X <- rmvnorm(n, mean = mu_X, sigma = Sigma)
  }
  
  # Compute sigma_i for each row
  log_sigma <- X %*% beta
  sigma_i <- as.vector(exp(log_sigma))
  if (any(sigma_i <= 0)) stop("Non-positive sigma encountered!")
  
  # Generate responses from GPD(mu, sigma_i, xi)
  y <- numeric(n)
  for (i in 1:n) {
    y[i] <- rgpd(1, loc = mu, scale = sigma_i[i], shape = xi)
  }
  
  return(list(y = y, X = X, sigma = sigma_i))
}

## Jacobi prior bases estimation
estimate_xi_traunc_cauchy <- function(xi_init,y, mu, sigma){
  neg_log_posterior_xi <- function(xi, y, mu, sigma) {
    # Inputs:
    # xi     : scalar shape parameter (numeric)
    # y      : observed response vector
    # mu     : scalar threshold (same for all y)
    # sigma  : scale parameters vector (same length as y)
    
    # Step 1: Check constraints
    z <- (y - mu) / sigma
    # Check for validity: support constraint
    if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
      return(1e10)  # Penalise invalid xi
    }
    
    # Step 2: Log-likelihood part
    loglik <- sum(-(1/xi + 1) * log(1 + xi * z))
    
    # Step 3: Log-prior part (Truncated Cauchy up to constant)
    logprior <- -log(1 + xi^2)
    
    # Step 4: Return log-posterior up to proportionality
    logpost <- loglik + logprior
    return(-logpost)
  }
  
  opt <- optimize(neg_log_posterior_xi, interval = c(-1e6, 1)
                  , y=y, mu=mu, sigma=sigma, tol = 1e-6)
  return(opt$minimum)
}

posterior_mode_sigma <- function(y, mu, xi) {
  # y  : observed vector (numeric)
  # mu : threshold (scalar)
  # xi : shape parameter (scalar)
  
  z <- y - mu
  term1 <- 1 + xi * z
  term2 <- sqrt(term1^2 + 4 * z)
  
  sigma_hat <- (-term1 + term2) / 2
  return(sigma_hat)
}


Jacobi_estimate <- function(y, X, mu, xi_init = 0, 
                              tol = 1e-6, max_iter = 1000,
                              verbose = TRUE) {
  
  n <- length(y)
  p <- ncol(X)
  
  # Precompute A and S
  A <- solve(t(X) %*% X) %*% t(X)
  
  # Initialise alpha, beta
  xi <- xi_init
  beta <- rep(0, p)
  
  beta_errors <- c()
  
  epsilon <- 1e-6  # lower bound for numerical stability
  
  for (iter in 1:max_iter) {
    eta <- X %*% beta
    sigma = exp(eta)
    # Optimise xi from posterior
    xi_new <- estimate_xi_traunc_cauchy(xi_init = xi,
                                        y=y, mu=mu, sigma=sigma)
    
    
    ## Estimate sigma
    sigma_new = posterior_mode_sigma(y=y,mu=mu,xi=xi_new)
    
    eta_new = log(sigma_new)
    
    # Check for NaNs and stop early
    if (any(is.nan(eta_new)) || is.nan(xi_new)) {
      warning("NaN encountered in eta_new or xi_new at iteration ", iter)
      if (verbose) {
        print(list(xi_new = xi_new))
      }
      break
    }
    
    beta_new <- A %*% eta_new
    
    if (any(is.nan(beta_new))) {
      warning("NaN encountered in beta_new at iteration ", iter)
      if (verbose) {
        print(list(eta_new = eta_new[1:5], beta_new = beta_new[1:5]))
      }
      break
    }
    
    # Compute alpha error (L2 norm)
    xi_diff = abs(xi_new-xi)
    # Compute beta error (L2 norm)
    beta_diff <- sqrt(sum((beta_new - beta)^2))
    beta_errors <- c(beta_errors, beta_diff)
    
    # Early stopping when beta converges within tolerance
    if (beta_diff < tol && xi_diff<tol) {
      if (verbose) message("Beta converged within tolerance at iteration ", iter)
      break
    }
    
    # Update for next iteration
    xi <- xi_new
    sigma <- sigma_new
    beta <- beta_new
    
    if(verbose){
      cat("iter:", iter,
          "| alpha_error:", round(xi_diff, 4),
          "| beta_error:", round(beta_diff, 4), "\n")
    }
    
  }
  if(verbose) cat("Numbe of Iteration :", iter,"\n")
  return(list(xi = xi, sigma = sigma, beta = beta))
}


## Full MLE

neg_log_lik_gpd <- function(par, y, X, mu) {
  # par = c(beta_1, ..., beta_p, xi)
  p <- ncol(X)
  beta <- par[1:p]
  xi <- par[p + 1]
  
  eta <- X %*% beta
  sigma <- exp(eta)
  z <- (y - mu) / sigma
  
  # check for valid support
  if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
    return(1e10)  # Penalise invalid xi or likelihood domain
  }
  
  term1 <- sum(log(sigma))
  term2 <- sum((1 / xi + 1) * log(1 + xi * z))
  
  return(term1 + term2)  # Negative log-likelihood
}

## Cauchy prior on beta and truncated Cauchy prior on xi

nlp_gpd_cauchy_prior <- function(par, y, X, mu) {
  p <- ncol(X)
  beta <- par[1:p]
  xi <- par[p + 1]
  
  eta <- X %*% beta
  sigma <- exp(eta)
  z <- (y - mu) / sigma
  
  # Check support
  if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
    return(1e10)  # assign very high neg-log-posterior
  }
  
  # Log-likelihood
  nll <- neg_log_lik_gpd(par = par, y = y, X=X, mu=mu)
  
  # Log prior for beta ~ Cauchy(0,1)
  nlp_beta <- sum(log(1 + beta^2))
  
  # Log prior for xi ~ Truncated Cauchy(0,1)
  nlp_xi <- log(1 + xi^2)
  
  return(nll + nlp_beta + nlp_xi)

}

## GPD regression MLE with L1 regularisation (LASSO) on beta
## using 5-fold cross-validation


neg_log_lik_gpd_lasso <- function(par, y, X, mu, lambda) {
  p <- ncol(X)
  beta <- par[1:p]
  xi <- par[p + 1]
  
  eta <- X %*% beta
  sigma <- exp(eta)
  z <- (y - mu) / sigma
  
  if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
    return(1e10)
  }
  
  # Negative log-likelihood
  nll <- sum(log(sigma) + (1 / xi + 1) * log(1 + xi * z))
  
  # L1 penalty (no penalty on xi)
  penalty <- lambda * sum(abs(beta))
  
  return(nll + penalty)
}

fit_gpd_lasso <- function(lambda, y, X, mu, init_par = NULL) {
  p <- ncol(X)
  if (is.null(init_par)) init_par <- c(rep(0, p), 0.1)
  
  fit <- optim(par = init_par,
               fn = neg_log_lik_gpd_lasso,
               y = y, X = X, mu = mu, lambda = lambda,
               method = "BFGS",
               control = list(maxit = 1000))
  
  return(list(par = fit$par, value = fit$value, converged = fit$convergence == 0))
}

cv_gpd_lasso <- function(y, X, mu, lambda_grid, k = 5) {
  n <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  
  cv_errors <- numeric(length(lambda_grid))
  
  for (l in seq_along(lambda_grid)) {
    lambda <- lambda_grid[l]
    fold_errors <- numeric(k)
    
    for (fold in 1:k) {
      test_idx <- which(folds == fold)
      train_idx <- setdiff(1:n, test_idx)
      
      y_train <- y[train_idx]
      X_train <- X[train_idx, , drop = FALSE]
      
      y_test <- y[test_idx]
      X_test <- X[test_idx, , drop = FALSE]
      
      fit <- fit_gpd_lasso(lambda, y_train, X_train, mu)
      if (!fit$converged) {
        fold_errors[fold] <- Inf
        next
      }
      
      beta_hat <- fit$par[1:ncol(X)]
      xi_hat <- fit$par[ncol(X) + 1]
      
      # Evaluate on test set (without penalty)
      eta_test <- X_test %*% beta_hat
      sigma_test <- exp(eta_test)
      z_test <- (y_test - mu) / sigma_test
      
      if (xi_hat >= 1 || any(1 + xi_hat * z_test <= 0) || !is.finite(xi_hat)) {
        fold_errors[fold] <- Inf
        next
      }
      
      nll_test <- sum(log(sigma_test) + (1 / xi_hat + 1) * log(1 + xi_hat * z_test))
      fold_errors[fold] <- nll_test
    }
    
    cv_errors[l] <- mean(fold_errors)
  }
  
  best_idx <- which.min(cv_errors)
  best_lambda <- lambda_grid[best_idx]
  
  return(list(best_lambda = best_lambda,
              cv_errors = cv_errors,
              lambda_grid = lambda_grid))
}

## GPD regression MLE with L2 regularisation (Ridge) on beta
## using 5-fold cross-validation

# Ridge-penalised negative log-likelihood function
neg_log_lik_gpd_ridge <- function(par, y, X, mu, lambda) {
  p <- ncol(X)
  beta <- par[1:p]
  xi <- par[p + 1]
  
  eta <- X %*% beta
  sigma <- exp(eta)
  z <- (y - mu) / sigma
  
  if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
    return(1e10)
  }
  
  # Negative log-likelihood
  nll <- sum(log(sigma) + (1 / xi + 1) * log(1 + xi * z))
  
  # L2 penalty (no penalty on xi)
  penalty <- lambda * sum(beta^2)
  
  return(nll + penalty)
}

# Fit model for a given lambda
fit_gpd_ridge <- function(lambda, y, X, mu, init_par = NULL) {
  p <- ncol(X)
  if (is.null(init_par)) init_par <- c(rep(0, p), 0.1)
  
  fit <- optim(par = init_par,
               fn = neg_log_lik_gpd_ridge,
               y = y, X = X, mu = mu, lambda = lambda,
               method = "BFGS",
               control = list(maxit = 1000))
  
  return(list(par = fit$par, value = fit$value, converged = fit$convergence == 0))
}

# Cross-validation over lambda grid
cv_gpd_ridge <- function(y, X, mu, lambda_grid, k = 5) {
  n <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  
  cv_errors <- numeric(length(lambda_grid))
  
  for (l in seq_along(lambda_grid)) {
    lambda <- lambda_grid[l]
    fold_errors <- numeric(k)
    
    for (fold in 1:k) {
      test_idx <- which(folds == fold)
      train_idx <- setdiff(1:n, test_idx)
      
      y_train <- y[train_idx]
      X_train <- X[train_idx, , drop = FALSE]
      
      y_test <- y[test_idx]
      X_test <- X[test_idx, , drop = FALSE]
      
      fit <- fit_gpd_ridge(lambda, y_train, X_train, mu)
      if (!fit$converged) {
        fold_errors[fold] <- Inf
        next
      }
      
      beta_hat <- fit$par[1:ncol(X)]
      xi_hat <- fit$par[ncol(X) + 1]
      
      eta_test <- X_test %*% beta_hat
      sigma_test <- exp(eta_test)
      z_test <- (y_test - mu) / sigma_test
      
      if (xi_hat >= 1 || any(1 + xi_hat * z_test <= 0) || !is.finite(xi_hat)) {
        fold_errors[fold] <- Inf
        next
      }
      
      nll_test <- sum(log(sigma_test) + (1 / xi_hat + 1) * log(1 + xi_hat * z_test))
      fold_errors[fold] <- nll_test
    }
    
    cv_errors[l] <- mean(fold_errors)
  }
  
  best_idx <- which.min(cv_errors)
  best_lambda <- lambda_grid[best_idx]
  
  return(list(best_lambda = best_lambda,
              cv_errors = cv_errors,
              lambda_grid = lambda_grid))
}


### Zellner's g-prior

## Negative log-posterior for GPD with Zellner's g-prior

neg_log_post_gpd_gprior <- function(par, y, X, mu, g_scale) {
  p <- ncol(X)
  beta <- par[1:p]
  xi <- par[p + 1]
  
  eta <- X %*% beta
  sigma <- exp(eta)
  z <- (y - mu) / sigma
  
  if (xi >= 1 || any(1 + xi * z <= 0) || !is.finite(xi)) {
    return(1e10)
  }
  
  # Log-likelihood
  nll <- neg_log_lik_gpd(par = par, y = y, X=X, mu=mu)
  
  # L2 penalty (no penalty on xi)
  M_inv <- solve(t(X) %*% X + diag(g_scale,nrow=p,ncol=p))
  Sigma <- M_inv*g_scale
  
  neglogprior <- -dmvnorm(beta, mean = rep(0, p)
                        ,sigma = Sigma,log=TRUE)
  
  return(nll + neglogprior)
}
## Fit GPD with Zellner's g-prior
fit_gpd_gprior <- function(y, X, mu, g_scale, init_par = NULL) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  p <- ncol(X)
  if (is.null(init_par)) init_par <- c(rep(0, p), 0.1)
  
  fit <- optim(par = init_par,
               fn = neg_log_post_gpd_gprior,
               y = y, X = X, mu = mu
               , g_scale = g_scale,  # <--- pass g here
               method = "BFGS",
               control = list(maxit = 1000))
  
  return(list(par = fit$par, value = fit$value, converged = fit$convergence == 0))
}

## Optional: Cross-validate g from a grid
cv_gpd_gprior <- function(y, X, mu, g_grid, k = 5) {
  n <- length(y)
  folds <- sample(rep(1:k, length.out = n))
  
  cv_errors <- numeric(length(g_grid))
  
  for (l in seq_along(g_grid)) {
    g_scale <- g_grid[l]
    fold_errors <- numeric(k)
    
    for (fold in 1:k) {
      test_idx <- which(folds == fold)
      train_idx <- setdiff(1:n, test_idx)
      
      y_train <- y[train_idx]
      X_train <- X[train_idx, , drop = FALSE]
      
      y_test <- y[test_idx]
      X_test <- X[test_idx, , drop = FALSE]
      
      fit <- fit_gpd_gprior(g_scale=g_scale,y= y_train
                            , X=X_train, mu=mu)
      if (!fit$converged) {
        fold_errors[fold] <- Inf
        next
      }
      
      beta_hat <- fit$par[1:ncol(X)]
      xi_hat <- fit$par[ncol(X) + 1]
      
      eta_test <- X_test %*% beta_hat
      sigma_test <- exp(eta_test)
      z_test <- (y_test - mu) / sigma_test
      
      if (xi_hat >= 1 || any(1 + xi_hat * z_test <= 0) || !is.finite(xi_hat)) {
        fold_errors[fold] <- Inf
        next
      }
      
      nll_test <- sum(log(sigma_test) + (1 / xi_hat + 1) * log(1 + xi_hat * z_test))
      fold_errors[fold] <- nll_test
    }
    
    cv_errors[l] <- mean(fold_errors)
  }
  
  best_idx <- which.min(cv_errors)
  best_g <- g_grid[best_idx]
  
  return(list(best_g = best_g,
              cv_errors = cv_errors,
              g_grid = g_grid))
}
