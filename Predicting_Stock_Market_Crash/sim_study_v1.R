# Load required library
library(evd)  # for rgpd
library(mvtnorm)
library(MASS)
source("functions_v1.R")

No.of.Datasets<-100
RNGkind(sample.kind = "Rounding")
set.seed(123)
seed_str<-sample(1:1000,No.of.Datasets,replace = T)

mu=2
p = 5
sigma = 1
rho = 0.5
n = 100

rmse_ridge = rmse_beta_ridge  = rmse_xi_ridge = 0 
rmse_lasso = rmse_beta_lasso = rmse_xi_lasso = 0 
rmse_cauchy = rmse_beta_cauchy =rmse_xi_cauchy = 0
rmse_gprior = rmse_beta_gprior =rmse_xi_gprior = 0

aic_cauchy = aic_lasso = aic_ridge = aic_gprior =  0
bic_cauchy = bic_lasso = bic_ridge = bic_gprior = 0

time_taken_cauchy = 0
time_taken_ridge = 0
time_taken_lasso = 0
time_taken_gprior = 0

B_hat = matrix(NA,nrow = No.of.Datasets, ncol = p)

train_size = n*0.80
k=1
for (k in 1:No.of.Datasets){
  set.seed(seed_str[k])
  
  beta_true <- rnorm(p)
  xi_true <- runif(1,min = -1/2,max=1/2)
  
  data <- simulate_gpd_regression(n=n, beta = beta_true
                                  , mu=mu, xi = xi_true
                                  ,rho = 0.0)
  y <- data$y
  X<- data$X
  
  X_train = X[1:train_size,]
  X_test = X[(train_size+1):n,]
  X_test = as.matrix(X_test)
  Data_test = data.frame(X_test)
  
  y_train = y[1:train_size]
  y_test = y[(train_size+1):n]
  Data = cbind.data.frame(y,X)
  colnames(Data_test) = colnames(Data)[2:length(colnames(Data))]
  
  ##---------------------------------------------
  
  ## Cauchy Prior--------------
  
  time_taken_cauchy_start = Sys.time()
      start_par <- c(rep(0,p),0.1)
      fit_cauchy <- optim(par = start_par,
                          fn = nlp_gpd_cauchy_prior,
                          y = y_train, X = X_train, mu = mu,
                          method = "BFGS",
                          control = list(maxit = 1000))
      
      
      # Results
      beta_hat_cauchy <- fit_cauchy$par[1:p]
      xi_hat_cauchy <- fit_cauchy$par[p + 1]
      eta_cauchy = X_test%*%beta_hat_cauchy
      sigma_cauchy = exp(eta_cauchy)
      y_pred_cauchy = mu+(sigma_cauchy/(1-xi_hat_cauchy))
      rmse_cauchy_scr = sqrt(mean((y_pred_cauchy-y_test)^2))
      rmse_beta_cauchy_scr = sqrt(mean((beta_hat_cauchy-beta_true)^2))
      rmse_xi_cauchy_scr = sqrt(mean((xi_hat_cauchy -xi_true)^2))
      metrics_cauchy <- compute_aic_bic(
                                par_hat = fit_cauchy$par,
                                y = y_train,
                                X = X_train,
                                mu = mu,
                                neg_log_lik_fn = neg_log_lik_gpd
                              )
      
  time_taken_cauchy_end = Sys.time()
  rmse_cauchy[k] = rmse_cauchy_scr
  rmse_beta_cauchy[k] = rmse_beta_cauchy_scr
  rmse_xi_cauchy[k] = rmse_xi_cauchy_scr
  aic_cauchy[k] = metrics_cauchy$AIC
  bic_cauchy[k] = metrics_cauchy$BIC
  time_taken_cauchy[k] = time_taken_cauchy_end - time_taken_cauchy_start
  
  ## Lasso--------------
  ## GPD regression MLE with L1 regularisation (LASSO) on beta
  ## using 5-fold cross-validation
  
  time_taken_lasso_start = Sys.time()
        lambda_grid <- 10^seq(-2, 2, length.out = 50)
        
        cv_result <- cv_gpd_lasso(y = y_train, X = X_train
                                  , mu = mu
                                  , lambda_grid = lambda_grid
                                  , k = 5)
        # Final fit
        fit_lasso <- fit_gpd_lasso(lambda = cv_result$best_lambda
                                   , y = y_train
                                   , X = X_train, mu = mu)
        beta_hat_lasso <- fit_lasso$par[1:p]
        xi_hat_lasso <- fit_lasso$par[p + 1]
        eta_lasso = X_test%*%beta_hat_lasso
        sigma_lasso = exp(eta_lasso)
        y_pred_lasso = mu+(sigma_lasso/(1-xi_hat_lasso))
        rmse_lasso_scr = sqrt(mean((y_pred_lasso-y_test)^2))
        rmse_beta_lasso_scr = sqrt(mean((beta_hat_lasso-beta_true)^2))
        rmse_xi_lasso_scr = sqrt(mean((xi_hat_lasso -xi_true)^2))
        metrics_lasso <- compute_aic_bic(
                              par_hat = fit_lasso$par,
                              y = y_train,
                              X = X_train,
                              mu = mu,
                              neg_log_lik_fn = neg_log_lik_gpd
                            )
        
  time_taken_lasso_end = Sys.time()
  rmse_lasso[k] = rmse_lasso_scr
  rmse_beta_lasso[k] = rmse_beta_lasso_scr
  rmse_xi_lasso[k] = rmse_xi_lasso_scr
  aic_lasso[k] = metrics_lasso$AIC
  bic_lasso[k] = metrics_lasso$BIC
  time_taken_lasso[k] = time_taken_lasso_end - time_taken_lasso_start
  
  ## Ridge--------------
  ## GPD regression MLE with L2 regularisation (Ridge) on beta
  ## using 5-fold cross-validation
  
  time_taken_ridge_start = Sys.time()
        lambda_grid <- 10^seq(-2, 2, length.out = 50)
        
        cv_result <- cv_gpd_ridge(y = y_train, X = X_train
                                  , mu = mu
                                  , lambda_grid = lambda_grid
                                  , k = 5)
        # Final fit
        fit_ridge <- fit_gpd_ridge(lambda = cv_result$best_lambda
                                   , y = y_train
                                   , X = X_train, mu = mu)
        
        beta_hat_ridge <- fit_ridge$par[1:p]
        xi_hat_ridge <- fit_ridge$par[p + 1]
        eta_ridge = X_test%*%beta_hat_ridge
        sigma_ridge = exp(eta_ridge)
        y_pred_ridge = mu+(sigma_ridge/(1-xi_hat_ridge))
        rmse_ridge_scr = sqrt(mean((y_pred_ridge-y_test)^2))
        rmse_beta_ridge_scr = sqrt(mean((beta_hat_ridge-beta_true)^2))
        rmse_xi_ridge_scr = sqrt(mean((xi_hat_ridge -xi_true)^2))
        metrics_ridge <- compute_aic_bic(
                              par_hat = fit_ridge$par,
                              y = y_train,
                              X = X_train,
                              mu = mu,
                              neg_log_lik_fn = neg_log_lik_gpd
                            )
  time_taken_ridge_end = Sys.time()
  rmse_ridge[k] = rmse_ridge_scr
  rmse_beta_ridge[k] = rmse_beta_ridge_scr
  rmse_xi_ridge[k] = rmse_xi_ridge_scr
  aic_ridge[k] = metrics_ridge$AIC
  bic_ridge[k] = metrics_ridge$BIC
  time_taken_ridge[k] = time_taken_ridge_end - time_taken_ridge_start
  
  ## GPD regression with Zellner's g-prior
  ## on beta using 5-fold cross-validation
  time_taken_gprior_start = Sys.time()
        g_grid <- 10^seq(-2, 2, length.out = 50)
        cv_result <- cv_gpd_gprior(y = y_train, X = X_train
                                   , mu = mu, g_grid = g_grid)
        
        # Final fit
        fit_gprior <- fit_gpd_gprior(
          g_scale = cv_result$best_g,
          y = y_train,
          X = X_train,
          mu = mu
        )
        beta_hat_gprior <- fit_gprior$par[1:p]
        xi_hat_gprior <- fit_gprior$par[p + 1]
        eta_gprior = X_test%*%beta_hat_gprior
        sigma_gprior = exp(eta_gprior)
        y_pred_gprior = mu+(sigma_gprior/(1-xi_hat_gprior))
        rmse_gprior_scr = sqrt(mean((y_pred_gprior-y_test)^2))
        rmse_beta_gprior_scr = sqrt(mean((beta_hat_gprior-beta_true)^2))
        rmse_xi_gprior_scr = sqrt(mean((xi_hat_gprior -xi_true)^2))
        metrics_gprior <- compute_aic_bic(
          par_hat = fit_gprior$par,
          y = y_train,
          X = X_train,
          mu = mu,
          neg_log_lik_fn = neg_log_lik_gpd
        )
  time_taken_gprior_end = Sys.time()
  rmse_gprior[k] = rmse_gprior_scr
  rmse_beta_gprior[k] = rmse_beta_gprior_scr
  rmse_xi_gprior[k] = rmse_xi_gprior_scr
  aic_gprior[k] = metrics_gprior$AIC
  bic_gprior[k] = metrics_gprior$BIC
  time_taken_gprior[k] = time_taken_gprior_end - time_taken_gprior_start
  if(k %% 10 ==0){cat("k = ",k,"\n")}      
}

rmse = c(rmse_cauchy = median(rmse_cauchy)
         ,rmse_lasso = median(rmse_lasso)
         ,rmse_ridge = median(rmse_ridge)
         ,rmse_gprior = median(rmse_gprior)
)

rmse_beta = c(rmse_beta_cauchy = median(rmse_beta_cauchy)
              ,rmse_beta_lasso = median(rmse_beta_lasso)
              ,rmse_beta_ridge = median(rmse_beta_ridge)
              ,rmse_beta_gprior = median(rmse_beta_gprior))
         
rmse_xi = c(rmse_xi_cauchy = median(rmse_xi_cauchy)
            ,rmse_xi_lasso = median(rmse_xi_lasso)
            ,rmse_xi_ridge = median(rmse_xi_ridge)
            ,rmse_xi_gprior = median(rmse_xi_gprior))

aic  = c(aic_cauchy = median(aic_cauchy)
         ,aic_lasso = median(aic_lasso)
         ,aic_ridge = median(aic_ridge)
         ,aic_gprior = median(aic_gprior))

bic  = c(bic_cauchy = median(bic_cauchy)
         ,bic_lasso = median(bic_lasso)
         ,bic_ridge = median(bic_ridge)
         ,bic_gprior = median(bic_gprior))



time_take = c(time_taken_cauchy = median(time_taken_cauchy)
              ,time_taken_lasso = median(time_taken_lasso)
              ,time_taken_ridge = median(time_taken_ridge)
              ,time_taken_gprior = median(time_taken_gprior))

metrics = cbind.data.frame(rmse,rmse_beta,rmse_xi,aic,bic,time_take)
rownames(metrics)<-c("Cauchy","Lasso","Ridge","g-prior")
colnames(metrics)<-c("RMSE_y","RMSE_beta","RMSE_xi","AIC","BIC","time_take")
metrics$time_multiples<-1
metrics$time_multiples[2]<-metrics$time_take[2]/metrics$time_take[1]
metrics$time_multiples[3]<-metrics$time_take[3]/metrics$time_take[1]
metrics$time_multiples[4]<-metrics$time_take[4]/metrics$time_take[1]
metrics <- t(metrics)

library(xtable)
xtable(metrics)