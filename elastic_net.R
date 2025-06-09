library(glmnet)   
library(Matrix)   
library(caret)    
library(dplyr)    
library(grpreg)
library(MASS)  
library(SGL)
library(sgs)
library(dfr)
library(cluster)
library(ggplot2)  
library(tidyr)  
library(RcppML)
library(splines)  

compute_metrics <- function(beta_true, beta_hat, y_true, y_pred) {
  specificity <- sum(beta_hat[beta_true == 0] == 0) / sum(beta_true == 0)
  sensitivity <- sum(beta_hat[beta_true != 0] != 0) / sum(beta_true != 0)
  mse_y <- mean((y_true - y_pred)^2)
  mse_beta <- mean((beta_true - beta_hat)^2)
  
  list(
    specificity = specificity, 
    sensitivity = sensitivity, 
    mse_y = mse_y, 
    mse_beta = mse_beta
  )
}

run_simulations <- function(p_values, n_simulations, 
                            n = 100, 
                            num_groups = 10, 
                            alpha_values = c(0.95), 
                            correlation = 0.5, 
                            train_prop = 0.8) {
  all_results <- data.frame()
  
  for (p in p_values) {
    groups <- rep(1:num_groups, each = p / num_groups)
    
    for (i in seq_len(n_simulations)) {            
      full_data <- gen_toy_data(
        p = p,
        n = n,
        rho = correlation,
        grouped = TRUE,
        groups = groups,
        group_sparsity = 0.5,
        var_sparsity = 0.5,
        noise_level = 1,
        orthogonal = FALSE,
        data_mean = 0,
        data_sd = 1,
        signal_mean = 0,
        signal_sd = sqrt(10),
        seed_id = i
      )
      
      X <- full_data$X
      y <- full_data$y
      true_beta <- full_data$true_beta

      set.seed(2000 * p + i)                         
      train_idx <- sample(seq_len(n), size = floor(train_prop * n))
      test_idx  <- setdiff(seq_len(n), train_idx)
      
      X_train <- X[train_idx, ]
      y_train <- y[train_idx]
      X_test  <- X[test_idx, ]
      y_test  <- y[test_idx]

      simulation_results <- fit_elastic_net(
        X_train = X_train,
        y_train = y_train,
        X_test  = X_test,
        y_test  = y_test,
        true_beta   = true_beta,
        alpha_values = alpha_values
      )

      simulation_results$features  <- p
      simulation_results$iteration <- i              

      all_results <- rbind(all_results, simulation_results)

      if (i %% 10 == 0) {
        write.csv(all_results,
                  sprintf("features_results_%d_partial_%03d.csv", p, i),
                  row.names = FALSE)
      }
    }
  }
  
  all_results
}

results <- run_simulations(
  p_values = c(500),
  n_simulations = 500
)

write.csv(results, "elastic_net.csv", row.names = FALSE)