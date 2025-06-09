# 导入库
library(Matrix)  
library(caret)  
library(dplyr)  
library(MASS)  
library(SGL)  
library(sgs)  
library(dfr)  
library(cluster)  
library(ggplot2)  
library(tidyr)  
library(RcppML)  
library(glmnet)  
  
compute_metrics <- function(beta_true, beta_hat, y_true, y_pred) {  
  specificity <- sum(beta_hat[beta_true == 0] == 0) / sum(beta_true == 0)  
  sensitivity <- sum(beta_hat[beta_true != 0] != 0) / sum(beta_true != 0)  
  mse_y <- mean((y_true - y_pred)^2)  
  mse_beta <- mean((beta_true - beta_hat)^2)  
  list(specificity = specificity, sensitivity = sensitivity, mse_y = mse_y, mse_beta = mse_beta)  
}  
 
reorder_group <- function(groups) {  
  unique_groups <- sort(unique(groups))  
  match(groups, unique_groups)  
}  

fit_elastic_net <- function(X_train, y_train, X_test, y_test, true_beta, alpha_values) {
  results <- data.frame()
  
  for (alpha_value in alpha_values) {
    tryCatch({
      cv_fit <- cv.glmnet(
        x = as.matrix(X_train),
        y = y_train,
        alpha = alpha_value,
        nfolds = 5,
        standardize = FALSE
      )

      best_lambda <- cv_fit$lambda.min
      final_model <- glmnet(
        x = as.matrix(X_train),
        y = y_train,
        alpha = alpha_value,
        lambda = best_lambda
      )

      beta_hat <- coef(final_model)[-1]  

      y_pred <- predict(final_model, newx = as.matrix(X_test), type = "response")

      metrics <- compute_metrics(
        beta_true = true_beta,
        beta_hat = beta_hat,
        y_true = y_test,
        y_pred = y_pred
      )

      result <- data.frame(
        alpha = alpha_value,
        mse_beta = metrics$mse_beta,
        mse_y = metrics$mse_y,
        sensitivity = metrics$sensitivity,
        specificity = metrics$specificity,
        num_selected = sum(beta_hat != 0)
      )
      
      results <- rbind(results, result)
      
    }, error = function(e) {
      warning(paste("Error in fit_elastic_net:", e$message))
    })
  }
  
  return(results)
}

p <- 500  
n <- 100  
num_groups <- 10  
alpha_values <- c(0.95)  
signal_strength_values <- c(0.1, 0.5, 1, 5, 10)  
n_simulations <- 100  
train_prop <- 0.8  
correlation <- 0.5  
 
all_results <- data.frame()  

for (signal_strength in signal_strength_values) { 
  
  for (sim_id in 1:n_simulations) {      
    cat("\n=== Testing signal strength:", signal_strength, "iteration:", sim_id, "===\n")   

    groups <- rep(1:num_groups, each = p / num_groups)  
 
    set.seed(as.integer(signal_strength * 10 + sim_id))    
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
      signal_sd = sqrt(signal_strength), 
      seed_id = as.integer(signal_strength * 10 + sim_id)   
    )    
    
    X <- full_data$X
    y <- full_data$y
    true_beta <- full_data$true_beta

    set.seed(2000 * p + sim_id)
    train_idx <- sample(1:n, size = floor(train_prop * n))
    test_idx <- setdiff(1:n, train_idx)

    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test <- X[test_idx, ]
    y_test <- y[test_idx]

    simulation_results <- fit_elastic_net(
      X_train = X_train,
      y_train = y_train,
      X_test = X_test,
      y_test = y_test,
      true_beta = true_beta,
      alpha_values = alpha_values
    )

    simulation_results$signal_strength <- signal_strength 
    simulation_results$iteration <- sim_id

    all_results <- rbind(all_results, simulation_results)

    if (sim_id %% 10 == 0) {
      temp_file <- paste0("features_results_signal_strength_partial_", signal_strength, "_", sim_id, ".csv")
      write.csv(all_results, temp_file, row.names = FALSE)
    }
  }
}

write.csv(all_results, "elastic_net_signal_strength.csv", row.names = FALSE)