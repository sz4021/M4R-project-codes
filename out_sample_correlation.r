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

p <- 500  
n <- 100  
num_groups <- 10  
alpha_values <- c(0, 0.95)  
correlation_values <- c(0.1, 0.3, 0.5, 0.7, 0.9)  
n_simulations <- 100  
train_prop <- 0.8  
 
all_results <- data.frame()  
 
for (correlation in correlation_values) {  
  
  for (sim_id in 1:n_simulations) {      
    cat("\n=== Testing correlation:", correlation, "iteration:", sim_id, "===\n") 

    groups <- rep(1:num_groups, each = p/num_groups)  

    set.seed(1000 * as.numeric(correlation) + sim_id)  
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
      seed_id = 1000 * as.numeric(correlation) + sim_id  
    )  
    
    X <- full_data$X  
    y <- full_data$y  
    true_beta <- full_data$true_beta  

    set.seed(2000 * as.numeric(correlation) + sim_id)  
    train_idx <- sample(1:n, size = floor(train_prop * n))  
    test_idx <- setdiff(1:n, train_idx)  

    X_train <- X[train_idx, ]  
    y_train <- y[train_idx]  
    X_test <- X[test_idx, ]  
    y_test <- y[test_idx]  
 
    simulation_results <- data.frame()  
    
    # 1. Lasso (alpha = 1)
    lasso_fit <- dfr_sgl.cv(
      X = X_train,
      y = y_train,
      groups = 1:ncol(X_train),
      type = "linear",
      alpha = 1,
      nfolds = 5,
      screen = FALSE,
      verbose = FALSE,
      intercept = TRUE,
      max_iter = 10000
    )
    
    lasso_beta <- as.numeric(lasso_fit$fit$beta[-1])  

    y_pred_lasso <- X_test %*% lasso_beta  

    lasso_metrics <- compute_metrics(  
      beta_true = true_beta,  
      beta_hat = lasso_beta,  
      y_true = y_test,  
      y_pred = y_pred_lasso  
    )  
    
    lasso_result <- data.frame(  
      correlation = correlation,  
      iteration = sim_id,  
      method = "Lasso",  
      mse_beta = lasso_metrics$mse_beta,  
      mse_y = lasso_metrics$mse_y,  
      sensitivity = lasso_metrics$sensitivity,  
      specificity = lasso_metrics$specificity,  
      num_selected = sum(lasso_beta != 0)  
    )  
    
    simulation_results <- rbind(simulation_results, lasso_result)  

    # SVD grouping  
    dec <- svd(X_train)  
    grp_svd <- apply(dec$v, 1, function(x) which.max(x * dec$d))  
    grp_svd <- reorder_group(grp_svd)  
    order_svd_grp <- order(grp_svd, decreasing = FALSE)  
    
    # QR grouping  
    qr_dec <- qr(X_train)  
    R <- qr.R(qr_dec)  
    grp_qr <- apply(R, 2, which.max)  
    grp_qr <- reorder_group(grp_qr)  
    order_qr_grp <- order(grp_qr, decreasing = FALSE)  
    
    # KMeans grouping  
    tryCatch({  
      cat("\nStarting KMeans analysis...")  
      
      kmeans_gap <- cluster::clusGap(t(X_train),  
                                     FUN = function(x, k) {  
                                       list(cluster = kmeans(x, centers = k, nstart = 25, iter.max = 500)$cluster)  
                                     },  
                                     K.max = 15,  
                                     B = 10)  
      
      cat("\nGap statistic completed")  
      print(summary(kmeans_gap))  
      
      optimal_k <- maxSE(kmeans_gap$Tab[,"gap"], kmeans_gap$Tab[,"SE.sim"], method = "firstSEmax")  
      cat("\nOptimal K determined:", optimal_k)
      
      optimal_k <- max(2, optimal_k)  
      cat("\nAdjusted optimal K:", optimal_k)  
      kmeans_result <- kmeans(t(X_train), centers = optimal_k, nstart = 10)  
      print(table(kmeans_result$cluster))  
      
      grp_kmeans <- reorder_group(kmeans_result$cluster)  
      order_kmeans_grp <- order(grp_kmeans, decreasing = FALSE)  
    }, error = function(e) {  
      print(traceback())  
      grp_kmeans <- rep(1:num_groups, length.out = ncol(X_train))  
      order_kmeans_grp <- 1:ncol(X_train)  
    })
    
    # NMF grouping  
    tryCatch({  
      nmf_result <- RcppML::nmf(X_train, k = num_groups, maxit = 500)  
      grp_nmf <- apply(nmf_result$h, 2, which.max)  
      grp_nmf <- reorder_group(grp_nmf)  
      order_nmf_grp <- order(grp_nmf, decreasing = FALSE)  
    }, error = function(e) {  
      cat("\nError in NMF for simulation", sim_id, ":", e$message)  
      grp_nmf <- rep(1:num_groups, length.out = ncol(X_train))  
      order_nmf_grp <- 1:ncol(X_train)
    })  
    
    groupings <- list(  
      KMeans = list(grp = grp_kmeans, order = order_kmeans_grp),  
      SVD = list(grp = grp_svd, order = order_svd_grp),  
      QR = list(grp = grp_qr, order = order_qr_grp),  
      NMF = list(grp = grp_nmf, order = order_nmf_grp)  
    )  
    
    
    for (method_name in names(groupings)) {  
      grouping_info <- groupings[[method_name]]  
      group_indices <- grouping_info$grp  
      order_grp <- grouping_info$order 
      
      for (alpha in alpha_values) {  
        model_type <- ifelse(alpha == 0, "Group Lasso", "Sparse Group Lasso")  
        full_method_name <- paste(method_name, model_type, sep = "_")  
        
        tryCatch({  
          fit <- dfr_sgl.cv(  
            X = X_train[, order_grp],  
            y = y_train,  
            groups = group_indices[order_grp],  
            type = "linear",  
            alpha = alpha,  
            nfolds = 5,  
            screen = FALSE,  
            verbose = FALSE,  
            intercept = TRUE,  
            max_iter = 10000
          )  
          
          beta_hat <- as.numeric(fit$fit$beta)[-1]  
          beta_original_order <- numeric(length(beta_hat))  
          beta_original_order[order_grp] <- beta_hat  

          y_pred <- X_test %*% beta_original_order  

          metrics <- compute_metrics(  
            beta_true = true_beta,  
            beta_hat = beta_original_order,  
            y_true = y_test,  
            y_pred = y_pred  
          )  
          
          method_result <- data.frame(  
            correlation = correlation,  
            iteration = sim_id,  
            method = full_method_name,  
            mse_beta = metrics$mse_beta,  
            mse_y = metrics$mse_y,  
            sensitivity = metrics$sensitivity,  
            specificity = metrics$specificity,  
            num_selected = sum(beta_original_order != 0)  
          )  
          
          simulation_results <- rbind(simulation_results, method_result)  
        }, error = function(e) {  
          cat("\nError in", full_method_name, "for simulation", sim_id, ":", e$message)  
        })  
      }  
    }  

    all_results <- rbind(all_results, simulation_results)  

    if (sim_id %% 10 == 0) {  
      temp_file <- paste0("correlation_results_", correlation, "_partial_", sim_id, ".csv")  
      write.csv(all_results, temp_file, row.names = FALSE)  
    }  
  }  
}  

write.csv(all_results, "out_sample_correlation.csv", row.names = FALSE)  

