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
library(dplyr) 
library(foreach)  
library(doParallel)  



# Compute metrics
compute_metrics <- function(beta_true, beta_hat, y_true, y_pred) {
  specificity <- sum(beta_hat[beta_true == 0] == 0) / sum(beta_true == 0)
  sensitivity <- sum(beta_hat[beta_true != 0] != 0) / sum(beta_true != 0)
  mse_y <- mean((y_true - y_pred)^2)
  mse_beta <- mean((beta_true - beta_hat)^2)
  list(specificity = specificity, sensitivity = sensitivity, mse_y = mse_y, mse_beta = mse_beta)
}


# Reorder group function
reorder_group <- function(groups) {
  unique_groups <- sort(unique(groups))
  match(groups, unique_groups)
}



# Main parameters  
p <- 500  
n <- 100   
num_groups <- 10  
group_size <- p / num_groups  
alpha_values <- c(0, 0.95)  # 0 for Group Lasso, 0.95 for Sparse Group Lasso  
num_iterations <- 100

all_data_sets <- list()  

for (i in 1:num_iterations) {  
  groups <- rep(1:num_groups, each = p/num_groups)  

  data <- gen_toy_data(  
    p = p,  
    n = n,  
    rho = 0.5,  
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

  set.seed(i*100)  
  fold_indices <- sample(rep(1:5, length.out = n))  
  
  all_data_sets[[i]] <- list(  
    X = data$X,  
    y = data$y,  
    beta = data$true_beta,  
    fold_indices = fold_indices  
  )  
}  

process_iteration <- function(i, all_data_sets, alpha_values) {  
  
  data_set <- all_data_sets[[i]]  
  X <- data_set$X  
  y <- data_set$y  
  beta <- data_set$beta  
  fold_indices <- data_set$fold_indices  
  
  performance_results <- list()  
  clustering_results <- list()  

  all_clustering <- list()  
  
  # SVD grouping  
  cat("\nPerforming SVD clustering...\n")  
  dec <- svd(X)  
  grp_svd <- apply(dec$v, 1, function(x) which.max(x * dec$d))  
  grp_svd <- reorder_group(grp_svd)  
  order_svd_grp <- order(grp_svd, decreasing = FALSE)  

  all_clustering$SVD <- data.frame(  
    feature_index = 1:ncol(X),  
    group = grp_svd,  
    method = "SVD",  
    iteration = i,  
    fold = 0  
  )  
  
  # QR grouping  
  cat("\nPerforming QR clustering...\n")  
  qr_dec <- qr(X)  
  R <- qr.R(qr_dec)  
  grp_qr <- apply(R, 2, which.max)  
  grp_qr <- reorder_group(grp_qr)  
  order_qr_grp <- order(grp_qr, decreasing = FALSE)  

  all_clustering$QR <- data.frame(  
    feature_index = 1:ncol(X),  
    group = grp_qr,  
    method = "QR",  
    iteration = i,  
    fold = 0  
  )  
  
  # KMeans grouping  
  cat("\nPerforming K-means clustering...\n")  
  kmeans_gap <- cluster::clusGap(t(X),  
                              FUN = function(x, k) {  
                                list(cluster = kmeans(x, centers = k, nstart = 25, iter.max = 500)$cluster)  
                              },  
                              K.max = 15,  
                              B = 50)  
  optimal_k <- maxSE(kmeans_gap$Tab[,"gap"], kmeans_gap$Tab[,"SE.sim"], method = "firstSEmax")  
  cat("\nOptimal number of clusters (K-means):", optimal_k, "\n")  
  kmeans_result <- kmeans(t(X), centers = optimal_k, nstart = 10)  
  grp_kmeans <- reorder_group(kmeans_result$cluster)  
  order_kmeans_grp <- order(grp_kmeans, decreasing = FALSE)  

  all_clustering$KMeans <- data.frame(  
    feature_index = 1:ncol(X),  
    group = grp_kmeans,  
    method = "KMeans",  
    iteration = i,  
    fold = 0  
  )  
  
  # NMF grouping  
  cat("\nPerforming NMF clustering...\n")  
  nmf_result <- RcppML::nmf(X, k = num_groups, maxit = 500)  
  grp_nmf <- apply(nmf_result$h, 2, which.max)  
  grp_nmf <- reorder_group(grp_nmf)  
  order_nmf_grp <- order(grp_nmf, decreasing = FALSE)  

  all_clustering$NMF <- data.frame(  
    feature_index = 1:ncol(X),  
    group = grp_nmf,  
    method = "NMF",  
    iteration = i,  
    fold = 0  
  )  
 
  for (method_name in names(all_clustering)) {  
    clustering_results[[length(clustering_results) + 1]] <- all_clustering[[method_name]]  
  }  

  groupings <- list(  
    SVD = list(grp = grp_svd, order = order_svd_grp),  
    QR = list(grp = grp_qr, order = order_qr_grp),  
    NMF = list(grp = grp_nmf, order = order_nmf_grp),  
    KMeans = list(grp = grp_kmeans, order = order_kmeans_grp)  
  )  
 
  for (fold in 1:5) {  
    train_indices <- which(fold_indices != fold)  
    test_indices <- which(fold_indices == fold)  
    
    X1 <- X[train_indices, ]  
    Y1 <- y[train_indices]    
    X2 <- X[test_indices, ]    
    Y2 <- y[test_indices]     
    
    # Lasso (alpha = 1)  
    cat("\nFitting Lasso (alpha = 1)...\n")  
    lasso_fit <- dfr_sgl.cv(  
      X = X1,  
      y = Y1,  
      groups = 1:ncol(X1),  
      type = "linear",  
      alpha = 1,  
      nfolds = 5,  
      screen = FALSE,  
      verbose = FALSE,  
      intercept = TRUE,  
      max_iter = 1000  
    )  
    
    lasso_beta <- as.numeric(lasso_fit$fit$beta[-1])  
    lasso_pred <- X2 %*% lasso_beta  # Predict on test set  
    num_selected_vars_lasso <- sum(lasso_beta != 0)  
    
    lasso_metrics <- compute_metrics(beta_true = beta, beta_hat = lasso_beta, y_true = Y2, y_pred = lasso_pred)  
    performance_results[[length(performance_results) + 1]] <- data.frame(  
      method = "Lasso",  
      specificity = lasso_metrics$specificity,  
      sensitivity = lasso_metrics$sensitivity,  
      mse_y = lasso_metrics$mse_y,  
      mse_beta = lasso_metrics$mse_beta,  
      num_selected_vars = num_selected_vars_lasso,  
      iteration = i,  
      fold = fold  
    )  

    for (method in names(groupings)) {  
      grouping_info <- groupings[[method]]  
      group_indices <- grouping_info$grp  
      order_grp <- grouping_info$order  
      
      cat("\n=== Processing method:", method, "===\n")  
      
      for (alpha in alpha_values) {  
        model_type <- ifelse(alpha == 0, "Group Lasso", "Sparse Group Lasso")  
        cat(sprintf("\nFitting %s (alpha = %.2f)...\n", model_type, alpha))  
        
        fit <- dfr_sgl.cv(  
          X = X1[, order_grp],  
          y = Y1,  
          groups = group_indices[order_grp],  
          type = "linear",  
          alpha = alpha,  
          nfolds = 5,  
          screen = FALSE,   
          verbose = FALSE,  
          intercept = TRUE,  
          max_iter = 1000  
        )  
        
        beta_hat <- as.numeric(fit$fit$beta)  
        beta_hat <- beta_hat[-1]  
        beta_original_order <- numeric(length(beta_hat))  
        beta_original_order[order_grp] <- beta_hat  
        
        y_pred <- X2 %*% beta_original_order  # Predict on test set  
        
        metrics <- compute_metrics(beta_true = beta, beta_hat = beta_original_order, y_true = Y2, y_pred = y_pred)  
        num_selected_vars <- sum(beta_original_order != 0)  
        
        performance_results[[length(performance_results) + 1]] <- data.frame(  
          method = paste(method, model_type, sep = "_"),  
          specificity = metrics$specificity,  
          sensitivity = metrics$sensitivity,  
          mse_y = metrics$mse_y,  
          mse_beta = metrics$mse_beta,  
          num_selected_vars = num_selected_vars,  
          iteration = i,  
          fold = fold  
        )  
      }  
    }  
  }  
  
  return(list(  
    performance_results = do.call(rbind, performance_results),  
    clustering_results = do.call(rbind, clustering_results)  
  ))  
}  

all_performance_results <- data.frame()  
all_clustering_results <- data.frame()  

for (i in 1:num_iterations) {  
  iteration_result <- process_iteration(i, all_data_sets, alpha_values)  

  all_performance_results <- rbind(all_performance_results, iteration_result$performance_results)  
  all_clustering_results <- rbind(all_clustering_results, iteration_result$clustering_results)  

  if (i %% 5 == 0 || i == num_iterations) {  

    write.csv(all_performance_results, file = paste0("intermediate_performance_", i, ".csv"), row.names = FALSE)  
    write.csv(all_clustering_results, file = paste0("intermediate_clustering_", i, ".csv"), row.names = FALSE)  
  }  
}  


write.csv(all_performance_results, file = "out_sample_cv_it.csv", row.names = FALSE)  