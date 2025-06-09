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
library(RcppML)
library(parallel)  
library(foreach)  
library(doParallel) 

# Reorder group function
reorder_group <- function(groups) {
  unique_groups <- sort(unique(groups))
  match(groups, unique_groups)
}


fit_real_data <- function(X, y, max_kmeans_k) {  
  # First, define helper functions  
  maxSE <- function(gap, SE.sim, method = c("firstSEmax", "Tibs2001SEmax", "globalSEmax",  
                                           "firstmax", "globalmax")) {  
    method <- match.arg(method)  
    f <- function(x, i) x[i] - x[i-1]  
    g <- function(x, i) x[i] - x[i+1]  
    G <- f(gap, 2:length(gap))  
    sigG <- f(SE.sim, 2:length(SE.sim))  
    
    # Select optimal k based on different methods  
    if (method == "firstSEmax") {  
      for (i in 1:(length(G))) {  
        if (G[i] >= -sigG[i]) {  
          return(i + 1)  # Return corresponding k value  
        }  
      }  
      return(length(gap))  # If no conditions are met  
    } else if (method == "globalSEmax") {  
      j <- which.max(G >= -sigG)  
      if (length(j) == 0) j <- length(gap)  
      else j <- j[1] + 1  
      return(j)  
    } else if (method == "Tibs2001SEmax") {  
      for (i in 1:(length(gap) - 2)) {  
        if (g(gap, i) <= 0 && g(gap, i+1) > 0) {  
          return(i + 1)  
        }  
      }  
      return(length(gap) - 1)  
    } else if (method == "firstmax") {  
      for (i in 1:(length(gap) - 1)) {  
        if (g(gap, i) <= 0) {  
          return(i)  
        }  
      }  
      return(length(gap))  
    } else if (method == "globalmax") {  
      return(which.max(gap))  
    }  
  }  
  
  # ===== Initial Checks =====  
  cat("\n===== Initial Checks =====\n")  
  cat("X dimensions:", dim(X)[1], "rows ×", dim(X)[2], "columns\n")  
  cat("y length:", length(y), "\n")  
  cat("colnames(X) length:", length(colnames(X)), "\n")  
  cat("max_kmeans_k:", max_kmeans_k, "\n\n")  
  
  # ===== SVD Grouping =====  
  cat("===== SVD Grouping =====\n")  
  dec <- svd(X)  
  cat("SVD decomposition completed\n")  
  cat("dec$u dimensions:", dim(dec$u)[1], "×", dim(dec$u)[2], "\n")  
  cat("dec$v dimensions:", dim(dec$v)[1], "×", dim(dec$v)[2], "\n")  
  cat("dec$d length:", length(dec$d), "\n")  
  
  grp_svd <- apply(dec$v, 1, function(x) which.max(x * dec$d))  
  cat("SVD grouping before: unique groups", length(unique(grp_svd)), "\n")  
  cat("SVD grouping before length:", length(grp_svd), "\n")  
  
  grp_svd <- reorder_group(grp_svd)  
  cat("SVD grouping after: unique groups", length(unique(grp_svd)), "\n")  
  cat("SVD grouping after length:", length(grp_svd), "\n")  
  
  order_svd_grp <- order(grp_svd, decreasing = FALSE)  
  cat("order_svd_grp length:", length(order_svd_grp), "\n\n")  
  
  # ===== QR Grouping =====  
  cat("===== QR Grouping =====\n")  
  qr_dec <- qr(X)  
  R <- qr.R(qr_dec)  
  cat("QR decomposition completed\n")  
  cat("R matrix dimensions:", dim(R)[1], "×", dim(R)[2], "\n")  
  
  grp_qr <- apply(R, 2, which.max)  
  cat("QR grouping before: unique groups", length(unique(grp_qr)), "\n")  
  cat("QR grouping before length:", length(grp_qr), "\n")  
  
  grp_qr <- reorder_group(grp_qr)  
  cat("QR grouping after: unique groups", length(unique(grp_qr)), "\n")  
  cat("QR grouping after length:", length(grp_qr), "\n")  
  
  order_qr_grp <- order(grp_qr, decreasing = FALSE)  
  cat("order_qr_grp length:", length(order_qr_grp), "\n\n")  
  
  tryCatch({  
    # Create grouping dataframes  
    grouping_info <- list(  
      SVD = data.frame(  
        variable = colnames(X),  
        group = grp_svd,  
        group_details = paste("SVD_Group", grp_svd)  
      ),  
      QR = data.frame(  
        variable = colnames(X),  
        group = grp_qr,  
        group_details = paste("QR_Group", grp_qr)  
      )  
    )  
    
    cat("Successfully created grouping dataframes\n")  
  }, error = function(e) {  
    cat("Error creating grouping dataframes:", conditionMessage(e), "\n")  
    return(NULL)  
  })  
  
  # Check the number of rows in each grouping dataframe  
  for (method in names(grouping_info)) {  
    cat(method, "grouping dataframe rows:", nrow(grouping_info[[method]]), "\n")  
  }  
  
  # ===== Grouping Information Summary =====  
  cat("\n===== Grouping Information Summary =====\n")  
  
  group_summary <- tryCatch({  
    lapply(grouping_info, function(method_groups) {  
      # Summarize variables by group  
      group_vars <- split(method_groups$variable, method_groups$group)  
      
      # Create detailed description for each group  
      group_summaries <- lapply(names(group_vars), function(group_num) {  
        vars <- group_vars[[group_num]]  
        cat("Group", group_num, "contains", length(vars), "variables\n")  
        list(  
          group_number = as.numeric(group_num),  
          variables = vars,  
          num_variables = length(vars)  
        )  
      })  
      
      return(group_summaries)  
    })  
  }, error = function(e) {  
    cat("Error creating group summary:", conditionMessage(e), "\n")  
    return(NULL)  
  })  
  
  # Save grouping information to CSV  
  cat("\n===== Save Grouping Information to CSV =====\n")  
  
  tryCatch({  
    lapply(names(group_summary), function(method) {  
      # Convert grouping information for each method to dataframe  
      group_df <- do.call(rbind, lapply(group_summary[[method]], function(group) {  
        data.frame(  
          method = method,  
          group_number = group$group_number,  
          variables = paste(group$variables, collapse = ", "),  
          num_variables = group$num_variables  
        )  
      }))  
      
      # Save as CSV  
      write.csv(  
        group_df,  
        file = paste0(method, "_grouping_summary.csv"),  
        row.names = FALSE  
      )  
      cat(method, "grouping information saved to CSV\n")  
    })  
  }, error = function(e) {  
    cat("Error saving CSV:", conditionMessage(e), "\n")  
  })  
  
  # ===== Lasso Regularization =====  
  cat("\n===== Lasso Regularization =====\n")  
  
  lasso_fit <- tryCatch({  
    # 开始日志记录  
    cat("Starting Lasso regularization...\n")  
    # 执行Lasso回归  
      fit <- dfr_sgl.cv(  
        X = X,  
        y = y,  
        groups = 1:ncol(X),  # 每个变量作为单独分组  
        type = "linear",  
        alpha = 1,  
        nfolds = 10,  
        screen = TRUE, 
        verbose = TRUE,  # 改为TRUE以获取更多详细信息  
        intercept = TRUE,  
        max_iter = 5000,
        lambda = seq(0.0001, 10, length.out = 50) 
      )  
    
    # 成功完成日志  
      cat("Lasso regularization complete\n")  
    # 打印一些基本结果信息  
      cat("Lambda min:", fit$lambda.min, "\n")  
      cat("Number of non-zero coefficients:", sum(fit$beta != 0), "\n")  
    
      fit  
  })
  # Lasso results  
  lasso_beta <- as.numeric(lasso_fit$fit$beta[-1])  # Remove intercept  
  lasso_y_pred <- X %*% lasso_beta  
  lasso_mse_y <- mean((y - lasso_y_pred)^2)  
  lasso_num_selected_vars <- sum(lasso_beta != 0)  
  
  # Get the specific features selected by Lasso  
  lasso_selected_features <- colnames(X)[lasso_beta != 0]  
  
  cat("Lasso MSE:", lasso_mse_y, "\n")  
  cat("Lasso selected variables:", lasso_num_selected_vars, "\n")  
  cat("Lasso selected features:", paste(lasso_selected_features, collapse=", "), "\n")  
  
  results <- list(  
    Lasso = list(  
      method = "Lasso",  
      mse_y = lasso_mse_y,  
      num_selected_vars = lasso_num_selected_vars,  
      selected_features = lasso_selected_features  
    )  
  )  
  
  # ===== Group Lasso and Sparse Group Lasso =====  
  cat("\n===== Group Lasso and Sparse Group Lasso =====\n")  
  
  for (method_name in names(grouping_info)) {  
    cat("\n----- Processing", method_name, "method -----\n")  
    
    current_groups <- grouping_info[[method_name]]$group  
    order_grp <- order(current_groups, decreasing = FALSE)  
    
    cat("current_groups length:", length(current_groups), "\n")  
    cat("order_grp length:", length(order_grp), "\n")  
    cat("unique(current_groups):", paste(sort(unique(current_groups)), collapse=", "), "\n")  
    cat("First 5 elements of order_grp:", paste(head(order_grp, 5), collapse=", "), "\n")  
    
    # Group Lasso and Sparse Group Lasso  
    for (alpha in c(0, 0.95)) {  
      model_type <- ifelse(alpha < 0.01, "Group Lasso", "Sparse Group Lasso")  
      
      cat("\nStarting", method_name, model_type, "regularization...\n")  
      cat("X[, order_grp] dimensions:", nrow(X), "×", length(order_grp), "\n")  
      cat("current_groups[order_grp] length:", length(current_groups[order_grp]), "\n")  
      
      # Check for duplicate or non-continuous indices  
      if (anyDuplicated(order_grp) > 0) {  
        cat("Warning: order_grp contains duplicate values! Position:", anyDuplicated(order_grp), "\n")  
      }  
      
      if (min(order_grp) < 1 || max(order_grp) > ncol(X)) {  
        cat("Warning: order_grp contains invalid indices! Range:", min(order_grp), "to", max(order_grp), "\n")  
      }  
      
      # Fit model  
      fit <- tryCatch({  
        dfr_sgl.cv(  
          X = X[, order_grp, drop=FALSE],  
          y = y,  
          groups = current_groups[order_grp],  
          type = "linear",  
          alpha = alpha,  
          nfolds = 10,  
          screen = FALSE,  
          verbose = FALSE,  
          intercept = TRUE,  
          max_iter = 5000,
          lambda = seq(0.001, 10, length.out = 50) 
        )  
      }, error = function(e) {  
        cat(method_name, model_type, "regularization error:", conditionMessage(e), "\n")  
        print(traceback())  
        return(NULL)  
      })  
      
      if (is.null(fit)) {  
        cat(method_name, model_type, "failed, skipping result calculation\n")  
        next  
      }  
      
      cat(method_name, model_type, "regularization complete\n")  
      
      # Extract beta and predictions  
      beta_hat <- as.numeric(fit$fit$beta)[-1]  # Remove intercept  
      cat("beta_hat length:", length(beta_hat), "\n")  
      
      beta_original_order <- numeric(length(beta_hat))  
      beta_original_order[order_grp] <- beta_hat  
      
      y_pred <- X %*% beta_original_order  
      
      # Calculate metrics  
      mse_y <- mean((y - y_pred)^2)  
      num_selected_vars <- sum(beta_original_order != 0)  
      
      # Get the specific features selected by this method  
      selected_features <- colnames(X)[beta_original_order != 0]  
      
      cat(method_name, model_type, "MSE:", mse_y, "\n")  
      cat(method_name, model_type, "selected variables:", num_selected_vars, "\n")   
      cat(method_name, model_type, "selected features:", paste(selected_features, collapse=", "), "\n")  
      
      # Store results  
      results[[length(results) + 1]] <- list(  
        method = paste(method_name, model_type, sep="_"),  
        mse_y = mse_y,  
        num_selected_vars = num_selected_vars,  
        selected_features = selected_features  
      )  
    }  
  }  
  
  cat("\n===== Function execution complete =====\n")  
  
  # Return grouping information and results  
  return(list(  
    grouping_info = grouping_info,  
    group_summary = group_summary,  
    results = results  
  ))  
}  


scheetz <- readRDS("/Users/zhou/Desktop/IC大四/M4R/real data/scheetz.rds")

results_scheetz <- fit_real_data(scheetz$X, scheetz$y, max_kmeans_k = 95)

# 首先将结果转换为数据框  
results_df_scheetz <- do.call(rbind, lapply(results_scheetz$results, as.data.frame))  

write.csv(results_df_scheetz, file = "/Users/zhou/Desktop/IC大四/M4R/real data/results_scheetz.csv", row.names = FALSE)

# 首先将结果转换为数据框  
grouping_scheetz <- do.call(rbind, lapply(results_scheetz$grouping_info, as.data.frame))  

write.csv(grouping_scheetz, file = "/Users/zhou/Desktop/IC大四/M4R/real data/grouping_scheetz.csv", row.names = FALSE)