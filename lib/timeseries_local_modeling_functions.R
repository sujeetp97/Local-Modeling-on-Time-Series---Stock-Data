to_Input_Matrices <- function(y, n, n_is_cols = T){
  N <-  NROW(y)
  # cat(paste0("\nInput -> \t|\ty : [", N, "]\t|\t n = ", n))
  start <- N
  if(n_is_cols == T){
    end <- n+1  
    counter_end <- n
  }else{
    end <- N-n+1
    counter_end <- N-n
  }
  
  Y <- y[start:end]
  
  X <- c()
  counter <- 1
  for(index in 1:counter_end){
    X <- c(X, y[ (start-counter) : (end - counter) ])
    counter <- counter + 1
  }
  
  if(n_is_cols == T){
    X <- matrix(data = X, ncol = n, byrow = FALSE)
  }
  else{
    X <- matrix(data = X, nrow = n, byrow = FALSE)
  }
  
  
  # cat(paste0("\nOutput\t|\tY : [", NROW(Y), ", 1]\t|\tX : [", NROW(X), ", ", NCOL(X), "]"))
  
  return(list(
    'X' = X,
    'Y' = Y
  ))
  
}

get_sorted_Xindices_by_xq <- function(X, xq){
  x_delta_mean <- lapply(abs(X-xq), FUN = mean)
  x_delta_dt <- as.data.table(data.frame("COL" = names(x_delta_mean), "DELTA" = as.vector(unlist(x_delta_mean))))
  x_delta_dt <- x_delta_dt[order(DELTA, decreasing = FALSE)]
  return(x_delta_dt)
}

get_constant_eLOO <- function(train_output){
  pred_errors <- c()
  for(index in c(1:NROW(train_output))){
    
    prediction <- colMeans(x = train_output[!(index)], na.rm = T)
    pred_errors <- c(pred_errors, mean(abs(as.numeric(prediction - train_output[index])), na.rm = T))
    
  }
  loo_error <- sum(pred_errors)/NROW(train_output)
  
  return(loo_error)
}

get_constant_model_LOO_error_results <- function(y, n, nearest_neighbors, h, offset = 0){
  
  if(n < h){
    return("n cannot be less than h")
  }
  
  input_matrix <- to_Input_Matrices(y = y, n = n, n_is_cols = F) # create input and output state-space resconstructed matrices
  
  LOO_error <- c()
  neighbors <- c()
  
  
  X <- as.data.table(data.frame(input_matrix$X))
  setnames(X, old = names(X), new = as.character(c(1:NCOL(X))))
  if(nearest_neighbors > (NCOL(X))){
    nearest_neighbors <- NCOL(X)
  }
  
  xq <- input_matrix$Y
  x_delta_dt <- get_sorted_Xindices_by_xq(X = X, xq = xq)
  X[, "0" := xq] # To get target values easily
  
  
  for(k in c(2:nearest_neighbors)){
    # cat(" ")
    # cat(k)
    closest_vectors <- X[, as.vector(x_delta_dt$COL[1:k]), with = F]
    closest_vectors_cols <- colnames(closest_vectors)
    target_cols <- as.character(as.numeric(closest_vectors_cols) - 1)
    train_output <- X[(1+offset):((1+(h-1)+offset)), target_cols, with = F]
    train_output <- transpose(train_output)
    LOO_error <- c(LOO_error, get_constant_eLOO(train_output = train_output))
    neighbors <- c(neighbors, k)
  }
  
  LOO_results_table <- as.data.table(data.frame(neighbors, LOO_error))
  
  return(LOO_results_table)
}


get_constant_model_LOO_error_results_for_n_and_k <- function(y, n, k, h, offset = 0, check_all_n = T){
  
  
  LOO_results <- get_constant_model_LOO_error_results(y = y, n = h, nearest_neighbors = k, h = h, offset = offset)
  LOO_results[, n := h]
  
  if(NROW(n) == 1 & n[1] == h){
    return(LOO_results)
  }else if (min(n) < 1)(
    return("Incorrect value for embedding dimension n")
  )
  
  if((check_all_n == T) & (NROW(n) == 1)){
    n <- c(h:n)
  }
  
  for(index in n){
    # cat(paste0("\n---------------------- Embedding Dimension : ", index, "-----------------------\n"))
    results <- get_constant_model_LOO_error_results(y = y, n = index, nearest_neighbors = k, h = h, offset = offset)
    results[, n := index]
    LOO_results <- rbindlist(l = list(LOO_results, results), use.names = T, fill = F, idcol = NULL)
    
  }
  
  return(LOO_results)
  
}

predict_time_series_with_nn_constant_model <- function(y, n, k, h, offset = 0){
  input_matrix <- to_Input_Matrices(y = y, n = n, n_is_cols = F) # create input and output state-space resconstructed matrices
  X <- as.data.table(data.frame(input_matrix$X))
  setnames(X, old = names(X), new = as.character(c(1:NCOL(X))))
  xq <- input_matrix$Y
  x_delta_dt <- get_sorted_Xindices_by_xq(X = X, xq = xq)
  X[, "0" := xq] # To get target values easily
  
  closest_vectors <- X[, as.vector(x_delta_dt$COL[1:k]), with = F]
  closest_vectors_cols <- colnames(closest_vectors)
  target_cols <- as.character(as.numeric(closest_vectors_cols) - 1)
  
  train_output <- X[(1+offset):((1+(h-1)+offset)), target_cols, with = F]
  train_output <- transpose(train_output)
  
  prediction <- colMeans(x = train_output, na.rm = T)
  
  
  return(prediction)
  
}



get_single_step_multi_output_predictions <- function(y, n, k, h, n_threshold, k_threshold, type = 'Recursive', update_k_threshold = F){
  cat(paste0("\nStarting '", type, "' Multi Step Single Output Process for Step size:", h, " ..."))
  predictions <- c()
  optimal_ns <-c()
  optimal_ks <- c()
  observed <- c()
  optimal_LOO_errors <- c()
  
  optimal_n <- 0
  offset <- 0
  check_all_n <- T
  for(h_index in c(1:h)){
    
    
    if(h_index == 1){
      cat(paste0("\nFinding the Optimal Embedding Dimension"))
    }
    
    results <- get_constant_model_LOO_error_results_for_n_and_k(y = y, n = n, k = k, h = 1, offset = offset, check_all_n = check_all_n)
    results <- results[order(n, neighbors, decreasing = F)]
    
    if(h_index == 1){
      results <- results[(n >= n_threshold) & (neighbors >= k_threshold)]
      optimal_n <- max(results[LOO_error == min(results$LOO_error)]$n)
      cat(paste0("\nOptimal Embedding Dimension Found - ", optimal_n))
      if(type == "DirRec" & (n < ((optimal_n + h) - 1))){
        return(paste0("Because type = DirRec \nFor Optimal n: ", optimal_n, "and h: ", h, " n has to be >= ", ((optimal_n + h) - 1)))
      }
    }else{
      results <- results[n == optimal_n]
    }
    
    
    if(h_index != 1 & update_k_threshold == T){
      results <- results[neighbors >= k_threshold]
    }
    optimal_k <- max(results[LOO_error == min(results$LOO_error)]$neighbors)
    
    optimal_LOO_error <- min(results[LOO_error == min(results$LOO_error)]$LOO_error)
    cat(paste0("\n", h_index, ".\t", "Optimal Embedding Dimension: ", optimal_n, "\t Optimal Neighbors: ", optimal_k))  
    
    
    prediction <- predict_time_series_with_nn_constant_model(y = y, n = optimal_n, k = optimal_k, h = 1, offset = offset)
    predictions <- c(predictions, prediction)
    optimal_ns <- c(optimal_ns, optimal_n)
    optimal_ks <- c(optimal_ks, optimal_k)
    optimal_LOO_errors <- c(optimal_LOO_errors, optimal_LOO_error)
    
    if(type == "Direct"){
      offset <- offset + 1
    }else{
      y <- c(y, prediction)  
    }
    
    if(h_index == 1 & update_k_threshold == T){
      k_threshold <- optimal_k
    }
    
    check_all_n <- F
    if(type == "DirRec"){
      optimal_n <- optimal_n + 1
    }
    n <- optimal_n
    
  }
  
  recursive_single_step_results <- as.data.table(data.frame(
    predictions, optimal_ns, optimal_ks, optimal_LOO_errors, forecast_step = c(1:h)
  ))
  
  return(recursive_single_step_results)
  
}

get_MIMO_predictions <- function(y = y, n = n, k = k, h = h, n_threshold, k_threshold){
  cat(paste0("\nStarting MIMO process for Step size:", h, " ..."))
  cat(paste0("\nFinding the Optimal Embedding Dimension"))
  MIMO_results <- get_constant_model_LOO_error_results_for_n_and_k(y = y, n = n, k = k, h = h, offset = 0)
  MIMO_results <- MIMO_results[n >= n_threshold & neighbors >= k_threshold]
  optimal_n <- min(MIMO_results[LOO_error == min(MIMO_results$LOO_error)]$n)
  cat(paste0("\nOptimal Embedding Dimension Found - ", optimal_n))
  optimal_k <- max(MIMO_results[LOO_error == min(MIMO_results$LOO_error)]$neighbors)
  cat(paste0("\nOptimal Num of Neighbors Found - ", optimal_k))
  LOO_error <- min(MIMO_results[LOO_error == min(MIMO_results$LOO_error)]$LOO_error)
  
  MIMO_predictions <- predict_time_series_with_nn_constant_model(y = y, n = optimal_n, k = optimal_k, h = h, offset = 0)
  
  return(list(
    "predictions" = MIMO_predictions,
    "optimal_n" = optimal_n,
    "optimal_k" = optimal_k,
    "LOO_error" = LOO_error
  ))
  
}

get_comparative_results <- function(training, testing, n, k, h, n_threshold, k_threshold, update_k_threshold = F){
  
  if(NROW(testing) != h){
    cat("'testing' should be a univariate time series of length h")
    return()
  }
  
  recursive_predictions <- get_single_step_multi_output_predictions(y = training, n = n, k = k, h = h, n_threshold = n_threshold, k_threshold = k_threshold, type = 'Recursive', update_k_threshold = update_k_threshold)
  
  direct_predictions <- get_single_step_multi_output_predictions(y = training, n = n, k = k, h = h, n_threshold = n_threshold, k_threshold = k_threshold, type = 'Direct', update_k_threshold = update_k_threshold)
  
  dirRec_predictions <- get_single_step_multi_output_predictions(y = training, n = n, k = k, h = h, n_threshold = n_threshold, k_threshold = k_threshold, type = 'DirRec', update_k_threshold = update_k_threshold)
  
  MIMO_predictions <- get_MIMO_predictions(y = training, n = n, k = k, h = h, n_threshold = n_threshold, k_threshold = k_threshold)
  
  
  predictions <- c(recursive_predictions$predictions, direct_predictions$predictions, dirRec_predictions$predictions, as.vector(unlist(MIMO_predictions$predictions)))
  prediction_err <- c(recursive_predictions$optimal_LOO_errors, direct_predictions$optimal_LOO_errors, dirRec_predictions$optimal_LOO_errors, rep(MIMO_predictions$LOO_error, h))
  prediction_type <- rep(c("Recursive", "Direct", "DirRec", "MIMO"), each = h)
  prediction_optimal_n <- c(recursive_predictions$optimal_ns, direct_predictions$optimal_ns, dirRec_predictions$optimal_ns, rep(MIMO_predictions$optimal_n, h))
  prediction_optimal_k <- c(recursive_predictions$optimal_ks, direct_predictions$optimal_ks, dirRec_predictions$optimal_ks, rep(MIMO_predictions$optimal_k, h))
  
  prediction_comparison = as.data.table(data.frame(
    "index" = rep(c(1:h), 4),
    "observed" = rep(testing, 4),
    predictions, prediction_err, prediction_type, prediction_optimal_n, prediction_optimal_k
  ))
  
  prediction_comparison[, delta := abs(predictions - observed)]
  prediction_comparison_summary <- prediction_comparison[, .(delta_mean = mean(delta), delta_sd = sd(delta)), by = prediction_type]
  
  pd <- position_dodge(0.05)
  comp_plot <- ggplot(data = prediction_comparison, aes(x = index, y = predictions, group = prediction_type, color = prediction_type)) + 
    geom_errorbar(aes(ymin = predictions-prediction_err, ymax = predictions+prediction_err), position = pd) + 
    geom_line(position = pd) + 
    geom_point(position = pd, size = 3) + 
    geom_line(data = prediction_comparison[prediction_type == "Recursive"], aes(x = index, y = observed), color = "black") 
  
  return(list(
    "comparison_dt" = prediction_comparison,
    "comparison_summary" = prediction_comparison_summary,
    "comparison_plot" = comp_plot
  ))
}


