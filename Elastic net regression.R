#Create dataset
source("clean scripts/Load and clean data with imputations.R")

elastic_net_regression_data <- PFAS_FLICA_covar_data %>% 
  select(3:10, 20:24, 26, 28:29, 32:41) %>% 
  mutate_all(as.numeric) %>% 
  filter(!is.na(IC1))

#Load additional libraries
library("glmnet")
library("janitor")
library("ggpubr")
library("ggbreak")

#Display small values in decimals
options(scipen=10000)

#Set seed for reproducibility
set.seed(1)

#Split dataset into 10 folds and set up for ridge regressions ------------------

#Create outcome (y) vectors to predict
yIC1 <- elastic_net_regression_data$IC1
yIC2 <- elastic_net_regression_data$IC2
yIC3 <- elastic_net_regression_data$IC3
yIC4 <- elastic_net_regression_data$IC4
yIC5 <- elastic_net_regression_data$IC5
yIC6 <- elastic_net_regression_data$IC6
yIC7 <- elastic_net_regression_data$IC7
yIC8 <- elastic_net_regression_data$IC8
yIC9 <- elastic_net_regression_data$IC9
yIC10 <- elastic_net_regression_data$IC10

#Create prediction (x) matrix:
x <- as.matrix(elastic_net_regression_data)[, -(17:26)]

# Create random partitions of the data. In this case we wish to to 10-fold CV, so partition the data into 10 folds.
folds = sample(rep(1:10, length = nrow(elastic_net_regression_data)))

# View the number of observations in each fold
table(folds)

# View / print to the console any specified folds. In the example below, 
# I print the 3rd fold (arbitrarily as an example):
x[folds==3,]
yIC1[folds==3]

set.seed(1)
#IC1 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC1_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC1_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC1[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC1[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC1[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC1[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC1_model_results <- rbind(IC1_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC1_coefficients <- rbind(IC1_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC1_coefficients_plot <- IC1_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC2 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC2_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC2_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC2[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC2[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC2[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC2[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC2_model_results <- rbind(IC2_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC2_coefficients <- rbind(IC2_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC2_coefficients_plot <- IC2_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC3 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC3_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC3_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC3[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC3[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC3[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC3[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC3_model_results <- rbind(IC3_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC3_coefficients <- rbind(IC3_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC3_coefficients_plot <- IC3_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC4 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC4_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC4_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC4[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC4[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC4[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC4[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC4_model_results <- rbind(IC4_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC4_coefficients <- rbind(IC4_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC4_coefficients_plot <- IC4_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC5 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC5_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC5_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC5[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC5[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC5[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC5[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC5_model_results <- rbind(IC5_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC5_coefficients <- rbind(IC5_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC5_coefficients_plot <- IC5_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC6 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC6_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC6_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC6[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC6[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC6[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC6[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC6_model_results <- rbind(IC6_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC6_coefficients <- rbind(IC6_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC6_coefficients_plot <- IC6_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC7 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC7_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC7_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC7[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC7[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC7[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC7[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC7_model_results <- rbind(IC7_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC7_coefficients <- rbind(IC7_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC7_coefficients_plot <- IC7_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC8 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC8_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC8_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC8[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC8[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC8[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC8[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC8_model_results <- rbind(IC8_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC8_coefficients <- rbind(IC8_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC8_coefficients_plot <- IC8_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC9 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC9_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC9_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC9[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC9[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC9[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC9[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC9_model_results <- rbind(IC9_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC9_coefficients <- rbind(IC9_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC9_coefficients_plot <- IC9_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

set.seed(1)
#IC10 ---------------------------------------------------------------------------
#Create an empty dataframe to store all results inside
IC10_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
IC10_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(x[folds!=k,], yIC10[folds!=k], alpha = a) #specify the model
      best_lambda <- cv_model$lambda.min #extract lambda value with lowest MSE
      cv_error <- min(cv_model$cvm) #extract the MSE of lambda.min
      lambda_and_error <- cbind(best_lambda, cv_error)
      lambdas_and_errors <- rbind(lambdas_and_errors, lambda_and_error)
    }
    
    print(str_c("Nested CV complete for alpha = ", a))
    alpha <- a
    lambda <- mean(lambdas_and_errors[,1])
    error <- mean(lambdas_and_errors[,2])
    alpha_lambda_a <- cbind(alpha, lambda, error)
    
    alpha_lambda <- rbind(alpha_lambda, alpha_lambda_a)
  }
  
  #We have identified the lambda with the lowest MSE for each alpha. 
  # Now we identify at which alpha we have the overall lowest MSE.
  optimal_parameters <- alpha_lambda %>% filter(error== min(error))
  
  print(str_c("Optimal alpha and lambda chosen by CV for fold k = ", k))
  
  #set seed for reproducibility
  set.seed(1) 
  
  # define models
  optimal_model_train <- glmnet(x[folds!=k,], yIC10[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(x[folds==k,], yIC10[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yIC10[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = x[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  IC10_model_results <- rbind(IC10_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(elastic_net_regression_data[1:16])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  IC10_coefficients <- rbind(IC10_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

#Extract the number of models and mean coefficient of each PFAS across all 10 folds

IC10_coefficients_plot <- IC10_coefficients %>% 
  mutate_all(as.numeric) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(Models = rowSums(.[1:10] != 0)) %>% #Number of models this variable included in (possible range 0:10)
  mutate(Coefficient = rowMeans(.[1:10])) %>% #Mean  of coefficient across all 10 models   
  mutate(SEM = apply(.[1:10], 1, sd)/sqrt(10)) %>% #SEM of coefficient across all 10 models
  head(8) %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname") %>% 
  select(PFAS, Models, Coefficient, SEM)

#Get normalised RMSE ---------------------------------------------------------
rangey <- max(elastic_net_regression_data$IC1) - min(elastic_net_regression_data$IC1)
IC1_model_results <- IC1_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC2) - min(elastic_net_regression_data$IC2)
IC2_model_results <- IC2_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC3) - min(elastic_net_regression_data$IC3)
IC3_model_results <- IC3_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC4) - min(elastic_net_regression_data$IC4)
IC4_model_results <- IC4_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC5) - min(elastic_net_regression_data$IC5)
IC5_model_results <- IC5_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC6) - min(elastic_net_regression_data$IC6)
IC6_model_results <- IC6_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC7) - min(elastic_net_regression_data$IC7)
IC7_model_results <- IC7_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC8) - min(elastic_net_regression_data$IC8)
IC8_model_results <- IC8_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC9) - min(elastic_net_regression_data$IC9)
IC9_model_results <- IC9_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

rangey <- max(elastic_net_regression_data$IC10) - min(elastic_net_regression_data$IC10)
IC10_model_results <- IC10_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)


#Plots -------------------------------------------------------------------------

#Create a common colour palette 
red_blue_colour_function <- colorRampPalette(c("#AB1C08", "white", "#0815AB"))
red_blue_gradient <- red_blue_colour_function(8)

#Get model parameter means and standard errors across all 10 folds
colMeans(IC1_model_results)
(sd(IC1_model_results$train_dev_ratio))/sqrt(10)
(sd(IC1_model_results$nrmse))/sqrt(10)
colMeans(IC2_model_results)
(sd(IC2_model_results$train_dev_ratio))/sqrt(10)
(sd(IC2_model_results$nrmse))/sqrt(10)
colMeans(IC3_model_results)
(sd(IC3_model_results$train_dev_ratio))/sqrt(10)
(sd(IC3_model_results$nrmse))/sqrt(10)
colMeans(IC4_model_results)
(sd(IC4_model_results$train_dev_ratio))/sqrt(10)
(sd(IC4_model_results$nrmse))/sqrt(10)
colMeans(IC5_model_results)
(sd(IC5_model_results$train_dev_ratio))/sqrt(10)
(sd(IC5_model_results$nrmse))/sqrt(10)
colMeans(IC6_model_results)
(sd(IC6_model_results$train_dev_ratio))/sqrt(10)
(sd(IC6_model_results$nrmse))/sqrt(10)
colMeans(IC7_model_results)
(sd(IC7_model_results$train_dev_ratio))/sqrt(10)
(sd(IC7_model_results$nrmse))/sqrt(10)
colMeans(IC8_model_results)
(sd(IC8_model_results$train_dev_ratio))/sqrt(10)
(sd(IC8_model_results$nrmse))/sqrt(10)
colMeans(IC9_model_results)
(sd(IC9_model_results$train_dev_ratio))/sqrt(10)
(sd(IC9_model_results$nrmse))/sqrt(10)
colMeans(IC10_model_results)
(sd(IC10_model_results$train_dev_ratio))/sqrt(10)
(sd(IC10_model_results$nrmse))/sqrt(10)

#Build plots

IC1_plot <- IC1_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC1") + 
  annotate("text", label = c("Dev ratio = 0.12 (±0.03)", "NRMSE = 0.24 (±0.023)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +   
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(size=14),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC2_plot <- IC2_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.605, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC2") + 
  annotate("text", label = c("Dev ratio = 0.06 (±0.017)", "NRMSE = 0.16 (±0.023)"), x=0.6, y=c(0.28, 0.18), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +   
  scale_y_continuous(breaks = -0.6, limits = c(-0.625, 0.3)) +
  scale_y_break(c(-0.27, -0.18), ticklabels = c(-0.2, -0.1, 0, 0.1, 0.2), , space = 0.25) +
  scale_y_break(c(-0.52, -0.29), ticklabels = c(-0.3), , space = 0.25) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"),
        axis.ticks.y.left = element_blank(),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())

IC2_plot <- IC2_plot +
  theme(plot.background = element_blank(), 
        panel.background = element_blank())

#IC2 y axis needs to be broken

IC3_plot <- IC3_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC3") + 
  annotate("text", label = c("Dev ratio = 0.03 (±0.014)", "NRMSE = 0.24 (±0.025)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC4_plot <- IC4_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.335, size = 4, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC4") + 
  annotate("text", label = c("Dev ratio = 0.31 (±0.04)", "NRMSE = 0.22 (±0.015)"), x=0.6, y=c(0.28, 0.21), size = 3.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.335, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 11.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 9), axis.text.y = element_text(size = 8),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC5_plot <- IC5_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC5") + 
  annotate("text", label = c("Dev ratio = 0.16 (±0.023)", "NRMSE = 0.18 (±0.021)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC6_plot <- IC6_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.325, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC6") + 
  annotate("text", label = c("Dev ratio = 0.046 (±0.022)", "NRMSE = 0.25 (±0.018)"), x=0.6, y=c(-0.11, -0.24), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +   
  ylim(-0.325, 3.15) +
  scale_y_break(c(0.2, 1.54), ticklabels = 1.6, space = 0.25) +
  scale_y_break(c(1.6, 3.1), ticklabels = 3.1, space = 0.25) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_text(angle=90, size=14),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank())

IC6_plot <- IC6_plot +
  theme(plot.background = element_blank(), 
        panel.background = element_blank())

#IC6 y-axis needs to be broken

IC7_plot <- IC7_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC7") + 
  annotate("text", label = c("Dev ratio = 0.021 (±0.018)", "NRMSE = 0.21 (±0.022)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC8_plot <- IC8_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC8") + 
  annotate("text", label = c("Dev ratio = 0.085 (±0.024)", "NRMSE = 0.14 (±0.028)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC9_plot <- IC9_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.335, size = 4, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC9") + 
  annotate("text", label = c("Dev ratio = 0.45 (±0.054)", "NRMSE = 0.23 (±0.030)"), x=0.6, y=c(0.28, 0.21), size = 3.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.335, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 11.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 9), axis.text.y = element_text(size = 8),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

IC10_plot <- IC10_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.225, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularised Coefficient", title = "IC10") + 
  annotate("text", label = c("Dev ratio = 0.026 (±0.018)", "NRMSE = 0.24 (±0.021)"), x=0.6, y=c(0.28, 0.21), size = 5.5, hjust = 0) +
  geom_hline(yintercept =0) +  
  ylim(-0.225, 0.295) +
  theme(plot.title = element_text(hjust = 0.5, size = 17.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 12),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

all_elastic_net_plots <- ggarrange(IC1_plot, print(IC2_plot), IC3_plot, IC4_plot, IC5_plot, 
                                   print(IC6_plot), IC7_plot, IC8_plot, IC9_plot, IC10_plot,
                                   nrow = 2, ncol = 5)

all_elastic_net_plots

ggsave("plots/vectors/elastic_net_regressions.svg", units = "px", width = 5000, height = 3500, bg="transparent")

IC4_IC9_elastic_net_plots <- ggarrange(IC4_plot, IC9_plot)

IC4_IC9_elastic_net_plots

ggsave("plots/vectors/IC4_IC9_elastic_net_plots.svg", units = "px", width = 2500, height = 1250, bg="transparent")

supplementary_elastic_net_plots <- ggarrange(IC1_plot, print(IC2_plot), IC3_plot, IC5_plot, 
                                   print(IC6_plot), IC7_plot, IC8_plot, IC10_plot,
                                   nrow = 2, ncol = 4)

supplementary_elastic_net_plots

ggsave("plots/vectors/supplementary_elastic_net_plots.svg", units = "px", width = 4000, height = 3000, bg="transparent")

