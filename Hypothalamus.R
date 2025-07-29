#Create dataset
source("Load and clean data with imputations.R")

#Load additional libraries
library("glmnet")
library("janitor")
library("ggtext")
library('sensemakr')
library('lm.beta')
library("ggpubr")

HTH_data <- PFAS_FLICA_covar_data[, c(3:10, 20:24, 26, 28:29, 49)] %>% 
  filter(!is.na(Hypothalamus)) %>% 
  select(Hypothalamus, everything())

#Correlation matrix and plot ---------------------------------------------------
cor_matrix <- as.matrix(cor(HTH_data, use='pairwise.complete.obs', method = 'spearman'))

#ordered bar plot of correlation coefficients
orange_blue_colour_function <- colorRampPalette(c("#EB984E", "white", "#2980B9"))
orange_blue_gradient <- orange_blue_colour_function(8)

HTH_PFAS_correlations <- cor_matrix[1, 2:9] %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename("PFAS" = "rowname", "coef" = ".")

HTH_PFAS_correlations %>% 
  ggplot(data = ., aes(x=reorder(PFAS, coef), y=coef)) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(fill = reorder(PFAS, coef))) + #Colour them in order
  
  scale_fill_manual(values = orange_blue_gradient) + #as colour scheme, use the blues_gradient object that I created
  
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + #remove legend
  labs(y = "Correlation Coefficient", title = "Hypothalamus MD") + #format titles
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) +
  geom_hline(yintercept =0) 

ggsave("plots/vectors/HTH_correlations.svg", units = "px", width = 1100, height = 1100)

#Elastic net regression with cross validation ----------------------------------

#Follow the ame general principle as the Elastic net regression.R script for the 10 brain components

#Create data
xHTH <- HTH_data[, -1] %>% as.matrix()
yHTH <- HTH_data$Hypothalamus

#Display small values in decimals
options(scipen=10000)

# Set seed for reproducibility
set.seed(1)

# Create random partitions of the data
folds = sample(rep(1:10, length = nrow(HTH_data)))

#Create an empty dataframe to store all results inside
HTH_model_results <- as.data.frame(matrix(ncol = 0, nrow = 0))
HTH_coefficients = as.data.frame(matrix(ncol = 0, nrow = 0))

for (k in 1:10){
  #Create an ampty dataframe to store the lambda with the lowest MSE at each alpha
  alpha_lambda <- as.data.frame(matrix(ncol = 0, nrow = 0))
  
  #Run an inner loop to optimize lambda 100 times and choose the mean
  for (a in seq(0, 1, length.out = 11)){
    
    #Create an ampty dataframe to store the lambda with the lowest MSE
    lambdas_and_errors = as.data.frame(matrix(ncol = 0, nrow = 0))
    
    for (i in 1:100){ 
      set.seed(i) #set the seed 100 times to make the model results reproducible
      cv_model <- cv.glmnet(xHTH[folds!=k,], yHTH[folds!=k], alpha = a) #specify the model
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
  optimal_model_train <- glmnet(xHTH[folds!=k,], yHTH[folds!=k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda) #model with training data
  optimal_model_test <- glmnet(xHTH[folds==k,], yHTH[folds==k], alpha = optimal_parameters$alpha, lambda = optimal_parameters$lambda)  #model with test data
  
  # extract deviation ratios
  train_dev_ratio <- optimal_model_train$dev.ratio
  test_dev_ratio <- optimal_model_test$dev.ratio
  
  # make predictions and calculate mse and r2
  y_test <- yHTH[folds==k] #define test y data
  y_predicted <- predict(optimal_model_train, newx = xHTH[folds==k,], s = optimal_parameters$lambda) #predict new y from model and new x
  prediction_rmse <- sqrt(mean((y_test - y_predicted)^2))
  prediction_r2 <- 1 - (sum((y_test - y_predicted)^2) / sum((y_test - mean(y_test))^2))
  
  #Save model results
  k_results <- cbind(k, alpha = optimal_parameters$alpha,lambda = optimal_parameters$lambda, train_dev_ratio, test_dev_ratio, prediction_rmse, prediction_r2)
  HTH_model_results <- rbind(HTH_model_results, k_results)
  
  k_coefficients <- coef(optimal_model_test)[-1,] %>% #We add [-1,] to remove the coefficient of the intercept
    matrix() %>% 
    as.data.frame() %>% 
    mutate(Variable = colnames(HTH_data[2:17])) %>% 
    rename("coef" = "V1") %>% 
    select(Variable, coef) %>% 
    t() %>% 
    row_to_names(1) 
  HTH_coefficients <- rbind(HTH_coefficients, k_coefficients)
  
  print(str_c("Model parameters saved with optimal alpha and lambda for fold k = ", k))
  
}

rangey <- max(HTH_data$Hypothalamus) - min(HTH_data$Hypothalamus)
HTH_model_results <- HTH_model_results %>% 
  mutate(nrmse = prediction_rmse/rangey)

HTH_coefficients_plot <- HTH_coefficients %>% 
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

#Create a common colour palette 
red_blue_colour_function <- colorRampPalette(c("#AB1C08", "white", "#0815AB"))
red_blue_gradient <- red_blue_colour_function(8)

#Get model parameter means and standard errors across all 10 folds
colMeans(HTH_model_results)
(sd(HTH_model_results$train_dev_ratio))/sqrt(10)
(sd(HTH_model_results$nrmse))/sqrt(10)



HTH_coefficients_plot %>% 
  ggplot(data = .) + 
  geom_bar(stat = "identity", colour = "black", alpha = 0.95, aes(x=reorder(PFAS, Coefficient), y=Coefficient, fill = reorder(PFAS, Coefficient))) + 
  geom_text(aes(x = PFAS, y = -0.0000125, size = 10, label = Models)) +
  geom_errorbar(aes(x=PFAS, ymin=(Coefficient-SEM), ymax=(Coefficient+SEM), width=0.5)) +
  scale_fill_manual(values = red_blue_gradient) + 
  labs(y = "Regularized Coefficient", title = "Hypothalamus MD") + 
  annotate("text", label = c("Dev ratio = 0.37 (±0.021)", "NRMSE = 0.19 (±0.019)"), x=1, y=c(0.0000065, 0.00000475), size = 4, hjust = 0) +
  
  geom_hline(yintercept =0) +   
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.position = "none", 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

ggsave("plots/vectors/HTH_elastic_net_regression.svg", units = "px", width = 1100, height = 1100)
#linear models with individual PFAS - PFOS-Br, PFHxS-Br -------------

#PFOS_Br

#Plot
p11 <- HTH_data %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "Hypothalamus MD ~ Maternal PFOS-Br", 
       subtitle = "R^2^=0.10, β=-0.29 (-0.57- -0.041), adjusted p=0.026", 
       x = "PFOS-Br", y = "MD (mm^2^/s, autoscaled)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_markdown(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 


#Minimally-adjusted model
PFOS_Br_vs_HTH_min_lm <- HTH_data %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_HTH_min_lm) 
lm.beta(PFOS_Br_vs_HTH_min_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFOS_Br_vs_HTH_full_lm <- HTH_data %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_HTH_full_lm) 
lm.beta(PFOS_Br_vs_HTH_full_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOS_Br)
sy <- sd(PFAS_FLICA_covar_data$Hypothalamus)
confint(PFOS_Br_vs_HTH_full_lm)[2,]
confint(PFOS_Br_vs_HTH_full_lm)[2,]*(sx/sy)

#PFHxS_Br

#Plot
p12 <- HTH_data %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "Hypothalamus MD ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.010, β=0.094 (-0.14-0.29), adjusted p=0.50", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_HTH_min_lm <- HTH_data %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_HTH_min_lm) 
lm.beta(PFHxS_Br_vs_HTH_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_HTH_full_lm <- HTH_data %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex  + Gravidity, data = .)
summary(PFHxS_Br_vs_HTH_full_lm)
lm.beta(PFHxS_Br_vs_HTH_full_lm) %>% coef()
sensemakr::partial_r2(PFHxS_Br_vs_HTH_full_lm)

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFHxS_Br)
sy <- sd(PFAS_FLICA_covar_data$Hypothalamus)
confint(PFHxS_Br_vs_HTH_full_lm)[2,]
confint(PFHxS_Br_vs_HTH_full_lm)[2,]*(sx/sy)

#Save plots -------------------------------------

ggarrange(p11, p12, nrow = 1, ncol = 2)
ggsave("plots/vectors/HTH_PFAS_scatterplots.svg", units = "px", width = 2650, height = 1250)
