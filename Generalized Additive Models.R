library("tidyverse")
library("mgcv")

PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% 
  rename("WQS" = "Total (WQS) PFAS", "Child_Age" = "Child Age", "Maternal_Age" = "Maternal Age")

outlier_dataset_visual <- outlier_dataset_visual %>% 
  rename("WQS" = "Total (WQS) PFAS", "Child_Age" = "Child Age", "Maternal_Age" = "Maternal Age")


#GAM for full dataset ----------------------------------------------------------------------------
WQS_GAM_results <- data.frame()
for (k in 32:41){
  set.seed(1)
  y <- PFAS_FLICA_covar_data[, k]
  y <- unlist(y)
  
  print(colnames(PFAS_FLICA_covar_data)[k])
  
  #give the general structure for the linear and gam models
  model_lm <- gam(y ~ WQS + Child_Age + Sex + Maternal_Age, data = PFAS_FLICA_covar_data)
  model_gam <- gam(y ~ s(WQS) + Child_Age + Sex + Maternal_Age, data = PFAS_FLICA_covar_data)
  
  #extract the  p-values
  pval_lm <- (summary(model_lm)$p.table %>% as.data.frame())[2, 4]
  pval_gam <- (summary(model_gam)$s.table %>% as.data.frame())[1, 4]
  
  #test lm v gam
  anova <- anova(model_lm, model_gam, test = "Chisq")
  test <- anova %>% as.data.frame()
  p_test <- test[2, 5]
  better_model <- if_else(test[2, 2] < test[1, 2], "GAM", "LM")
  
  #variable name
  component <- colnames(PFAS_FLICA_covar_data)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(component) %>%  
    mutate(pval_lm = pval_lm) %>% 
    mutate(pval_gam = pval_gam) %>% 
    mutate(p_test = p_test) %>% 
    mutate(better_model = better_model)
  
  #for each iteration, append the results to our results dataframe
  WQS_GAM_results <- rbind(WQS_GAM_results, results_k)
  
}

WQS_GAM_results <- WQS_GAM_results %>% 
  arrange(pval_gam) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval_gam*nrow(.)/pval_rank > 1, 1.0, pval_gam*nrow(.)/pval_rank))

#GAM for 2 datapoints removed ----------------------------------------------------------

WQS_outlier_GAM_results <- data.frame()
for (k in 32:41){
  set.seed(1)
  y <- outlier_dataset_visual[, k]
  y <- unlist(y)
  
  print(colnames(outlier_dataset_visual)[k])
  
  #give the general structure for the linear and gam models
  model_lm <- gam(y ~ WQS + Child_Age + Sex + Maternal_Age, data = outlier_dataset_visual)
  model_gam <- gam(y ~ s(WQS) + Child_Age + Sex + Maternal_Age, data = outlier_dataset_visual)
  
  #extract the  p-values
  pval_lm <- (summary(model_lm)$p.table %>% as.data.frame())[2, 4]
  pval_gam <- (summary(model_gam)$s.table %>% as.data.frame())[1, 4]
  
  #test lm v gam
  anova <- anova(model_lm, model_gam, test = "Chisq")
  test <- anova %>% as.data.frame()
  p_test <- test[2, 5]
  better_model <- if_else(test[2, 2] < test[1, 2], "GAM", "LM")
  
  #variable name
  component <- colnames(outlier_dataset_visual)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(component) %>%  
    mutate(pval_lm = pval_lm) %>% 
    mutate(pval_gam = pval_gam) %>% 
    mutate(p_test = p_test) %>% 
    mutate(better_model = better_model)
  
  #for each iteration, append the results to our results dataframe
  WQS_outlier_GAM_results <- rbind(WQS_outlier_GAM_results, results_k)
  
}

WQS_outlier_GAM_results <- WQS_outlier_GAM_results %>% 
  arrange(pval_gam) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval_gam*nrow(.)/pval_rank > 1, 1.0, pval_gam*nrow(.)/pval_rank))

#Plots for GAM full dataset ----------------------------------------------------

p1 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC1)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC1", 
        subtitle = "p=0.082, pFDR=0.41",
       x = "Total (WQS) PFAS", y = "IC") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p2 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC2)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC2", 
       subtitle = "p=0.72, pFDR=0.72",
       x = "Total (WQS) PFAS", y = "IC2") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p3 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC3)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC3", 
       subtitle = "p=0.23, pFDR=0.57",
       x = "Total (WQS) PFAS", y = "IC3") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p4 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC4", 
       subtitle = "p=0.32, pFDR=0.45",
       x = "Total (WQS) PFAS", y = "IC4") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p5 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC5", 
       subtitle = "p=0.15, pFDR=0.49",
       x = "Total (WQS) PFAS", y = "IC5") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p6 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC6)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC6", 
       subtitle = "p=0.35, pFDR=0.44",
       x = "Total (WQS) PFAS", y = "IC") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p7 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC7)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC7", 
       subtitle = "p=0.31, pFDR=0.51",
       x = "Total (WQS) PFAS", y = "IC7") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p8 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC8)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC8", 
       subtitle = "p=0.016, pFDR=0.16",
       x = "Total (WQS) PFAS", y = "IC8") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p9 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC9", 
       subtitle = "p=0.25, pFDR=0.50",
       x = "Total (WQS) PFAS", y = "IC9") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p10 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = WQS, y = IC10)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "#EB984E") +
  labs(title = "IC10", 
       subtitle = "p=0.42, pFDR=0.47",
       x = "Total (WQS) PFAS", y = "IC10") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 5)
ggsave("plots/vectors/GAM_WQS_all_data.svg", units = "px", width = 5000, height = 2000)

#Plots for GAM 2 datapoints removed ----------------------------------------------------


p1 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC1)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC1", 
       subtitle = "p=0.12, pFDR=0.31",
       x = "Total (WQS) PFAS", y = "IC") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p2 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC2)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC2", 
       subtitle = "p=0.31, pFDR=0.39",
       x = "Total (WQS) PFAS", y = "IC2") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p3 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC3)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC3", 
       subtitle = "p=0.12, pFDR=0.40",
       x = "Total (WQS) PFAS", y = "IC3") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p4 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC4", 
       subtitle = "p=0.17, pFDR=0.33",
       x = "Total (WQS) PFAS", y = "IC4") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p5 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC5", 
       subtitle = "p=0.020, pFDR=0.20",
       x = "Total (WQS) PFAS", y = "IC5") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p6 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC6)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC6", 
       subtitle = "p=0.31, pFDR=0.45",
       x = "Total (WQS) PFAS", y = "IC") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p7 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC7)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC7", 
       subtitle = "p=0.095, pFDR=0.48",
       x = "Total (WQS) PFAS", y = "IC7") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p8 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC8)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC8", 
       subtitle = "p=0.52, pFDR=0.52",
       x = "Total (WQS) PFAS", y = "IC8") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p9 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC9", 
       subtitle = "p=0.17, pFDR=0.28",
       x = "Total (WQS) PFAS", y = "IC9") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

p10 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = WQS, y = IC10)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "gam", size = 2, colour = "Maroon") +
  labs(title = "IC10", 
       subtitle = "p=0.43, pFDR=0.48",
       x = "Total (WQS) PFAS", y = "IC10") +
  theme(plot.title = element_text(hjust = 0.5, size=16),
        plot.subtitle = element_text(hjust = 0.5, size=14), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 5)
ggsave("plots/vectors/GAM_WQS_outliers_removed.svg", units = "px", width = 5000, height = 2000)
