#Create dataset
library("tidyverse")
library("sensemakr")
library("lm.beta")
library("ggtext")
library("pscl")
library("ggpubr")

sdq_5y <- read_sav("data/Data_5y_MRI_SDQ.sav")[, c(1, 35:40)]

colnames(sdq_5y) <- c("FB_code", "emotional_problems", "conduct_problems",
                   "hyperactivity_inattention", "peer_relationship_problems", 
                   "prosocial_behaviour", "total_difficulties")

flica_data <- read_csv("data/FLICA_age_sex_bmi.csv")

flica_and_sdq <- left_join(flica_data, sdq_5y, by = "FB_code")

flica_and_sdq <- flica_and_sdq %>% 
  mutate(emotional_problems_binary = if_else(emotional_problems < 4, 0, 1)) %>% 
  mutate(conduct_problems_binary = if_else(conduct_problems < 3, 0, 1)) %>% 
  mutate(hyperactivity_inattention_binary = if_else(hyperactivity_inattention < 6, 0, 1)) %>% 
  mutate(peer_relationship_problems_binary = if_else(peer_relationship_problems < 3, 0, 1)) %>% 
  mutate(prosocial_behaviour_binary = if_else(prosocial_behaviour > 7, 0, 1)) %>% 
  mutate(total_difficulties_binary = if_else(total_difficulties < 14, 0, 1))

pfas_and_sdq <- read.csv("data/PFAS_FLICA_covar_data_FB_code.csv", check.names = F)[, 1:19] %>% 
  left_join(., (sdq_5y %>%  mutate(FB_code = str_sub(FB_code, end = -2))), by = "FB_code") %>% 
  mutate(emotional_problems_binary = if_else(emotional_problems < 4, 0, 1)) %>% 
  mutate(conduct_problems_binary = if_else(conduct_problems < 3, 0, 1)) %>% 
  mutate(hyperactivity_inattention_binary = if_else(hyperactivity_inattention < 6, 0, 1)) %>% 
  mutate(peer_relationship_problems_binary = if_else(peer_relationship_problems < 3, 0, 1)) %>% 
  mutate(prosocial_behaviour_binary = if_else(prosocial_behaviour > 7, 0, 1)) %>% 
  mutate(total_difficulties_binary = if_else(total_difficulties < 14, 0, 1)) %>% 
  filter(!is.na(total_difficulties_binary))

HTH_AC_codes <- read.csv("data/HTH_5Y_AC_codes.txt", header = FALSE)
FB_AC_code <- read.delim("data/5V-AC-2-FB-5V-MRI.txt") %>% 
  select(-AC) %>% 
  rename("FB_code" = ID_C, "AC_code" = X) %>% 
  mutate(AC_code = as.character(AC_code))

HTH_MD_5_values <- read.csv("data/HTH_5Y_MD_values.txt", header = FALSE)

HTH_MD_5 <- data.frame(c("AC_code" = HTH_AC_codes, 
                         "Hypothalamus" = HTH_MD_5_values)) %>% 
  rename("AC_code" = AC_code.V1, "Hypothalamus" = Hypothalamus.V1) %>% 
  mutate(AC_code = as.character(AC_code)) %>% 
  left_join(FB_AC_code, HTH_MD_5, by = "AC_code") %>% 
  select(FB_code, Hypothalamus) 

HTH_sdq <- left_join(HTH_MD_5, sdq_5y, by = "FB_code") %>% 
  mutate(emotional_problems_binary = if_else(emotional_problems < 4, 0, 1)) %>% 
  mutate(conduct_problems_binary = if_else(conduct_problems < 3, 0, 1)) %>% 
  mutate(hyperactivity_inattention_binary = if_else(hyperactivity_inattention < 6, 0, 1)) %>% 
  mutate(peer_relationship_problems_binary = if_else(peer_relationship_problems < 3, 0, 1)) %>% 
  mutate(prosocial_behaviour_binary = if_else(prosocial_behaviour > 7, 0, 1)) %>% 
  mutate(total_difficulties_binary = if_else(total_difficulties < 14, 0, 1)) %>% 
  filter(!is.na(total_difficulties_binary))


#IC5 linear models -------------------------------------------------------------
IC5_results <- data.frame()
for (k in 16:21){
  set.seed(1)
  y <- flica_and_sdq[, k]
  y <- unlist(y)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- lm(y ~ IC5, data = flica_and_sdq)
  
  #extract coefficient for the metabolite
  coef  <- as.data.frame(lm.beta(model)$standardized.coefficients)[2,1]
  
  #extract partial r2 for the metabolite
  rsq <- as.data.frame(partial_r2(model))[2,1]
  
  #extract the  p-value 
  pval <- (coef(summary(model)) %>% as.data.frame())[2, 4]
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>% 
    mutate(coef = coef) %>% 
    mutate(rsq = rsq)
  
  #for each iteration, append the results to our results dataframe
  IC5_results <- rbind(IC5_results, results_k)
  
}


#IC9 linear models -------------------------------------------------------------
IC9_results <- data.frame()
for (k in 16:21){
  set.seed(1)
  y <- flica_and_sdq[, k]
  y <- unlist(y)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- lm(y ~ IC9, data = flica_and_sdq)
  
  #extract coefficient for the metabolite
  coef  <- as.data.frame(lm.beta(model)$standardized.coefficients)[2,1]
  
  #extract partial r2 for the metabolite
  rsq <- as.data.frame(partial_r2(model))[2,1]
  
  #extract the  p-value 
  pval <- (coef(summary(model)) %>% as.data.frame())[2, 4]
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>% 
    mutate(coef = coef) %>% 
    mutate(rsq = rsq)
  
  #for each iteration, append the results to our results dataframe
  IC9_results <- rbind(IC9_results, results_k)
  
}

#IC5 logistic regression models -------------------------------------------------------------
IC5_logreg_results <- data.frame()
for (k in 22:27){
  set.seed(1)
  y <- flica_and_sdq[, k]
  y <- unlist(y)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- glm(y ~ IC5 + Sex + Age_at_scan, data = flica_and_sdq, family = binomial)
  
  #extract coefficient for the metabolite
  coef  <- as.data.frame(lm.beta(model)$standardized.coefficients)[2,1]
  
  #extract partial r2 for the model
  pseudo_rsq <- pR2(model)[4]
  
  #extract the  p-value 
  pval <- (coef(summary(model)) %>% as.data.frame())[2, 4]
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>% 
    mutate(coef = coef) %>% 
    mutate(pseudo_rsq = pseudo_rsq)
  
  #for each iteration, append the results to our results dataframe
  IC5_logreg_results <- rbind(IC5_logreg_results, results_k)
  
}


#IC9 logistic regression models -------------------------------------------------------------
IC9_logreg_results <- data.frame()
for (k in 22:27){
  set.seed(1)
  y <- flica_and_sdq[, k]
  y <- unlist(y)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- glm(y ~ IC9, data = flica_and_sdq, family = binomial)
  
  #extract coefficient for the metabolite
  coef  <- as.data.frame(lm.beta(model)$standardized.coefficients)[2,1]
  
  #extract partial r2 for the model
  pseudo_rsq <- pR2(model)[4]
  
  #extract the  p-value 
  pval <- (coef(summary(model)) %>% as.data.frame())[2, 4]
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>% 
    mutate(coef = coef) %>% 
    mutate(pseudo_rsq = pseudo_rsq)
  
  #for each iteration, append the results to our results dataframe
  IC9_logreg_results <- rbind(IC9_logreg_results, results_k)
  
}


#IC5 t-test --------------------------------------------------------------------

IC5_ttest_results <- data.frame()
for (k in 22:27){
  set.seed(1)
  x <- flica_and_sdq[, k]
  x <- unlist(x)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- t.test(IC5 ~ x, data = flica_and_sdq)
  
  #extract the t-statistic
  tstat <- model$statistic
  
  #extract the  p-value 
  pval <- model$p.value
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>%  
    mutate(tstat = tstat)
  
  #for each iteration, append the results to our results dataframe
  IC5_ttest_results <- rbind(IC5_ttest_results, results_k)
  
}

IC5_ttest_results <- IC5_ttest_results %>%  
  mutate(p_FWER = if_else(pval*nrow(.) > 1, 1.0, pval*nrow(.))) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank))

#IC9 t-test --------------------------------------------------------------------

IC9_ttest_results <- data.frame()
for (k in 22:27){
  set.seed(1)
  x <- flica_and_sdq[, k]
  x <- unlist(x)
  
  print(colnames(flica_and_sdq)[k])
  
  #give the general structure for the linear model
  model <- t.test(IC9 ~ x, data = flica_and_sdq)
  
  #extract the t-statistic
  tstat <- model$statistic
  
  #extract the  p-value 
  pval <- model$p.value
  
  #variable name
  sdq <- colnames(flica_and_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>%  
    mutate(tstat = tstat)
  
  #for each iteration, append the results to our results dataframe
  IC9_ttest_results <- rbind(IC9_ttest_results, results_k)
  
}

IC9_ttest_results <- IC9_ttest_results %>%  
  mutate(p_FWER = if_else(pval*nrow(.) > 1, 1.0, pval*nrow(.))) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank))

#HTH t-test --------------------------------------------------------------------

HTH_ttest_results <- data.frame()
for (k in 9:14){
  set.seed(1)
  x <- HTH_sdq[, k]
  x <- unlist(x)
  
  print(colnames(HTH_sdq)[k])
  
  #give the general structure for the linear model
  model <- t.test(Hypothalamus ~ x, data = HTH_sdq)
  
  #extract the t-statistic
  tstat <- model$statistic
  
  #extract the  p-value 
  pval <- model$p.value
  
  #variable name
  sdq <- colnames(HTH_sdq)[k]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pval = pval) %>%  
    mutate(tstat = tstat)
  
  #for each iteration, append the results to our results dataframe
  HTH_ttest_results <- rbind(HTH_ttest_results, results_k)
  
}

HTH_ttest_results <- HTH_ttest_results %>%  
  mutate(p_FWER = if_else(pval*nrow(.) > 1, 1.0, pval*nrow(.))) %>% 
  arrange(pval) %>% 
  mutate(pval_rank = 1:nrow(.)) %>% 
  mutate(p_FDR = if_else(pval*nrow(.)/pval_rank > 1, 1.0, pval*nrow(.)/pval_rank))

#All PFAS t-test ----------------------------------------------------------------------------

PFAS_ttest_results <- data.frame()
for (k in 26:31){
  set.seed(1)
  x <- pfas_and_sdq[, k]
  x <- unlist(x)
  
  print(colnames(pfas_and_sdq)[k])
  
  for (i in 2:9) {
    y <- pfas_and_sdq[, i]
    y <- unlist(y)
    
    print(colnames(pfas_and_sdq)[i])
  
  #give the general structure for the linear model
  model <- t.test(y ~ x, data = pfas_and_sdq)
  
  #extract the t-statistic
  tstat <- model$statistic
  
  #extract the  p-value 
  pval <- model$p.value
  
  #variable name
  sdq <- colnames(pfas_and_sdq)[k]
  
  #pfas name
  pfas <- colnames(pfas_and_sdq)[i]
  
  #save results of each t-test
  results_k <- as.data.frame(sdq) %>%  
    mutate(pfas = pfas) %>% 
    mutate(pval = pval) %>%  
    mutate(tstat = tstat)
  
  #for each iteration, append the results to our results dataframe
  PFAS_ttest_results <- rbind(PFAS_ttest_results, results_k)
  
}}

#Box plots for IC5 --------------------------------------------------
p1 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(emotional_problems_binary)) %>% 
  mutate(sdq = recode(emotional_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Emotional problems", subtitle = "t=1.2, p=0.27, pFDR=0.54") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p2 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(conduct_problems_binary)) %>% 
  mutate(sdq = recode(conduct_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Conduct problems", subtitle = "t=-1.6, p=0.12, pFDR=0.36") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p3 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(hyperactivity_inattention_binary)) %>% 
  mutate(sdq = recode(hyperactivity_inattention_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Hyperactivity/inattention", subtitle = "t=1.0, p=0.31, pFDR=0.46") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p4 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(peer_relationship_problems_binary)) %>% 
  mutate(sdq = recode(peer_relationship_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Peer relationship problems", subtitle = "t=-2.4, p=0.017, pFDR=0.10") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p5 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(prosocial_behaviour_binary)) %>% 
  mutate(sdq = recode(prosocial_behaviour_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Abnormal prosocial behaviour", subtitle = "t=0.17, p=0.86, pFDR=0.86") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p6 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(total_difficulties_binary)) %>% 
  mutate(sdq = recode(total_difficulties_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC5), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2980B9", "#F8E92E")) +
  labs(title = "Total difficulties", subtitle = "t=-0.74, p=0.46, pFDR=0.55") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
ggsave("plots/vectors/IC5_SDQ.svg", units = "px", width = 3000, height = 2000)

#Box plots for IC9 --------------------------------------------------
p1 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(emotional_problems_binary)) %>% 
  mutate(sdq = recode(emotional_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Emotional problems", subtitle = "t=0.03, p=0.98, pFDR=0.98") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p2 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(conduct_problems_binary)) %>% 
  mutate(sdq = recode(conduct_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Conduct problems", subtitle = "t=-0.04, p=0.96, pFDR=0.99") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p3 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(hyperactivity_inattention_binary)) %>% 
  mutate(sdq = recode(hyperactivity_inattention_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Hyperactivity/inattention", subtitle = "t=-0.94, p=0.36, pFDR=0.99") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p4 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(peer_relationship_problems_binary)) %>% 
  mutate(sdq = recode(peer_relationship_problems_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Peer relationship problems", subtitle = "t=-1.6, p=0.11, pFDR=0.68") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p5 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(prosocial_behaviour_binary)) %>% 
  mutate(sdq = recode(prosocial_behaviour_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Abnormal prosocial behaviour", subtitle = "t=0.88, p=0.38, pFDR=0.76") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

p6 <- flica_and_sdq %>% 
  mutate(sdq = as.factor(total_difficulties_binary)) %>% 
  mutate(sdq = recode(total_difficulties_binary, "0" = "No", "1" = "Yes")) %>% 
  filter(!is.na(sdq)) %>% 
  ggplot(., aes(x = sdq, y = IC9), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#2d76f5", "#e8dc27")) +
  labs(title = "Total difficulties", subtitle = "t=-0.14, p=0.88, pFDR=0.99") +
  geom_boxplot(aes(fill = sdq)) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        plot.title = element_text(hjust = 0.5, size = 16), 
        plot.subtitle = element_text(hjust = 0.5, size = 14))

ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
ggsave("plots/vectors/IC9_SDQ.svg", units = "px", width = 3000, height = 2000)
