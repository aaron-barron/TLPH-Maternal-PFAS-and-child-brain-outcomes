source("clean scripts/Load and clean data with imputations.R")
library("ggpubr")
library("svglite")
#Estimate outliers BEFORE transforming and autoscaling the data.

#Check outliers in final dataset for each PFAS ----------------------------------------
PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% filter(!is.na(IC1))

quartiles <- data.frame()
for (k in 3:10){
  set.seed(1)
  y <- PFAS_FLICA_covar_data[, k]
  y <- unlist(y)
  print(colnames(PFAS_FLICA_covar_data)[k])
  quartile_k <- quantile(y, probs = c(0, 0.25, 0.5, 0.75, 1))
  quartile_k_df <- data.frame(quartile_k) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(PFAS = colnames(PFAS_FLICA_covar_data)[k])
  
  quartiles <- rbind(quartiles, quartile_k_df)
}

quartiles <- quartiles %>% 
  select(PFAS, everything()) %>% 
  mutate(IQR = `75%` - `25%`) %>% 
  mutate(lower = `25%` - (1.5*IQR)) %>% 
  mutate(upper = `75%` + (1.5*IQR)) %>%
  mutate(extreme_lower = `25%` - (3*IQR)) %>% 
  mutate(extreme_upper = `75%` + (3*IQR))


#PFHpS n=1, extreme=0
#PFHxS_Br n=2, extreme=2
#PFNA n=2, extreme=0
#PFOA_Br n=4, extreme=3
#PFOA_L n=2, extreme=2
#PFOS_Br n=2, extreme=0
#PFOS_L n=5, extreme=2
#WQS n=0


#Outlier removal ----------------------------------------

outlier_dataset_IQR3 <- PFAS_FLICA_covar_data %>% 
  mutate(PFHpS = if_else(PFHpS < quartiles$extreme_upper[1], PFHpS, NA)) %>% 
  mutate(PFHxS_Br = if_else(PFHxS_Br < quartiles$extreme_upper[2], PFHxS_Br, NA)) %>%
  mutate(PFNA = if_else(PFNA < quartiles$extreme_upper[3], PFNA, NA)) %>%
  mutate(PFOA_Br = if_else(PFOA_Br < quartiles$extreme_upper[4], PFOA_Br, NA)) %>%
  mutate(PFOA_L = if_else(PFOA_L < quartiles$extreme_upper[5], PFOA_L, NA)) %>%
  mutate(PFOS_Br = if_else(PFOS_Br < quartiles$extreme_upper[6], PFOS_Br, NA)) %>%
  mutate(PFOS_L = if_else(PFOS_L < quartiles$extreme_upper[7], PFOS_L, NA)) %>%
  mutate(`Total (WQS) PFAS` = if_else(`Total (WQS) PFAS` < quartiles$extreme_upper[8], `Total (WQS) PFAS`, NA)) %>%
  mutate_at(3:19, log2) %>%  #Log2 transform all analytes
  mutate_at(c(3:19, 21:23, 30:31), autoscale_function) #autoscale all analytes and numeric covariates

outlier_dataset_IQR1.5 <- PFAS_FLICA_covar_data %>% 
  mutate(PFHpS = if_else(PFHpS < quartiles$upper[1], PFHpS, NA)) %>% 
  mutate(PFHxS_Br = if_else(PFHxS_Br < quartiles$upper[2], PFHxS_Br, NA)) %>%
  mutate(PFNA = if_else(PFNA < quartiles$upper[3], PFNA, NA)) %>%
  mutate(PFOA_Br = if_else(PFOA_Br < quartiles$upper[4], PFOA_Br, NA)) %>%
  mutate(PFOA_L = if_else(PFOA_L < quartiles$upper[5], PFOA_L, NA)) %>%
  mutate(PFOS_Br = if_else(PFOS_Br < quartiles$upper[6], PFOS_Br, NA)) %>%
  mutate(PFOS_L = if_else(PFOS_L < quartiles$upper[7], PFOS_L, NA)) %>%
  mutate(`Total (WQS) PFAS` = if_else(`Total (WQS) PFAS` < quartiles$upper[8], `Total (WQS) PFAS`, NA)) %>%
  mutate_at(3:19, log2) %>%  #Log2 transform all analytes
  mutate_at(c(3:19, 21:23, 30:31), autoscale_function) #autoscale all analytes and numeric covariates

outlier_dataset_visual <- PFAS_FLICA_covar_data %>% 
  mutate(PFHpS = if_else(PFHpS < quartiles$upper[1], PFHpS, NA)) %>% 
  mutate(PFHxS_Br = if_else(PFHxS_Br < quartiles$upper[2], PFHxS_Br, NA)) %>%
  mutate(PFNA = if_else(PFNA < quartiles$upper[3], PFNA, NA)) %>%
  mutate(PFOA_Br = if_else(PFOA_Br < quartiles$upper[4], PFOA_Br, NA)) %>%
  mutate(PFOA_L = if_else(PFOA_L < quartiles$upper[5], PFOA_L, NA)) %>%
  mutate(PFOS_Br = if_else(PFOS_Br < quartiles$upper[6], PFOS_Br, NA)) %>%
  mutate(PFOS_L = if_else(PFOS_L < quartiles$upper[7], PFOS_L, NA)) %>%
  mutate(`Total (WQS) PFAS` = if_else(`Total (WQS) PFAS` < quartiles$upper[8], `Total (WQS) PFAS`, NA)) %>%
  mutate(PFHxS_Br = if_else(PFHxS_Br > 0.01, PFHxS_Br, NA)) %>%
  mutate(PFNA = if_else(PFNA > 0.005, PFNA, NA)) %>%
  mutate(`Total (WQS) PFAS` = if_else(`Total (WQS) PFAS` > 0.025, `Total (WQS) PFAS`, NA)) %>%
  mutate_at(3:19, log2) %>%  #Log2 transform all analytes
  mutate_at(c(3:19, 21:23, 30:31), autoscale_function) #autoscale all analytes and numeric covariates


#IC5 IQR3---------------------------------------------------------------------------
#Total WQS PFAS

#Plot
p1 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = `Total (WQS) PFAS`, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Total (WQS) Maternal PFAS", 
       subtitle = "R^2^= 0.034, β = 0.18, p = 0.42 after full adjustment", 
       x = "Total (WQS) PFAS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
Total_PFAS_vs_IC5_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex, data = .)
summary(Total_PFAS_vs_IC5_min_lm) #p = 0.418
lm.beta(Total_PFAS_vs_IC5_min_lm) %>% coef() #B = 0.18
sensemakr::partial_r2(Total_PFAS_vs_IC5_min_lm) #par R2 = 0.03

#Fully-adjusted model
Total_PFAS_vs_IC5_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex + `Maternal Age`
     , data = .)
summary(Total_PFAS_vs_IC5_full_lm) 
lm.beta(Total_PFAS_vs_IC5_full_lm) %>% coef() 
sensemakr::partial_r2(Total_PFAS_vs_IC5_full_lm) 

#PFNA

#Plot
p2 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFNA, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFNA", 
       subtitle = "R^2^= 0.13, β = 0.39, p = 0.024 after full adjustment") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFNA_vs_IC5_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_IC5_min_lm)
lm.beta(PFNA_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_min_lm)

#Fully-adjusted model
PFNA_vs_IC5_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_IC5_full_lm)
lm.beta(PFNA_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_full_lm)

#PFOA_L

#Plot
p3 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOA_L, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.11, β = 0.33, p = 0.042 after full adjustment", x = "PFOA-L") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_L_vs_IC5_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_IC5_min_lm)
lm.beta(PFOA_L_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_min_lm)

#Fully-adjusted model
PFOA_L_vs_IC5_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_IC5_full_lm)
lm.beta(PFOA_L_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_full_lm)


#IC5 IQR1.5---------------------------------------------------------------------------
#Total WQS PFAS

#Plot
p1 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = `Total (WQS) PFAS`, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Total (WQS) Maternal PFAS", 
       subtitle = "R^2^= 0.034, β = 0.18, p = 0.42 after full adjustment", 
       x = "Total (WQS) PFAS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
Total_PFAS_vs_IC5_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex, data = .)
summary(Total_PFAS_vs_IC5_min_lm) #p = 0.418
lm.beta(Total_PFAS_vs_IC5_min_lm) %>% coef() #B = 0.18
sensemakr::partial_r2(Total_PFAS_vs_IC5_min_lm) #par R2 = 0.03

#Fully-adjusted model
Total_PFAS_vs_IC5_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex + `Maternal Age`
     , data = .)
summary(Total_PFAS_vs_IC5_full_lm) 
lm.beta(Total_PFAS_vs_IC5_full_lm) %>% coef() 
sensemakr::partial_r2(Total_PFAS_vs_IC5_full_lm) 

#PFNA

#Plot
p2 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFNA, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFNA", 
       subtitle = "R^2^= 0.13, β = 0.35, p = 0.032 after full adjustment") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFNA_vs_IC5_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_IC5_min_lm)
lm.beta(PFNA_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_min_lm)

#Fully-adjusted model
PFNA_vs_IC5_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_IC5_full_lm)
lm.beta(PFNA_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_full_lm)

#PFOA_L

#Plot
p3 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOA_L, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.11, β = 0.33, p = 0.042 after full adjustment", x = "PFOA-L") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_L_vs_IC5_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_IC5_min_lm)
lm.beta(PFOA_L_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_min_lm)

#Fully-adjusted model
PFOA_L_vs_IC5_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_IC5_full_lm)
lm.beta(PFOA_L_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_full_lm)


#IC5 visual---------------------------------------------------------------------------
#Total WQS PFAS

#Plot
p1 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = `Total (WQS) PFAS`, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Total (WQS) Maternal PFAS", 
       subtitle = "R^2^= 0.12, β = 0.33, p = 0.04 after full adjustment", 
       x = "Total (WQS) PFAS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
Total_PFAS_vs_IC5_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex, data = .)
summary(Total_PFAS_vs_IC5_min_lm) #p = 0.418
lm.beta(Total_PFAS_vs_IC5_min_lm) %>% coef() #B = 0.18
sensemakr::partial_r2(Total_PFAS_vs_IC5_min_lm) #par R2 = 0.03

#Fully-adjusted model
Total_PFAS_vs_IC5_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex + `Maternal Age`
     , data = .)
summary(Total_PFAS_vs_IC5_full_lm) 
lm.beta(Total_PFAS_vs_IC5_full_lm) %>% coef() 
sensemakr::partial_r2(Total_PFAS_vs_IC5_full_lm) 

#PFNA

#Plot
p2 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFNA, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFNA", 
       subtitle = "R^2^= 0.08, β = 0.29, p = 0.12 after full adjustment") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFNA_vs_IC5_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_IC5_min_lm)
lm.beta(PFNA_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_min_lm)

#Fully-adjusted model
PFNA_vs_IC5_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_IC5_full_lm)
lm.beta(PFNA_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_full_lm)

#PFOA_L

#Plot
p3 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOA_L, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.11, β = 0.33, p = 0.042 after full adjustment", x = "PFOA-L") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_L_vs_IC5_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_IC5_min_lm)
lm.beta(PFOA_L_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_min_lm)

#Fully-adjusted model
PFOA_L_vs_IC5_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_IC5_full_lm)
lm.beta(PFOA_L_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_full_lm)


#CC body IQR3 --------------------------------------------------------------------------------------

#PFNA ~ cc body FA 
p4 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFNA, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFNA", 
       subtitle = "R^2^= 0.090, β = 0.32, p = 0.037 after full adjustment",
       x = "PFNA",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFNA_vs_cc_body_FA_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_cc_body_FA_min_lm)
lm.beta(PFNA_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_min_lm) 

PFNA_vs_cc_body_FA_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_cc_body_FA_full_lm) 
lm.beta(PFNA_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_full_lm) 


#PFOA-L ~ cc body FA 
p5 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOA_L, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.05, β = 0.22, p = 0.13 after full adjustment",
       x = "PFOA-L",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFOA_L_vs_cc_body_FA_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_cc_body_FA_min_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_min_lm) 

PFOA_L_vs_cc_body_FA_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_cc_body_FA_full_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_full_lm) 


#CC body IQR1.5 --------------------------------------------------------------------------------------

#PFNA ~ cc body FA 
p4 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFNA, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFNA", 
       subtitle = "R^2^= 0.092, β = 0.31, p = 0.043 after full adjustment",
       x = "PFNA",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFNA_vs_cc_body_FA_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_cc_body_FA_min_lm)
lm.beta(PFNA_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_min_lm) 

PFNA_vs_cc_body_FA_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_cc_body_FA_full_lm) 
lm.beta(PFNA_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_full_lm) 


#PFOA-L ~ cc body FA 
p5 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOA_L, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.050, β = 0.22, p = 0.13 after full adjustment",
       x = "PFOA-L",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFOA_L_vs_cc_body_FA_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_cc_body_FA_min_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_min_lm) 

PFOA_L_vs_cc_body_FA_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_cc_body_FA_full_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_full_lm) 


#CC body visual --------------------------------------------------------------------------------------

#PFNA ~ cc body FA 
p4 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFNA, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFNA", 
       subtitle = "R^2^= 0.079, β = 0.33, p = 0.064 after full adjustment",
       x = "PFNA",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFNA_vs_cc_body_FA_min_lm <- outlier_dataset_visual %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_cc_body_FA_min_lm)
lm.beta(PFNA_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_min_lm) 

PFNA_vs_cc_body_FA_full_lm <- outlier_dataset_visual %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_cc_body_FA_full_lm) 
lm.beta(PFNA_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_full_lm) 


#PFOA-L ~ cc body FA 
p5 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOA_L, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.050, β = 0.22, p = 0.13 after full adjustment",
       x = "PFOA-L",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFOA_L_vs_cc_body_FA_min_lm <- outlier_dataset_visual %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_cc_body_FA_min_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_min_lm) 

PFOA_L_vs_cc_body_FA_full_lm <- outlier_dataset_visual %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_cc_body_FA_full_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_full_lm) 


#IC4 IQR3 --------------------------------------------------------------------------------------------------
#PFOA-Br

#Plot
p6 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DB9D95") +
  labs(title = "IC4 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.018, β=-0.12, p=0.38 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC4_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC4_min_lm) 
lm.beta(PFOA_Br_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC4_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC4_full_lm) 
lm.beta(PFOA_Br_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_full_lm) 

#PFHpS

#Plot
p7 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFHpS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC4 ~ Maternal PFHpS", 
       subtitle = "R^2^=0.042, β=0.18, p=0.16 after full adjustment", 
       x = "PFHpS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHpS_vs_IC4_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex, data = .)
summary(PFHpS_vs_IC4_min_lm) 
lm.beta(PFHpS_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_min_lm) 

#Fully-adjusted model
PFHpS_vs_IC4_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFHpS_vs_IC4_full_lm) #p = 0.174
lm.beta(PFHpS_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_full_lm) 

#IC4 IQR1.5 --------------------------------------------------------------------------------------------------
#PFOA-Br

#Plot
p6 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DB9D95") +
  labs(title = "IC4 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.011, β=-0.095, p=0.50 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC4_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC4_min_lm) 
lm.beta(PFOA_Br_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC4_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC4_full_lm) 
lm.beta(PFOA_Br_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_full_lm) 

#PFHpS

#Plot
p7 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFHpS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC4 ~ Maternal PFHpS", 
       subtitle = "R^2^=0.035, β=0.18, p=0.21 after full adjustment", 
       x = "PFHpS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHpS_vs_IC4_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex, data = .)
summary(PFHpS_vs_IC4_min_lm) 
lm.beta(PFHpS_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_min_lm) 

#Fully-adjusted model
PFHpS_vs_IC4_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFHpS_vs_IC4_full_lm) #p = 0.174
lm.beta(PFHpS_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_full_lm) 


#IC4 visual --------------------------------------------------------------------------------------------------
#PFOA-Br

#Plot
outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DB9D95") +
  labs(title = "IC4 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.016, β=-0.12, p=0.39 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) +
  xlim(-1.8, 2.05) 

#Minimally-adjusted model
PFOA_Br_vs_IC4_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC4_min_lm) 
lm.beta(PFOA_Br_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC4_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC4_full_lm) 
lm.beta(PFOA_Br_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_full_lm) 

#PFHpS

#Plot
outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFHpS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC4 ~ Maternal PFHpS", 
       subtitle = "R^2^=0.042, β=0.18, p=0.16 after full adjustment", 
       x = "PFHpS") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) +
  xlim(-1.8, 2.05) 

#Minimally-adjusted model
PFHpS_vs_IC4_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex, data = .)
summary(PFHpS_vs_IC4_min_lm) 
lm.beta(PFHpS_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_min_lm) 

#Fully-adjusted model
PFHpS_vs_IC4_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFHpS_vs_IC4_full_lm) #p = 0.174
lm.beta(PFHpS_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_full_lm) 

#IC9 IQR3 ----------------------------------------------------------------------------------------------------

#PFOA-Br

#Plot
p8 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#09037B") +
  labs(title = "IC9 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.13, β=0.34, p=0.035 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC9_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC9_min_lm) 
lm.beta(PFOA_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC9_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC9_full_lm)
lm.beta(PFOA_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_full_lm) 

#PFHxS_Br 

#Plot
p9 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "IC9 ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.24, β=-0.46, p<0.001 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_IC9_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_IC9_min_lm) 
lm.beta(PFHxS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_IC9_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex + Gravidity, data = .)
summary(PFHxS_Br_vs_IC9_full_lm) 
lm.beta(PFHxS_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_full_lm)

#PFOS_Br 

#Plot
p10 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC9 ~ Maternal PFOS-Br", 
       subtitle = "R^2^<0.001, β=0.020, p>0.99 after full adjustment", 
       x = "PFOS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOS_Br_vs_IC9_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_IC9_min_lm) 
lm.beta(PFOS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFOS_Br_vs_IC9_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_IC9_full_lm)
lm.beta(PFOS_Br_vs_IC9_full_lm) %>% coef()
sensemakr::partial_r2(PFOS_Br_vs_IC9_full_lm)

#IC9 IQR1.5 ----------------------------------------------------------------------------------------------------

#PFOA-Br

#Plot
p8 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#09037B") +
  labs(title = "IC9 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.12, β=0.33, p=0.044 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC9_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC9_min_lm) 
lm.beta(PFOA_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC9_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC9_full_lm)
lm.beta(PFOA_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_full_lm) 

#PFHxS_Br 

#Plot
p9 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "IC9 ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.24, β=-0.46, p<0.001 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_IC9_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_IC9_min_lm) 
lm.beta(PFHxS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_IC9_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex + Gravidity, data = .)
summary(PFHxS_Br_vs_IC9_full_lm) 
lm.beta(PFHxS_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_full_lm)

#PFOS_Br 

#Plot
p10 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC9 ~ Maternal PFOS-Br", 
       subtitle = "R^2^<0.001, β=0.026, p>0.99 after full adjustment", 
       x = "PFOS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOS_Br_vs_IC9_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_IC9_min_lm) 
lm.beta(PFOS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFOS_Br_vs_IC9_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_IC9_full_lm)
lm.beta(PFOS_Br_vs_IC9_full_lm) %>% coef()
sensemakr::partial_r2(PFOS_Br_vs_IC9_full_lm)

#IC9 visual ----------------------------------------------------------------------------------------------------

#PFOA-Br

#Plot
p8 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#09037B") +
  labs(title = "IC9 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.11, β=0.32, p=0.044 after full adjustment", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC9_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC9_min_lm) 
lm.beta(PFOA_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC9_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC9_full_lm)
lm.beta(PFOA_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_full_lm) 

#PFHxS_Br 

#Plot
p9 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "IC9 ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.21, β=-0.44, p=0.0028 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_IC9_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_IC9_min_lm) 
lm.beta(PFHxS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_IC9_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex + Gravidity, data = .)
summary(PFHxS_Br_vs_IC9_full_lm) 
lm.beta(PFHxS_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_full_lm)

#PFOS_Br 

#Plot
p10 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC9 ~ Maternal PFOS-Br", 
       subtitle = "R^2^<0.001, β=0.020, p>0.99 after full adjustment", 
       x = "PFOS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOS_Br_vs_IC9_min_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_IC9_min_lm) 
lm.beta(PFOS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFOS_Br_vs_IC9_full_lm <- outlier_dataset_visual %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_IC9_full_lm)
lm.beta(PFOS_Br_vs_IC9_full_lm) %>% coef()
sensemakr::partial_r2(PFOS_Br_vs_IC9_full_lm)

#Hypothalamus IQR3 -------------

#PFOS_Br

#Plot
p11 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "Hypothalamus MD ~ Maternal PFOS-Br", 
       subtitle = "R^2^=0.11, β=-0.30, p=0.025 after full adjustment", 
       x = "PFOS-Br", y = "MD (mm^2^/s, autoscaled)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_markdown(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 


#Minimally-adjusted model
PFOS_Br_vs_HTH_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_HTH_min_lm) 
lm.beta(PFOS_Br_vs_HTH_min_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFOS_Br_vs_HTH_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_HTH_full_lm) 
lm.beta(PFOS_Br_vs_HTH_full_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_full_lm) 

#PFHxS_Br

#Plot
p12 <- outlier_dataset_IQR3 %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "Hypothalamus MD ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.028, β=0.16, p=0.26 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_HTH_min_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_HTH_min_lm) 
lm.beta(PFHxS_Br_vs_HTH_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_HTH_full_lm <- outlier_dataset_IQR3 %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex  + Gravidity, data = .)
summary(PFHxS_Br_vs_HTH_full_lm)
lm.beta(PFHxS_Br_vs_HTH_full_lm) %>% coef()
sensemakr::partial_r2(PFHxS_Br_vs_HTH_full_lm)
#Hypothalamus IQR1.5 -------------

#PFOS_Br

#Plot
p11 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "Hypothalamus MD ~ Maternal PFOS-Br", 
       subtitle = "R^2^=0.05, β=-0.20, p=0.14 after full adjustment", 
       x = "PFOS-Br", y = "MD (mm^2^/s, autoscaled)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_markdown(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 


#Minimally-adjusted model
PFOS_Br_vs_HTH_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_HTH_min_lm) 
lm.beta(PFOS_Br_vs_HTH_min_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFOS_Br_vs_HTH_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_HTH_full_lm) 
lm.beta(PFOS_Br_vs_HTH_full_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_full_lm) 

#PFHxS_Br

#Plot
p12 <- outlier_dataset_IQR1.5 %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "Hypothalamus MD ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.028, β=0.16, p=0.26 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_HTH_min_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_HTH_min_lm) 
lm.beta(PFHxS_Br_vs_HTH_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_HTH_full_lm <- outlier_dataset_IQR1.5 %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex  + Gravidity, data = .)
summary(PFHxS_Br_vs_HTH_full_lm)
lm.beta(PFHxS_Br_vs_HTH_full_lm) %>% coef()
sensemakr::partial_r2(PFHxS_Br_vs_HTH_full_lm)
#Hypothalamus visual -------------

#PFOS_Br

#Plot
p11 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "Hypothalamus MD ~ Maternal PFOS-Br", 
       subtitle = "R^2^=0.050, β=-0.20, p=0.14 after full adjustment", 
       x = "PFOS-Br", y = "MD (mm^2^/s, autoscaled)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_markdown(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 


#Minimally-adjusted model
PFOS_Br_vs_HTH_min_lm <- outlier_dataset_visual %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_HTH_min_lm) 
lm.beta(PFOS_Br_vs_HTH_min_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFOS_Br_vs_HTH_full_lm <- outlier_dataset_visual %>% 
  lm(formula = Hypothalamus ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_HTH_full_lm) 
lm.beta(PFOS_Br_vs_HTH_full_lm) %>% coef() 
partial_r2(PFOS_Br_vs_HTH_full_lm) 

#PFHxS_Br

#Plot
p12 <- outlier_dataset_visual %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = Hypothalamus)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "Hypothalamus MD ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.021, β=0.14, p=0.34 after full adjustment", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_HTH_min_lm <- outlier_dataset_visual %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_HTH_min_lm) 
lm.beta(PFHxS_Br_vs_HTH_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_HTH_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_HTH_full_lm <- outlier_dataset_visual %>% 
  lm(formula = Hypothalamus ~ PFHxS_Br + `Child Age` + Sex  + Gravidity, data = .)
summary(PFHxS_Br_vs_HTH_full_lm)
lm.beta(PFHxS_Br_vs_HTH_full_lm) %>% coef()
sensemakr::partial_r2(PFHxS_Br_vs_HTH_full_lm)

#Save figures --------------------------

#IC5
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
ggsave("plots/vectors/outlier_visual_IC5_lm.svg", units = "px", width = 4000, height = 1250)

#CC body FA
ggarrange(p4, p5, nrow = 1, ncol = 2)
ggsave("plots/vectors/outlier_visual_cc_body_FA_lm.svg", units = "px", width = 2650, height = 1250)

#IC4
ggarrange(p6, p7, nrow = 1, ncol = 2)
ggsave("plots/vectors/outlier_visual_IC4_lm.svg", units = "px", width = 2650, height = 1250)

#IC9
ggarrange(p9, p10, p8, nrow = 1, ncol = 3)
ggsave("plots/vectors/outlier_visual_IC9_lm.svg", units = "px", width = 4000, height = 1250)

#Hypothalamus
ggarrange(p11, p12, nrow = 1, ncol = 2)
ggsave("plots/vectors/outlier_visual_HTH_lm.svg", units = "px", width = 2650, height = 1250)
