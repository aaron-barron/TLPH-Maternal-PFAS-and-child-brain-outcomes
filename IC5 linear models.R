#Create dataset
source("Load and clean data.R")

#Load additional libraries
library("sensemakr")
library("lm.beta")
library("ggtext")
library("poolr")

#Calculate meff for adjusting p-values -----------------------------------------

cor_mat_IC5 <- cor(PFAS_FLICA_covar_data[,c(5, 7, 10)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_IC5) <- 1

meff(R = cor_mat_IC5, method = "galwey") # = 2, multiply p-values by 2

PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% filter(!is.na(IC1))

#Total WQS PFAS ----------------------------------------------------------------

#Plot
p1 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = `Total (WQS) PFAS`, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Total (WQS) Maternal PFAS", 
       subtitle = "R^2^= 0.034, β=0.18 (-0.10-0.46), adjusted p=0.42", 
       x = "Total (WQS) PFAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
Total_PFAS_vs_IC5_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex, data = .)
summary(Total_PFAS_vs_IC5_min_lm) #p = 0.418
lm.beta(Total_PFAS_vs_IC5_min_lm) %>% coef() #B = 0.18
sensemakr::partial_r2(Total_PFAS_vs_IC5_min_lm) #par R2 = 0.03

#Fully-adjusted model
Total_PFAS_vs_IC5_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ `Total (WQS) PFAS` + `Child Age` + Sex + `Maternal Age`
     , data = .)
summary(Total_PFAS_vs_IC5_full_lm) 
lm.beta(Total_PFAS_vs_IC5_full_lm) %>% coef() 
sensemakr::partial_r2(Total_PFAS_vs_IC5_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$`Total (WQS) PFAS`)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(Total_PFAS_vs_IC5_full_lm)[2,]
confint(Total_PFAS_vs_IC5_full_lm)[2,]*(sx/sy)

#PFNA --------------------------------------------------------------------------

#Plot
p2 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = PFNA, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFNA", 
       subtitle = "R^2^= 0.13, β=0.39 (0.09-0.69), adjusted p = 0.024") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFNA_vs_IC5_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_IC5_min_lm)
lm.beta(PFNA_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_min_lm)

#Fully-adjusted model
PFNA_vs_IC5_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_IC5_full_lm)
lm.beta(PFNA_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFNA_vs_IC5_full_lm)

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFNA)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(PFNA_vs_IC5_full_lm)[2,]
confint(PFNA_vs_IC5_full_lm)[2,]*(sx/sy)

#PFOA_L-------------------------------------------------------------------------

#Plot
p3 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC5)) %>% 
  ggplot(data = ., aes(x = PFOA_L, y = IC5)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "IC5 ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.13, β=0.36 (0.09-0.63), adjusted p = 0.022", x = "PFOA-L") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_L_vs_IC5_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_IC5_min_lm)
lm.beta(PFOA_L_vs_IC5_min_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_min_lm)

#Fully-adjusted model
PFOA_L_vs_IC5_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_IC5_full_lm)
lm.beta(PFOA_L_vs_IC5_full_lm) %>% coef()
sensemakr::partial_r2(PFOA_L_vs_IC5_full_lm)

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOA_L)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(PFOA_L_vs_IC5_full_lm)[2,]
confint(PFOA_L_vs_IC5_full_lm)[2,]*(sx/sy)

#Post hoc hypothesis tests: corpus callosum FA ~ PFNA and PFOA-L ---------------

#Calculate meff for adjusting p-values
cor_mat_BCC <- cor(PFAS_FLICA_covar_data[,c(5, 7)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_BCC) <- 1

meff(R = cor_mat_BCC, method = "galwey") # = 1, multiply p-values by 1

#PFOA-L ~ cc body FA 
p5 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = PFOA_L, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFOA-L", 
       subtitle = "R^2^= 0.072, β=0.27 (-0.01-0-58), adjusted p=0.065",
       x = "PFOA-L",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFOA_L_vs_cc_body_FA_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex, data = .)
summary(PFOA_L_vs_cc_body_FA_min_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_min_lm) 

PFOA_L_vs_cc_body_FA_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = BCC_FA ~ PFOA_L + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_L_vs_cc_body_FA_full_lm) 
lm.beta(PFOA_L_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_L_vs_cc_body_FA_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOA_L)
sy <- sd(PFAS_FLICA_covar_data$BCC_FA)
confint(PFOA_L_vs_cc_body_FA_full_lm)[2,]
confint(PFOA_L_vs_cc_body_FA_full_lm)[2,]*(sx/sy)

#PFNA ~ cc body FA 
p4 <- PFAS_FLICA_covar_data %>% 
  ggplot(data = ., aes(x = PFNA, y = BCC_FA)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Dark Blue") +
  labs(title = "Corpus callosum body FA ~ Maternal PFNA", 
       subtitle = "R^2^= 0.090, β=0.32 (0.016-0.64), adjusted p=0.037",
       x = "PFNA",
       y = "FA") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

PFNA_vs_cc_body_FA_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex, data = .)
summary(PFNA_vs_cc_body_FA_min_lm)
lm.beta(PFNA_vs_cc_body_FA_min_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_min_lm) 

PFNA_vs_cc_body_FA_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = BCC_FA ~ PFNA + `Child Age` + Sex + `Maternal Age` + `Maternal Smoking`, data = .)
summary(PFNA_vs_cc_body_FA_full_lm) 
lm.beta(PFNA_vs_cc_body_FA_full_lm) %>% coef() 
sensemakr::partial_r2(PFNA_vs_cc_body_FA_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFNA)
sy <- sd(PFAS_FLICA_covar_data$BCC_FA)
confint(PFNA_vs_cc_body_FA_full_lm)[2,]
confint(PFNA_vs_cc_body_FA_full_lm)[2,]*(sx/sy)

#Save plots --------------------

library("ggpubr")
ggarrange(p1, p2, p3, nrow = 1, ncol = 3)
ggsave("plots/vectors/IC5_PFAS_scatterplots.svg", units = "px", width = 4000, height = 1250)

ggarrange(p4, p5, nrow = 1, ncol = 2)
ggsave("plots/vectors/CC_body_FA_PFAS_scatterplots.svg", units = "px", width = 2650, height = 1250)
