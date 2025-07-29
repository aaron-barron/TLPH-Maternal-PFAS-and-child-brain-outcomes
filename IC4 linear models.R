#Create dataset
source("Load and clean data.R")

#Load additional libraries
library("sensemakr")
library("lm.beta")
library("ggtext")
library("poolr")

#Calculate meff for adjusting p-values -----------------------------------------

cor_mat_IC4 <- cor(PFAS_FLICA_covar_data[,c(3, 6)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_IC4) <- 1

meff(R = cor_mat_IC4, method = "galwey") # = 1, multiply p-values by 1

#PFOA-Br -----------------------------------------------------------------------

#Plot
p6 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC4)) %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DB9D95") +
  labs(title = "IC4 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.016, β=-0.12 (-0.37-0.14), adjusted p=0.38", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC4_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC4_min_lm) 
lm.beta(PFOA_Br_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC4_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC4 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC4_full_lm) 
lm.beta(PFOA_Br_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC4_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOA_Br)
sy <- sd(PFAS_FLICA_covar_data$IC4)
confint(PFOA_Br_vs_IC4_full_lm)[2,]
confint(PFOA_Br_vs_IC4_full_lm)[2,]*(sx/sy)

#PFHpS -------------------------------------------------------------------------

#Plot
p7 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC4)) %>% 
  ggplot(data = ., aes(x = PFHpS, y = IC4)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC4 ~ Maternal PFHpS", 
       subtitle = "R^2^=0.042, β=0.18 (-0.06-0.35), adjusted p=0.16", 
       x = "PFHpS") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHpS_vs_IC4_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex, data = .)
summary(PFHpS_vs_IC4_min_lm) 
lm.beta(PFHpS_vs_IC4_min_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_min_lm) 

#Fully-adjusted model
PFHpS_vs_IC4_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC4 ~ PFHpS + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFHpS_vs_IC4_full_lm) #p = 0.174
lm.beta(PFHpS_vs_IC4_full_lm) %>% coef() 
sensemakr::partial_r2(PFHpS_vs_IC4_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFHpS)
sy <- sd(PFAS_FLICA_covar_data$IC4)
confint(PFHpS_vs_IC4_full_lm)[2,]
confint(PFHpS_vs_IC4_full_lm)[2,]*(sx/sy)

#Save plots ----------------------------

#IC4
ggarrange(p6, p7, nrow = 1, ncol = 2)
ggsave("plots/vectors/IC4_PFAS_scatterplots.svg", units = "px", width = 2650, height = 1250)
