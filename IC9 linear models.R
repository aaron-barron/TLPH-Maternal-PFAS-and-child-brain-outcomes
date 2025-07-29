#Create dataset
source("Load and clean data.R")

#Load additional libraries
library("sensemakr")
library("lm.beta")
library("ggtext")
library("poolr")

#Calculate meff for adjusting p-values -----------------------------------------

cor_mat_IC9 <- cor(PFAS_FLICA_covar_data[,c(4, 6, 7)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_IC9) <- 1

meff(R = cor_mat_IC9, method = "galwey") # = 2, multiply p-values by 2

#PFOA-Br -----------------------------------------------------------------------

#Plot
p8 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC9)) %>% 
  ggplot(data = ., aes(x = PFOA_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#09037B") +
  labs(title = "IC9 ~ Maternal PFOA-Br", 
       subtitle = "R^2^=0.14, β=0.36 (0.10-0.62), adjusted p=0.016", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOA_Br_vs_IC9_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex, data = .)
summary(PFOA_Br_vs_IC9_min_lm) 
lm.beta(PFOA_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_min_lm) 

#Fully-adjusted model
PFOA_Br_vs_IC9_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFOA_Br + `Child Age` + Sex + `Maternal BMI`, data = .)
summary(PFOA_Br_vs_IC9_full_lm)
lm.beta(PFOA_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFOA_Br_vs_IC9_full_lm) 

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOA_Br)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFOA_Br_vs_IC9_full_lm)[2,]
confint(PFOA_Br_vs_IC9_full_lm)[2,]*(sx/sy)

#PFHxS_Br ----------------------------------------------------------------------

#Plot
p9 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC9)) %>% 
  ggplot(data = ., aes(x = PFHxS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "Maroon") +
  labs(title = "IC9 ~ Maternal PFHxS-Br", 
       subtitle = "R^2^=0.12, β=-0.32 (-0.58- -0.06), adjusted p=0.036", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFHxS_Br_vs_IC9_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex, data = .)
summary(PFHxS_Br_vs_IC9_min_lm) 
lm.beta(PFHxS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFHxS_Br_vs_IC9_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFHxS_Br + `Child Age` + Sex + Gravidity, data = .)
summary(PFHxS_Br_vs_IC9_full_lm) 
lm.beta(PFHxS_Br_vs_IC9_full_lm) %>% coef() 
sensemakr::partial_r2(PFHxS_Br_vs_IC9_full_lm)

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFHxS_Br)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFHxS_Br_vs_IC9_full_lm)[2,]
confint(PFHxS_Br_vs_IC9_full_lm)[2,]*(sx/sy)


#PFOS_Br -----------------------------------------------------------------------

#Plot
p10 <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(IC9)) %>% 
  ggplot(data = ., aes(x = PFOS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", size = 2, colour = "#DBDDF3") +
  labs(title = "IC9 ~ Maternal PFOS-Br", 
       subtitle = "R^2^<0.001, β=0.020 (-0.26-0.30), adjusted p>0.99", 
       x = "PFOS-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90")) 

#Minimally-adjusted model
PFOS_Br_vs_IC9_min_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex, data = .)
summary(PFOS_Br_vs_IC9_min_lm) 
lm.beta(PFOS_Br_vs_IC9_min_lm) %>% coef() 
sensemakr::partial_r2(PFOS_Br_vs_IC9_min_lm)

#Fully-adjusted model
PFOS_Br_vs_IC9_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~ PFOS_Br + `Child Age` + Sex + `Maternal BMI` + `Maternal Alcohol` + Gravidity, data = .)
summary(PFOS_Br_vs_IC9_full_lm)
lm.beta(PFOS_Br_vs_IC9_full_lm) %>% coef()
sensemakr::partial_r2(PFOS_Br_vs_IC9_full_lm)

#standardized beta 95% CI
sx <- sd(PFAS_FLICA_covar_data$PFOS_Br)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFOS_Br_vs_IC9_full_lm)[2,]
confint(PFOS_Br_vs_IC9_full_lm)[2,]*(sx/sy)
#save plots -----------------------

ggarrange(p9, p10, p8, nrow = 1, ncol = 3)
ggsave("plots/vectors/IC9_PFAS_scatterplots.svg", units = "px", width = 4000, height = 1250)
