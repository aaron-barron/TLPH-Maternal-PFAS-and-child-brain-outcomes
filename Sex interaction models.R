#Create dataset
source("Load and clean data.R")

#Load additional libraries
library("sensemakr")
library("lm.beta")
library("ggtext")
library("poolr")

PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% filter(!is.na(IC5))

#Calculate meff for adjusting p-values -----------------------------------------

#For IC5
cor_mat_IC5 <- cor(PFAS_FLICA_covar_data[,c(5, 7, 10)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_IC5) <- 1

meff(R = cor_mat_IC5, method = "galwey") # = 2, multiply p-values by 2

#For IC9
cor_mat_IC9 <- cor(PFAS_FLICA_covar_data[,c(4, 6, 7)], use='pairwise.complete.obs', method = 'spearman') %>%  
  as.matrix() 

diag(cor_mat_IC9) <- 1

meff(R = cor_mat_IC9, method = "galwey") # = 2, multiply p-values by 2

#Sex-specific analyses for IC5 -------------------------------------------------

#Just compare IC5 value by sex
PFAS_FLICA_covar_data %>% 
  mutate(Sex = as.factor(Sex)) %>% 
  mutate(Sex = recode(Sex, "1" = "Male", "2" = "Female")) %>% 
  ggplot(., aes(x=Sex, y=IC5, fill=Sex)) +
  scale_fill_manual(values = c("#C4193B", "#1998C4")) +
  geom_boxplot(linewidth = 1.25, alpha = 0.8) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5), axis.title.y = element_text(size = 14)) +
  labs(title = "IC5 ~ Sex", subtitle = "t = 1.74, p = 0.089") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

t.test(IC5 ~ as.factor(Sex), PFAS_FLICA_covar_data)

ggsave("plots/vectors/IC5_sex.svg", units = "px", width = 1500, height = 1500)

#Interaction: Total (WQS) PFAS and Sex
sex1 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = `Total (WQS) PFAS`, y = IC5)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  
  labs(title = "IC5 ~ Maternal Total (WQS) PFAS x Sex", 
       subtitle = "R^2^=0.020, β=-0.42 (-1.3-0.49), adjusted p=0.70", 
       x = "Total (WQS) PFAS") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        legend.position = "none")

WQS_vs_IC5_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~  `Child Age`  + `Maternal Age` + `Total (WQS) PFAS`*Sex, data = .)
summary(WQS_vs_IC5_sex_full_lm) 
lm.beta(WQS_vs_IC5_sex_full_lm) %>% coef() 
partial_r2(WQS_vs_IC5_sex_full_lm) 

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$`Total (WQS) PFAS`*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(WQS_vs_IC5_sex_full_lm)[6,]
confint(WQS_vs_IC5_sex_full_lm)[6,]*(sx/sy)

#Interaction: PFNA and Sex
sex2 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = PFNA, y = IC5)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  labs(title = "IC5 ~ Maternal PFNA x Sex", 
       subtitle = "R^2^ < 0.01, β=0.18 (-0.64-1.02), adjusted p>0.99", 
       x = "PFNA") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        legend.position = "none") 


PFNA_vs_IC5_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~  `Child Age`  + `Maternal Age` + `Maternal Smoking` + PFNA*Sex, data = .)
summary(PFNA_vs_IC5_sex_full_lm) 
lm.beta(PFNA_vs_IC5_sex_full_lm) %>% coef() 
partial_r2(PFNA_vs_IC5_sex_full_lm) 

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$PFNA*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(PFNA_vs_IC5_sex_full_lm)[7,]
confint(PFNA_vs_IC5_sex_full_lm)[7,]*(sx/sy)

#Interaction: PFOA-L and Sex
sex3 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = PFOA_L, y = IC5)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  labs(title = "IC5 ~ Maternal PFOA-L x Sex", 
       subtitle = "R^2^=0.02, β=-0.40 (-0.56-0.89), adjusted p=0.68", 
       x = "PFOA-L") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        legend.position = "none") 

PFOA_L_vs_IC5_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC5 ~  `Child Age`  + `Maternal BMI`  + PFOA_L*Sex, data = .)
summary(PFOA_L_vs_IC5_sex_full_lm)
lm.beta(PFOA_L_vs_IC5_sex_full_lm) %>% coef() 
partial_r2(PFOA_L_vs_IC5_sex_full_lm)

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$PFOA_L*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC5)
confint(PFNA_vs_IC5_sex_full_lm)[7,]
confint(PFNA_vs_IC5_sex_full_lm)[7,]*(sx/sy)

ggarrange(sex1, sex2, sex3, ncol = 3, nrow = 1)
ggsave("plots/vectors/IC5_sex_interaction_scatterplots.svg", units = "px", width = 4000, height = 1250)

#Sex-specific analyses for IC9 ------------------------------------------------------

#Just compare IC9 value by sex
PFAS_FLICA_covar_data %>% 
  mutate(Sex = as.factor(Sex)) %>% 
  mutate(Sex = recode(Sex, "1" = "Male", "2" = "Female")) %>% 
  ggplot(., aes(x=Sex, y=IC9, fill=Sex)) +
  scale_fill_manual(values = c("#C4193B", "#1998C4")) +
  geom_boxplot(linewidth = 1.25, alpha = 0.8) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 15), plot.subtitle = element_text(hjust = 0.5), axis.title.y = element_text(size = 14)) +
  labs(title = "IC9 ~ Sex", subtitle = "t = 2, p = 0.003") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

t.test(IC9 ~ as.factor(Sex), PFAS_FLICA_covar_data)

ggsave("plots/vectors/IC9_sex.svg", units = "px", width = 1500, height = 1500)

#Interaction: PFOA-Br and Sex
sex4 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = PFOA_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  labs(title = "IC9 ~ Maternal PFOA-Br x Sex", 
       subtitle = "R^2^<0.01, β=-0.05 (-0.82-0.71), adjusted p>0.99", 
       x = "PFOA-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        legend.position = "none")

PFOA_Br_vs_IC9_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~  `Child Age`  + `Maternal BMI` + PFOA_Br*Sex, data = .)
summary(PFOA_Br_vs_IC9_sex_full_lm) 
lm.beta(PFOA_Br_vs_IC9_sex_full_lm) %>% coef() 
partial_r2(PFOA_Br_vs_IC9_sex_full_lm) 

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$PFOA_Br*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFOA_Br_vs_IC9_sex_full_lm)[6,]
confint(PFOA_Br_vs_IC9_sex_full_lm)[6,]*(sx/sy)

#Interaction: PFOS-Br and Sex
sex5 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = PFOS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  labs(title = "IC9 ~ Maternal PFOS-Br x Sex", 
       subtitle = "R^2^<0.01, β=-0.19 (-1.04-0.67), p>0.99", 
       x = "PFOS-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        legend.position = "none") 

PFOS_Br_vs_IC9_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~  `Child Age`  + `Maternal BMI` + `Maternal Alcohol` +Gravidity +PFOS_Br*Sex, data = .)
summary(PFOS_Br_vs_IC9_sex_full_lm) 
lm.beta(PFOS_Br_vs_IC9_sex_full_lm) %>% coef() 
partial_r2(PFOS_Br_vs_IC9_sex_full_lm) 

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$PFOS_Br*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFOS_Br_vs_IC9_sex_full_lm)[8,]
confint(PFOS_Br_vs_IC9_sex_full_lm)[8,]*(sx/sy)

#Interaction: PFHxS-Br and Sex
sex6 <- PFAS_FLICA_covar_data %>% 
  mutate(Sex = recode(Sex,  "1" = "Male", "2" = "Female")) %>%
  ggplot(data = ., aes(x = PFHxS_Br, y = IC9)) +
  geom_point(size = 2, alpha = 0.75, aes(colour = Sex)) +
  geom_smooth(method = "lm", size = 2, aes(colour = Sex)) +
  
  scale_colour_manual(values = c("#1998C4", "#C4193B")) +
  labs(title = "IC9 ~ Maternal PFHxS-Br x Sex", 
       subtitle = "R^2^<0.01, β=-0.28 (-1.14-0.57), adjusted p>0.99", 
       x = "PFHxS-Br") +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_markdown(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"), 
        legend.position = "none") 

PFHxS_Br_vs_IC9_sex_full_lm <- PFAS_FLICA_covar_data %>% 
  lm(formula = IC9 ~  `Child Age`  + Gravidity + PFHxS_Br*Sex, data = .)
summary(PFHxS_Br_vs_IC9_sex_full_lm) 
lm.beta(PFHxS_Br_vs_IC9_sex_full_lm) %>% coef() 
partial_r2(PFHxS_Br_vs_IC9_sex_full_lm)

#standardized beta 95% CI
PFAS_FLICA_covar_data$x <- PFAS_FLICA_covar_data$PFHxS_Br*PFAS_FLICA_covar_data$Sex
sx <- sd(PFAS_FLICA_covar_data$x)
sy <- sd(PFAS_FLICA_covar_data$IC9)
confint(PFHxS_Br_vs_IC9_sex_full_lm)[6,]
confint(PFHxS_Br_vs_IC9_sex_full_lm)[6,]*(sx/sy)

ggarrange(sex6, sex4, sex5,  ncol = 3, nrow = 1)
ggsave("plots/vectors/IC9_sex_interaction_scatterplots.svg", units = "px", width = 4000, height = 1250)
