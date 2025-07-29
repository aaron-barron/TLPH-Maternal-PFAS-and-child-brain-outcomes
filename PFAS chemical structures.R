#Create dataset
library(ggtext)
#In practice, to save time, instead of re-running these scripts each time, the 
#user could save the results tables and read them back into this script.

source("Load and clean data.R")
source("Correlation network analysis.R")
source("Elastic net regression.R")

#Elastic net coefficients ------------------------------------------------------

#Create datasets
all_regularized_coefs <- rbind((IC1_coefficients_plot %>% mutate(IC = "IC1")),
                               (IC2_coefficients_plot %>% mutate(IC = "IC2")),
                               (IC3_coefficients_plot %>% mutate(IC = "IC3")),
                               (IC4_coefficients_plot %>% mutate(IC = "IC4")),
                               (IC5_coefficients_plot %>% mutate(IC = "IC5")),
                               (IC6_coefficients_plot %>% mutate(IC = "IC6")),
                               (IC7_coefficients_plot %>% mutate(IC = "IC7")),
                               (IC8_coefficients_plot %>% mutate(IC = "IC8")),
                               (IC9_coefficients_plot %>% mutate(IC = "IC9")),
                               (IC10_coefficients_plot %>% mutate(IC = "I10")))

all_regularized_coefs_structure <- all_regularized_coefs %>% 
  mutate(coef_magnitude = if_else(Coefficient < 0, Coefficient*(-1), Coefficient)) %>%
  mutate(carbons = case_when(PFAS == "PFHpS" ~ 7, PFAS == "PFHxS_Br" ~ 6, PFAS == "PFNA" ~ 9, PFAS == "PFOA_Br" ~ 8,
                             PFAS == "PFOA_L" ~ 8, PFAS == "PFOS_Br" ~ 8, PFAS == "PFOS_L" ~ 8, TRUE ~ NA)) %>% 
  
  mutate(structure = case_when(PFAS == "PFHpS" ~ "Linear", PFAS == "PFHxS_Br" ~ "Branched", PFAS == "PFNA" ~ "Linear", PFAS == "PFOA_Br" ~ "Branched",
                               PFAS == "PFOA_L" ~ "Linear", PFAS == "PFOS_Br" ~ "Branched", PFAS == "PFOS_L" ~ "Linear", TRUE ~ NA)) %>% 
  
  mutate(functional_group = case_when(PFAS == "PFHpS" ~ "Sulfonic Acid", PFAS == "PFHxS_Br" ~ "Sulfonic Acid", PFAS == "PFNA" ~ "Carboxylic Acid", PFAS == "PFOA_Br" ~ "Carboxylic Acid",
                                      PFAS == "PFOA_L" ~ "Carboxylic Acid", PFAS == "PFOS_Br" ~ "Sulfonic Acid", PFAS == "PFOS_L" ~ "Sulfonic Acid", TRUE ~ NA)) %>% 
  
  mutate(length = case_when(PFAS == "PFHxS_Br" ~ "PFHxS_Br", TRUE ~ "Long-chain")) %>% 
  
  mutate(mw = case_when(PFAS == "PFHpS" ~ 450, PFAS == "PFHxS_Br" ~ 400, PFAS == "PFNA" ~ 464, PFAS == "PFOA_Br" ~ 415,
                        PFAS == "PFOA_L" ~ 415, PFAS == "PFOS_Br" ~ 500, PFAS == "PFOS_L" ~ 500, TRUE ~ NA))



IC9_coef_chemistry <- all_regularized_coefs_structure %>% filter(IC == "IC9")

#Test linear vs branched
ggplot(data = IC9_coef_chemistry %>% filter(!is.na(structure)), aes(x=structure, y=coef_magnitude)) +
  
  geom_boxplot(aes(fill = structure), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "PFAS Structure", y = "Magnitude of IC9 Regularized Coefficient", title = "PFAS structure", 
       subtitle = "t = 3.35, p = 0.066") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        plot.subtitle = element_text(hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

t.test(coef_magnitude ~ structure, data = IC9_coef_chemistry)

ggsave("plots/vectors/IC9_chemistry.svg", units = "px", width = 1250, height = 1250)

all_regularized_coefs_structure %>% 
  filter(!is.na(structure)) %>% 
  ggplot(.) +
  geom_density(aes(x=Coefficient, fill = as.factor(structure), colour = as.factor(structure)), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  scale_colour_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "Regularized Coefficients", y = "Density", colour = "PFAS structure", fill = "PFAS structure") +
  theme_pubclean() +
  xlim(-0.35, 0.35) +
  annotate("text", label = "t = 1.54, p = 0.13", x=-0.25, y=12.5) + 
  ylim(0, 17.5) +
  theme(legend.position = "none", legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

t.test(coef_magnitude ~ structure, data = all_regularized_coefs_structure)
ggsave("plots/vectors/elasticnet_branching.svg", units = "px", width = 1250, height = 1250)

#Test sulfonates vs carboxylates
ggplot(data = IC9_coef_chemistry %>% filter(!is.na(functional_group)), aes(x=functional_group, y=coef_magnitude)) +
  
  geom_boxplot(aes(fill = functional_group), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "Functional Group", y = "Magnitude of IC9 Regularized Coefficient", title = "Functional group", 
       subtitle = "t = -0.20, p = 0.85") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90"))

t.test(coef_magnitude ~ functional_group, data = IC9_coef_chemistry)


all_regularized_coefs_structure %>% 
  filter(!is.na(functional_group)) %>% 
  ggplot(.) +
  geom_density(aes(x=Coefficient, fill = as.factor(functional_group), colour = as.factor(functional_group)), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  scale_colour_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "Regularized Coefficients", y = "Density", colour = "Functional group", fill = "Functional group") +
  theme_pubclean() +
  xlim(-0.35, 0.35) +
  annotate("text", label = "t = 0.97, p = 0.34", x=-0.25, y=12.5) + 
  ylim(0, 17.5) +
  theme(legend.position = "none", legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))

t.test(coef_magnitude ~ functional_group, data = all_regularized_coefs_structure)

ggsave("plots/vectors/elasticnet_functionalgroup.svg", units = "px", width = 1250, height = 1250)
#Correlation coefficients ------------------------------------------------------

#Create datasets
all_correlation_coefs_structure <- heatmap_data %>% 
  rename("PFAS" = "rowname", "IC" = "name", "Coefficient" = "value") %>% 
  mutate(coef_magnitude = if_else(Coefficient < 0, Coefficient*(-1), Coefficient)) %>%
  mutate(carbons = case_when(PFAS == "PFHpS" ~ 7, PFAS == "PFHxS_Br" ~ 6, PFAS == "PFNA" ~ 9, PFAS == "PFOA_Br" ~ 8,
                             PFAS == "PFOA_L" ~ 8, PFAS == "PFOS_Br" ~ 8, PFAS == "PFOS_L" ~ 8, TRUE ~ NA)) %>% 
  
  mutate(structure = case_when(PFAS == "PFHpS" ~ "Linear", PFAS == "PFHxS_Br" ~ "Branched", PFAS == "PFNA" ~ "Linear", PFAS == "PFOA_Br" ~ "Branched",
                               PFAS == "PFOA_L" ~ "Linear", PFAS == "PFOS_Br" ~ "Branched", PFAS == "PFOS_L" ~ "Linear", TRUE ~ NA)) %>% 
  
  mutate(functional_group = case_when(PFAS == "PFHpS" ~ "Sulfonic Acid", PFAS == "PFHxS_Br" ~ "Sulfonic Acid", PFAS == "PFNA" ~ "Carboxylic Acid", PFAS == "PFOA_Br" ~ "Carboxylic Acid",
                                      PFAS == "PFOA_L" ~ "Carboxylic Acid", PFAS == "PFOS_Br" ~ "Sulfonic Acid", PFAS == "PFOS_L" ~ "Sulfonic Acid", TRUE ~ NA)) %>% 
  
  mutate(length = case_when(PFAS == "PFHxS_Br" ~ "PFHxS_Br", TRUE ~ "Long-chain")) %>% 
  
  mutate(mw = case_when(PFAS == "PFHpS" ~ 450, PFAS == "PFHxS_Br" ~ 400, PFAS == "PFNA" ~ 464, PFAS == "PFOA_Br" ~ 415,
                        PFAS == "PFOA_L" ~ 415, PFAS == "PFOS_Br" ~ 500, PFAS == "PFOS_L" ~ 500, TRUE ~ NA))



IC5_coef_chemistry <- all_correlation_coefs_structure %>% filter(IC == "IC5")

#Test linear vs branched
ggplot(data = IC5_coef_chemistry %>% filter(!is.na(structure)), aes(x=structure, y=coef_magnitude)) +
  
  geom_boxplot(aes(fill = structure), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "PFAS Structure", y = " IC5 Correlation Coefficient", title = "PFAS structure", 
       subtitle = "t = -1.16, p = 0.295") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

t.test(coef_magnitude ~ structure, data = IC5_coef_chemistry) 


all_correlation_coefs_structure %>% 
  filter(!is.na(structure)) %>% 
  ggplot(.) +
  geom_density(aes(x=Coefficient, fill = as.factor(structure), colour = as.factor(structure)), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  scale_colour_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "Correlation Coefficients", y = "Density", colour = "PFAS structure", fill = "PFAS structure") +
  theme_pubclean() +
  xlim(-0.5, 0.6) +
  annotate("text", label = "t = -0.49, p = 0.623", x=-0.3, y=3.5) +
  theme(legend.position = "right", legend.background = element_blank()) +
  ylim(0, 4.5) +
  theme(axis.title.y = element_blank(), panel.background = element_blank(), plot.background = element_blank())

t.test(coef_magnitude ~ structure, data = all_correlation_coefs_structure) 
ggsave("plots/vectors/correlations_branching_legend.svg", units = "px", width = 1500, height = 1500)


#Test sulfonates vs carboxylates 
ggplot(data = IC5_coef_chemistry %>% filter(!is.na(functional_group)), aes(x=functional_group, y=coef_magnitude)) +
  
  geom_boxplot(aes(fill = functional_group), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "Functional Group", y = "Magnitude of IC5 Spearman Coefficient", title = "Functional Group", 
       subtitle = "t = 0.76, p = 0.50") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90"))

t.test(coef_magnitude ~ functional_group, data = IC5_coef_chemistry) 


all_correlation_coefs_structure %>% 
  filter(!is.na(functional_group)) %>% 
  ggplot(.) +
  geom_density(aes(x=Coefficient, fill = as.factor(functional_group), colour = as.factor(functional_group)), alpha = 0.3, linewidth = 1.5) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  scale_colour_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "Correlation Coefficients", y = "Density", colour = "Functional group", fill = "Functional group") +
  theme_pubclean() +
  xlim(-0.5, 0.6) +
  annotate("text", label = "t = 3.47, p = 0.001", x=-0.3, y=3.5) +
  theme(legend.position = "none", legend.background = element_blank(), panel.background = element_blank(), plot.background = element_blank()) +
  ylim(0, 4.5)

t.test(coef_magnitude ~ functional_group, data = all_correlation_coefs_structure) 
ggsave("plots/vectors/correlations_functionalgroup.svg", units = "px", width = 1500, height = 1500)

#Hypothalamus ------------------------------------------------------------------

#Create datasets
HTH_chemistry_data <- left_join(HTH_coefficients_plot, HTH_PFAS_correlations, by = "PFAS") %>% 
  rename("Correlation" = "coef", "Regularized" = "Coefficient") %>% 
  mutate(correlation_magnitude = if_else(Correlation < 0, Correlation*(-1), Correlation)) %>% 
  mutate(regularized_magnitude = if_else(Regularized < 0, Regularized*(-1), Regularized)) %>% 
  mutate(carbons = c(7, 6, 9, 8, 8, 8, 8, NA)) %>% 
  mutate(structure = c("Linear", "Branched", "Linear", "Branched", "Linear", "Branched", "Linear", NA)) %>% 
  mutate(mw = c(450, 400, 464, 415, 415, 500, 500, NA)) %>% 
  mutate(functional_group = c("Sulfonic Acid", "Sulfonic Acid", "Carboxylic Acid", "Carboxylic Acid",
                              "Carboxylic Acid", "Sulfonic Acid", "Sulfonic Acid", NA)) %>% 
  mutate(length = case_when(PFAS == "PFHxS_Br" ~ "PFHxS_Br", TRUE ~ "Long-chain"))


#Test linear vs branched
t.test(correlation_magnitude ~ structure, data = HTH_chemistry_data)
t.test(regularized_magnitude ~ structure, data = HTH_chemistry_data)

branch_correlation <- ggplot(data = HTH_chemistry_data %>% filter(!is.na(structure)), aes(x=structure, y=correlation_magnitude)) +
  geom_boxplot(aes(fill = structure), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "PFAS Structure", y = "Magnitude of correlation coefficient", title = "PFAS structure", 
       subtitle = "t = 1.79, p = 0.14") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90"))


branch_elasticnet <- ggplot(data = HTH_chemistry_data %>% filter(!is.na(structure)), aes(x=structure, y=regularized_magnitude)) +
  geom_boxplot(aes(fill = structure), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#8b3bbc", "#3bbc85")) +
  labs(x = "PFAS Structure", y = "Magnitude of Regularized coefficient", title = "PFAS structure", 
       subtitle = "t = 1.19, p = 0.33") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90")) +
  theme(axis.title.y = element_blank())


#Test sulfonates vs carboxylates
t.test(correlation_magnitude ~ functional_group, data = HTH_chemistry_data)
t.test(regularized_magnitude ~ functional_group, data = HTH_chemistry_data)

fg_correlation <-ggplot(data = HTH_chemistry_data %>% filter(!is.na(functional_group)), aes(x=functional_group, y=correlation_magnitude)) +
  geom_boxplot(aes(fill = functional_group), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "PFAS Structure", y = "Magnitude of correlation coefficient", title = "Functional group", 
       subtitle = "t = 0.20, p = 0.85") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90"))


fg_elasticnet <- ggplot(data = HTH_chemistry_data %>% filter(!is.na(functional_group)), aes(x=functional_group, y=regularized_magnitude)) +
  geom_boxplot(aes(fill = functional_group), colour = "black", alpha = 0.8, linewidth = 1.25) +
  scale_fill_manual(values = c("#ded51c", "#cb2512")) +
  labs(x = "PFAS Structure", y = "Magnitude of regularized coefficient", title = "Functional group", 
       subtitle = "t = -4.05, p = 0.015") +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_markdown(hjust = 0.5), 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white", linetype = 1),
        panel.grid = element_line(color = "gray90"))

ggarrange(fg_elasticnet, branch_elasticnet)
