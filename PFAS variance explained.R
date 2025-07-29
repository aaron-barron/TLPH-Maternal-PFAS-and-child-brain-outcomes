#Create dataset
source("Load and clean data.R")

#Load additional libraries
library("scater")
library("ggpubr")

#Get variance explained --------------------------------------------------------
data_PFAS <- PFAS_FLICA_covar_data[3:10] %>% t() %>% as.matrix()
data_variables <- PFAS_FLICA_covar_data[22:29] 

variance <- getVarianceExplained(x=data_PFAS, variables=data_variables)

plotExplanatoryVariables(variance, theme_size = 10) 

variance_dataframe <- variance %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("Variable") %>%  
  mutate(number = c(1:8))

#Plots -------------------------------------------------------------------------

#Colours
colours <- c("#ddadef", "#DAF7A6", "#FFC300", "#bd5674", "#EB984E", "#3bbc78", "#AED6F1", "#3b9dbc")

#PFHpS
p1 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFHpS), y=PFHpS)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFHpS") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") 

#PFHxS_Br
p2 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFHxS_Br), y=PFHxS_Br)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +   
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFHxS_Br") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

#PFNA
p3 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFNA), y=PFNA)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFNA") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

#PFOA_Br
p4 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFOA_Br), y=PFOA_Br)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFOA_Br") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

#PFOA_L
p5 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFOA_L), y=PFOA_L)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFOA_L") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed")


#PFOS_Br
p6 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFOS_Br), y=PFOS_Br)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFOS_Br") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

#PFOS_L
p7 <- ggplot(variance_dataframe, aes(x=reorder(Variable, PFOS_L), y=PFOS_L)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "PFOS_L") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

#Total (WQS) PFAS
p8 <- ggplot(variance_dataframe, aes(x=reorder(Variable, `Total (WQS) PFAS`), y=`Total (WQS) PFAS`)) +
  geom_bar(stat = "identity", colour = "black", aes(fill = reorder(Variable, number))) +
  scale_fill_manual(values = colours) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  labs(y = "% Variance Explained", title = "Total (WQS) PFAS") + #format titles
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, size = 10), plot.background = element_blank(), panel.background = element_blank(), panel.grid = element_line(color = "gray90") ) +
  geom_hline(yintercept =0.5, linetype="dashed") + theme(axis.title.y = element_blank())

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 2, ncol = 4)
ggsave("plots/vectors/PFAS_variance_explained.svg", units = "px", width = 3000, height = 1500)
