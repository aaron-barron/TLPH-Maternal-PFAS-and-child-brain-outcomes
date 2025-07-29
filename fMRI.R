#Create dataset
source("clean scripts/Load and clean data.R")

fMRI_clusters <- PFAS_FLICA_covar_data %>% 
  filter(!is.na(PFNA_cluster))

#Scatterplots of significant ReHo clusters from SPM voxel-wise analyses --------

fmri1 <- ggplot(data = fMRI_clusters, aes(x = PFNA, y = PFNA_cluster)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linewidth = 2, colour = "#CB4335") + 
  labs(y = "ReHo") +
  ylim(0.6, 2.1) +
  theme(plot.background = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

fmri2 <- ggplot(data = fMRI_clusters, aes(x = PFOA_L, y = PFOA_L_cluster)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linewidth = 2, colour = "#F1C40F") + 
  labs(y = "ReHo", x = "PFOA-L") +
  ylim(0.6, 2.1) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

fmri4 <- ggplot(data = fMRI_clusters, aes(x = PFOA_Br, y = PFOA_Br_cluster)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linewidth = 2, colour = "#28B463") +
  labs(y = "ReHo", x = "PFOA-Br") +
  ylim(0.6, 2.1) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

fmri3 <- ggplot(data = fMRI_clusters, aes(x = PFHxS_Br, y = PFHxS_Br_cluster)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", linewidth = 2, colour = "Dark Blue") + 
  labs(y = "ReHo", x = "PFHxS") +
  ylim(0.6, 2.1) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

ggarrange(fmri1, fmri2, fmri3, fmri4, ncol = 4, nrow = 1)

ggsave("plots/vectors/fmri_scatterplots.svg", units = "px", width = 4000, height = 1000)

#Compare ReHo cluster peaks to distribution of intensity of IC9 GMV ------------

#Load the data, a NIfTI file which has been converted to ASCII table format
IC9_mi1 <- read.table("data/IC9_mi1_ascii00000") %>% 
  pivot_longer(1:182)

mean(IC9_mi1$value) #0.3332902
sd(IC9_mi1$value) #1.127273

#The distance of the ReHo peaks (lowest one is IC9 = 4.56) from the mean:
(4.56 - mean(IC9_mi1$value))/sd(IC9_mi1$value) #3.7495

#Plot histogram (of non-zero voxels)
IC9_mi1 %>%  
  filter(value != 0) %>% 
  ggplot(data = ., aes(x = value, y = after_stat(count / sum(count))*100)) + 
  geom_histogram(fill = "#C4193B", alpha = 0.8, binwidth = 0.5) +
  geom_vline(xintercept = c(4.5, 4.6, 4.8, 6.7), linetype="dashed") + 
  geom_line(data=tibble(x=c(1.25, 4.25), y = 28), aes(x=x, y=y)) +
  annotate("text", label = "Coordinates
of ReHo 
cluster peaks", x=9.5, y=20, size = 5) +
  annotate("text", label = "3.75 SD", x=2.65, y=30, size = 5) +
  labs(x = "IC9 Grey Matter Volume Intensity", y = "%") +
  scale_x_continuous(n.breaks = 10) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"),
        axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15))

ggsave("plots/vectors/fmri_IC9_histogram.svg", units = "px", width = 2000, height = 2000)
