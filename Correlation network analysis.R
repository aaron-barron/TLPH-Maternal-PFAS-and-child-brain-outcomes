#Create dataset
source("Load and clean data with imputations.R")

#Load additional libraries
library("corrplot")
library("qpgraph")
library("circlize")
library("ComplexHeatmap")

#Compute correlations ----------------------------------------------------------

#Remove superfluous variables: ID and Batch, hypothalamus MD, ReHo clusters
correlation_network_data <- PFAS_FLICA_covar_data %>% 
  select(3:47) %>% 
  select(-Parity, -Income)

#Compute matrix
cor_matrix <- as.matrix(cor(correlation_network_data, use='pairwise.complete.obs', method = 'spearman'))

#Corrplot
svglite("plots/vectors/corrplot.svg", width=1250, height=1250, bg = "transparent")
corrplot(corr = cor(correlation_network_data, use='pairwise.complete.obs', method = 'spearman'), 
                     type = "lower", method = "color", number.cex = 0.2, tl.col = "black", bg = "transparent")

dev.off()


#Compute non-rejection rates (NRRs) --------------------------------------------
set.seed(12)

NRR_matrix <- qpNrr(cor_matrix) %>% as.matrix()

NRR_dataframe <- as.data.frame(NRR_matrix)

#Plot distribution of all NRRs
NRR_dataframe %>% 
  pivot_longer(cols = 1:43) %>% 
  ggplot(data = ., aes(x = value, y = after_stat(count / sum(count))*100)) + 
  geom_histogram(fill = "#2980B9", alpha = 0.9) +
  geom_vline(xintercept = 0.15, linetype="dashed") + 
  annotate("text", label = c("Projected", "Discarded"), x=c(0.04, 0.26), y=12) +
  labs(x = "Non-Rejection Rate", y = "%") +
  scale_x_continuous(n.breaks = 10) +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

#Plot distribution of all NRRs below 1.0
NRR_dataframe %>% 
  pivot_longer(cols = 1:43) %>% 
  filter(value < 1) %>% 
  ggplot(data = ., aes(x = value, y = after_stat(count / sum(count))*100)) + 
  geom_histogram(fill = "#2980B9", alpha = 0.9) +
  geom_vline(xintercept = 0.15, linetype="dashed") + 
  annotate("text", label = c("Projected", "Discarded"), x=c(0.075, 0.225), y=12) +
  labs(x = "Non-Rejection Rate", y = "%") +
  scale_x_continuous(n.breaks = 10)

#Based on this distribution, NRR < 0.15 was chosen, and corresponding annotations added to the histograms
#Save a table of binary T/F values of whether each correlation passes this threshold
NRR_0.15 <- as.data.frame(as(NRR_matrix, "matrix") < 0.15)

#Prepare data for Circos plot --------------------------------------------------

#Replace all intra-dataset associations with zeros

cor_matrix_projection <- as.data.frame(cor_matrix) 
cor_matrix_projection[1:8, 1:8] <- 0 #removes PFAS x PFAS
cor_matrix_projection[9:14, 9:14] <- 0 #removes bile acids x bile acids
cor_matrix_projection[15, 15] <- 0 #removes CRP x CRP
cor_matrix_projection[16:17, 16:17] <- 0 #removes thyroid x thyroid
cor_matrix_projection[18:27, 18:27] <- 0 #removes confounders x confounders
cor_matrix_projection[29:37, 29:37] <- 0 #removes ICs x ICs
cor_matrix_projection[38:43, 38:43] <- 0 #removes SDQ x SDQ

cor_matrix_projection_long <- cor_matrix_projection %>% 
  rownames_to_column() %>% 
  pivot_longer(2:44)

#Replace correlations with NRR > 0.15 with zeros

NRR_0.15_long <- rownames_to_column(NRR_0.15) %>% 
  pivot_longer(2:44)

cor_matrix_projection_NRR_0.15 <- cbind(cor_matrix_projection_long, NRR_0.15_long)

colnames(cor_matrix_projection_NRR_0.15) <- c("rowname", "name", "value", "NRR_rowname", "NRR_name", "NRR_value")

cor_matrix_projection_NRR_0.15 <- cor_matrix_projection_NRR_0.15 %>% 
  mutate(new_data = if_else(NRR_value == F, 0, value)) %>% 
  replace(is.na(.), 0) %>% 
  select(rowname, name, new_data) %>% 
  pivot_wider(values_from = "new_data", names_from = "name") %>% 
  column_to_rownames()

# Create a colour list vector for all variables
grid_col <- c(PFHpS= "lightblue", PFHxS_Br= "lightblue", PFNA= "lightblue", PFOA_Br= "lightblue",
              PFOA_L= "lightblue", PFOS_Br= "lightblue", PFOS_L= "lightblue",  `Total (WQS) PFAS`= "lightblue",
              `Maternal Age`= "maroon", `Maternal Education` = "maroon", Gravidity = "maroon", `Maternal Smoking` = "maroon", 
              `Maternal BMI` = "maroon", Sex= "maroon", `Child Age` = "maroon", `Birth Weight` = "maroon", `Gestational Weeks` = "maroon",
              `Maternal Alcohol`= "maroon", IC1= "lightgrey", IC2= "lightgrey", IC3= "lightgrey",
              IC4= "lightgrey", IC5= "lightgrey", IC6= "lightgrey", IC7= "lightgrey",
              IC8= "lightgrey", IC9= "lightgrey", IC10= "lightgrey", `Emotional Problems` = "#EB984E", 
              `Conduct Problems` = "#EB984E", `Hyperactivity/Inattention` = "#EB984E", 
              `Peer Relationship Problems` = "#EB984E", `Prosocial Behaviour` = "#EB984E",
              `Total Difficulties`= "#EB984E", CRP= "#F8E92E", TSH = "#B87EF2", T4 = "#B87EF2",
              `Unconjugated-Primary` = "#72E16D", `Unconjugated-Secondary` = "#72E16D", `Glycine-Primary` = "#72E16D",
              `Glycine-Secondary` = "#72E16D", `Taurine-Primary` = "#72E16D", `Taurine-Secondary` = "#72E16D")

# Create a legend for all 

variable_types <- as.data.frame(x = c("Maternal PFAS", "Maternal CRP", "Maternal Thyroid Hormones",
                                      "Maternal Bile Acids", "Child Brain Structural IC", 
                                      "Clinical/Demographic Variable"
                                      ,"SDQ"
))
variable_colours <- as.data.frame(x = c("lightblue", "#F8E92E", "#B87EF2", "#72E16D","lightgrey","maroon"
                                        ,"#EB984E"
))

variable_types <- cbind(variable_types, variable_colours)

colnames(variable_types) <- c("type", "colour")

legend_variables = Legend(at = c(variable_types$type), type = "points", 
                          legend_gp = gpar(col = variable_types$colour), title_position = "topleft", 
                          title = "Variable Type", size = unit(12.5, "points"), grid_height = unit(12.5, "points"),
                          grid_width = unit(12.5, "points"),
                          labels_gp = gpar(fontsize = 12.5),
                          title_gp = gpar(fontsize = 12.5, fontface = "bold"), nrow = 2)

legend_links = Legend(at = c("Positive", "Negative"), type = "lines", 
                      legend_gp = gpar(col = c("#2980B9", "#EB984E"), lwd = 7.5), title_position = "topleft", 
                      title = "Association", grid_height = unit(12.5, "points"),
                      grid_width = unit(12.5, "points"),
                      labels_gp = gpar(fontsize = 12.5),
                      title_gp = gpar(fontsize = 12.5, fontface = "bold"), nrow = 1)

legend_full = packLegend(legend_links, legend_variables, direction = "horizontal")

#Circos plot for NRR 0.15 ------------------------------------------------------

#Create the data matrix 
circos_matrix_NRR_0.15 <- cor_matrix_projection_NRR_0.15 %>% 
  rownames_to_column() %>% 
  pivot_longer(2:44) %>% 
  filter(value == 0 | !duplicated(value)) %>% 
  pivot_wider(values_from = "value", names_from = "name") %>% 
  column_to_rownames() %>%  as.matrix()


#Create a matrix for the colours of the links - blue for positive values and orange for negative values
link_colours <- cor_matrix_projection_NRR_0.15 %>%  
  rownames_to_column() %>% 
  pivot_longer(2:44) %>% 
  filter(value == 0 | !duplicated(value)) %>% 
  mutate(colour = case_when(value > 0 ~ "#2980B9", value < 0 ~ "#EB984E",)) %>% 
  
  mutate(colour = case_when(value > 0 & value < 0.15 ~ adjustcolor(colour, alpha = 0.3),
                            value > 0.15 & value < 0.2 ~ adjustcolor(colour, alpha = 0.4),
                            value > 0.2 & value < 0.25 ~ adjustcolor(colour, alpha = 0.5),
                            value > 0.25 & value < 0.3 ~ adjustcolor(colour, alpha = 0.6),
                            value > 0.3 & value < 0.35 ~ adjustcolor(colour, alpha = 0.7),
                            value > 0.35 & value < 0.4 ~ adjustcolor(colour, alpha = 0.8),
                            value > 0.4  ~ adjustcolor(colour, alpha = 0.9),
                            value < 0 & value > -0.15 ~ adjustcolor(colour, alpha = 0.3),
                            value < -0.15 & value > -0.2 ~ adjustcolor(colour, alpha = 0.4),
                            value < -0.2 & value > -0.25 ~ adjustcolor(colour, alpha = 0.5),
                            value < -0.25 & value > -0.3 ~ adjustcolor(colour, alpha = 0.6),
                            value < -0.3 & value > -0.35 ~ adjustcolor(colour, alpha = 0.7),
                            value < -0.35 & value > -0.4 ~ adjustcolor(colour, alpha = 0.8),
                            value < -0.4  ~ adjustcolor(colour, alpha = 0.9))) %>% 
  select(-value) %>% 
  pivot_wider(values_from = "colour", names_from = "name") %>% 
  column_to_rownames() %>%  as.matrix()

#Run this line each time I want to edit the diagram, as opposed to build upon it
svg("plots/vectors/circos_plot.svg", width=12, height=12, bg = "transparent")

par(bg="transparent")

circos.clear()

#Rotate the plot as you wish. This needs to be run after circos.clear, but before making and calling the plot. So lots of back and forth.
circos.par(start.degree = 220)

#Create the chord diagram
chordDiagram(circos_matrix_NRR_0.15, 
             grid.col = grid_col, 
             annotationTrack = "grid",
             preAllocateTracks = 1, #This creates an empty "track", which we will later populate with labels/names
             col = link_colours,
             link.border = "white",
             order = colnames(circos_matrix_NRR_0.15)
             
)



#Create sector labels
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 2.5, sector.name, 
              facing = "clockwise", 
              niceFacing = TRUE, #adjust for human eyes (T/F)
              adj = c(0.1, 0.1), #label rotation
              cex=1.125) #fontsize
  
  
}, bg.border = NA)

# Add the legend to the plot
draw(legend_full, x = unit(0.1, "npc"), y = unit(0.01, "npc"), just = c("left", "bottom"))

#Save
dev.off()
circos.clear()

#Heatmap  ----------------------------------------------------------------------

#Create data
heatmap_data <- cor_matrix[1:8, 28:37] %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(2:11) %>% 
  arrange(desc(name)) %>% 
  mutate(name = factor(name, ordered = T, levels = 
                         c("IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8", "IC9", "IC10"))) %>% 
  mutate(rowname = factor(rowname, ordered = T, levels = 
                            c("Total (WQS) PFAS", "PFOS_L","PFOS_Br", "PFOA_L", "PFOA_Br", "PFNA", "PFHxS_Br", "PFHpS")))

#Create colour scale
heatmap_colours_function <- colorRampPalette(c("#EB984E", "white", "#2980B9"))
heatmap_colours <- heatmap_colours_function(100)

#Plot
ggplot(data = heatmap_data, aes(x = name, y = rowname)) +
  geom_tile(aes(fill =  value)) +
  labs(fill = "Correlation
Coefficient") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "left") +
  theme(panel.background = element_blank(), plot.background = element_blank(), legend.text = element_blank(), legend.title = element_blank()) +
  scale_fill_gradient2(low = "#EB984E", mid  = "white", high = "#2980B9") +
  scale_y_discrete(position = "right")

ggsave("plots/individual panels/2B unannotated.svg", units = "px", width = 3000, height = 1000)
