#Set seed for reproducibility
set.seed(10)

#Load libraries
library("conflicted")
library("tidyverse")
library("missForest")
library("batchtma")
library("wqs")

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("cbind", "base")
conflict_prefer("rbind", "base")

#Read dataset into environment
PFAS_FLICA_covar_data <- read.csv("PFAS_FLICA_full_dataset.csv", check.names = F)

#Impute missing covariate data with mice (BMI and birthweight) -----------------

PFAS_FLICA_covar_data_imputations <- PFAS_FLICA_covar_data[, 2:46] %>% 
  select(-`Child Age`) #we don't wish to impute "age at time of scan" for children not scanned

colnames(PFAS_FLICA_covar_data_imputations) <- make.names(colnames(PFAS_FLICA_covar_data_imputations)) 

PFAS_FLICA_covar_data_imputations <- missForest(PFAS_FLICA_covar_data_imputations)$ximp



PFAS_FLICA_covar_data$`Maternal BMI` <- PFAS_FLICA_covar_data_imputations$Maternal.BMI
PFAS_FLICA_covar_data$`Birth Weight` <- PFAS_FLICA_covar_data_imputations$Birth.Weight

#Impute missing analyte values with half minima (variables 3:34) ---------------
na_to_zero_function <- function(x) {replace(x, is.na(x), 0)}
replace_function <- function(x) {replace(x, x == 0, min(x[x>0])/2)}

PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% 
  mutate_at(3:34, na_to_zero_function) %>% 
  mutate_at(3:34, replace_function)



#PCA and check for batch effects among PFAS ------------------------------------
PFAS_PCA <- prcomp(PFAS_FLICA_covar_data[3:9], center = T, scale. = T)

summary(PFAS_PCA)

PFAS_PCA_dataframe <- as.data.frame(PFAS_PCA$x) %>% 
  select(PC1, PC2) %>% 
  mutate(Batch = PFAS_FLICA_covar_data$Batch) 

ggplot(PFAS_PCA_dataframe, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = as.factor(Batch)), size = 2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Batch)) , colour = "black", alpha = 0.25, level = 0.995) +
  labs(title = "PCA for PFAS", x = "PC1 (39.3%)", y = "PC2 (18%)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

ggsave("plots/vectors/pca_pfas.svg", units = "px", width = 1250, height = 1250)

#PCA and check for batch effects among bile acids ------------------------------
BA_PCA <- prcomp(PFAS_FLICA_covar_data[10:31], center = T, scale. = T)

summary(BA_PCA)

BA_PCA_dataframe <- as.data.frame(BA_PCA$x) %>% 
  select(PC1, PC2) %>% 
  mutate(Batch = PFAS_FLICA_covar_data$Batch) 

ggplot(BA_PCA_dataframe, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = as.factor(Batch)), size = 2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Batch)) , colour = "black", alpha = 0.25, level = 0.995) +
  labs(title = "PCA for Bile Acids", x = "PC1 (25.1%)", y = "PC2 (12.2%)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

ggsave("plots/vectors/pca_ba.svg", units = "px", width = 1250, height = 1250)

BA_corrected  <- PFAS_FLICA_covar_data[c(1:2, 10:31)] %>% 
  adjust_batch(markers = 3:24, batch = Batch, 
               method = quantnorm)

BA_corrected_PCA <- prcomp(BA_corrected[25:46], center = T, scale. = T)

summary(BA_corrected_PCA)

BA_corrected_PCA_dataframe <- as.data.frame(BA_corrected_PCA$x) %>% 
  select(PC1, PC2) %>% 
  mutate(Batch = PFAS_FLICA_covar_data$Batch) 

ggplot(BA_corrected_PCA_dataframe, aes(x=PC1, y=PC2)) +
  geom_point(aes(colour = as.factor(Batch)), size = 2) +
  geom_point(shape = 1,size = 2,colour = "black") +
  stat_ellipse(geom = "polygon", aes(fill = as.factor(Batch)) , colour = "black", alpha = 0.25, level = 0.995) +
  labs(title = "PCA for Bile Acids (Corrected)", x = "PC1 (24.3%)", y = "PC2 (12.0%)") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  theme(plot.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_line(color = "gray90"))

ggsave("plots/vectors/pca_ba_corrected.svg", units = "px", width = 1250, height = 1250)

BA_corrected <- BA_corrected[c(1:2, 25:46)]

colnames(BA_corrected) <- colnames(PFAS_FLICA_covar_data)[c(1:2, 10:31)]

#Total PFAS variable with WQS regression ---------------------------------------
set.seed(5)

PFAS <- PFAS_FLICA_covar_data[, 3:9] %>%
  mutate(sum_PFAS = rowSums(.[1:7])) %>% 
  select(sum_PFAS, everything())

# Outcome variable = simple sum-total PFAS
y.train <- PFAS$sum_PFAS

# Dependent variables = individual PFAS
x.train <- PFAS %>% select(-sum_PFAS) 

# Perform WQS. Output will be a vector corresponding to the number of rows in PFAS dataframe.
PFAS_WQS <- wqs.est(y.train, x.train, B = 100)

#The output I need is save in PFAS_WQS$WQS

#Add the new variable to my main PFAS dataframe
PFAS <- PFAS %>% 
  mutate("Total (WQS) PFAS" = PFAS_WQS$WQS)

#WQS created a few zero values, so replace these with half-minima 
PFAS <- PFAS %>% 
  mutate_at(9, replace_function)

#Make donut plot of weights
WQS_weights <- colnames(PFAS[, 2:8]) %>% 
  as.data.frame() %>% 
  mutate(weights = PFAS_WQS$weights*100) %>% 
  arrange(weights) %>% 
  mutate(ymax = cumsum(weights)) %>% 
  mutate(ymin = c(0, head(ymax, n=-1))) %>% 
  mutate(weights = as.numeric(str_sub(weights, end = 4))) %>% 
  rename("PFAS" = ".") %>% 
  mutate(PFAS2 = PFAS) %>% 
  mutate(weights2 = weights) %>% 
  mutate(percent = "%") %>% 
  mutate(colon = ": ") %>% 
  unite(col = "label", sep = "", c(PFAS2, colon, weights2, percent)) %>% 
  mutate(labelPosition = (ymax + ymin) / 2) 


ggplot(data = WQS_weights, aes(ymax=ymax, ymin=ymin, xmax =8, xmin=7, fill= reorder(label, weights))) +
  geom_rect() +
  scale_fill_brewer(palette= 16) +
  coord_polar(theta="y", start = 6) + # Try to remove that to understand how the chart is built initially
  xlim(c(6, 8)) + # Try to remove that to see how to make a pie chart
  theme_void() +
  labs(title = "WQS-PFAS Weighting") +
  theme(legend.title=element_blank(), plot.title = element_text(hjust = 1), legend.position = "bottom", 
        legend.direction = "vertical", title = element_text(size = 10))

ggsave("plots/individual panels/1C.svg", units = "px", width = 1500, height = 1250)

#Group bile acids by class -----------------------------------------------------
BA_groups <- BA_corrected %>% 
  select(-Batch) %>% 
  mutate(`Unconjugated-Primary` = (CA + CDCA)/2) %>% 
  mutate(`Unconjugated-Secondary` = (DCA + HDCA + LCA + UDCA)/4) %>%
  mutate(`Glycine-Primary` = (GCA + GCDCA + GHCA)/3) %>%
  mutate(`Glycine-Secondary` = (GDCA + GHDCA + GHDCA_X2 + GLCA)/4) %>%
  mutate(`Taurine-Primary` = (`tauro-beta-Muricholicacid` + TCA + TCDCA + TMCA + THCA)/5) %>%
  mutate(`Taurine-Secondary` = (TLCA + TUDCA + TDCA + THDCA)/4) %>%
  select(-CA, -CDCA, -DCA, -GCA, -GCDCA, -GDCA, -GHCA, -GHDCA, -GHDCA_X2, -GLCA, -HDCA, -LCA, -`tauro-beta-Muricholicacid`, -TCA, -TCDCA, -TDCA, -THCA, -THDCA, -TLCA, -TMCA, -TUDCA, -UDCA)


#Add WQS PFAS and bile acid groups to main dataset -----------------------------
PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>% 
  mutate(`Total (WQS) PFAS` = PFAS$`Total (WQS) PFAS`) %>% 
  left_join(., BA_groups, by = "ID") %>% 
  select(1:9, 69, 70:75, 32:68)

#Transform and autoscale -------------------------------------------------------
autoscale_function <- function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T)}

PFAS_FLICA_covar_data <- PFAS_FLICA_covar_data %>%
  mutate_at(3:19, log2) %>%  #Log2 transform all analytes
  mutate_at(c(3:19, 21:23, 30:31), autoscale_function) #autoscale all analytes and numeric covariates


