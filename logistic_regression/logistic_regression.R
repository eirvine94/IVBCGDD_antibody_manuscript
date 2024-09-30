#####
# Author: Edward B. Irvine, PhD
#
# Description: 
# This script processes antibody feature data to compute area under the curve (AUC) for 
# both plasma and bronchoalveolar lavage (BAL) samples across different time points. It includes data 
# preprocessing, variance filtering, and logistic regression analysis to identify relationships between 
# antibody features and infection status. The script generates a volcano plot to visualize significant 
# antibody features identified through regression analysis.
#
# Created: 25 September 2021
#
# Modified: 30 September 2024
#####

##################
##### Housekeeping ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
##################

# Load required libraries
library(dplyr)
library(MESS)
library(stringr)
library(ggplot2)
library(ggrepel)

# Import data
plasma <- read.csv("eth_plasma_FC_plsda.csv")
BAL <- read.csv("eth_BAL_FC_plsda.csv")
meta <- read.csv("IVDD_meta_data.csv")





######################
##### Define functions ---------------------------------------------------------------------------------------------------------------------------------------------------------------
######################

# Remove low variance features
var_threshold <- function(dat, threshold) {
  #####
  # This function removes low variance features from a dataset.
  # Inputs:
  #   - dat: A data frame with the first column containing labels (e.g., "Animal") and subsequent columns containing numeric features.
  #   - threshold: A numeric value specifying the minimum variance threshold. Features with variance below this value will be removed.
  # Output:
  #   - A data frame containing only the "Animal" column and features with variance above the specified threshold.
  #####
  
  # Remove rows with missing values
  dat_2 <- na.omit(dat)
  
  # Calculate variance for each feature (excluding the first column)
  var_frame <- data.frame(sapply(dat_2[ , 2:ncol(dat_2)], var))
  colnames(var_frame)[1] <- "var_auc"
  
  # Filter features with variance greater than the threshold
  var_frame$feature <- rownames(var_frame)
  var_frame <- var_frame[var_frame$var_auc > threshold, ]
  
  # Keep the label column ("Animal") and selected features
  keep <- c("Animal", rownames(var_frame))
  
  return(dat[ , keep])
}





######################
##### Pre-process data ------------------------------------------------------------------------------------------------------------------------------------------------------------------
######################

animals <- meta$Animal
features <- colnames(plasma)[4:ncol(plasma)]

###
### BAL AUC
###
# Make time series data
BAL$Timepoint[BAL$Timepoint == "Pre-BCG"] <- 0
BAL$Timepoint[BAL$Timepoint == "Week 4 Post-BCG"] <- 4
BAL$Timepoint[BAL$Timepoint == "Week 9 Post-BCG"] <- 9
BAL$Timepoint[BAL$Timepoint == "Week 12 Post-BCG"] <- 12
BAL$Timepoint <- as.numeric(BAL$Timepoint)

# Create list. Names are animal names. Elements are time series data for all antibody features
BAL_list <- vector("list", length(animals))
names(BAL_list) <- animals
for (animal in animals) {
  BAL_list[[animal]] <- BAL %>%
    filter(Animal == animal) %>%
    arrange(Timepoint)
}

# Initialize data frame to store AUC data
BAL_auc <- data.frame(matrix(ncol = 91, nrow = 34))
colnames(BAL_auc) <- c("Animal", features)
BAL_auc$Animal <- animals
rownames(BAL_auc) <- animals

# Compute AUC for each antibody feature in each animal
for (animal in animals) {
  if (nrow(BAL_list[[animal]]) == 4) {
    for (feature in features) {
      if (!any(is.na(BAL_list[[animal]][ , feature]))) {
        BAL_auc[animal, feature] <- auc(BAL_list[[animal]][ , "Timepoint"], BAL_list[[animal]][ , feature])
      }
    }
  }
}

# Remove low variance features
BAL_auc <- var_threshold(dat = BAL_auc, threshold = 1)

# Combine meta and BAL auc data
BAL_auc <- data.frame(cbind(meta$Animal, meta$dose, meta$log_total_cfu, meta$mtb_dose, meta$cohort, meta$gender, meta$holding, meta$age, BAL_auc[ , 2:ncol(BAL_auc)]))
colnames(BAL_auc) <- str_remove(colnames(BAL_auc), "meta.")

# Prep covariates for analysis
BAL_auc$cohort <- as.factor(BAL_auc$cohort)
BAL_auc$gender <- as.factor(BAL_auc$gender)
BAL_auc$holding <- as.factor(BAL_auc$holding)
BAL_auc$log_total_cfu[BAL_auc$log_total_cfu == 0] <- 0
BAL_auc$log_total_cfu[BAL_auc$log_total_cfu > 0] <- 1
BAL_auc$log_total_cfu <- as.factor(BAL_auc$log_total_cfu)



###
### Plasma AUC
###
# Make time series data
plasma$Timepoint[plasma$Timepoint == "Pre-BCG"] <- 0
plasma$Timepoint[plasma$Timepoint == "Week 4 Post-BCG"] <- 4
plasma$Timepoint[plasma$Timepoint == "Week 9 Post-BCG"] <- 9
plasma$Timepoint[plasma$Timepoint == "Week 12 Post-BCG"] <- 12
plasma$Timepoint[plasma$Timepoint == "2 Weeks Pre-Mtb"] <- 22
plasma$Timepoint[plasma$Timepoint == "Week 2 Post-Mtb"] <- 26
plasma$Timepoint[plasma$Timepoint == "Week 4 Post-Mtb"] <- 28
plasma$Timepoint[plasma$Timepoint == "Week 12 Post-Mtb"] <- 36
plasma$Timepoint <- as.numeric(plasma$Timepoint)


### Pre-Mtb challenge
# Make pre-Mtb challenge subset
plasma_preTB <- plasma[plasma$Timepoint < 24, ]

# Create list. Names are animal names. Elements are time series data for all antibody features
plasma_preTB_list <- vector("list", length(animals))
names(plasma_preTB_list) <- animals
for (animal in animals) {
  plasma_preTB_list[[animal]] <- plasma_preTB %>%
    filter(Animal == animal) %>%
    arrange(Timepoint)
}

# Initialize data frame to store AUC data
plasma_preTB_auc <- data.frame(matrix(ncol = 91, nrow = 34))
colnames(plasma_preTB_auc) <- c("Animal", features)
plasma_preTB_auc$Animal <- animals
rownames(plasma_preTB_auc) <- animals

# Compute AUC for each antibody feature in each animal
for (animal in animals) {
  if (nrow(plasma_preTB_list[[animal]]) == 5) {
    for (feature in features) {
      if (!any(is.na(plasma_preTB_list[[animal]][ , feature]))) {
        plasma_preTB_auc[animal, feature] <- auc(plasma_preTB_list[[animal]][ , "Timepoint"], plasma_preTB_list[[animal]][ , feature])
      }
    }
  }
}

# Remove low variance features
plasma_preTB_auc <- var_threshold(dat = plasma_preTB_auc, threshold = 30)



### Post-Mtb challenge
# Make post-Mtb challenge subset
plasma_postTB <- plasma[plasma$Timepoint > 24, ]

# Create list. Names are animal names. Elements are time series data for all antibody features
plasma_postTB_list <- vector("list", length(animals))
names(plasma_postTB_list) <- animals
for (animal in animals) {
  plasma_postTB_list[[animal]] <- plasma_postTB %>%
    filter(Animal == animal) %>%
    arrange(Timepoint)
}

# Initialize data frame to store AUC data
plasma_postTB_auc <- data.frame(matrix(ncol = 91, nrow = 34))
colnames(plasma_postTB_auc) <- c("Animal", features)
plasma_postTB_auc$Animal <- animals
rownames(plasma_postTB_auc) <- animals

# Compute AUC for each antibody feature in each animal
for (animal in animals) {
  if (nrow(plasma_postTB_list[[animal]]) == 3) {
    for (feature in features) {
      if (!any(is.na(plasma_postTB_list[[animal]][ , feature]))) {
        plasma_postTB_auc[animal, feature] <- auc(plasma_postTB_list[[animal]][ , "Timepoint"], plasma_postTB_list[[animal]][ , feature])
      }
    }
  }
}

# Remove low variance features
plasma_postTB_auc <- var_threshold(dat = plasma_postTB_auc, threshold = 15)



###
### Combine BAL + Plasma
###
colnames(BAL_auc)[9:ncol(BAL_auc)] <- str_c("BAL.Pre.", colnames(BAL_auc)[9:ncol(BAL_auc)])
colnames(plasma_preTB_auc)[2:ncol(plasma_preTB_auc)] <- str_c("Plasma.Pre.", colnames(plasma_preTB_auc)[2:ncol(plasma_preTB_auc)])
colnames(plasma_postTB_auc)[2:ncol(plasma_postTB_auc)] <- str_c("Plasma.Post.", colnames(plasma_postTB_auc)[2:ncol(plasma_postTB_auc)])
combined_auc <- data.frame(cbind(BAL_auc, plasma_preTB_auc[ , 2:ncol(plasma_preTB_auc)], plasma_postTB_auc[ , 2:ncol(plasma_postTB_auc)]))
dat <- combined_auc





#########################
##### Logistic regression ---------------------------------------------------------------------------------------------------------------------------------------------------------------
#########################

# Set plot aesthetic parameters
options(ggrepel.max.overlaps = Inf)
plasma.col <- "#FA8334"
BAL.col <- "#048A81"
norm.size <- 14
small.size <- 12
shapes <- c(16, 15)
    
# Initialize data frame to hold regression output
lm_out_doseCTR <- setNames(data.frame(matrix(ncol = 3, nrow = ncol(dat) - 7)), c("feature", "p_val", "coeff"))

# Fit logistic regression models for each antibody feature
count <- 1
for (i in colnames(dat[9:ncol(dat)])) {
  
  # Fit model
  fit <- glm(log_total_cfu ~ age + gender + dose + eval(parse(text = i)), data = dat, family = "binomial") 
  fit.summary <- data.frame(coef(summary(fit)))
  
  # Pull out logistic regression coefficients for antibody feature
  coeff <- data.frame(fit.summary$Estimate)
  feature_coeff <- coeff[5, 1] 
  
  # Pull out p-value associated with given antibody feature
  p_vals <- data.frame(fit.summary$Pr...z..)
  feature_p_val <- p_vals[5, 1] 
  
  # Store feature name, coefficient and p-value
  lm_out_doseCTR$feature[count] <- i
  lm_out_doseCTR$coeff[count] <- feature_coeff
  lm_out_doseCTR$p_val[count] <- feature_p_val
  
  count <- count + 1
}

# FDR p-value adjustment 
lm_out_doseCTR$fdr <- p.adjust(lm_out_doseCTR$p_val, method = "fdr")

# Create column storing the -log10 p-value
lm_out_doseCTR$neglog <- -log10(lm_out_doseCTR$p_val)

# Add metadata
lm_out_doseCTR$Compartment <- rep("Plasma", nrow(lm_out_doseCTR))
lm_out_doseCTR$Compartment[grepl("BAL", lm_out_doseCTR$feature)] <- "BAL"
lm_out_doseCTR$cCompartment <- factor(lm_out_doseCTR$Compartment, levels = c("Plasma", "BAL"))
lm_out_doseCTR$Timepoint <- rep("Pre-Mtb", nrow(lm_out_doseCTR))
lm_out_doseCTR$Timepoint[grepl("Post", lm_out_doseCTR$feature)] <- "Post-Mtb"
lm_out_doseCTR$Timepoint <- factor(lm_out_doseCTR$Timepoint, levels = c("Pre-Mtb", "Post-Mtb"))

# Set point labels
lm_out_doseCTR$label <- lm_out_doseCTR$feature
lm_out_doseCTR$label[lm_out_doseCTR$p_val > 0.05] <- ""
lm_out_doseCTR$label <- str_remove(lm_out_doseCTR$label, "Plasma.Pre.")
lm_out_doseCTR$label <- str_remove(lm_out_doseCTR$label, "Plasma.Post.")
lm_out_doseCTR$label <- str_remove(lm_out_doseCTR$label, "BAL.Pre.")

# Volcano plot
#png("lm_doseCTR_volcano.png", width = 5, height = 5, units = 'in', res = 300)
lm_plot_doseCTR <- lm_out_doseCTR %>%
  ggplot(aes(x = coeff, y = neglog, label = label, color = Compartment, shape = Timepoint)) +
  theme_linedraw() +
  geom_point(alpha = 0.7, size = 2) +
  xlab("Regression coefficient") +
  ylab("-log10(p-value)") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", size = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.35) +
  geom_text_repel(size = 2, color = "black", segment.size = 0.2) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = c(BAL.col, plasma.col)) +
  theme(plot.title = element_text(size = norm.size, face = "bold"), 
        axis.text = element_text(size = small.size), 
        axis.title = element_text(size = norm.size),
        legend.title = element_text(size = norm.size), 
        legend.text = element_text(size = small.size), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        aspect.ratio = 5/5) 
lm_plot_doseCTR
#dev.off()
    