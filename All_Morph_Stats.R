## Stats for ALL Morphology Traits




# Install Packages --------------------------------------------------------
#install.packages("viridis")

# Load required libraries
library(readr) # for read_csv
library(dplyr)     # for data wrangling
library(tidyr)     # for pivot_longer
library(car)       # for leveneTest (homogeneity of variance)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(emmeans)
library(candisc)
library(viridis)


# Define file path

file_data <- "Merged_Ancient_Maize_Cob_Data_Sheet.csv"

# Read the CSV into a data frame
maize_data <- read_csv(file_data)



# General Morph Stats -----------------------------------------------------


# Select only numeric morphological measurements + Site
morphology_data <- maize_data %>%
  select(Site, Length_mm, Diameter_mm, Cupule_number, "Mean_Cupule_Width(mm)", 
         "Mean_cupule_height(mm)",
         Mean_kernel_row, Total_Wt_g)


# Reshape the data into long format
morphology_long <- morphology_data %>%
  pivot_longer(cols = -Site, names_to = "Trait", values_to = "Value")

# Check normality assumption for ANOVA
shapiro_results <- morphology_long %>%
  group_by(Trait) %>%
  summarise(p_value = shapiro.test(Value)$p.value)

# Run one-way ANOVA to test if traits differ significantly
anova_model <- aov(Value ~ Trait, data = morphology_long)
anova_summary <- summary(anova_model)

# Check homogeneity of variances
levene_test <- leveneTest(Value ~ Trait, data = morphology_long)

# Run post hoc Tukey's test if ANOVA is significant
tukey_results <- TukeyHSD(anova_model)
list(tukey_results)

# Print results
list(ANOVA = anova_summary, LeveneTest = levene_test, TukeyHSD = tukey_results)



# Which traits differ sig?

# Welch's ANOVA for unequal variance
oneway.test(Value ~ Trait, data = morphology_long, var.equal = FALSE)

#install.packages("rstatix")


# Run Games-Howell test
games_howell_test(morphology_long, Value ~ Trait)

print(games_howell_test(morphology_long, Value ~ Trait), n = Inf)


ggplot(morphology_long, aes(x = Trait, y = Value, fill = Trait)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distribution of Traits")





# MANOVA for general differences across sites -----------------------------

  
# Make sure 'Site' is a factor
maize_data$Site <- as.factor(maize_data$Site)

# 4. Define which morphological traits you want in the MANOVA
#    Adjust as needed to match your actual column names
traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)

# 5. Subset the data to keep only Site and these traits
morph_data <- maize_data %>%
  select(Site, all_of(traits))

# 6. Create a matrix of just the numeric traits
morph_matrix <- as.matrix(morph_data[, traits])

# 7. Fit the MANOVA
#    - 'Site' is the grouping factor
#    - the morphological traits (in 'morph_matrix') are the response variables
manova_model <- manova(morph_matrix ~ Site, data = morph_data)

# 8. Summarize the MANOVA results
#    Pillai’s test is often the most robust, but you can also use "Wilks", "Hotelling-Lawley", or "Roy"
manova_summary <- summary(manova_model, test = "Pillai")
print(manova_summary)


# 9. Look at univariate ANOVAs for each trait (part of the MANOVA output)
univariate_results <- summary.aov(manova_model)
univariate_results



# # MH SITES ONLY ---------------------------------------------------------


# Filter data for Middle Horizon only
mh_data <- maize_data %>%
  filter(Period == "MH")

# Make sure 'Site' is a factor
mh_data$Site <- as.factor(mh_data$Site)

# Define morphological traits (adjust names if needed)
traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)

# Subset data to keep Site and traits
morph_data_mh <- mh_data %>%
  select(Site, all_of(traits))

# Create matrix of numeric traits
morph_matrix_mh <- as.matrix(morph_data_mh[, traits])

# Fit MANOVA: traits ~ Site (within MH period)
manova_model_mh <- manova(morph_matrix_mh ~ Site, data = morph_data_mh)

# Summarize MANOVA results using Pillai's trace
manova_summary_mh <- summary(manova_model_mh, test = "Pillai")
print(manova_summary_mh)

# Look at univariate ANOVA results for each trait
univariate_results_mh <- summary.aov(manova_model_mh)
print(univariate_results_mh)


# Load package for Tukey's HSD if not already
if (!require(multcomp)) install.packages("multcomp")
library(multcomp)

# For each significant trait, run ANOVA and then Tukey HSD post-hoc test

significant_traits <- c("Cupule_number", "Mean_Cupule_Width(mm)", "Mean_kernel_row")

posthoc_results <- list()

# for (trait in significant_traits) {
  # Fit an ANOVA model with Site as factor
  # formula_str <- as.formula(paste(trait, "~ Site"))
  # anova_model <- aov(formula_str, data = morph_data_mh)
  
  for(trait in traits){
    formula_str <- as.formula(paste0("`", trait, "` ~ Site"))
    anova_model <- aov(formula_str, data = morph_data_mh)
    tukey_res <- TukeyHSD(anova_model)
    posthoc_results[[trait]] <- summary(tukey_res)
  }


# Print post-hoc test results
posthoc_results

for (trait in traits) {
  cat("\nPost-hoc results for trait:", trait, "\n")
  print(posthoc_results[[trait]])
}



# TRAITS AMONG MH SITES
# Filter data for MH sites only
mh_sites <- c("M1-95-Chen-Chen", "M10-11-Cemetery", "M10-11-Templete", "M12", "RM-07-M43", "RM-08-M43")
maize_data_mh <- subset(maize_data, Site %in% mh_sites)


traits <- c("Length_mm", "Diameter_mm", "Cupule_number", "Mean_Cupule_Width(mm)", 
            "Mean_cupule_height(mm)", "Mean_kernel_row", "Total_Wt_g")

for (trait in traits) {
  # Wrap trait with backticks in formula to handle special characters
  formula_str <- as.formula(paste0("`", trait, "` ~ Site"))
  aov_model <- aov(formula_str, data = maize_data_mh)
  tukey_res <- TukeyHSD(aov_model)
  posthoc_results_mh[[trait]] <- tukey_res
}

# Print results
for (trait in traits) {
  cat("\nPost-hoc results for trait:", trait, "\n")
  print(posthoc_results_mh[[trait]])
}





# ALL SITES Post-Hoc test ALL TRAITS --------------------------------------

# Example for Length_mm
aov_length <- aov(Length_mm ~ Site, data = maize_data)
TukeyHSD(aov_length)

# Example for Diameter_mm
aov_diameter <- aov(Diameter_mm ~ Site, data = maize_data)
TukeyHSD(aov_diameter)

# Example for Cupule_number
aov_cupule_number <- aov(Cupule_number ~ Site, data = maize_data)
TukeyHSD(aov_cupule_number)

# Example for Mean_Cupule_Width(mm)
aov_cupule_width <- aov(`Mean_Cupule_Width(mm)` ~ Site, data = maize_data)
TukeyHSD(aov_cupule_width)

# Example for Mean_cupule_height(mm)
aov_cupule_height <- aov(`Mean_cupule_height(mm)` ~ Site, data = maize_data)
TukeyHSD(aov_cupule_height)

# Example for Mean_kernel_row
aov_kernel <- aov(Mean_kernel_row ~ Site, data = maize_data)
TukeyHSD(aov_kernel)

# Example for Total_Wt_g
aov_weight <- aov(Total_Wt_g ~ Site, data = maize_data)
TukeyHSD(aov_weight)





# Plot this
  # for stat_compare_means

# Cob Length
ggplot(maize_data, aes(x = Site, y = Length_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Length_mm, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cob Length (mm) by Site")

# Cob Diameter
ggplot(maize_data, aes(x = Site, y = Diameter_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Diameter_mm, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cob Diameter (mm) by Site")

# Cupule Number
ggplot(maize_data, aes(x = Site, y = Cupule_number)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Cupule_number, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cupule Number (circumference) by Site")

# Cupule Width
ggplot(maize_data, aes(x = Site, y = `Mean_Cupule_Width(mm)`)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$`Mean_Cupule_Width(mm)`, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Mean Cupule Width (mm) by Site")

# Cupule Height
ggplot(maize_data, aes(x = Site, y = `Mean_cupule_height(mm)`)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$`Mean_cupule_height(mm)`, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Mean Cupule Height (mm) by Site")


# cupule surface area
maize_data$Surface_Area <- maize_data$`Mean_Cupule_Width(mm)` * maize_data$`Mean_cupule_height(mm)`


# Surface Area by Site
ggplot(maize_data, aes(x = Site, y = Surface_Area)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Surface_Area, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cupule Surface Area (mm^2) by Site") +
  ylab("Surface Area (mm²)")


# Calculate counts per site
site_counts <- maize_data %>%
  group_by(Site) %>%
  summarise(n = n())

# Merge counts back into the main data for plotting
maize_data_with_counts <- maize_data %>%
  left_join(site_counts, by = "Site")

# Plot with counts
ggplot(maize_data_with_counts, aes(x = Site, y = Length_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data_with_counts$Length_mm, na.rm = TRUE)) +
  geom_text(data = site_counts, aes(x = Site, y = 0, label = paste0("n=", n)), vjust = 1.5, size = 3) +
  theme_bw() +
  ggtitle("Cob Length (mm) by Site") +
  ylab("Length (mm)")




# Mean Kernel Row Number
ggplot(maize_data, aes(x = Site, y = Mean_kernel_row)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Mean_kernel_row, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Mean Kernel Row Number by Site")

# Weight
ggplot(maize_data, aes(x = Site, y = Total_Wt_g)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Total_Wt_g, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cob Weight (g) by Site")






# Plot CDA ----------------------------------------------------------------


# Kind of like a PCA



# Fit MANOVA model with all traits
manova_model <- manova(cbind(
  Length_mm,
  Diameter_mm,
  Cupule_number,
  `Mean_Cupule_Width(mm)`,
  `Mean_cupule_height(mm)`,
  Mean_kernel_row,
  Total_Wt_g
) ~ Site, data = maize_data)

# Canonical Discriminant Analysis
cda <- candisc(manova_model)

# View canonical scores (for plotting, etc.)
head(cda$scores)

# Plot canonical variates
plot(cda)
title("CDA for Significant Traits Across Sites")







# FOCUS on cob size, kernel over site ------------------------------------

  
# Define the custom order
site_order <- c("M73", "M103", "M7", "RM-07-M43", "RM-08-M43", 
                "M10-11-Cemetery", "M10-11-Templete", "M1-95-Chen-Chen", 
                "M12", "M11", "M44")

# Set the factor levels for the Site column
maize_data$Site <- factor(maize_data$Site, levels = site_order)




# Define a function to assign time periods based on the Site
maize_data$Time_Period <- case_when(
  maize_data$Site %in% c("M73", "M103") ~ "LF",            # Late Formative
  maize_data$Site == "M7" ~ "LF/MH",                       # Late Formative / Middle Horizon
  maize_data$Site %in% c("M11", "M44") ~ "LIP",            # Late Intermediate Period
  TRUE ~ "MH"                                              # Middle Horizon
)

# Boxplot for Length with Time Period as fill
ggplot(maize_data, aes(x = Site, y = Length_mm, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("LF" = "lightblue", "LF/MH" = "orange", "LIP" = "green", "MH" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Cob Length by Site and Time Period",
       x = "Site",
       y = "Length (mm)",
       fill = "Time Period")


# Boxplot for Weight with Time Period as fill
ggplot(maize_data, aes(x = Site, y = Total_Wt_g, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("LF" = "lightblue", "LF/MH" = "orange", "LIP" = "green", "MH" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Cob Weight by Site and Time Period",
       x = "Site",
       y = "Weight (g)",
       fill = "Time Period")



# Boxplot for Diameter_mm with Time Period as fill
ggplot(maize_data, aes(x = Site, y = Diameter_mm, fill = Time_Period)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  scale_fill_manual(values = c("LF" = "lightblue", "LF/MH" = "orange", "LIP" = "green", "MH" = "purple")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of Cob Diameter by Site and Time Period",
       x = "Site",
       y = "Diameter (mm)",
       fill = "Time Period")




# MANOVA for general differences across time periods ----------------------

  
# Make sure 'Period' is a factor
maize_data$Period <- as.factor(maize_data$Period)

# 4. Define which morphological traits you want in the MANOVA
#    Adjust as needed to match your actual column names
traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)

# 5. Subset the data to keep only Site and these traits

morph_data <- maize_data %>%
  dplyr::select(Period, all_of(traits))



# 6. Create a matrix of just the numeric traits
morph_matrix <- as.matrix(morph_data[, traits])

# 7. Fit the MANOVA
#    - 'Site' is the grouping factor
#    - the morphological traits (in 'morph_matrix') are the response variables
manova_model <- manova(morph_matrix ~ Period, data = morph_data)

# 8. Summarize the MANOVA results
#    Pillai’s test is often the most robust, but you can also use "Wilks", "Hotelling-Lawley", or "Roy"
manova_summary <- summary(manova_model, test = "Pillai")
print(manova_summary)


# 9. Look at univariate ANOVAs for each trait (part of the MANOVA output)
univariate_results <- summary.aov(manova_model)
univariate_results



# Post-Hoc test to look at length, diameter, cupule number, 
# cupule width, cupule height, kernel row number, and weight by Period

# Example for Length_mm
aov_length <- aov(Length_mm ~ Period, data = maize_data)
TukeyHSD(aov_length)

# Example for Diameter_mm
aov_diameter <- aov(Diameter_mm ~ Period, data = maize_data)
TukeyHSD(aov_diameter)

# Example for Cupule_number
aov_cupule_number <- aov(Cupule_number ~ Period, data = maize_data)
TukeyHSD(aov_cupule_number)

# Example for Mean_Cupule_Width(mm)
aov_cupule_width <- aov(`Mean_Cupule_Width(mm)` ~ Period, data = maize_data)
TukeyHSD(aov_cupule_width)

# Example for Mean_cupule_height(mm)
aov_cupule_height <- aov(`Mean_cupule_height(mm)` ~ Period, data = maize_data)
TukeyHSD(aov_cupule_height)

# Example for Mean_kernel_row
aov_kernel <- aov(Mean_kernel_row ~ Period, data = maize_data)
TukeyHSD(aov_kernel)

# Example for Total_Wt_g
aov_weight <- aov(Total_Wt_g ~ Period, data = maize_data)
TukeyHSD(aov_weight)




# Between Time Periods ----------------------------------------------------


# Get estimated marginal means for Period
pairwise_results <- emmeans(manova_model, ~ Period)

# Late Formative (LF) to Middle Horizon (MH)
contrast_LF_MH <- contrast(pairwise_results, 
                           method = "trt.vs.ctrl", 
                           ref = "LF")  # LF is the reference
print(contrast_LF_MH)

# Within Middle Horizon period — since you have MH and LF/MH (transitional)
# You can test pairwise differences including LF/MH and MH
contrast_within_MH <- contrast(pairwise_results, 
                               method = "pairwise", 
                               adjust = "BH")
print(contrast_within_MH)

# Focus on MH related groups only
mh_periods <- c("LF", "LF/MH", "MH")

pairwise_mh <- emmeans(manova_model, pairwise ~ Period, 
                       adjust = "BH")$contrasts %>%
  as.data.frame() %>%
  filter(contrast %in% c("LF - (LF/MH)", "(LF/MH) - MH", "LF - MH"))

print(pairwise_mh)

# Middle Horizon (MH) to Late Intermediate Period (LIP)
contrast_MH_LIP <- contrast(pairwise_results, 
                            method = "trt.vs.ctrl", 
                            ref = "MH")  # MH is the reference
print(contrast_MH_LIP)







# Focus on Cupule_number, Mean_Cupule_Width(mm), Mean_cupule_heigh --------


# Plot this

# For Cupule_number
ggplot(maize_data, aes(x = Period, y = Cupule_number)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 25.5) +
  theme_bw() +
  ggtitle("Cupule Number (circumference) by Time")


# For `Mean_Cupule_Width(mm)`
ggplot(maize_data, aes(x = Period, y = `Mean_Cupule_Width(mm)`)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 8.5) +
  theme_bw() +
  ggtitle("Mean_Cupule_Width(mm) by Time")

# For `Length(mm)`
ggplot(maize_data, aes(x = Period, y = Length_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 125.5) +
  theme_bw() +
  ggtitle("Length (mm) by Time")


# For Mean_cupule_height(mm)
ggplot(maize_data, aes(x = Period, y = `Mean_cupule_height(mm)`)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 5) +
  theme_bw() +
  ggtitle("Mean_cupule_height(mm) by Time")



# Create a new column for Surface Area (Width * Height)
maize_data$Surface_Area <- maize_data$`Mean_Cupule_Width(mm)` * maize_data$`Mean_cupule_height(mm)`

# Plot for Surface Area (Width * Height)
ggplot(maize_data, aes(x = Period, y = Surface_Area)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 35.5) +
  theme_bw() +
  ggtitle("Cupule Surface Area (Width x Height) by Time")



# For Diameter_mm
ggplot(maize_data, aes(x = Period, y = Diameter_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 75) +
  theme_bw() +
  ggtitle("Diameter_mm by Time")



# For Diameter_mm
ggplot(maize_data, aes(x = Period, y = Diameter_mm)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 75) +
  theme_bw() +
  ggtitle("Diameter_mm by Time")




# CDA Plot ----------------------------------------------------------------


# Fit MANOVA model with all traits
manova_model <- manova(cbind(
  Cupule_number,
  `Mean_Cupule_Width(mm)`,
  `Mean_cupule_height(mm)`,
  Diameter_mm
) ~ Period, data = maize_data)

# Canonical Discriminant Analysis
cda <- candisc(manova_model)

# View canonical scores (for plotting, etc.)
head(cda$scores)

# Plot canonical variates
plot(cda)
title("CDA for Significant Traits Across Time")



# Fit MANOVA model with the new Surface_Area variable
manova_model <- manova(cbind(
  Cupule_number,
  Cupule_Surface_Area,  # Using the Surface Area instead of Mean_Cupule_Width(mm) and Mean_cupule_height(mm)
  Diameter_mm
) ~ Period, data = maize_data)

# Canonical Discriminant Analysis
cda <- candisc(manova_model)

# View canonical scores (for plotting, etc.)
head(cda$scores)

# Plot canonical variates
plot(cda)
title("CDA for Significant Traits Across Time")




# Look into change over time w Individual Traits --------------------------


## FOCUS ON LENGTH


# Reorder the 'Period' factor levels
maize_data$Period <- factor(maize_data$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Create the boxplot with the ordered periods
ggplot(maize_data, aes(x = Period, y = Length_mm, fill = Period)) +
  geom_boxplot() +
  labs(title = "Distribution of Maize Cob Length Over Time Periods",
       x = "Time Period",
       y = "Length (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Create the histogram with the ordered periods
ggplot(maize_data, aes(x = Length_mm, fill = Period)) +
  geom_histogram(position = "dodge", binwidth = 2, color = "black") +
  facet_wrap(~ Period, scales = "free_y") +  # Create separate histograms for each period
  labs(title = "Histogram of Maize Cob Length Over Time Periods",
       x = "Length (mm)",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## FOCUS ON DIAMETER


# Reorder the 'Period' factor levels
maize_data$Period <- factor(maize_data$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Create the boxplot with the ordered periods
ggplot(maize_data, aes(x = Period, y = Diameter_mm, fill = Period)) +
  geom_boxplot() +
  labs(title = "Distribution of Maize Cob Diameter Over Time Periods",
       x = "Time Period",
       y = "Diameter (mm)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# FOCUS ON WEIGHT


ggplot(maize_data, aes(x = Period, y = Total_Wt_g, fill = Period)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2, fill = "black") +
  scale_fill_viridis_d() +
  labs(
    title = "Distribution of Cob Weight Over Time",
    x = "Period",
    y = "Cob Weight (g)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")





# FOCUS ON CUPULE CHARACTERISTICS

# Reorder Period levels
maize_data$Period <- factor(maize_data$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Select and pivot the relevant columns
maize_long <- maize_data %>%
  select(Period, 
         Cupule_number,
         `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`,
         Mean_kernel_row) %>%
  pivot_longer(cols = -Period,
               names_to = "Trait",
               values_to = "Value")

# Optional: Clean trait labels for prettier plotting
maize_long$Trait <- recode(maize_long$Trait,
                           "Cupule_number" = "Cupule Number",
                           "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
                           "Mean_cupule_height(mm)" = "Cupule Height (mm)",
                           "Mean_kernel_row" = "Kernel Rows")

# Make faceted boxplots
ggplot(maize_long, aes(x = Period, y = Value, fill = Period)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Time Period", y = "Value", title = "Distribution of Morphological Traits Over Time") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))




# PCA ---------------------------------------------------------------------


# By Site
#Step 1: Select and rename numeric columns (optional but cleaner)
morpho_wide <- maize_data %>%
  dplyr::select(Site,
                Length = Length_mm,
                Diameter = Diameter_mm,
                CupuleNumber = Cupule_number,
                CupuleWidth = `Mean_Cupule_Width(mm)`,
                CupuleHeight = `Mean_cupule_height(mm)`,
                KernelRows = Mean_kernel_row,
                Weight = Total_Wt_g) %>%
    drop_na()  # PCA can't handle NA values
  

# Step 2: Run PCA on just the numeric columns
morpho_pca <- prcomp(dplyr::select(morpho_wide, -Site), scale. = TRUE)


# Step 3: View PCA summary
summary(morpho_pca)

# Step 4: Create a scores data frame and add Site info back
scores <- as.data.frame(morpho_pca$x)
scores$Site <- morpho_wide$Site

# Step 5 (optional): Plot the first two principal components

ggplot(scores, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Ancient Maize Cob Morphology by Site",
       x = paste0("PC1 (", round(summary(morpho_pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(morpho_pca)$importance[2,2]*100, 1), "%)"))



# By Trait

# Step 1: Select and rename numeric traits only
trait_matrix <- maize_data %>%
  dplyr::select(Length = Length_mm,
         Diameter = Diameter_mm,
         CupuleNumber = Cupule_number,
         CupuleWidth = `Mean_Cupule_Width(mm)`,
         CupuleHeight = `Mean_cupule_height(mm)`,
         KernelRows = Mean_kernel_row,
         Weight = Total_Wt_g) %>%
  drop_na()

# Step 2: Transpose the data so traits are rows and cobs are columns
trait_matrix_t <- t(trait_matrix)

# Step 3: Run PCA (transpose makes traits the "observations")
pca_traits <- prcomp(trait_matrix_t, scale. = TRUE)

# Step 4: View PCA loadings (how traits relate to each PC)
summary(pca_traits)
pca_traits$rotation  # Shows how cobs contribute to each trait-PC

# Step 5: Plot the traits in PC space

trait_scores <- as.data.frame(pca_traits$x)
trait_scores$Trait <- rownames(trait_scores)


ggplot(trait_scores, aes(x = PC1, y = PC2, label = Trait, color = Trait)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, size = 5) +
  scale_color_viridis_d() +   # discrete viridis colors
  theme_minimal() +
  labs(title = "PCA of Traits Across Maize Cobs",
       x = paste0("PC1 (", round(summary(pca_traits)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_traits)$importance[2,2]*100, 1), "%)"),
       color = "Trait")  # legend title


# IMPROVED PCA for TRAITS
# Step 1: Select and rename numeric traits only
trait_matrix <- maize_data %>%
  dplyr::select(Length = Length_mm,
         Diameter = Diameter_mm,
         CupuleNumber = Cupule_number,
         CupuleWidth = `Mean_Cupule_Width(mm)`,
         CupuleHeight = `Mean_cupule_height(mm)`,
         KernelRows = Mean_kernel_row,
         Weight = Total_Wt_g) %>%
  drop_na()

# Step 2: Transpose the data so traits are rows and cobs are columns
trait_matrix_t <- t(trait_matrix)

# Step 3: Run PCA (transpose makes traits the "observations")
pca_traits <- prcomp(trait_matrix_t, scale. = TRUE)

# Step 4: View PCA loadings (how traits relate to each PC)
summary(pca_traits)
pca_traits$rotation  # Shows how cobs contribute to each trait-PC

# Step 5: Plot the traits in PC space with distinct colors for each trait

trait_scores <- as.data.frame(pca_traits$x)
trait_scores$Trait <- rownames(trait_scores)


ggplot(trait_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Trait), size = 4) +
  scale_color_viridis_d() +
  geom_text(aes(label = Trait), vjust = -0.8, size = 5) +
  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "PCA of Traits Across Maize Cobs",
       x = paste0("PC1 (", round(summary(pca_traits)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_traits)$importance[2,2]*100, 1), "%)")
)

# To view proportion of variance
summary(pca_traits)





# PCA for PC1 Period ------------------------------------------------------

# Step 1: Drop rows with NA in any trait and keep index
complete_rows <- maize_data %>%
  dplyr::select(Length = Length_mm,
         Diameter = Diameter_mm,
         CupuleNumber = Cupule_number,
         CupuleWidth = `Mean_Cupule_Width(mm)`,
         CupuleHeight = `Mean_cupule_height(mm)`,
         KernelRows = Mean_kernel_row,
         Weight = Total_Wt_g) %>%
  drop_na()

# Step 2: Subset maize_data to match only complete rows
# This assumes row order has not been altered!
maize_metadata <- maize_data[complete.cases(
  maize_data[, c("Length_mm", "Diameter_mm", "Cupule_number", 
                 "Mean_Cupule_Width(mm)", "Mean_cupule_height(mm)", 
                 "Mean_kernel_row", "Total_Wt_g")]
), ]

# Step 3: Run PCA
pca_cobs <- prcomp(complete_rows, scale. = TRUE)

# Step 4: Create score data frame and merge metadata
pc_scores <- as.data.frame(pca_cobs$x)
pc_scores$Period <- maize_metadata$Period
pc_scores$Site <- maize_metadata$Site
pc_scores$Sample_ID <- maize_metadata$Sample_ID

# Step 5: Plot PC1 by Period

# Relevel 'Period' to desired order: LF, LF/MH, LIP, MH
pc_scores$Period <- factor(pc_scores$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Now re-plot with the updated factor levels

ggplot(pc_scores, aes(x = Period, y = PC1, fill = Period)) +
  geom_boxplot(outlier.color = NA) +       # boxplots with no outlier points
  geom_jitter(width = 0.2, alpha = 0.4, color = "black") +  # black jitter points
  labs(title = "PC1 (Cob Size) Across Time Periods",
       y = "PC1 Score", x = "Period") +
  theme_minimal() +
  theme(panel.grid = element_blank())

ggplot(pc_scores, aes(x = Period, y = PC1, fill = Period)) +
  geom_boxplot(outlier.color = NA) +       
  geom_jitter(width = 0.2, alpha = 0.4, color = "black") +  
  labs(title = "PC1 (Cob Size) Across Time Periods",
       y = "PC1 Score", x = "Period") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha('white', 0.6), color = NA))



## Focus on PC1 for Site
# Reorder Site levels
pc_scores$Site <- factor(pc_scores$Site, levels = c(
  "M73", "M103", "M7", "RM-07-M43", "RM-08-M43", 
  "M10-11-Cemetery", "M10-11-Templete", "M1-95-Chen-Chen", 
  "M12", "M11", "M44"
))

# Plot PC1 across Sites
ggplot(pc_scores, aes(x = Site, y = PC1, fill = Period)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  labs(title = "PC1 (Cob Size) Across Archaeological Sites",
       x = "Site", y = "PC1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")




# With sample size

# Count number of samples per site
site_counts <- pc_scores %>%
  group_by(Site) %>%
  summarise(n = n())

# Merge counts into main PC score data
pc_scores_labeled <- pc_scores %>%
  left_join(site_counts, by = "Site") %>%
  mutate(Site_label = paste0(Site, "\n(n=", n, ")"))

# Reorder Site_label using your preferred order
site_order <- c("M73", "M103", "M7", "RM-07-M43", "RM-08-M43", 
                "M10-11-Cemetery", "M10-11-Templete", "M1-95-Chen-Chen", 
                "M12", "M11", "M44")

# Match original Site order to the new Site_label
pc_scores_labeled$Site_label <- factor(pc_scores_labeled$Site_label,
                                       levels = pc_scores_labeled %>%
                                         filter(Site %in% site_order) %>%
                                         distinct(Site, Site_label) %>%
                                         slice(match(site_order, Site)) %>%
                                         pull(Site_label)
)

# Plot
ggplot(pc_scores_labeled, aes(x = Site_label, y = PC1, fill = Period)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  labs(title = "PC1 (Cob Size) Across Archaeological Sites",
       x = "Site", y = "PC1 Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")






# Chen chen v Omo style ---------------------------------------------------


  # Define Chen Chen and Omo style groups
  pc_scores_labeled <- pc_scores_labeled %>%
  mutate(Style_Group = case_when(
    Site %in% c("M1-95-Chen-Chen", "RM-07-M43", "RM-08-M43") ~ "Chen Chen",
    Site == "M12" ~ "Omo",
    TRUE ~ NA_character_
  ))

# Drop rows with NA in Style_Group
pc_scores_grouped <- pc_scores_labeled %>%
  filter(!is.na(Style_Group))


library(ggplot2)
library(dplyr)

# QQ Plots for each group
ggplot(pc_scores_grouped, aes(sample = PC1, color = Style_Group)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Style_Group) +
  theme_minimal() +
  labs(title = "QQ Plots for Chen Chen and Omo Style Groups")

# Test for normality in each group
shapiro_test <- pc_scores_grouped %>%
  group_by(Style_Group) %>%
  summarise(p_value = shapiro.test(PC1)$p.value)

print(shapiro_test)

# Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(PC1 ~ Style_Group, data = pc_scores_grouped)

print(wilcox_test_result)


# The Wilcoxon rank-sum test result shows:

#W = 613, which is the test statistic.

#p-value = 0.3807, which is greater than 0.05.

#Since the p-value is greater than 0.05, we fail to reject the null hypothesis, 
#meaning there is no statistically significant difference in PC1 (cob size) 
# scores between the Chen Chen style sites and the Omo style site.





# Chen chen v Omo style v M10 Temple
# Define Chen Chen, Omo, and Ritual Complex style groups
pc_scores_labeled <- pc_scores_labeled %>%
  mutate(Style_Group = case_when(
    Site %in% c("M1-95-Chen-Chen", "RM-07-M43", "RM-08-M43") ~ "Chen Chen",
    Site == "M12" ~ "Omo",
    Site %in% c("M10-11-Cemetery", "M10-11-Templete") ~ "Ritual Complex",  # Add Ritual Complex group
    TRUE ~ NA_character_
  ))

# Drop rows with NA in Style_Group
pc_scores_grouped <- pc_scores_labeled %>%
  filter(!is.na(Style_Group))


# QQ Plots for each group
ggplot(pc_scores_grouped, aes(sample = PC1, color = Style_Group)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Style_Group) +
  theme_minimal() +
  labs(title = "QQ Plots for Chen Chen, Omo, and Ritual Complex Style Groups")

# Test for normality in each group
shapiro_test <- pc_scores_grouped %>%
  group_by(Style_Group) %>%
  summarise(p_value = shapiro.test(PC1)$p.value)

print(shapiro_test)

# Kruskal-Wallis test for non-normal data
kruskal_test_result <- kruskal.test(PC1 ~ Style_Group, data = pc_scores_grouped)
print(kruskal_test_result)

# Pairwise Wilcoxon tests with p-value adjustment
pairwise_wilcoxon <- pairwise.wilcox.test(pc_scores_grouped$PC1, pc_scores_grouped$Style_Group, p.adjust.method = "BH")
print(pairwise_wilcoxon)

# No sig difference between groups based on PC1