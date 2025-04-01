## ANOVA for Morphology

-------------------------------------------------------------------------
# install.packages("readxl")
  
  # Load required libraries
library(readr)     # for read_csv
library(dplyr)     # for data wrangling
library(tidyr)     # for pivot_longer
library(car)       # for leveneTest (homogeneity of variance)

# Define file path

file_data <- "cleaned_Ancient_Maize_Cob_Data_Sheet.csv"


# Read the CSV into a data frame
maize_data <- read_csv(file_data)

# Select only numeric morphological measurements + Site
morphology_data <- maize_data %>%
  select(Site, Length_mm, Diameter_mm, Cupule_number, "Mean_Cupule_Width(mm)", 
         "Mean_cupule_height(mm)", Mean_cupule_depth,
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

# Print results
list(ANOVA = anova_summary, LeveneTest = levene_test, TukeyHSD = tukey_results)


# Which traits differ sig?

# Welch's ANOVA for unequal variance
oneway.test(Value ~ Trait, data = morphology_long, var.equal = FALSE)

#install.packages("rstatix")
library(rstatix)

# Run Games-Howell test
games_howell_test(morphology_long, Value ~ Trait)

library(ggplot2)
ggplot(morphology_long, aes(x = Trait, y = Value, fill = Trait)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distribution of Traits")


library(ggpubr)
ggpaired(games_howell_test(morphology_long, Value ~ Trait), 
         x = "group1", y = "estimate", 
         color = "group2") +
  theme_minimal() +
  ggtitle("Pairwise Differences in Trait Means")



-------------------------------------------------------------------------
  ## MANOVA for general differences across sites
  
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
  "Mean_cupule_depth",
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





# Post-Hoc test to look at cupule depth and KRN

# Example for Mean_cupule_depth
aov_depth <- aov(Mean_cupule_depth ~ Site, data = maize_data)
TukeyHSD(aov_depth)

# Example for Mean_kernel_row
aov_kernel <- aov(Mean_kernel_row ~ Site, data = maize_data)
TukeyHSD(aov_kernel)




# Plot this

library(ggplot2)

# For Mean_cupule_depth
ggplot(maize_data, aes(x = Site, y = Mean_cupule_depth)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 7.5) +
  theme_bw() +
  ggtitle("Cupule Depth by Site")

# For Mean_kernel_row
ggplot(maize_data, aes(x = Site, y = Mean_kernel_row)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = 32.5) +
  theme_bw() +
  ggtitle("Kernel Row by Site")




# Kind of like a PCA
# install.packages("candisc")
library(candisc)

# Suppose manova_model is your MANOVA:
# manova_model <- manova(cbind(Length_mm, Diameter_mm, ...) ~ Site, data = maize_data)

# Perform canonical discriminant analysis
cda <- candisc(manova_model)

# cda$scores contains the canonical variates
# Plot the first two canonical variates
plot(cda)

