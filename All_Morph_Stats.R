## ANOVA for ALL Morphology

-------------------------------------------------------------------------
  # install.packages("readxl")
  
  # Load required libraries
library(readr)     # for read_csv
library(dplyr)     # for data wrangling
library(tidyr)     # for pivot_longer
library(car)       # for leveneTest (homogeneity of variance)

# Define file path

file_data <- "Merged_Ancient_Maize_Cob_Data_Sheet.csv"

# Read the CSV into a data frame
maize_data <- read_csv(file_data)

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

# Print results
list(ANOVA = anova_summary, LeveneTest = levene_test, TukeyHSD = tukey_results)



# Which traits differ sig?

# Welch's ANOVA for unequal variance
oneway.test(Value ~ Trait, data = morphology_long, var.equal = FALSE)

#install.packages("rstatix")
library(rstatix)

# Run Games-Howell test
games_howell_test(morphology_long, Value ~ Trait)

print(games_howell_test(morphology_long, Value ~ Trait), n = Inf)



library(ggplot2)
ggplot(morphology_long, aes(x = Trait, y = Value, fill = Trait)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Distribution of Traits")





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






# Post-Hoc test to look at length, diameter, cupule number, 
# cupule width, cupule height, mean kernel row, weight


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

library(ggplot2)
library(ggpubr)  # for stat_compare_means

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

library(ggplot2)
library(ggpubr)  # for stat_compare_means

# Surface Area by Site
ggplot(maize_data, aes(x = Site, y = Surface_Area)) +
  geom_boxplot() +
  stat_compare_means(method = "anova", label.y = max(maize_data$Surface_Area, na.rm = TRUE)) +
  theme_bw() +
  ggtitle("Cupule Surface Area (mm^2) by Site") +
  ylab("Surface Area (mm²)")





library(dplyr)
library(ggplot2)
library(ggpubr)

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





# Kind of like a PCA
# Load required package
# install.packages("candisc") # Uncomment if not yet installed
library(candisc)

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





-------------------------------------------------------------------------
  ## FOCUS on cob size, kernel over site

  
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



-------------------------------------------------------------------------
  ## MANOVA for general differences across time periods
  
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
  select(Period, all_of(traits))

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





# Focus on Cupule_number, Mean_Cupule_Width(mm), 
# Mean_cupule_height(mm), and Diameter_mm
# Plot this

library(ggplot2)

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



# Kind of like a PCA
# Load required package
# install.packages("candisc") # Uncomment if not yet installed
library(candisc)

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








-------------------------------------------------------------------------
# Look into change over time

## FOCUS ON LENGTH


# Load necessary libraries
library(ggplot2)






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


library(ggplot2)

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

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

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


