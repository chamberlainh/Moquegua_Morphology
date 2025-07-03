## Stats for ALL Morphology Traits




# Install Packages --------------------------------------------------------


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
library(rstatix)
library(broom)
library(writexl)
library(RColorBrewer)
library(multcomp)



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
  dplyr::select(Site, all_of(traits))

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

# ANALYSIS INDICATES ALL TRAITS EXCEPT WEIGHT SIGNIFICIANT



# 9. Look at univariate ANOVAs for each trait (part of the MANOVA output)
univariate_results <- summary.aov(manova_model)
univariate_results


# SAVE general trait results

# Welch's ANOVA
welch_anova <- oneway.test(Value ~ Trait, data = morphology_long, var.equal = FALSE)
welch_df <- data.frame(
  Df = welch_anova$parameter[1],
  F = welch_anova$statistic,
  p_value = welch_anova$p.value,
  Method = "Welch's ANOVA"
)

# Games-Howell test
games_df <- morphology_long %>%
  games_howell_test(Value ~ Trait)

# MANOVA summary
# Extract the test results (coefficients) from the MANOVA summary
manova_df <- as.data.frame(manova_summary$stats)

# Rename columns for clarity
colnames(manova_df) <- c("Df", "Pillai", "approx_F", "num_Df", "den_Df", "Pr(>F)")

# Add rownames as a column
manova_df$Effect <- rownames(manova_df)

# Univariate ANOVAs
univariate_list <- summary.aov(manova_model)
# Convert each univariate result to a tidy table and combine
univariate_df <- do.call(rbind, lapply(names(univariate_list), function(trait) {
  broom::tidy(univariate_list[[trait]]) %>%
    mutate(Trait = trait)
}))

# Summary Stats
# Start fresh from full maize_data
morphology_data <- maize_data %>%
  dplyr::select(Site, Length_mm, Diameter_mm, Cupule_number,
         `Mean_Cupule_Width(mm)`, `Mean_cupule_height(mm)`,
         Mean_kernel_row, Total_Wt_g)

# Then pivot
morphology_long <- morphology_data %>%
  pivot_longer(cols = -Site, names_to = "Trait", values_to = "Value")


morph_summary_stats <- morphology_long %>%
  group_by(Trait) %>%
  summarise(
    N = sum(!is.na(Value)),
    Missing = sum(is.na(Value)),
    Mean = mean(Value, na.rm = TRUE),
    SD = sd(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE)
  ) %>%
  ungroup()




# Save
write_xlsx(
  list(
    "Welch_ANOVA" = welch_df,
    "Games_Howell" = games_df,
    "MANOVA_Pillai" = manova_df,
    "Univariate_ANOVAs" = univariate_df,
    "Summary_Stats" = morph_summary_stats
  ),
  path = "Supplementary_Morphology_Stats.xlsx"
)



# M73 and M12 from CDA -----------------------------


## Noticed in CDA M73 (cupule number) and M12 (cob dimensions) seem unique
maize_data$Site <- as.factor(maize_data$Site)
mod_cupule <- aov(Cupule_number ~ Site, data = maize_data)
tukey_cupule <- TukeyHSD(mod_cupule)

# View M73 comparisons
tukey_cupule$Site[grep("M73", rownames(tukey_cupule$Site)), ]
mod_weight <- aov(Total_Wt_g ~ Site, data = maize_data)
tukey_weight <- TukeyHSD(mod_weight)

# View M12 comparisons
tukey_weight$Site[grep("M12", rownames(tukey_weight$Site)), ]




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
  dplyr::select(Site, all_of(traits))

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

significant_traits <- c("Cupule_number", "Mean_Cupule_Width(mm)", "Mean_kernel_row", "Mean_cupule_height(mm)")

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
# Calculate and store Tukey HSD post-hoc tests
posthoc_results <- list()

for (trait in significant_traits) {
  cat("\nPost-hoc results for trait:", trait, "\n")
  
  # Use backticks around variable names to handle special characters
  formula_str <- paste0("`", trait, "` ~ Site")
  aov_model <- aov(as.formula(formula_str), data = morph_data_mh)
  
  # Tukey HSD
  tukey <- TukeyHSD(aov_model)
  
  # Store result
  posthoc_results[[trait]] <- tukey
  
  # Print result
  print(tukey)
}


# # LIP SITES ONLY ---------------------------------------------------------


# Filter data for LIP only
lip_data <- maize_data %>%
  filter(Period == "LIP")

# Make sure 'Site' is a factor
lip_data$Site <- as.factor(lip_data$Site)

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
morph_data_lip <- lip_data %>%
  dplyr::select(Site, all_of(traits))

# Create matrix of numeric traits
morph_matrix_lip <- as.matrix(morph_data_lip[, traits])

# Fit MANOVA: traits ~ Site (within LIP period)
manova_model_lip <- manova(morph_matrix_lip ~ Site, data = morph_data_lip)

# Summarize MANOVA results using Pillai's trace
manova_summary_lip <- summary(manova_model_lip, test = "Pillai")
print(manova_summary_lip)

# Look at univariate ANOVA results for each trait
univariate_results_lip <- summary.aov(manova_model_lip)
print(univariate_results_lip)


# For significant trait total weight, run ANOVA and then Tukey HSD post-hoc test

significant_traits <- c("Total_Wt_g")

posthoc_results <- list()

# for (trait in significant_traits) {
# Fit an ANOVA model with Site as factor
# formula_str <- as.formula(paste(trait, "~ Site"))
# anova_model <- aov(formula_str, data = morph_data_mh)

for(trait in traits){
  formula_str <- as.formula(paste0("`", trait, "` ~ Site"))
  anova_model <- aov(formula_str, data = morph_data_lip)
  tukey_res <- TukeyHSD(anova_model)
  posthoc_results[[trait]] <- summary(tukey_res)
}


# Print post-hoc test results
# Calculate and store Tukey HSD post-hoc tests
posthoc_results <- list()

for (trait in significant_traits) {
  cat("\nPost-hoc results for trait:", trait, "\n")
  
  # Use backticks around variable names to handle special characters
  formula_str <- paste0("`", trait, "` ~ Site")
  aov_model <- aov(as.formula(formula_str), data = morph_data_lip)
  
  # Tukey HSD
  tukey <- TukeyHSD(aov_model)
  
  # Store result
  posthoc_results[[trait]] <- tukey
  
  # Print result
  print(tukey)
}


# SAVE AS A SPREADSHEET
# Make MANOVA df for lip
manova_df_lip <- as.data.frame(manova_summary_mh$stats)
colnames(manova_df_lip) <- c("Df", "Pillai", "approx_F", "num_Df", "den_Df", "Pr(>F)")
manova_df_lip$Effect <- rownames(manova_df_lip)
manova_df_lip <- manova_df_lip[, c("Effect", "Df", "Pillai", "approx_F", "num_Df", "den_Df", "Pr(>F)")]

# Make univariate ANOVA df for lip
univariate_df_lip <- do.call(rbind, lapply(names(univariate_results_lip), function(trait) {
  broom::tidy(univariate_results_lip[[trait]]) %>%
    mutate(Trait = trait)
}))

# Make Tukey from ANOVA df for MH

posthoc_results <- list()

for (trait in traits) {
  formula_str <- as.formula(paste0("`", trait, "` ~ Site"))
  anova_model <- aov(formula_str, data = morph_data_lip)
  tukey_res <- TukeyHSD(anova_model)
  posthoc_results[[trait]] <- tukey_res  # store raw Tukey result
}

tukey_lip_df <- do.call(rbind, lapply(names(posthoc_results), function(trait) {
  # Extract just the site comparison table
  tukey_result <- posthoc_results[[trait]]$Site
  df <- as.data.frame(tukey_result)
  df$Comparison <- rownames(tukey_result)
  df$Trait <- trait
  df
})) %>%
  dplyr::select(Trait, Comparison, diff, lwr, upr, `p adj`)

write_xlsx(
  list(
    "MANOVA_LIP" = manova_df_lip,
    "ANOVA_LIP" = univariate_df_lip,
    "TukeyHSD_LIP" = tukey_lip_df
  ),
  path = "Supplementary_LIP_Traits_with_AllSites.xlsx"
)



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



# SAVE AS A SPREADSHEET
# Make MANOVA df for MH
manova_df_mh <- as.data.frame(manova_summary_mh$stats)
colnames(manova_df_mh) <- c("Df", "Pillai", "approx_F", "num_Df", "den_Df", "Pr(>F)")
manova_df_mh$Effect <- rownames(manova_df_mh)
manova_df_mh <- manova_df_mh[, c("Effect", "Df", "Pillai", "approx_F", "num_Df", "den_Df", "Pr(>F)")]

# Make univariate ANOVA df for MH
univariate_df_mh <- do.call(rbind, lapply(names(univariate_results_mh), function(trait) {
  broom::tidy(univariate_results_mh[[trait]]) %>%
    mutate(Trait = trait)
}))

# Make Tukey from ANOVA df for MH

posthoc_results <- list()

for (trait in traits) {
  formula_str <- as.formula(paste0("`", trait, "` ~ Site"))
  anova_model <- aov(formula_str, data = morph_data_mh)
  tukey_res <- TukeyHSD(anova_model)
  posthoc_results[[trait]] <- tukey_res  # store raw Tukey result
}

tukey_mh_df <- do.call(rbind, lapply(names(posthoc_results), function(trait) {
  # Extract just the site comparison table
  tukey_result <- posthoc_results[[trait]]$Site
  df <- as.data.frame(tukey_result)
  df$Comparison <- rownames(tukey_result)
  df$Trait <- trait
  df
})) %>%
  dplyr::select(Trait, Comparison, diff, lwr, upr, `p adj`)


# Make Tukey for All Traits df for MH
traits <- c(
  "Length_mm", "Diameter_mm", "Cupule_number", "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)", "Mean_kernel_row", "Total_Wt_g"
)


# Combine all
tukey_all_df <- do.call(rbind, tukey_all_sites) %>%
  dplyr::select(Trait, Comparison, diff, lwr, upr, `p adj`)

write_xlsx(
  list(
    "MANOVA_MH" = manova_df_mh,
    "ANOVA_MH" = univariate_df_mh,
    "TukeyHSD_MH" = tukey_mh_df
  ),
  path = "Supplementary_MH_Traits_with_AllSites.xlsx"
)




# Plot MH Traits Across Sites NEED TO FIX

# 1. Define label replacements
trait_labels <- c(
  "Cupule_number" = "Cupule Number",
  "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
  "Mean_cupule_height(mm)" = "Cupule Height (mm)",
  "Mean_kernel_row" = "Kernel Rows",
  "Length_mm" = "Length (mm)",
  "Diameter_mm" = "Diameter (mm)",
  "Total_Wt_g" = "Total Weight (g)"
)

# 2. Reshape data
mh_long <- morph_data_mh %>%
  pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") %>%
  mutate(
    Site = factor(Site),
    Trait = dplyr::recode(Trait, !!!trait_labels)
  )

# 3. Prepare p-values
pval_df <- tukey_mh_df %>%
  separate(Comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(Trait = dplyr::recode(Trait, !!!trait_labels)) %>%
  rename(p = `p adj`) %>%
  left_join(
    mh_long %>%
      group_by(Trait) %>%
      summarise(y.position = max(Value, na.rm = TRUE) * 1.1),
    by = "Trait"
  )

pval_df <- pval_df %>%
  mutate(p = formatC(p, format = "e", digits = 2))  # scientific notation, or use `round(p, 3)

# Plot with P's
mh_long <- mh_long %>%
  mutate(Trait = factor(
    Trait,
    levels = c(
      "Diameter (mm)",
      "Length (mm)",
      "Total Weight (g)",
      "Cupule Height (mm)",
      "Cupule Number",
      "Cupule Width (mm)",
      "Kernel Rows"
    )
  ))


ggplot(mh_long, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +  # Adds extra headroom for p-values
  facet_wrap(~ Trait, scales = "free_y") +
  stat_compare_means(method = "anova", label = "p.format", size = 3) +  # Smaller p-value text
  labs(
    title = "Morphological Trait Variation Across MH Sites",
    x = "Site", y = "Trait Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )


# Split MH Plots into Cobs and Cupule Characteristics
# Define trait groupings
cob_traits <- c("Length (mm)", "Diameter (mm)", "Total Weight (g)")
cupule_traits <- c("Cupule Number", "Cupule Width (mm)", "Cupule Height (mm)", "Kernel Rows")

ggplot(filter(mh_long, Trait %in% cob_traits), aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  facet_wrap(~ Trait, scales = "free_y") +
  stat_compare_means(method = "anova", label = "p.format", size = 3) +
  labs(
    title = "Cob Dimensions Across MH Sites",
    x = "Site", y = "Trait Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )

ggplot(filter(mh_long, Trait %in% cupule_traits), aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  facet_wrap(~ Trait, scales = "free_y") +
  stat_compare_means(method = "anova", label = "p.format", size = 3) +
  labs(
    title = "Cupule Characteristics Across MH Sites",
    x = "Site", y = "Trait Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "none"
  )







# Plot individual anovas for stat_compare_means

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




# Look at M7 for LF/MH Wari transition ----------------------------------------------------------------

maize_data <- maize_data %>%
  filter(Site != "M44") %>%  # Exclude M44
  mutate(
    Cultural_Period = case_when(
      Site %in% c("M73", "M103") ~ "Early_Huaracane",
      Site == "M7" ~ "Transitional_Wari",
      TRUE ~ "Later_Tiwanaku"
    )
  )


# Ensure all traits are numeric and clean
traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)

morph_data <- maize_data %>%
  dplyr::select(Cultural_Period, all_of(traits)) %>%
  filter(complete.cases(.))  # Remove rows with NA

# MANOVA model
manova_model <- manova(as.matrix(morph_data[, traits]) ~ Cultural_Period, data = morph_data)

# Summary
summary(manova_model, test = "Wilks")


results <- lapply(traits, function(trait) {
  formula <- as.formula(paste0("`", trait, "` ~ Cultural_Period"))
  aov_model <- aov(formula, data = morph_data)
  summary(aov_model)
})
names(results) <- traits

results

tukey_tests <- lapply(traits, function(trait) {
  model <- aov(as.formula(paste0("`", trait, "` ~ Cultural_Period")), data = morph_data)
  TukeyHSD(model)
})
names(tukey_tests) <- traits



# Example: Cupule number comparison
#tukey_tests$Cupule_number$Cultural_Period[grep("Transitional_Wari", rownames(tukey_tests$Cupule_number$Cultural_Period)), ]

# Plot
# ggplot(morph_data, aes(x = Cultural_Period, y = Cupule_number, fill = Cultural_Period)) +
#  geom_boxplot() +
#  theme_minimal() +
#  labs(title = "Cupule Number by Cultural Period")

# ALL TRAITS for Wari

m7_tukey_results <- lapply(traits, function(trait) {
  model <- aov(as.formula(paste0("`", trait, "` ~ Cultural_Period")), data = morph_data)
  tukey <- TukeyHSD(model)
  df <- as.data.frame(tukey$Cultural_Period)
  df$Comparison <- rownames(df)
  df$Trait <- trait
  df %>% filter(grepl("Transitional_Wari", Comparison))
})

# Combine all results into one data frame
m7_all <- do.call(rbind, m7_tukey_results)

# Optional: tidy the column order
m7_all <- m7_all %>%
  dplyr::select(Trait, Comparison, diff, lwr, upr, `p adj`) %>%
  arrange(Trait)

# View it
View(m7_all)  # or just print(m7_all) for console

# SAVE
write_xlsx(
  list(
    "TukeyHSD_M7_Transition" = m7_all
  ),
  path = "TukeyHSD_M7_Transition.xlsx"
)


# Plot CDA ----------------------------------------------------------------


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


# Plot canonical variates with custom axis limits
plot(cda,
     main = "CDA for Traits Across Sites",
     xlim = c(-8, 14),
     ylim = c(-8, 8),
     cex = 0.3  # Smaller text for trait labels
)

# Add legend
legend(x = 15, y = 0,  # Adjust these values based on your plot
       legend = levels(maize_data_no_outlier$Site),
       pch = 1:length(levels(maize_data_no_outlier$Site)),
       col = 1:length(levels(maize_data_no_outlier$Site)),
       title = "Site", 
       cex = 0.7, 
       bty = "n")



# TRY WITHOUT WEIGHT BC Not Stat Sig

# Fit MANOVA model with all traits
maize_data_no_outlier <- maize_data[-185, ]


manova_model <- manova(cbind(
  Length_mm,
  Diameter_mm,
  Cupule_number,
  `Mean_Cupule_Width(mm)`,
  `Mean_cupule_height(mm)`,
  Mean_kernel_row
) ~ Site, data = maize_data_no_outlier)

# Canonical Discriminant Analysis
cda <- candisc(manova_model)

# View canonical scores (for plotting, etc.)
head(cda$scores)
summary(scores$Can2)

# Where is the M10 outlier?
outlier_row <- scores[which.min(scores$Can2), ]
print(outlier_row)


# Plot canonical variates with custom axis limits
plot(cda,
     main = "CDA for Significant Traits Across Sites (No Weight)",
     xlim = c(-8, 14),
     ylim = c(-8, 8)
)

# Add legend
legend(
  "bottomright",
  legend = levels(maize_data_no_outlier$Site),
  pch = 1:length(levels(maize_data_no_outlier$Site)),
  col = 1:length(levels(maize_data_no_outlier$Site)),
  title = "Site", 
  cex = 0.6, 
  bty = "n"
)


# Check M12 and M73 again

# M73 and M12 from CDA -----------------------------


## Noticed in CDA M73 (cupule number) and M12 (cob dimensions) seem unique
maize_data$Site <- as.factor(maize_data$Site)
mod_cupule <- aov(Cupule_number ~ Site, data = maize_data)
tukey_cupule <- TukeyHSD(mod_cupule)

# View M73 comparisons SIG
tukey_cupule$Site[grep("M73", rownames(tukey_cupule$Site)), ]

# M12 NOT SIG
mod_KRN <- aov(Mean_kernel_row ~ Site, data = maize_data)
tukey_KRN <- TukeyHSD(mod_KRN)

tukey_KRN$Site[grep("M12", rownames(tukey_KRN$Site)), ]


# SAVE

# Extract M73 from cupule number Tukey results
m73_cupule_df <- as.data.frame(tukey_cupule$Site)
m73_cupule_df$Comparison <- rownames(tukey_cupule$Site)
m73_cupule_df <- m73_cupule_df[grep("M73", m73_cupule_df$Comparison), ]
m73_cupule_df$Trait <- "Cupule_number"

# Extract M12 from kernel row number Tukey results
m12_krn_df <- as.data.frame(tukey_KRN$Site)
m12_krn_df$Comparison <- rownames(tukey_KRN$Site)
m12_krn_df <- m12_krn_df[grep("M12", m12_krn_df$Comparison), ]
m12_krn_df$Trait <- "Mean_kernel_row"

cda_traits_df <- rbind(m73_cupule_df, m12_krn_df)

# Optional: clean column order
cda_traits_df <- cda_traits_df[, c("Trait", "Comparison", "diff", "lwr", "upr", "p adj")]

write_xlsx(
  list(
    "CDA_M12_M73_Traits" = cda_traits_df
  ),
  path = "CDA_M12_M73_Traits.xlsx"
)


# OLD BOXPLOTS  ------------------------------------


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


# Focus on LIP related groups only
# Get estimated marginal means for Period
pairwise_results <- emmeans(manova_model, ~ Period)

# Late Formative (LF) to Middle Horizon (MH)
contrast_LIP <- contrast(pairwise_results, 
                           method = "trt.vs.ctrl", 
                           ref = "LIP")  # LIP is the reference
print(contrast_LIP)

# Get site counts
site_counts_by_period <- maize_data %>%
  dplyr::select(Site, Period) %>%
  distinct() %>%                     # Remove duplicate site-period pairs
  group_by(Period) %>%
  summarise(Site_Count = n_distinct(Site)) %>%
  arrange(Period)

print(site_counts_by_period)

n_per_site_period <- maize_data %>%
  group_by(Site, Period) %>%
  summarise(N = n(), .groups = "drop") %>%
  arrange(Period, Site)

print(n_per_site_period)

# Save
lf_mh_df <- as.data.frame(contrast_LF_MH)
within_mh_df <- as.data.frame(contrast_within_MH)
lip_contrast_df <- as.data.frame(contrast_LIP)

write_xlsx(
  list(
    "LF_vs_MH" = lf_mh_df,
    "Within_MH" = within_mh_df,
    "LIP_vs_Others" = lip_contrast_df
  ),
  path = "Supplmentary_Period_Contrasts.xlsx"
)



# Summary Stats Across Period  --------

# Create summary statistics for each trait by Period
trait_summary <- maize_data %>%
  group_by(Period) %>%
  summarise(
    N_Length = sum(!is.na(Length_mm)),
    Mean_Length = mean(Length_mm, na.rm = TRUE),
    SD_Length = sd(Length_mm, na.rm = TRUE),
    Min_Length = min(Length_mm, na.rm = TRUE),
    Max_Length = max(Length_mm, na.rm = TRUE),
    
    N_Diameter = sum(!is.na(Diameter_mm)),
    Mean_Diameter = mean(Diameter_mm, na.rm = TRUE),
    SD_Diameter = sd(Diameter_mm, na.rm = TRUE),
    Min_Diameter = min(Diameter_mm, na.rm = TRUE),
    Max_Diameter = max(Diameter_mm, na.rm = TRUE),
    
    N_Cupules = sum(!is.na(Cupule_number)),
    Mean_Cupules = mean(Cupule_number, na.rm = TRUE),
    SD_Cupules = sd(Cupule_number, na.rm = TRUE),
    Min_Cupules = min(Cupule_number, na.rm = TRUE),
    Max_Cupules = max(Cupule_number, na.rm = TRUE),
    
    N_CupuleWidth = sum(!is.na(`Mean_Cupule_Width(mm)`)),
    Mean_CupuleWidth = mean(`Mean_Cupule_Width(mm)`, na.rm = TRUE),
    SD_CupuleWidth = sd(`Mean_Cupule_Width(mm)`, na.rm = TRUE),
    Min_CupuleWidth = min(`Mean_Cupule_Width(mm)`, na.rm = TRUE),
    Max_CupuleWidth = max(`Mean_Cupule_Width(mm)`, na.rm = TRUE),
    
    N_CupuleHeight = sum(!is.na(`Mean_cupule_height(mm)`)),
    Mean_CupuleHeight = mean(`Mean_cupule_height(mm)`, na.rm = TRUE),
    SD_CupuleHeight = sd(`Mean_cupule_height(mm)`, na.rm = TRUE),
    Min_CupuleHeight = min(`Mean_cupule_height(mm)`, na.rm = TRUE),
    Max_CupuleHeight = max(`Mean_cupule_height(mm)`, na.rm = TRUE),
    
    N_KernelRow = sum(!is.na(Mean_kernel_row)),
    Mean_KernelRow = mean(Mean_kernel_row, na.rm = TRUE),
    SD_KernelRow = sd(Mean_kernel_row, na.rm = TRUE),
    Min_KernelRow = min(Mean_kernel_row, na.rm = TRUE),
    Max_KernelRow = max(Mean_kernel_row, na.rm = TRUE),
    
    N_Weight = sum(!is.na(Total_Wt_g)),
    Mean_Weight = mean(Total_Wt_g, na.rm = TRUE),
    SD_Weight = sd(Total_Wt_g, na.rm = TRUE),
    Min_Weight = min(Total_Wt_g, na.rm = TRUE),
    Max_Weight = max(Total_Wt_g, na.rm = TRUE)
  )

# View the table
print(trait_summary)

# Optional: save it
writexl::write_xlsx(trait_summary, path = "Summary_Traits_By_Period.xlsx")


# Cupule depth for MH  --------
# Read your CSV file
maize_depth_data <- read_csv("cleaned_Ancient_Maize_Cob_Data_Sheet.csv")

depth_summary <- maize_depth_data %>%
  filter(!is.na(Mean_cupule_depth)) %>%
  # group_by(Site) %>%
  summarise(
    N = n(),
    Mean = mean(Mean_cupule_depth, na.rm = TRUE),
    SD = sd(Mean_cupule_depth, na.rm = TRUE),
    Min = min(Mean_cupule_depth, na.rm = TRUE),
    Max = max(Mean_cupule_depth, na.rm = TRUE)
  )

print(depth_summary)

# Combine your existing summaries and add depth_summary as a new sheet
write_xlsx(
  list(
    "Summary_Stats" = trait_summary,
    "Cupule_Depth" = depth_summary
  ),
  path = "Summary_Traits_By_Period.xlsx"
)


# Linear Regression across Period DOESNT WORK  --------
# Ensure Period is a factor (you've already done this)
maize_data$Period <- as.factor(maize_data$Period)

# Define traits of interest
traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)

traits <- c(
  "Length_mm",
  "Diameter_mm",
  "Cupule_number",
  "Mean_Cupule_Width(mm)",
  "Mean_cupule_height(mm)",
  "Mean_kernel_row",
  "Total_Wt_g"
)
lm_results$Length_mm


# Run linear regression for each trait
lm_results <- lapply(traits, function(trait) {
  formula <- as.formula(paste(trait, "~ Period"))
  model <- lm(formula, data = maize_data)
  summary(model)
})
names(lm_results) <- traits  # Label each result by trait

trait_pvals <- sapply(lm_results, function(model_summary) {
  anova(model_summary)$`Pr(>F)`[1]  # p-value for Period
})

# Order from most to least significant
trait_pvals[order(trait_pvals)]









# Focus on Cupule_number, Mean_Cupule_Width(mm), Mean_cupule_height --------


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





# FOCUS ON CUPULE CHARACTERISTICS ---------------------------------------------------------------------

# Reorder Period levels
maize_data$Period <- factor(maize_data$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Select and pivot the relevant columns
maize_long <- maize_data %>%
  dplyr::select(Period, 
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
                           "Kernel_row_number" = "Kernel Rows")

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




# CREATE FACETED BOXPLOTS FOR ALL RAW TRAITS OVER TIME
#install.packages("ggpubr")
library(ggpubr)
library(tidyr)


# Reorder Period levels
maize_data$Period <- factor(maize_data$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Select and pivot the relevant columns
maize_long <- maize_data %>%
  dplyr::select(Period, 
                Cupule_number,
                `Mean_Cupule_Width(mm)`,
                `Mean_cupule_height(mm)`,
                Mean_kernel_row,
                Length_mm,
                Diameter_mm,
                Total_Wt_g) %>%
  pivot_longer(cols = -Period,
               names_to = "Trait",
               values_to = "Value")

# Optional: Clean trait labels for prettier plotting3
maize_long$Trait <- as.character(maize_long$Trait)
maize_long$Trait <- dplyr::recode(maize_long$Trait,
                                  "Cupule_number" = "Cupule Number",
                                  "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
                                  "Mean_cupule_height(mm)" = "Cupule Height (mm)",
                                  "Mean_kernel_row" = "Kernel Rows",
                                  "Length_mm" = "Length (mm)",
                                  "Diameter_mm" = "Diameter (mm)",
                                  "Total_Wt_g" = "Total Weight (g)"
)


# Reorder Trait levels for plotting
maize_long$Trait <- factor(maize_long$Trait, levels = c(
  "Length (mm)",
  "Diameter (mm)",
  "Total Weight (g)",
  "Kernel Rows",
  "Cupule Number",
  "Cupule Width (mm)",
  "Cupule Height (mm)"
))



period_colors <- c(
  "LF" = "#fc8d62",       # greenish
  "LF/MH" = "#66c2a5",   # orangish
  "MH" = "#8da0cb",      # bluish
  "LIP" = "#e78ac3"      # pinkish
)

scale_fill_manual(values = period_colors)


# Make faceted boxplots
ggplot(maize_long, aes(x = Period, y = Value, fill = Period)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_fill_manual(values = period_colors) +
  stat_compare_means(method = "anova", label = "p.format") +  # Global ANOVA p-value per facet
  labs(x = "Time Period", y = "Value", title = "Distribution of Morphological Traits Over Time") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.4, 0.2),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha('white', 0.6), color = NA))




# Linear Regression ---------------------------------------------------------------------


maize_data$Period_numeric <- dplyr::recode(maize_data$Period,
                                           "LF" = 1,
                                           "LF/MH" = 2,
                                           "MH" = 3,
                                           "LIP" = 4
)

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

maize_long <- maize_data %>%
  dplyr::select(Period_numeric, all_of(traits)) %>%
  pivot_longer(cols = -Period_numeric, names_to = "Trait", values_to = "Value")

maize_long$Trait <- dplyr::recode(maize_long$Trait,
                                  "Cupule_number" = "Cupule Number",
                                  "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
                                  "Mean_cupule_height(mm)" = "Cupule Height (mm)",
                                  "Mean_kernel_row" = "Kernel Rows",
                                  "Length_mm" = "Length (mm)",
                                  "Diameter_mm" = "Diameter (mm)",
                                  "Total_Wt_g" = "Total Weight (g)")


#maize_long$Trait <- recode(maize_long$Trait,
#                           "Cupule_number" = "Cupule Number",
#                           "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
#                           "Mean_cupule_height(mm)" = "Cupule Height (mm)",
#                           "Mean_kernel_row" = "Kernel Rows",
#                           "Length_mm" = "Length (mm)",
#                           "Diameter_mm" = "Diameter (mm)",
#                           "Total_Wt_g" = "Total Weight (g)"
#)



#maize_long$Trait <- factor(maize_long$Trait, levels = c(
#  "Length (mm)",
#  "Diameter (mm)",
#  "Total Weight (g)",
#  "Kernel Rows",
#  "Cupule Number",
#  "Cupule Width (mm)",
#  "Cupule Height (mm)"
#))


ggplot(maize_long, aes(x = Period_numeric, y = Value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  stat_cor(method = "pearson", label.x = 1, label.y.npc = "bottom", size = 3) +
  stat_regline_equation(label.y.npc = "top", size = 3) +
  facet_wrap(~ Trait, scales = "free_y") +
  labs(
    x = "Chronological Period (numeric)",
    y = "Trait Value",
    title = "Linear Regression of Morphological Traits Over Time"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

# Updated Axis ---------------------------------------------------------------------

# Reshape the data
maize_long <- maize_data %>%
  dplyr::select(Period, Period_numeric,
         Length_mm,
         Diameter_mm,
         Total_Wt_g,
         Mean_kernel_row,
         Cupule_number,
         `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`) %>%
  pivot_longer(cols = -c(Period, Period_numeric),
               names_to = "Trait",
               values_to = "Value")

# Recode Trait labels for display
maize_long$Trait <- dplyr::recode(maize_long$Trait,
                                  "Length_mm" = "Length (mm)",
                                  "Diameter_mm" = "Diameter (mm)",
                                  "Total_Wt_g" = "Total Weight (g)",
                                  "Mean_kernel_row" = "Kernel Rows",
                                  "Cupule_number" = "Cupule Number",
                                  "Mean_Cupule_Width(mm)" = "Cupule Width (mm)",
                                  "Mean_cupule_height(mm)" = "Cupule Height (mm)"
)

# Reorder facets
#maize_long$Trait <- factor(maize_long$Trait, levels = c(
#  "Length (mm)",
#  "Diameter (mm)",
#  "Total Weight (g)",
#  "Kernel Rows",
#  "Cupule Number",
#  "Cupule Width (mm)",
#  "Cupule Height (mm)"
#))


ggplot(maize_long, aes(x = Period_numeric, y = Value)) +
  geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
  labs(x = "Time Period", y = "Trait Value",
       title = "Linear Regression of Morphological Traits Across Time Periods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.4, 0.2),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha('white', 0.6), color = NA))

library(ggpubr)

# ggplot(maize_long, aes(x = Period_numeric, y = Value)) +
#  geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
#  geom_smooth(method = "lm", se = TRUE, color = "black") +
#  stat_regline_equation(
#    aes(label =  paste(..eq.label.., ..p.label.., sep = "~~~")),
#    label.x.npc = "left", label.y.npc = 0.95,
#    size = 3,
#    color = "black"
#  ) +
#  facet_wrap(~ Trait, scales = "free_y") +
#  scale_color_manual(values = period_colors) +
#  scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
#  labs(x = "Time Period", y = "Trait Value",
#       title = "Linear Regression of Morphological Traits Across Time Periods") +
#  theme_minimal() +
#  theme(strip.text = element_text(size = 12, face = "bold"),
#        axis.text.x = element_text(angle = 45, hjust = 1)) +
#  theme(panel.grid = element_blank(),
#        legend.position = c(0.4, 0.2),
#        legend.justification = c("right", "top"),
#        legend.background = element_rect(fill = alpha('white', 0.6), color = NA))

# ADD P Values
#install.packages("ggpmisc")
library(ggpmisc)

ggplot(maize_long, aes(x = Period_numeric, y = Value)) +
  geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, 
    parse = TRUE,
    label.x.npc = "center",
    label.y.npc = 1.0,
    size = 4
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
  labs(x = "Time Period", y = "Trait Value",
       title = "Linear Regression of Morphological Traits Across Time Periods") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid = element_blank(),
        legend.position = c(0.4, 0.2),
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha('white', 0.6), color = NA))




# Linear w 2 plots ---------------------------------------------------------------------

# Create a list of traits to filter
cob_traits <- c("Length (mm)", "Diameter (mm)", "Total Weight (g)")
cupule_traits <- c("Kernel Rows", "Cupule Number", "Cupule Width (mm)", "Cupule Height (mm)")

# Plot 1: Cob Dimensions --------------------------------------------------
ggplot(filter(maize_long, Trait %in% cob_traits),
       aes(x = Period_numeric, y = Value)) +
  geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE,
    label.x.npc = "center", label.y.npc = 1.0, size = 3
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
  labs(
    x = "Time Period", y = "Trait Value",
    title = "Linear Regression of Cob Dimensions Across Time Periods"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

# Plot 2: Cupule Characteristics ------------------------------------------
ggplot(filter(maize_long, Trait %in% cupule_traits),
       aes(x = Period_numeric, y = Value)) +
  geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  stat_poly_eq(
    aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x, parse = TRUE,
    label.x.npc = "center", label.y.npc = 1.0, size = 3
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  scale_color_manual(values = period_colors) +
  scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
  labs(
    x = "Time Period", y = "Trait Value",
    title = "Linear Regression of Cupule Characteristics Across Time Periods"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )



  


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


# UPDATE FOR COLOR SCHEME 

ggplot(trait_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Trait), size = 4) +
  scale_color_manual(values = c(
    "#FFB3BA", "#FFDFBA", "#A6E3E9", "#BAFFC9", "#BAE1FF", "#D5BAFF", "#FFCCE5"
  )) +
  geom_text(aes(label = Trait), vjust = -0.8, size = 5) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "PCA of All Traits Across All Sites",
       x = paste0("PC1 (", round(summary(pca_traits)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_traits)$importance[2,2]*100, 1), "%)"))


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
  "M73", "M103", "M7", "M43", "M43", 
  "M10", "M10", "M1", 
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
site_order <- c("M73", "M103", "M7", "M43", "M43", 
                "M10", "M10", "M1", 
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
  Site %in% c("M1", "M43", "M43") ~ "Chen Chen",
  Site == "M12" ~ "Omo",
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
  labs(title = "QQ Plots for Chen Chen and Omo Style Groups")

# Test for normality in each group TOO CONSERVATIVE
#shapiro_test <- pc_scores_grouped %>%
#  group_by(Style_Group) %>%
#  summarise(p_value = shapiro.test(PC1)$p.value)

#print(shapiro_test)

# Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(PC1 ~ Style_Group, data = pc_scores_grouped)

print(wilcox_test_result)


# Small sample size of Omo 

set.seed(42)

# Observed difference in medians or means
obs_diff <- with(pc_scores_grouped, 
                 median(PC1[Style_Group == "Chen Chen"]) -
                   median(PC1[Style_Group == "Omo"]))

# Permute group labels
perm_diffs <- replicate(10000, {
  shuffled <- sample(pc_scores_grouped$Style_Group)
  median(pc_scores_grouped$PC1[shuffled == "Chen Chen"]) -
    median(pc_scores_grouped$PC1[shuffled == "Omo"])
})

# Two-sided p-value
perm_p <- mean(abs(perm_diffs) >= abs(obs_diff))
cat("Permutation test p-value:", perm_p, "\n")






# Chen chen v Omo style v M10 Temple
# Define Chen Chen, Omo, and Ritual Complex style groups
pc_scores_labeled <- pc_scores_labeled %>%
  mutate(Style_Group = case_when(
    Site %in% c("M1-95-Chen-Chen", "M43", "M43") ~ "Chen Chen",
    Site == "M12" ~ "Omo",
    Site %in% c("M10") ~ "M10",  # Add Ritual Complex group
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
#shapiro_test <- pc_scores_grouped %>%
#  group_by(Style_Group) %>%
#  summarise(p_value = shapiro.test(PC1)$p.value)

#print(shapiro_test)

# Kruskal-Wallis test for non-normal data
#kruskal_test_result <- kruskal.test(PC1 ~ Style_Group, data = pc_scores_grouped)
#print(kruskal_test_result)

# Pairwise Wilcoxon tests with p-value adjustment
pairwise_wilcoxon <- pairwise.wilcox.test(pc_scores_grouped$PC1, pc_scores_grouped$Style_Group, p.adjust.method = "BH")
print(pairwise_wilcoxon)

# SAVE


# ---- 1. Wilcoxon rank-sum test summary ----
wilcox_df <- data.frame(
  W = wilcox_test_result$statistic,
  p_value = wilcox_test_result$p.value
)

# ---- 2. Permutation test result ----
perm_df <- data.frame(
  Observed_Diff = obs_diff,
  Permutation_p_value = perm_p
)

# ---- 3. Pairwise Wilcoxon tests ----
pairwise_df <- as.data.frame(pairwise_wilcoxon$p.value)
pairwise_df$Comparison <- rownames(pairwise_df)

# Optional: reshape if you want long-format
# pairwise_long <- tidyr::pivot_longer(pairwise_df, -Comparison, names_to = "vs", values_to = "adjusted_p")

# ---- Write all to an Excel file ----
write_xlsx(
  list(
    "OmoChenChen_Wilcoxon_Test" = wilcox_df,
    "OmoChenChen_Permutation_Test" = perm_df,
    "OmoChenChenRitual_Pairwise" = pairwise_df
  ),
  path = "Supplmentary_Omo_ChenChen_M10.xlsx"
)


# AMS Dates ---------------------------------------------------

ams_dates <- read_csv("AMS_Dates_Final.csv")


maize_dated <- maize_data %>%
  inner_join(ams_dates %>% select(Sample_ID, `calibrated_2_sigma`), by = "Sample_ID")


morphology_long <- maize_dated %>%
  dplyr::select(Sample_ID, Period, calibrated_2_sigma, 
                Length_mm,
                Diameter_mm,
                Cupule_number,
                `Mean_Cupule_Width(mm)`,
                `Mean_cupule_height(mm)`,
                Mean_kernel_row,
                Total_Wt_g) %>%
  pivot_longer(cols = -c(Sample_ID, Period, calibrated_2_sigma), 
               names_to = "Trait", values_to = "Value")%>%
  filter(!is.na(Value))


# Sample size too small for real stats, but can look for patterns
# Summary Stats
trait_summary <- morphology_long %>%
  group_by(Trait, Period) %>%
  summarise(
    n = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(morphology_long, aes(x = Period, y = Value, fill = Period)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_bw() +
  ggtitle("Trait Distributions by Period (n = 20)")


# Ensure Period is a factor for correct order (adjust levels as needed)
trait_summary$Period <- factor(trait_summary$Period, levels = c("LF", "LF/MH", "MH", "LIP"))

# Line plot of mean trait values by period
ggplot(trait_summary, aes(x = Period, y = mean, group = Trait, color = Trait)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean Trait Values by Period",
    x = "Period",
    y = "Mean Trait Value",
    color = "Trait"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# MH Stats only
mh_data <- maize_dated %>%
  filter(Period == "MH")

mh_data %>%
  select(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) %>%
  summary()


mh_long <- mh_data %>%
  pivot_longer(
    cols = c(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
             `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g),
    names_to = "Trait", values_to = "Value"
  )

# Plot MH
ggplot(mh_long, aes(x = Trait, y = Value)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Trait Distributions Within MH", y = "Value", x = "Trait") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Look at pairwise correlations
mh_numeric <- mh_data %>%
  select(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g)

cor_matrix <- cor(mh_numeric, use = "pairwise.complete.obs")
round(cor_matrix, 2)

# Prepare data to save
mh_summary <- mh_data %>%
  select(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) %>%
  summarise(across(everything(),
                   list(n = ~sum(!is.na(.)),
                        mean = ~mean(., na.rm = TRUE),
                        sd = ~sd(., na.rm = TRUE),
                        min = ~min(., na.rm = TRUE),
                        max = ~max(., na.rm = TRUE)))) %>%
  pivot_longer(everything(),
               names_to = c("Trait", ".value"),
               names_sep = "_")

cor_matrix <- mh_data %>%
  select(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) %>%
  cor(use = "pairwise.complete.obs") %>%
  round(2) %>%
  as.data.frame()

cor_matrix <- tibble::rownames_to_column(cor_matrix, "Trait")

# Save
write_xlsx(
  list(
    "Trait_Summary" = trait_summary,
    "MH_Summary" = mh_summary,
    "MH_Correlations" = cor_matrix
  ),
  path = "Trait_Summary_AMS.xlsx"
)

# Trends within sorted MH by AMS dates

mh_data <- maize_dated %>%
  filter(Period == "MH" & !is.na(calibrated_2_sigma))

unique(mh_data$calibrated_2_sigma)

mh_data <- mh_data %>%
  mutate(cal_2sigma_clean = gsub("\\s+", "", calibrated_2_sigma))  # remove whitespace


mh_data$cal_2sigma_clean <- factor(mh_data$cal_2sigma_clean, levels = sort(unique(mh_data$cal_2sigma_clean)))

mh_long <- mh_data %>%
  select(Sample_ID, cal_2sigma_clean,
         Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) %>%
  pivot_longer(
    cols = -c(Sample_ID, cal_2sigma_clean),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))


# Define custom order (adjust as needed)
date_lookup <- c(
    "M43=6385.01"    = "AD 782 - 880",
    "M43=6261.06"    = "AD 798 - 885",
    "M1=124039.06"   = "AD 800 - 890",
    "M43=6385.02"    = "AD 820 - 894",
    "M43=6385.10"    = "AD 825 - 896",
    "M43=6174.07"    = "AD 826 - 903",
    "M43=6503.06"    = "AD 881 - 979",
    "M43=6174.02"    = "AD 884 - 980",
    "M43=6174.05"    = "AD 892 - 995",
    "M43=6321.09"    = "AD 916 - 973",
    "M43=6174.09"    = "AD 942 - 994",
    "M10=12744.02"   = "AD 920 - 957",
    "M10=7942.17"    = "AD  1063 - 1175"
  )

mh_long <- mh_long %>%
  mutate(calibrated_2_sigma = date_lookup[Sample_ID])

cal_2sigma_order <- unique(date_lookup)  # use the correct order from your lookup
mh_long$calibrated_2_sigma <- factor(mh_long$calibrated_2_sigma, levels = cal_2sigma_order)


ggplot(mh_long, aes(x = calibrated_2_sigma, y = Value, group = Trait, color = Trait)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "MH Trait Trends by Sample Calibrated 2-Sigma Date",
    x = "Calibrated 2-Sigma Date",
    y = "Trait Value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Step 1: Calculate summary stats for each trait and date
mh_trait_stats <- mh_long %>%
  group_by(Trait) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("n = ", n, "\nM = ", round(mean, 1), "\nSD = ", round(sd, 1)))


# Step 2: Plot with stat_summary + summary labels
ggplot(mh_long, aes(x = calibrated_2_sigma, y = Value, group = Trait, color = Trait)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_text(
    data = mh_trait_stats,
    aes(x = 1, y = Inf, label = label),
    inherit.aes = FALSE,
    size = 3,
    vjust = 1.1,
    hjust = 0,
    color = "black"
  ) +
  facet_wrap(~ Trait, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "AMS MH Trait Trends Over Time (Per-Trait Summary Stats)",
    x = "Calibrated 2-Sigma Date",
    y = "Trait Value"
  ) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


