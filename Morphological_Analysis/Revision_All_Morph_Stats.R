# =============================================================================
# DEBUGGED VERSION, 22 Oct 2026. Every change is marked "# FIX 22-Oct" in the code.
# Inputs: Merged_Ancient_Maize_Cob_Data_Sheet_with_IDs.csv (358 cobs, Specimen_ID),
#         cleaned_Ancient_Maize_Cob_Data_Sheet.csv (cupule depth only),
#         Revision_Supplementary_AMS_full.csv (SHCal20 dates).
# Changes that alter results:
#  1. M7 transition section no longer overwrites maize_data. It had dropped M44
#     (Cerro Chamos) from every later analysis, including the published CDA figure,
#     and lumped M11 (LIP, Tumilaca) into "Later_Tiwanaku".
#  2. site_order used old site names, which turned M43, M10 and M1 into NA.
#  3. Chen Chen style group used "M1-95-Chen-Chen", so Chen Chen M1 was dropped.
#  4. AMS section now uses the SHCal20 table; the old file had IntCal20 ranges and
#     wrong Omo values. Plot labels are built from the table, not hard-coded.
#  5. CDA figure axis labels were hard-coded ("Can1 (74.2%)", "Can2 (12.6%)");
#     now computed from the current run.
# Changes that fix errors:
#  6. select()/recode() masked by MASS (via multcomp) and car -> pinned to dplyr.
#  7. Duplicated factor levels ("M43","M43","M10","M10") -> error in factor().
#  8. Trait-by-Period regressions: missing backticks around trait names; anova()
#     was applied to summary objects.
#  9. Cupule_Surface_Area was undefined -> Surface_Area.
# 10. wilcox_test_result / obs_diff were never created -> the Omo/Chen Chen/M10
#     supplement now reports Kruskal-Wallis + BH-adjusted pairwise Wilcoxon tests.
# 11. tukey_all_sites undefined (unused output) -> disabled.
# 12. MH summary table split trait names at "_" -> trait names were garbled.
# REVISION 22-Oct (statistics, for R2.4); marked "# REVISION 22-Oct":
# 13. Removed the cross-trait ANOVA/Tukey/Welch/Games-Howell (and their supplement sheets).
# 14. Removed the transposed PCA ("PC1 95%"); Figure 3 is now a biplot of the scaled
#     PCA of cobs (Figure3_PCA_biplot.png, Supplementary_PCA.xlsx).
# 15. Removed ordinal period regressions and the emmeans contrasts averaged across
#     traits. New section at the end: Welch ANOVA per trait + Games-Howell pairwise
#     tests with 95% CIs and Hedges' g, all on the same 344 complete cases
#     (Supplementary_Period_Comparisons.xlsx, Figure4_traits_by_period.png).
# 16. Pillai's test throughout.
# 17. Chen Chen (M1, M43) vs Omo (M12), M10 excluded: Wilcoxon + effect size r
#     (Supplementary_ChenChen_vs_Omo.xlsx).
# 18. Saved outputs that were console-only: Supplementary_Period_MANOVA.xlsx,
#     MANOVA + univariate sheets in TukeyHSD_M7_Transition.xlsx, Supplementary_CDA.xlsx.
# 19. CDA no longer drops maize_data[-185, ] (a row dropped by position; after the data
#     correction it would have removed a valid cob). Fixed axis limits removed.
# Data corrections are in the datasheet itself (M43_6261_03 length/diameter swap;
# M103_1161_30 diameter 0.05 -> 5.0, to verify against in-person measurement).
# =============================================================================

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

# FIX 22-Oct: multcomp loads MASS, whose select() masks dplyr::select(), and car masks
# dplyr::recode(). Pin the dplyr versions so select()/recode() calls below work.
select <- dplyr::select
filter <- dplyr::filter
recode <- dplyr::recode

# Set working directory
setwd("C:/Users/chamb/OneDrive/Documents/Moquegua_Morphology/Morphological_Analysis")

# Define file path

# FIX 22-Oct: new datasheet with Specimen_ID, Lot and UCIAMS columns
file_data <- "Merged_Ancient_Maize_Cob_Data_Sheet_with_IDs.csv"

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

# REVISION 22-Oct: disabled. Cross-trait ANOVA/Tukey/Welch/Games-Howell compared raw values of traits with different units (R2.4: remove).
if (FALSE) {
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
  
  
  
  
  
}  # end disabled block

# MANOVA for general differences across sites -----------------------------


# Make sure 'Site' is a factor
mh_data$Site <- as.factor(mh_data$Site)

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

# REVISION 22-Oct: disabled. Welch/Games-Howell across traits removed (R2.4).
if (FALSE) {
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
  
}  # end disabled block

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
  dplyr::select(Site, all_of(traits)) %>%
  filter(complete.cases(.))

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
  dplyr::select(Site, all_of(traits)) %>%
  filter(complete.cases(.))

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
manova_df_lip <- as.data.frame(manova_summary_lip$stats)
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
# FIX 22-Oct: tukey_all_sites is never created and tukey_all_df is not saved; disabled
# tukey_all_df <- do.call(rbind, tukey_all_sites) %>%
#   dplyr::select(Trait, Comparison, diff, lwr, upr, `p adj`)

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


# Split MH Plots into Cobs and Cupule Characteristics FIG FOR PUB
# Define trait groupings
cob_traits <- c("Length (mm)", "Diameter (mm)", "Total Weight (g)")
cupule_traits <- c("Cupule Number", "Cupule Width (mm)", "Cupule Height (mm)", "Kernel Rows")

ggplot(filter(mh_long, Trait %in% cob_traits), aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  facet_wrap(~ Trait, scales = "free_y") +
  stat_compare_means(method = "anova", label = "p.format", size = 4) +
  labs(
    title = "A. Cob Dimensions Across MH Sites",
    x = "Site", y = "Trait Value"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16),                       # Axis titles
    axis.text = element_text(size = 14),                        # Axis tick labels
    legend.title = element_text(size = 16),                     # Legend title
    legend.text = element_text(size = 14),   
    panel.grid = element_blank()
  )

#Save
ggsave("Cob_Dimensions_MH.jpeg", dpi = 300, width = 10, height = 7, units = "in")

ggplot(filter(mh_long, Trait %in% cupule_traits), aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "black") +
  geom_blank(aes(y = Value * 1.15)) +
  facet_wrap(~ Trait, scales = "free_y") +
  stat_compare_means(method = "anova", label = "p.format", size = 4) +
  labs(
    title = "B. Cupule Characteristics Across MH Sites",
    x = "Site", y = "Trait Value"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1), 
    axis.title = element_text(size = 16),                       # Axis titles
    axis.text = element_text(size = 14),                        # Axis tick labels
    legend.title = element_text(size = 16),                     # Legend title
    legend.text = element_text(size = 14),   
    panel.grid = element_blank()
  )

#Save
ggsave("Cupule_Characteristics_MH.jpeg", dpi = 300, width = 10, height = 7, units = "in")





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

# FIX 22-Oct: this used to overwrite maize_data, which silently dropped M44 (Cerro Chamos)
# from every later analysis (CDA figure, PCA, plots). It also lumped M11 (LIP, Tumilaca)
# into "Later_Tiwanaku". Now a separate object, excluding both LIP sites.
maize_transition <- maize_data %>%
  filter(!Site %in% c("M44", "M11")) %>%  # Exclude LIP sites
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

morph_data <- maize_transition %>%  # FIX 22-Oct
  dplyr::select(Cultural_Period, all_of(traits)) %>%
  filter(complete.cases(.))  # Remove rows with NA

# MANOVA model
manova_model <- manova(as.matrix(morph_data[, traits]) ~ Cultural_Period, data = morph_data)

# Summary
m7_manova_summary <- summary(manova_model, test = "Pillai")  # REVISION 22-Oct: Pillai throughout
print(m7_manova_summary)
m7_manova_df <- as.data.frame(m7_manova_summary$stats); m7_manova_df$Effect <- rownames(m7_manova_df)


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
print(m7_all)  # REVISION 22-Oct: View() fails outside RStudio

# SAVE
write_xlsx(
  list(
    "TukeyHSD_M7_Transition" = m7_all,
    "MANOVA_Pillai" = m7_manova_df,    # REVISION 22-Oct: saved
    "Univariate_ANOVAs" = do.call(rbind, lapply(names(results), function(tr)
      broom::tidy(results[[tr]][[1]]) %>% mutate(Trait = tr))),
    "n_by_group" = as.data.frame(table(morph_data$Cultural_Period))
  ),
  path = "TukeyHSD_M7_Transition.xlsx"
)


# Plot CDA FIG FOR PUB ----------------------------------------------------------------


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


# TEST FOR PUB
while (dev.cur() > 1) dev.off()

# Extract scores and manually add the Site info
# Extract the model frame used in MANOVA
model_data <- model.frame(manova_model)

# Add Site to the CDA scores
scores <- cda$scores
scores$Site <- model_data$Site  # This will match row-for-row
scores$Site <- as.factor(scores$Site)


# Now safely get the levels
site_levels <- levels(scores$Site)
n_sites <- length(site_levels)

# Colors and symbols
colors <- 1:n_sites
pch_vals <- 1:n_sites

# Extract eigenvalues and compute proportion of variance
eig_vals <- cda$eigenvalues
var_explained <- eig_vals / sum(eig_vals) * 100

# Round for display
can1_percent <- round(var_explained[1], 1)
can2_percent <- round(var_explained[2], 1)

xlab <- paste0("Can1 (", can1_percent, "%)")
ylab <- paste0("Can2 (", can2_percent, "%)")



# FIX 22-Oct: title() was called before plot() existed; labels now passed to plot() below
plot(
  scores$Can1, scores$Can2,
  col = colors[scores$Site],
  pch = pch_vals[scores$Site],
  xlab = xlab,  # FIX 22-Oct: was hard-coded "Can1 (74.2%)" from an older run
  ylab = ylab,  # FIX 22-Oct: was hard-coded "Can2 (12.6%)"
  main = "CDA Trait Differences Across Sites"
)

# Legend
legend("topright",
       legend = site_levels,
       col = colors,
       pch = pch_vals,
       title = "Site",
       cex = 0.8)

# Extract trait loadings for the first two canonical dimensions
trait_loadings <- cda$structure[, 1:2]
# Scale loadings for visualization
arrow_scale <- 5  # Adjust as needed
trait_loadings_scaled <- trait_loadings * arrow_scale
# Add arrows
arrows(
  x0 = 0, y0 = 0,
  x1 = trait_loadings_scaled[, 1],
  y1 = trait_loadings_scaled[, 2],
  length = 0.1, col = "darkgray", lwd = 1.5
)

# Add trait labels at end of arrows
text(
  x = trait_loadings_scaled[, 1],
  y = trait_loadings_scaled[, 2],
  labels = rownames(trait_loadings),
  col = "black",
  cex = 0.9,
  pos = 4
)



# ORIGNAL CDA FIG FOR PUB
# REVISION 22-Oct: no exclusion. maize_data[-185, ] dropped a row BY POSITION; the "outlier"
# (M10_12478_11) was a data-entry error, now corrected, and row 185 is now a different,
# valid cob (M10_8004_04). The CDA is refit on all complete cases (n = 344).
maize_data_no_outlier <- maize_data   # object name kept so the code below runs unchanged
maize_data_no_outlier$Site <- factor(maize_data_no_outlier$Site)
manova_model <- manova(cbind(Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
                             `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) ~ Site,
                       data = maize_data_no_outlier)
cda <- candisc(manova_model)

# Save CDA statistics (Figure 5)
cda_pct <- cda$eigenvalues / sum(cda$eigenvalues) * 100
write_xlsx(list(
  Canonical_dimensions = data.frame(Dimension = paste0("Can", seq_along(cda$eigenvalues)),
                                    Eigenvalue = cda$eigenvalues,
                                    Canonical_R2 = cda$canrsq,
                                    Percent = cda_pct,
                                    Cumulative = cumsum(cda_pct)),
  Structure_loadings = data.frame(Trait = rownames(cda$structure), cda$structure),
  n_by_site = as.data.frame(table(model.frame(manova_model)$Site))),
  path = "Supplementary_CDA.xlsx")

# Open JPEG device FIG FOR PUB
jpeg("CDA_Traits_Sites.jpg", width = 10, height = 5, units = "in", res = 300)


# Plot canonical variates with custom axis limits WORKS
plot(cda,
     main = "CDA for Traits Across Sites",
     cex = 0.3        # REVISION 22-Oct: fixed xlim/ylim removed so no points are cut off
)


# Add legend 
legend("topright",  # REVISION 22-Oct: was fixed coordinates (x = 21, y = 8)
       legend = levels(maize_data_no_outlier$Site),
       pch = 1:length(levels(maize_data_no_outlier$Site)),
       col = 1:length(levels(maize_data_no_outlier$Site)),
       title = "Site", 
       cex = 0.7, 
       bty = "n")

# Close device if saving
dev.off()

# TRY WITHOUT WEIGHT BC Not Stat Sig 

# Fit MANOVA model with all traits
maize_data_no_outlier <- maize_data  # REVISION 22-Oct: no positional exclusion


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
# summary(scores$Can2)  # REVISION 22-Oct: 'scores' not yet defined here

# Where is the M10 outlier?
#outlier_row <- scores[which.min(scores$Can2), ]
#print(outlier_row)

scores <- cda$scores

scores <- cda$scores
plot(scores$Can1, scores$Can2, col = scores$Site,
     pch = as.numeric(scores$Site),
     xlab = "Can1", ylab = "Can2", main = "CDA Scores by Site")

# Set legend as factor
maize_data_no_outlier$Site <- factor(maize_data_no_outlier$Site)

# Plot canonical variates with custom axis limits
plot(cda,
     main = "CDA for Significant Traits Across Sites (No Weight)",
     xlim = c(-8, 14),
     ylim = c(-8, 8),
     cex = 5.0  # Increase this for larger point symbols
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
# FIX 22-Oct: old site names turned M43, M10 and M1 into NA
site_order <- c("M73", "M103", "M7", "M43", "M10", "M1", "M12", "M11", "M44")

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

# REVISION 22-Oct: save the period MANOVA so the Results numbers trace to a file
period_manova_df <- as.data.frame(manova_summary$stats)
period_manova_df$Effect <- rownames(period_manova_df)
period_univ_df <- do.call(rbind, lapply(names(univariate_results), function(tr) {
  broom::tidy(univariate_results[[tr]]) %>% mutate(Trait = trimws(sub("Response", "", tr)))
}))
write_xlsx(list(MANOVA_Pillai_Period = period_manova_df,
                Univariate_ANOVAs_Period = period_univ_df,
                n = data.frame(n_complete_cases = nrow(residuals(manova_model)))),
           path = "Supplementary_Period_MANOVA.xlsx")
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




# REVISION 22-Oct: disabled. emmeans on the MANOVA averaged raw differences across traits with different units (R2.4). Replaced by 'Period comparisons per trait' at the end of the script.
if (FALSE) {
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
  
  
  
}  # end disabled block

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


# REVISION 22-Oct: disabled. Ordinal regression across periods (R2.4: implies a continuous trajectory).
if (FALSE) {
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
  # lm_results$Length_mm   # FIX 22-Oct: was called before lm_results existed
  
  
  # Run linear regression for each trait
  lm_results <- lapply(traits, function(trait) {
    # FIX 22-Oct: backticks needed for names like Mean_Cupule_Width(mm); keep the lm object
    formula <- as.formula(paste0("`", trait, "` ~ Period"))
    lm(formula, data = maize_data)
  })
  names(lm_results) <- traits  # Label each result by trait
  
  trait_pvals <- sapply(lm_results, function(model) {
    anova(model)$`Pr(>F)`[1]  # p-value for Period
  })
  
  # Order from most to least significant
  trait_pvals[order(trait_pvals)]
  
  
  
  
  
  
  
  
  
}  # end disabled block

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
  Surface_Area,  # FIX 22-Oct: was Cupule_Surface_Area (undefined). Using the Surface Area instead of Mean_Cupule_Width(mm) and Mean_cupule_height(mm)
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




# REVISION 22-Oct: disabled. Ordinal regressions and Figure 4 regression panels (R2.4). Figure 4 is now boxplots by period (see end of script).
if (FALSE) {
  # Linear Regression ---------------------------------------------------------------------
  period_colors <- c(
    "LF" = "#fc8d62",       # greenish
    "LF/MH" = "#66c2a5",   # orangish
    "MH" = "#8da0cb",      # bluish
    "LIP" = "#e78ac3"      # pinkish
  )
  
  scale_fill_manual(values = period_colors)
  
  
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
  
  
  
  
  # Linear w 2 plots FIG FOR PUB ---------------------------------------------------------------------
  
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
      label.x.npc = "center", label.y.npc = 1.0, size = 4
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    scale_color_manual(values = period_colors) +
    scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
    labs(
      x = "Time Period", y = "Trait Value",
      title = "A. Linear Regression of Cob Dimensions Across Time Periods"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1), 
      axis.title = element_text(size = 16),                       # Axis titles
      axis.text = element_text(size = 14),                        # Axis tick labels
      legend.title = element_text(size = 16),                     # Legend title
      legend.text = element_text(size = 14),   
      panel.grid = element_blank()
    )
  
  #Save
  ggsave("Linear_Regression_Cob_Dimesions.jpeg", dpi = 300, width = 10, height = 7, units = "in")
  
  # Plot 2: Cupule Characteristics ------------------------------------------
  ggplot(filter(maize_long, Trait %in% cupule_traits),
         aes(x = Period_numeric, y = Value)) +
    geom_point(aes(color = Period), position = position_jitter(width = 0.1), alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    geom_blank(aes(y = Value * 1.15)) +
    stat_poly_eq(
      aes(label = paste(..eq.label.., ..p.value.label.., sep = "~~~")),
      formula = y ~ x, parse = TRUE,
      label.x.npc = "center", label.y.npc = 1.0, size = 4
    ) +
    facet_wrap(~ Trait, scales = "free_y") +
    scale_color_manual(values = period_colors) +
    scale_x_continuous(breaks = 1:4, labels = c("LF", "LF/MH", "MH", "LIP")) +
    labs(
      x = "Time Period", y = "Trait Value",
      title = "B. Linear Regression of Cupule Characteristics Across Time Periods"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1), 
      axis.title = element_text(size = 16),                       # Axis titles
      axis.text = element_text(size = 14),                        # Axis tick labels
      legend.title = element_text(size = 16),                     # Legend title
      legend.text = element_text(size = 14),   
      panel.grid = element_blank()
    )
  
  
  #Save
  ggsave("Linear_Regression_Cupule_Characteristics.jpeg", dpi = 300, width = 10, height = 7, units = "in")
  
  
  
  
}  # end disabled block

# PCA FIGURE FOR PUB ---------------------------------------------------------------------


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

# REVISION 22-Oct: Figure 3 = this PCA (7 traits, scaled, complete cases), not the transposed PCA
pca_var <- round(summary(morpho_pca)$importance[2, 1:2] * 100, 1)
pca_load <- as.data.frame(morpho_pca$rotation[, 1:2]); pca_load$Trait <- rownames(pca_load)
pca_scores_fig <- as.data.frame(morpho_pca$x[, 1:2])
pca_scores_fig$Period <- factor(maize_data$Period[complete.cases(maize_data[, c("Length_mm","Diameter_mm","Cupule_number","Mean_Cupule_Width(mm)","Mean_cupule_height(mm)","Mean_kernel_row","Total_Wt_g")])],
                                levels = c("LF", "LF/MH", "MH", "LIP"))
arrow_scale <- 3
fig3 <- ggplot(pca_scores_fig, aes(PC1, PC2, colour = Period)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_segment(data = pca_load, aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "cm"))) +
  geom_text(data = pca_load, aes(x = PC1 * arrow_scale * 1.12, y = PC2 * arrow_scale * 1.12, label = Trait),
            inherit.aes = FALSE, size = 3) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"), y = paste0("PC2 (", pca_var[2], "%)"),
       title = paste0("PCA of seven cob traits (n = ", nrow(pca_scores_fig), ")")) +
  theme_minimal()
print(fig3)
ggsave("Figure3_PCA_biplot.png", fig3, width = 7, height = 5.5, dpi = 300)
write_xlsx(list(Variance = data.frame(PC = colnames(morpho_pca$x),
                                      Proportion = summary(morpho_pca)$importance[2, ],
                                      Cumulative = summary(morpho_pca)$importance[3, ]),
                Loadings = data.frame(Trait = rownames(morpho_pca$rotation), morpho_pca$rotation)),
           path = "Supplementary_PCA.xlsx")

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



# REVISION 22-Oct: disabled. Transposed PCA (traits as observations) is not a valid PCA; its '95% PC1' must not be reported.
if (FALSE) {
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
  
  # FIG FOR PUB - PCA AS SYMBOLS
  ggplot(trait_scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(shape = Trait), size = 4) +
    geom_text(aes(label = Trait), vjust = -0.8, size = 5) +
    scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8)) +  # Customize as needed
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 18),        # Main title size
      axis.title = element_text(size = 16),                       # Axis titles
      axis.text = element_text(size = 14),                        # Axis tick labels
      legend.title = element_text(size = 16),                     # Legend title
      legend.text = element_text(size = 14)                       # Legend labels
    ) +
    labs(
      title = "PCA of All Traits Across All Sites",
      x = paste0("PC1 (", round(summary(pca_traits)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_traits)$importance[2,2]*100, 1), "%)")
    )
  
  
  #Save
  ggsave("PCA_plot.jpeg", dpi = 300, width = 10, height = 7, units = "in")
  
  # To view proportion of variance
  summary(pca_traits)
  
  
  
  
  
}  # end disabled block

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
pc_scores$Specimen_ID <- maize_metadata$Specimen_ID

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
pc_scores$Site <- factor(pc_scores$Site, levels = c("M73", "M103", "M7", "M43", "M10", "M1", "M12", "M11", "M44"))  # FIX 22-Oct: duplicated levels



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
site_order <- c("M73", "M103", "M7", "M43", "M10", "M1", "M12", "M11", "M44")  # FIX 22-Oct: duplicated levels

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

# Chen chen v Omo style v M10 Temple
# Define Chen Chen, Omo, and Ritual Complex style groups
pc_scores_labeled <- pc_scores_labeled %>%
  mutate(Style_Group = case_when(
    # REVISION 22-Oct: two groups; M10 excluded from this comparison (decision 22 Oct)
    Site %in% c("M1", "M43") ~ "Chen Chen",
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
  labs(title = "QQ Plots for Chen Chen, Omo, and Ritual Complex Style Groups")

# Test for normality in each group
#shapiro_test <- pc_scores_grouped %>%
#  group_by(Style_Group) %>%
#  summarise(p_value = shapiro.test(PC1)$p.value)

#print(shapiro_test)

# Kruskal-Wallis test for non-normal data
# REVISION 22-Oct: two groups -> Wilcoxon rank-sum test with effect size r
wilcox_test_result <- wilcox.test(PC1 ~ Style_Group, data = pc_scores_grouped, exact = FALSE)
print(wilcox_test_result)
z_w <- qnorm(wilcox_test_result$p.value / 2)
wilcox_r <- abs(z_w) / sqrt(nrow(pc_scores_grouped))

# Pairwise Wilcoxon tests with p-value adjustment
pairwise_wilcoxon <- pairwise.wilcox.test(pc_scores_grouped$PC1, pc_scores_grouped$Style_Group, p.adjust.method = "BH")
print(pairwise_wilcoxon)

# SAVE


# FIX 22-Oct: wilcox_test_result and obs_diff were never created in this script, so the
# Wilcoxon and permutation sheets could not be regenerated. With three groups, report the
# Kruskal-Wallis test plus BH-adjusted pairwise Wilcoxon tests instead.
kw_df <- data.frame(
  W = unname(wilcox_test_result$statistic),
  p_value = wilcox_test_result$p.value,
  effect_size_r = wilcox_r,
  n_ChenChen = sum(pc_scores_grouped$Style_Group == "Chen Chen"),
  n_Omo = sum(pc_scores_grouped$Style_Group == "Omo")
)

pairwise_df <- as.data.frame(pairwise_wilcoxon$p.value)
pairwise_df$Comparison <- rownames(pairwise_df)

# Optional: reshape if you want long-format
# pairwise_long <- tidyr::pivot_longer(pairwise_df, -Comparison, names_to = "vs", values_to = "adjusted_p")

# ---- Write all to an Excel file ----
write_xlsx(
  list(
    "ChenChen_vs_Omo_Wilcoxon" = kw_df
  ),
  path = "Supplementary_ChenChen_vs_Omo.xlsx"
)


# AMS Dates ---------------------------------------------------

# FIX 22-Oct: SHCal20 table from the OxCal log. The old AMS_Dates_Final.csv had IntCal20 ranges,
# old sample IDs, and wrong Omo (M10) values.
ams_dates <- read_csv("Revision_Supplementary_AMS_full.csv") %>%
  filter(!is.na(`Specimen ID`)) %>%
  rename(Specimen_ID = `Specimen ID`,
         calibrated_2_sigma = `Calibrated 95.4% (intervals, probability)`) %>%
  mutate(age_BP = as.numeric(sub(" .*", "", `14C age (BP)`)),
         # overall 95.4% span (earliest to latest limit), used as a plot label
         span_95 = sapply(calibrated_2_sigma, function(x) {
           yrs <- as.numeric(unlist(regmatches(gsub("\\([^)]*\\)", "", x),
                                               gregexpr("[0-9]+", gsub("\\([^)]*\\)", "", x)))))
           paste0("AD ", min(yrs), "\u2013", max(yrs))
         }))


maize_dated <- maize_data %>%
  inner_join(ams_dates %>% select(Specimen_ID, `calibrated_2_sigma`), by = "Specimen_ID")


morphology_long <- maize_dated %>%
  dplyr::select(Specimen_ID, Period, calibrated_2_sigma, 
                Length_mm,
                Diameter_mm,
                Cupule_number,
                `Mean_Cupule_Width(mm)`,
                `Mean_cupule_height(mm)`,
                Mean_kernel_row,
                Total_Wt_g) %>%
  pivot_longer(cols = -c(Specimen_ID, Period, calibrated_2_sigma), 
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
               # FIX 22-Oct: names_sep = "_" split trait names like Mean_kernel_row_mean wrongly
               names_pattern = "(.*)_(n|mean|sd|min|max)$")

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
  select(Specimen_ID, cal_2sigma_clean,
         Length_mm, Diameter_mm, Cupule_number, `Mean_Cupule_Width(mm)`,
         `Mean_cupule_height(mm)`, Mean_kernel_row, Total_Wt_g) %>%
  pivot_longer(
    cols = -c(Specimen_ID, cal_2sigma_clean),
    names_to = "Trait",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value))


# FIX 22-Oct: the hard-coded date_lookup held old IntCal20 ranges (incl. the wrong Omo values).
# Labels now come from the SHCal20 table, ordered oldest to youngest by 14C age.
mh_dates <- ams_dates %>% filter(Specimen_ID %in% mh_long$Specimen_ID) %>% arrange(desc(age_BP))
date_lookup <- setNames(paste0(mh_dates$Specimen_ID, ": ", mh_dates$span_95), mh_dates$Specimen_ID)

mh_long <- mh_long %>%
  mutate(calibrated_2_sigma = date_lookup[Specimen_ID])

cal_2sigma_order <- unname(date_lookup)
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




# Period comparisons per trait (REVISION 22-Oct, replaces regressions and averaged contrasts) ----
# Periods are treated as unordered groups. All tests use the same complete-case dataset
# as the MANOVA (n = 344), so every n in the Results matches.
traits <- c("Length_mm", "Diameter_mm", "Cupule_number", "Mean_Cupule_Width(mm)",
            "Mean_cupule_height(mm)", "Mean_kernel_row", "Total_Wt_g")
cc_data <- maize_data[complete.cases(maize_data[, traits]), ] %>%
  mutate(Period = factor(as.character(Period), levels = c("LF", "LF/MH", "MH", "LIP")))

# Sample-size accounting (R2.5): all cobs vs complete cases, by period and site
n_accounting <- maize_data %>%
  mutate(Complete = complete.cases(maize_data[, traits])) %>%
  group_by(Period, Site) %>%
  summarise(N_all = n(), N_complete = sum(Complete), .groups = "drop") %>%
  arrange(factor(as.character(Period), levels = c("LF", "LF/MH", "MH", "LIP")), Site)
print(n_accounting, n = Inf)

hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp * (1 - 3 / (4 * (nx + ny) - 9))
}

omnibus <- list(); pairwise <- list()
for (tr in traits) {
  f <- as.formula(paste0("`", tr, "` ~ Period"))
  w <- oneway.test(f, data = cc_data, var.equal = FALSE)          # Welch ANOVA
  a <- anova(lm(f, data = cc_data))
  eta2 <- a$`Sum Sq`[1] / sum(a$`Sum Sq`)
  omnibus[[tr]] <- data.frame(Trait = tr, F_Welch = unname(w$statistic),
                              df1 = unname(w$parameter[1]), df2 = unname(w$parameter[2]),
                              p = w$p.value, eta_squared = eta2, n = nrow(cc_data))
  d <- data.frame(Value = cc_data[[tr]], Period = cc_data$Period)
  gh <- rstatix::games_howell_test(d, Value ~ Period, conf.level = 0.95)
  gh$Trait <- tr
  gh$n1 <- as.vector(table(d$Period)[gh$group1]); gh$n2 <- as.vector(table(d$Period)[gh$group2])
  gh$hedges_g <- mapply(function(g1, g2) hedges_g(d$Value[d$Period == g2], d$Value[d$Period == g1]),
                        gh$group1, gh$group2)
  pairwise[[tr]] <- gh
}
omnibus_df <- bind_rows(omnibus) %>% mutate(p_BH = p.adjust(p, method = "BH"))
pairwise_df <- bind_rows(pairwise) %>%
  dplyr::select(Trait, group1, group2, n1, n2, estimate, conf.low, conf.high, p.adj, hedges_g)
print(omnibus_df)

write_xlsx(list(Sample_sizes = n_accounting,
                Welch_ANOVA_by_trait = omnibus_df,
                Pairwise_GamesHowell = pairwise_df,
                Notes = data.frame(Note = c(
                  "Periods treated as unordered groups; complete cases only (same cobs as the MANOVA).",
                  "Welch ANOVA per trait; p_BH = Benjamini-Hochberg adjusted across the seven traits.",
                  "Pairwise: Games-Howell between periods within each trait; estimate = mean(group2) - mean(group1) with 95% CI.",
                  "hedges_g = standardized difference, group2 minus group1.",
                  "eta_squared from the standard one-way ANOVA of the same data."))),
           path = "Supplementary_Period_Comparisons.xlsx")

# Figure 4 (replaces regression panels): traits by period, unordered groups
fig4 <- cc_data %>%
  dplyr::select(Period, all_of(traits)) %>%
  pivot_longer(-Period, names_to = "Trait", values_to = "Value") %>%
  ggplot(aes(Period, Value)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
  facet_wrap(~ Trait, scales = "free_y") +
  labs(x = NULL, y = NULL, title = paste0("Cob traits by period (n = ", nrow(cc_data), ")")) +
  theme_bw()
print(fig4)
ggsave("Figure4_traits_by_period.png", fig4, width = 9, height = 7, dpi = 300)