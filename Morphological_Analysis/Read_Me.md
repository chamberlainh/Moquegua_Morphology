# Morphological Analysis

## Overview

This folder contains R project, R script, and data files for morphological analysis. For 3D models, see the repository https://webserver.anthro.univie.ac.at/heascloud/index.php/f/31564.

## File Contents

### Data Files

- `Merged_Ancient_Maize_Cob_Data_Sheet.csv` – Raw trait measurements for all specimens from Permit 1 and Permit 2.

- `cleaned_Ancient_Maize_Cob_Data_Sheet.csv` – Raw trait measurements for speicmens from Permit 1.

- `cleaned_Ancient_Maize_Cob_Data_Sheet.xlsx` – Raw trait measurements for speicmens from Permit 1.

### Generated Exploration Files
  
- `CDA_M12_M73_Traits.xlsx` - Formatted data file for exploration of CDA results for sites M12 and M73.

- `Maize_Period_Contrasts.xlsx` - Formatted data file for exploration of changes over time periods.
  
- `Style_PC1_Comparison_Results.xlsx` - Formatted data file for exploration of PC1 results.

- `Summary_Traits_By_Period.xlsx` - Formatted data file for summary traits by chronology.

- `Traits_Summary_AMS.xlsx` - Formatted data file for summary of traits for only AMS radiocarbon dated samples.

- `TukeyHSD_M7_Transition.xlsx` - Formatted data file for exploration of trait changes at M7.

### R Scripts

- `Morph_Stats.R` – R script for Permit 1 specimens statistical analysis and visualization of morphological measurements.
  
- `All_Morph_Stats.R` – R script for Permit 1 and Permit 2 specimens statistical analysis and visualization of morphological measurements.

### Plots & Visualizations

- See generated visuals in Visuals folder

### Project Files

- `Moquegua_Morphology.Rproj` – R project file to organize and manage the analysis.

- `.gitignore` – Specifies files and directories to be ignored by Git.

- `README.md` – This document.

## Usage
### Permit 1 specimens only

- Load the dataset into R using `read.csv("cleaned_Ancient_Maize_Cob_Data_Sheet.csv")`.

- Run `Morph_Stats.R` to generate descriptive statistics and visualizations.

- Check the output .png files for graphical summaries of the data.
  
### All specimens from Permit 1 and Permit 2
- Load the dataset into R using `read.csv("Merged_Ancient_Maize_Cob_Data_Sheet.csv")`.

- Run `All_Morph_Stats.R` to generate descriptive statistics and visualizations.

- Check the output .png and .xlsx files for summaries of the data.

