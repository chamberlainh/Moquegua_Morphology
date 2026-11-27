# Moquegua Valley Ancient Maize: Cob Morphology and AMS Radiocarbon Dates

Data and analysis code for:

> Chamberlain, H., et al. (in review). Shape Shifting: Morphological Variation in Ancient Maize (*Zea mays* ssp. *mays*) from the Moquegua Valley, Peru (758 BC–AD 1210). *Journal of Ethnobiology*. [DOI when published]

The study compares cob traits of 344 archaeological maize cobs from nine sites in the Moquegua Valley, southern Peru, spanning the Late Formative (LF), the LF/Middle Horizon transition (LF/MH), the Middle Horizon (MH), and the early Late Intermediate Period (LIP). Twenty cobs were directly dated by AMS radiocarbon. A subset of 36 cobs was analyzed with geometric morphometrics of 3D-scan outlines.

## Repository contents

| Folder | Contents |
|---|---|
| `Morphological_Analysis/` | Datasheets, the analysis script, the R project file, and the geometric morphometrics landmark file |
| `Supplementary/` | Supplementary tables and images cited in the paper |
| `Figures/` | Final figures as they appear in the paper |
| `AMS/` | Lab report, OxCal calibration and models for the 20 radiocarbon dates |
| `Visuals/` | Earlier exploratory plots, kept for transparency. **Not used in the paper**; several show analyses that were replaced during revision |

Each folder has its own README.

## Reproducing the analysis

1. Open `Morphological_Analysis/Moquegua_Morphology.Rproj` in RStudio (this sets the working directory).
2. Run `Morphological_Analysis/Maize_analysis_v3_2026-11-27.R` from top to bottom.
3. The script writes the supplementary tables (`.xlsx`) and Figures 3, 4, 5, 7 and 8 (`.png`) to the same folder.

Analyses were run in R 4.4.2. Packages: readr, dplyr, tidyr, ggplot2, ggpubr, ggrepel, rstatix, car, emmeans, candisc, broom, writexl, multcomp, viridis, RColorBrewer.

Figure 6 (outline shape, geometric morphometrics) was made in the EVAN Toolbox from `Morphological_Analysis/Maize_Morphologika.txt`; see the Methods of the paper.

## Sample sizes

- 344 cobs are analyzed: all cobs with complete data for the seven traits (288 measured directly, 56 from photographs and 3D models). The datasheet holds 358 cobs; 14 second-phase cobs lack a cupule count and are excluded from all trait analyses.
- 36 cobs (9 per period) were used for geometric morphometrics.
- 20 cobs were AMS dated.

## 3D models

3D surface models of the second-phase cobs: https://webserver.anthro.univie.ac.at/heascloud/index.php/f/31564. A video of the outline-plane procedure: https://doi.org/10.5281/zenodo.23500772.

## Citation

Please cite the paper above. To cite this repository, use its Zenodo DOI: [add after creating the Zenodo release].

## License

Code: MIT License. Data, tables and figures: CC BY 4.0. See `LICENSE`.

## Contact

Heather Chamberlain-Irwin, University of Vienna. heather.chamberlain.irwin@univie.ac.at
