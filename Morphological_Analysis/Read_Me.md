# Morphological Analysis

Datasheets, analysis script and R project for the cob-trait analyses in the paper.

## How to run

Open `Moquegua_Morphology.Rproj`, then run `Maize_analysis_v3_2026-11-27.R` from top to bottom. The script reads the three data files below and writes every supplementary table and Figures 3, 4, 5, 7 and 8 to this folder. Plots that were explored but are not in the paper are switched off; set `SHOW_EXPLORATORY <- TRUE` near the top of the script to draw them.

## Data files

| File | Contents |
|---|---|
| `Merged_Ancient_Maize_Cob_Data_Sheet_with_IDs.csv` | **Main datasheet.** One row per cob (358 cobs; 344 with complete data). Columns are described below |
| `cleaned_Ancient_Maize_Cob_Data_Sheet.csv` / `.xlsx` | First-phase cobs (288, measured directly) with the three separate cupule measurements per cob and cupule depth |
| `Revision_Supplementary_AMS_full.csv` | The 20 AMS dates, calibrated with SHCal20, joined to specimen IDs |
| `Maize_Morphologika.txt` | 20 pseudo-landmarks per cob for the 36 cobs in the geometric morphometric analysis (Morphologika format). Cobs are listed by their original codes (ZMDNA24_xx); `Supplementary/Supplementary_3D_Model_Key.xlsx` links these to the specimen IDs |

## Datasheet columns (`Merged_Ancient_Maize_Cob_Data_Sheet_with_IDs.csv`)

| Column | Meaning |
|---|---|
| `Specimen_ID` | Specimen ID used in the paper: Site_Lot_Number (e.g. M43_6261_03) |
| `Lot` | Excavation lot (provenience unit) |
| `Old_Sample_ID` | Original museum or paleogenomic code (e.g. ZMDNA24_30) |
| `UCIAMS` | Lab number for AMS-dated cobs |
| `Site` | Site code (M73 Tres Quebradas, M103 Montalvo, M7 Cerro Trapiche, M43 Rio Muerto, M10 Omo, M1 Chen Chen, M12 Omo, M11 Omo, M44 Cerro Chamos) |
| `Period` | LF, LF/MH, MH or LIP (assigned from archaeological context) |
| `Phase/Project`, `Phase` | Cultural affiliation and collection phase as recorded |
| `Cob_shape` | Qualitative outline shape |
| `Length_mm`, `Diameter_mm` | Cob length and midpoint diameter, caliper (mm) |
| `Cupule_number` | Number of cupules counted along the length of the cob |
| `Mean_Cupule_Width(mm)`, `Mean_cupule_height(mm)` | Mean of three measurements (top, middle, bottom), mm |
| `Mean_kernel_row` | **Kernel row number**: mean number of cupule (kernel) rows around the cob at top, middle and bottom. Not kernels per row |
| `Total_Wt_g` | Cob weight (g) |

First-phase cobs were measured directly; second-phase cupule traits were measured from images and converted with 72.4 pixels per cm (see `Supplementary/Supplementary_Batch2_Cob_Data_and_Conversion_Check.xlsx`). Three data-entry errors were corrected during revision (M43_6261_03, M103_1161_30, M10_12478_11).

## Outputs written by the script

Supplementary tables cited in the paper (copies in `Supplementary/`): `Supplementary_Morphology_Stats.xlsx`, `Supplementary_PCA.xlsx`, `Supplementary_Period_MANOVA.xlsx`, `Supplementary_Period_Comparisons.xlsx`, `Supplementary_CDA.xlsx`, `Supplementary_MH_Traits_with_AllSites.xlsx`, `Supplementary_LIP_Traits_with_AllSites.xlsx`, `Supplementary_ChenChen_vs_Omo.xlsx`, `Table2_Summary_by_Period_344.xlsx`.

Figures (copies in `Figures/`): `Figure3_PCA_biplot.png`, `Figure4_traits_by_period.png`, `Figure5_CDA_Traits_Sites.png`, `Figure7_Cob_Dimensions_MH.png`, `Figure8_Cupule_Characteristics_MH.png`.

Additional outputs, not cited in the paper: `TukeyHSD_M7_Transition.xlsx`, `CDA_M12_M73_Traits.xlsx`, `Summary_Traits_By_Period.xlsx`, `Trait_Summary_AMS.xlsx`.
