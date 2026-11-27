# Repository cleanup before submission

Work through this once, then commit everything in one go (e.g. "Clean repository for submission").
Paths are relative to the repository root. "Archive" = move to `Visuals/Exploratory_archive/`.

## 1. Add / replace

- [ ] `README.md`, `LICENSE` (root): replace/add with the new versions
- [ ] `Morphological_Analysis/Read_Me.md`, `Supplementary/Read_Me.md`, `AMS/Read_Me.md`, `Visuals/Read_Me.md`: replace with the new versions
- [ ] New folder `Figures/` with its `Read_Me.md`
- [ ] `Morphological_Analysis/Maize_analysis_v3_2026-11-27.R`: add (the current script; no `setwd()`)
- [ ] Run the script once from a fresh R session (Session → Restart R) so every output below is from the same run

## 2. Figures/ (copy the final versions here)

- [ ] `Figure3_PCA_biplot.png`, `Figure4_traits_by_period.png`, `Figure5_CDA_Traits_Sites.png`, `Figure7_Cob_Dimensions_MH.png`, `Figure8_Cupule_Characteristics_MH.png` (from `Morphological_Analysis/`)
- [ ] Figure 6 from Petra (name it `Figure6_GM_outline_PCA.png` or similar and update `Figures/Read_Me.md`)

## 3. Supplementary/ (one copy of each, from the final run)

- [ ] Replace with the new run's versions: `Supplementary_Morphology_Stats.xlsx`, `Supplementary_MH_Traits_with_AllSites.xlsx`, `Supplementary_LIP_Traits_with_AllSites.xlsx`
- [ ] Add from the new run: `Supplementary_PCA.xlsx`, `Supplementary_Period_MANOVA.xlsx`, `Supplementary_Period_Comparisons.xlsx`, `Supplementary_CDA.xlsx`, `Supplementary_ChenChen_vs_Omo.xlsx`
- [ ] Add: `Revision_Supplementary_AMS_full.csv` (copy; keep the original in `Morphological_Analysis/`, the script reads it there)
- [ ] Keep: `Supplementary_3D_Model_Key.xlsx`, `Supplementary_Batch2_Cob_Data_and_Conversion_Check.xlsx`, `Supplementary_Plane_outline.png`
- [ ] Delete (replaced during revision): `Supplementary_Period_Contrasts.xlsx` (Dunnett contrasts), `Supplmentary_Omo_ChenChen_M10.xlsx` (not produced by the script)

## 4. Morphological_Analysis/

Keep:
- [ ] `Merged_Ancient_Maize_Cob_Data_Sheet_with_IDs.csv`, `cleaned_Ancient_Maize_Cob_Data_Sheet.csv`, `cleaned_Ancient_Maize_Cob_Data_Sheet.xlsx`, `Revision_Supplementary_AMS_full.csv`, `Maize_Morphologika.txt`, `Moquegua_Morphology.Rproj`, `Maize_analysis_v3_2026-11-27.R`
- [ ] Outputs of the final run (the `.xlsx` and Figure `.png` files the script writes)

Delete:
- [ ] `statistical-shape-analysis-with-applications-in-r-2ed.pdf` (broken 2-byte file; copyrighted book)
- [ ] `Merged_Ancient_Maize_Cob_Data_Sheet.csv` (old datasheet without IDs or corrections)
- [ ] Old figure copies: `Figure7_Cob_Dimensions_MH.jpeg`, `Figure8_Cupule_Characteristics_MH.jpeg`, `Figure_5_CDA_Traits_Sites.jpg`, `Figure5_CDA.png`, `CDA_Traits_Sites.jpg`, `CDA_Traits_Across_Sites.jpeg`, `Cob_Dimensions_MH.jpeg`, `Cupule_Characteristics_MH.jpeg`, `AMS_MH_Line` (no extension; duplicate of `AMS_MH_Line.png`)
- [ ] Replaced analyses: `PCA_plot.jpeg` (traits-as-observations PCA), `Linear_Regression_Cob_Dimesions.jpeg`, `Linear_Regression_Cupule_Characteristics.jpeg`, `Maize_Period_Contrasts.xlsx`, `Supplmentary_Period_Contrasts.xlsx`, `Supplmentary_Omo_ChenChen_M10.xlsx`
- [ ] Old scripts: `All_Morph_Stats.R`, `Morph_Stats.R`, `Revision_All_Morph_Stats.R` (or keep only if you rename the new script to this name instead; one script in the end)

Optional (exploratory outputs, not cited): `AMS_MH_Line.png`, `Style_PC1_Comparison_Results.xlsx` can be archived or deleted. `TukeyHSD_M7_Transition.xlsx`, `CDA_M12_M73_Traits.xlsx`, `Summary_Traits_By_Period.xlsx`, `Trait_Summary_AMS.xlsx` are still written by the script; keep them (listed in the README as "not cited").

## 5. Visuals/

- [ ] Move everything except `Read_Me.md` into `Visuals/Exploratory_archive/` (or delete the folder's contents if you prefer a lean repository)
- [ ] Files without an extension (`Biplot`, `CDA_All_Traits`, `CDA_No_Weight`, `Linear_Cob`, `Linear_Cupule`, `MH_Trait_Sites`, `Omo_v_Chen-Chen`, `Omo_v_Chen-Chen_v_Ritual`, `PCA_Site`, `PCA_Site_Period_Trait`, `PCA_Trait`, `PC_Cob_Size_*`, `Boxplot_Traits_MH_Time`) and `PCA.bin`: delete, or add `.png` if they are images you want to keep

## 6. AMS/

- [ ] Rerun the Rio Muerto model with specimen IDs (open item), or keep the note in `AMS/Read_Me.md`

## 7. After committing

- [ ] Open the repository in a private browser window and check the README renders and links work
- [ ] Create a GitHub release, then archive it on Zenodo (GitHub → Zenodo integration) to get a DOI; add the DOI to the README and to the paper's Data Availability statement
