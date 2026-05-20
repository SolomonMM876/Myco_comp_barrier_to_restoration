# Mycorrhizal Community Composition as a Barrier to Ecological Restoration

Reproducible analysis pipeline comparing mycorrhizal fungal community
recovery at two Australian mine rehabilitation case studies: **Hunter Mining**
(Hunter Valley, NSW) and **Snowy Hydro** (Snowy Mountains, NSW).

## Repository Structure

```
Myco_comp_barrier_to_restoration/
├── Myco_Comp_Restoration_Barrier.Rmd   # Master RMarkdown — knit to reproduce all outputs
├── Myco_comp_barrier_to_restoration.Rproj
├── R_scripts/
│   ├── Site_maps.R                     # Site map figures
│   └── hmsc/                           # HMSC model preparation & post-processing
│       ├── dataPrep_hmsc_hunter_reference.R
│       ├── dataPrep_hmsc_hunter_rehab.R
│       ├── dataPrep_hmsc_snowies_reference.R
│       ├── dataPrep_hmsc_snowies_rehab.R
│       └── functions_hmsc.R
├── rawdata/
│   ├── ITS/                            # ITS amplicon sequencing data (~40–50 MB each)
│   │   ├── hunter_ITS_long_20210604.txt
│   │   └── snowies_ITS_long_20210604.txt
│   └── metadata/
│       ├── hunter/                     # Hunter Mining site metadata (CSVs + XLSX)
│       └── snowies/                    # Snowy Hydro site metadata (CSVs)
├── models/                             # Pre-fitted HMSC & CCA model objects (.RData)
│   ├── hunter/
│   │   ├── cca/                        # ordistep_reference.RData, ordistep_rehab.RData
│   │   ├── hmsc_reference/             # mcmc_output.RData, mcmc_convergence.RData
│   │   └── hmsc_rehab/
│   └── snowies/
│       ├── cca/
│       ├── hmsc_reference/
│       └── hmsc_rehab/
├── deriveddata/                        # Cached processed data (.rds)
├── submission_files/                   # Bundled data archive for reproducibility
│   ├── rawdata/                        #   (mirrors rawdata/, deriveddata/, models/)
│   ├── deriveddata/
│   └── models/
└── output/
    ├── plots/                          # All publication figures (auto-generated)
    └── tables/                         # All statistics tables (auto-generated)
```

## Reproducing from `submission_files`

The `submission_files/` folder (or its zipped archive) contains **all data
required to knit the RMarkdown**, bundled into three subfolders:

| Folder | Contents | Size |
|--------|----------|------|
| `rawdata/` | ITS amplicon tables and site metadata CSVs | ~85 MB |
| `deriveddata/` | Pre-processed RDS caches (genus-level OTU tables, predictors, etc.) | ~4 MB |
| `models/` | Pre-fitted HMSC MCMC outputs and CCA ordistep `.RData` files | ~1.5 GB |

 **Unzip** the `submission_files` archive (or locate the `submission_files/`
   folder) so you have the three subfolders: `rawdata/`, `deriveddata/`, and
   `models/`.

 **Copy the three folders into the project root** — so your directory looks like:
   ```
   Myco_comp_barrier_to_restoration/
   ├── Myco_Comp_Restoration_Barrier.Rmd
   ├── Myco_comp_barrier_to_restoration.Rproj
   ├── R_scripts/
   ├── rawdata/          ← from submission_files
   ├── deriveddata/      ← from submission_files
   └── models/           ← from submission_files
   ```

**Install R package dependencies** 

**Set the reproducibility flags** at the top of
   `Myco_Comp_Restoration_Barrier.Rmd`:
   ```r
   RELOAD_FROM_RAW  <- TRUE   # Process raw data from rawdata/ (~2 min)
   RERUN_INDICATORS <- TRUE   # Re-run indicator taxa analysis
   RECOMPUTE_HMSC   <- TRUE   # Compute HMSC associations from model objects
   RERUN_CCA        <- FALSE  # Load pre-fitted CCA models (set TRUE to refit)
   ```

**Knit** `Myco_Comp_Restoration_Barrier.Rmd` — all figures and tables are generated
 in `output/plots/` and `output/tables/`.


After the first run, processed data is cached in `deriveddata/`. Switch
`RELOAD_FROM_RAW <- FALSE` for fast subsequent knits.

## Dependencies

### R Packages

```r
install.packages(c(
  "tidyverse", "sf", "patchwork", "lme4", "car", "broom.mixed",
  "vegan", "fossil", "Hmsc", "coda", "circlize",
  "flextable", "officer", "tmap", "maptiles",
  "indicspecies", "permute", "here",
  "htmlwidgets", "base64enc", "ggrepel"
))

# chorddiag is not on CRAN — install from GitHub:
# remotes::install_github("mattflor/chorddiag")
```

## Analyses Included

| Section | Description | Key Output |
|---------|-------------|------------|
| Study Site Maps | Interactive & static satellite maps of sampling locations | `Figure_Map_*.png` |
| Soil Characteristics | Soil abiotic comparisons (TOC, P, PLFA) | `Figure_Soil_*.png` |
| Environmental PCA | PCA biplots of environmental predictors | `Figure_PCA_Combined.png` |
| Mycorrhizal Communities | EcM/AM richness & abundance boxplots | `Figure_Mycorrhizal_*.png` |
| Guild Relationships | EcM vs AM scatter plots | `Figure_Mycorrhizal_Relationships_*.png` |
| PLFA–DNA Comparison | ITS vs PLFA validation | `plfa_dna_comparison_*.png` |
| Indicator Taxa | IndVal analysis (multipatt) | `Table_Indicator_Taxa_*.docx` |
| CCA | Constrained correspondence analysis | `publication_cca_combined.png` |
| HMSC Model Fitting | Optional: refit HMSC models by sourcing `R_scripts/hmsc/dataPrep_*.R` | model `.RData` files |
| HMSC Diagnostics | Model convergence & fit | `Diagnostics_*.pdf` |
| HMSC Evaluation | Summary metrics table | `Table_HMSC_Model_Evaluation.docx` |
| HMSC Chord Diagrams | Static co-occurrence chord figures | `publication_hmsc_chord_*.png` |
| HMSC Interactive Chords | Interactive HTML chord diagrams (positive/negative side-by-side) | `publication_hmsc_chord_interactive.html` |
| Genera vs Site Recovery | Crosshair recovery plots | `Figure_EcM_Master_SiteRecovery.png` |

## Statistical Methods

- **Linear Mixed-Effects Models** (`lme4`): Treatment effects tested with
  Type II Wald F-tests (`car::Anova`). Random effects:
  - Hunter Mining: `(1 | Plot_ID)`
  - Snowy Hydro: `(1 | Plot_ID / Subplot)`
- **HMSC**: Hierarchical Models of Species Communities for joint species
  distribution modelling (presence-absence: probit; abundance: normal).
- **CCA**: Constrained Correspondence Analysis with bidirectional stepwise variable selection (`vegan::ordistep`).
- **Indicator Taxa**: IndVal analysis via `indicspecies::multipatt`.

## Notes

- `Myco_Comp_Restoration_Barrier.Rmd` is the primary reproducible entry point.
  Most analyses are coded inline. One folder is sourced externally:
  - `R_scripts/hmsc/dataPrep_hmsc_*.R` — HMSC data prep & model fitting
    (sourced only if the `hmsc-run-models` chunk is uncommented)
- `R_scripts/Site_maps.R` standalone script retained to generate manuscript plots
- All paths in the `R_scripts/hmsc/` scripts resolve relative to `proj_root`
  (set automatically via `here::here()` when sourced from the Rmd).
- HMSC models are **never re-fit** during a normal knit — only pre-saved
  model objects are loaded and post-processed.
- Raw ITS sequencing data (`rawdata/ITS/`) is **not tracked by git** due to
  file size (~40–50 MB per file).

## License


