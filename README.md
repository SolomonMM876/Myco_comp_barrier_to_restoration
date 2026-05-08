# Mycorrhizal Community Composition as a Barrier to Ecological Restoration

Reproducible analysis pipeline comparing mycorrhizal fungal community
recovery at two Australian mine rehabilitation case studies: **Hunter Mining**
(Hunter Valley, NSW) and **Snowy Hydro** (Snowy Mountains, NSW).

## Repository Structure

```
Myco_comp_barrier_to_restoration/
├── analysis.Rmd              # Master RMarkdown — knit to reproduce all outputs
├── R_scripts/                
│   ├── Site_maps.R           # Site map figures
│   └── hmsc/                 # HMSC model preparation & post-processing
│       ├── dataPrep_hmsc_hunter_reference.R
│       ├── dataPrep_hmsc_hunter_rehab.R
│       ├── dataPrep_hmsc_snowies_reference.R
│       ├── dataPrep_hmsc_snowies_rehab.R
│       └── functions_hmsc.R
├── rawdata/
│   ├── ITS/                  # ITS amplicon sequencing data — NOT tracked by git (~40–50 MB each)
│   ├── metadata/
│   │   ├── hunter/           # Hunter Mining site metadata
│   │   └── snowies/          # Snowy Hydro site metadata
├── models/                   # Pre-fitted HMSC & CCA model objects (.RData) — NOT tracked by git
│   ├── hunter/
│   │   ├── cca/              # ordistep_reference.RData, ordistep_rehab.RData
│   │   ├── hmsc_reference/   # mcmc_input/output/convergence.RData
│   │   └── hmsc_rehab/
│   └── snowies/
│       ├── cca/
│       ├── hmsc_reference/
│       └── hmsc_rehab/
├── deriveddata/              # Cached processed data (.rds) — auto-generated, not tracked
└── output/
    ├── plots/                # All publication figures (auto-generated, not tracked)
    └── tables/               # All statistics tables (auto-generated, not tracked)
```

## Quick Start

1. **Open** the R project: `Myco_comp_barrier_to_restoration.Rproj`
2. **Install dependencies** (see below)
3. **Obtain model files**: Place pre-fitted HMSC and CCA model `.RData`
   files in the `models/` directory (see structure above). These are not
   tracked by git due to size. Alternatively, uncomment the relevant
   `source()` calls in the `hmsc-run-models` chunk to refit from scratch
   (takes several hours per model).
4. **Knit** `analysis.Rmd` — this generates all figures and tables in `output/`

### First Run

On first run, set these flags at the top of `analysis.Rmd`:

```r
RELOAD_FROM_RAW  <- TRUE   # Process raw data (slow, ~2 min)
RERUN_INDICATORS <- FALSE  # Use pre-generated indicator taxa table
RECOMPUTE_HMSC   <- TRUE   # Compute HMSC associations from model objects
```

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

- `analysis.Rmd` is the primary reproducible entry point. Most analyses
  are coded inline. One folder is sourced externally:
  - `R_scripts/hmsc/dataPrep_hmsc_*.R` — HMSC data prep & model fitting
    (sourced only if the `hmsc-run-models` chunk is uncommented)
- `R_scripts/Site_maps.R` standalone script retained for to generate manuscript plots
- All paths in the `R_scripts/hmsc/` scripts resolve relative to `proj_root`
  (set automatically via `here::here()` when sourced from the Rmd).
- HMSC models are **never re-fit** during a normal knit — only pre-saved
  model objects are loaded and post-processed.
- Raw ITS sequencing data (`rawdata/ITS/`) is **not tracked by git** due to
  file size (~40–50 MB per file).

## License


