# PARKS_CES

This repository contains scripts and data used to model, extrapolate, map, and analyze cultural ecosystem services (CES) supplied by urban parks in São Paulo, Brazil.

The project focuses on two cultural ecosystem services:

1. **Potential Mental Restoration (PMR)** — the potential of urban parks to promote perceived psychological restoration.
2. **Potential Active Recreation (PAR)** — the potential of urban parks to support active recreational use.

The workflow combines park-user interview data, park landscape metrics, park-quality indicators, distance-decay functions, and spatial demand surfaces to evaluate the spatial distribution of CES supply, flow, demand, and mismatch across São Paulo.

## Project overview

Urban parks are key spaces for everyday human–nature contact in cities. However, cultural ecosystem services are difficult to quantify because they depend not only on the amount of green space, but also on park quality, landscape structure, accessibility, human perception, and actual use.

This project uses empirical information from park users and park-level attributes to:

- estimate perceived mental restoration and active recreation in sampled parks;
- extrapolate these services to a broader set of urban parks in São Paulo;
- apply distance-decay functions to estimate spatial service flow;
- build demand surfaces for each service;
- identify areas where demand exceeds modeled supply;
- detect spatial clusters of CES mismatch using Moran’s I and LISA analyses.

## Repository contents

### Scripts

| File | Description |
|---|---|
| `PAR_PMR_selecting_variables.R` | Screens candidate predictor variables for PMR and PAR models and compares alternative model structures. |
| `Model_extrapol_val_PMR.R` | Models Potential Mental Restoration, extrapolates PMR values to the full park dataset, and validates predictions. |
| `Model_extrapol_val_PAR.R` | Models Potential Active Recreation, extrapolates PAR values to the full park dataset, and validates predictions. |
| `mismatch_analysis.R` | Computes supply–demand mismatch surfaces for PMR and PAR using CES supply-flow and demand variables. |
| `LISA+Morans_mismatch.R` | Runs spatial autocorrelation analyses and identifies local clusters of CES mismatch using Moran’s I and LISA. |
| `biplot_script.R` | Produces bivariate supply–demand maps for PMR and PAR. |

### Data files

| File | Description |
|---|---|
| `CES_parks_final.csv` | Main interview-level dataset combining park-user information, perceived restoration, activity classification, park identifiers, park-quality indicators, and landscape metrics. |
| `all_parks_new.csv` | Park-level dataset used for extrapolating PMR and PAR to the broader set of São Paulo urban parks. |
| `parks_CES_FINAL.shp` and associated files | Shapefile containing park geometries and CES-related attributes. |
| `grid_data_CES.shp` and associated files | Spatial grid used for mapping CES supply-flow, demand, and mismatch surfaces. Some associated files may be too large for regular GitHub upload. |
| `CES_parks_final.qmd` | Quarto document associated with the CES parks analysis. |
| `parks_CES_FINAL.qmd` | Quarto document associated with the final parks CES dataset. |

## Analytical workflow

### 1. Select predictor variables

Run:

```r
source("PAR_PMR_selecting_variables.R")
```

This script evaluates candidate park-level predictors for the PMR and PAR models. It helps identify which landscape and park-quality variables are most relevant for explaining the park-specific effects extracted from the empirical interview models.

### 2. Model and extrapolate Potential Mental Restoration

Run:

```r
source("Model_extrapol_val_PMR.R")
```

This script models perceived restoration based on park-user data. It controls for individual-level characteristics, extracts park-level effects, relates these effects to park landscape and quality attributes, and extrapolates PMR values to the broader set of parks.

The script also includes a validation step to evaluate the predictive performance of the extrapolation.

### 3. Model and extrapolate Potential Active Recreation

Run:

```r
source("Model_extrapol_val_PAR.R")
```

This script models the probability of active recreational use in parks. Activities reported by users are classified into active recreation versus passive or non-recreational use.

The model is then used to extrapolate PAR values to the broader park dataset, followed by validation of the extrapolated predictions.

### 4. Compute supply–demand mismatch

Run:

```r
source("mismatch_analysis.R")
```

This script uses the spatial grid to combine CES supply-flow surfaces with demand surfaces.

Mismatch is calculated as:

```text
Mismatch = standardized demand - standardized supply
```

Positive mismatch values indicate areas where demand exceeds modeled supply. Negative values indicate areas where modeled supply is higher than demand.

### 5. Run Moran’s I and LISA analyses

Run:

```r
source("LISA+Morans_mismatch.R")
```

This script evaluates spatial autocorrelation in CES mismatch and identifies local clusters of supply–demand imbalance.

The analysis helps detect areas where high demand and low supply are spatially concentrated.

### 6. Produce bivariate supply–demand maps

Run:

```r
source("biplot_script.R")
```

This script generates bivariate maps combining supply and demand classes for each CES. These maps help visualize areas with different combinations of low/high supply and low/high demand.

## Required R packages

The scripts use several R packages for spatial analysis, modeling, validation, and visualization.

Install the main required packages with:

```r
install.packages(c(
  "sf",
  "terra",
  "spdep",
  "ggplot2",
  "dplyr",
  "scales",
  "lme4",
  "lmerTest",
  "car",
  "cvTools",
  "ggrepel",
  "gridExtra"
))
```

Additional packages may be required depending on the specific script and local setup.

## How to use this repository

1. Download or clone this repository.

```bash
git clone https://github.com/DougCirino/PARKS_CES.git
```

2. Open R or RStudio in the repository folder.

3. Make sure the data files are in the same folder expected by the scripts.

4. Install the required R packages.

5. Run the scripts in the following order:

```r
source("PAR_PMR_selecting_variables.R")
source("Model_extrapol_val_PMR.R")
source("Model_extrapol_val_PAR.R")
source("mismatch_analysis.R")
source("LISA+Morans_mismatch.R")
source("biplot_script.R")
```

Depending on the machine and on the availability of large spatial files, some scripts may require local adjustment of file paths.

## Data availability and large files

Some spatial files associated with this project may be too large for regular GitHub upload, especially shapefile components and raster/grid outputs.

If a file is missing from the GitHub repository because of size limits, users should request access to the large spatial files or obtain them from the external data repository indicated by the author.

Large files may include:

- `grid_data_CES.dbf`
- large shapefile components;
- raster outputs;
- intermediate spatial products.

The main scripts and smaller processed datasets are provided in this repository to document and reproduce the analytical workflow.

## Methodological summary

The project follows a supply–flow–demand logic for cultural ecosystem services.

First, empirical park-user data are used to estimate two CES indicators: Potential Mental Restoration and Potential Active Recreation. These indicators are modeled using interview data and then related to park-level landscape and quality variables.

Second, the fitted relationships are used to extrapolate CES values to a broader set of urban parks in São Paulo.

Third, distance-decay functions are applied to represent the spatial flow of CES from parks to surrounding urban areas. Larger or more attractive parks may influence broader areas, while smaller parks have more spatially restricted effects.

Fourth, demand surfaces are constructed to represent areas where the need for CES is higher, considering population and socioeconomic vulnerability.

Finally, supply-flow and demand are compared to identify mismatch areas, followed by spatial autocorrelation analyses to detect significant clusters of imbalance.

## Suggested citation

Cirino, D. W., Felappi, J. F., Lupinetti-Cunha, A., Tarragô, G. M., Carrasco, L. R., & Metzger, J. P. . *Bridging perception and landscape structure: mapping urban streetscape aesthetics as nature’s contributions to people*.

## Contact

Douglas W. Cirino  
University of São Paulo  
Institute of Biosciences  
São Paulo, Brazil
