# shinyForce3D: 3D Visualisation of Ground Reaction Forces (GRF)

This repository contains an R/Shiny application and supporting R package for visualising ground reaction force (GRF) vector data in dynamic 3D, static 3D, and 2D formats. The tool supports single and dual dataset comparisons, customisable sphere/arrow sizing, peak detection, and export of animated GIFs.

## Features

- **Interactive 3D**: Real‑time rotation and exploration of GRF vectors as colored spheres or arrows.
- **GIF Export**: Generate rotating 3D animations with all annotations.
- **Static 3D**: Publication‑quality snapshots from orthogonal views, with optional butterfly lines and peak vector highlighting.
- **2D Time-Series (Cartesian)**: Plot Fx, Fy, Fz, and resultant force over stance phase, with mean ± SD ribbons and peak markers.
- **2D Polar Coordinates**: Visualise magnitude (r), inclination (θ), and azimuth (φ) to distinguish directional vs. magnitude changes.
- **Pedotti Plot**: Animated or static depiction of two selected GRF components (e.g., Fy vs. Fz), ideal for coordination pattern analysis.

## Example data

There are two seperate example CSVs files in the `data/` directory, one with a single trial and one with multiple (10) trials:
- `data/single_trial.csv`  
- `data/multiple_trials.csv`  

Upload these directly via the app’s file‐picker to see how multi-trial and single-trial inputs work.

## Installation

```bash
# Clone the repository
git clone https://github.com/The-DanKadlec/shinyForce3D
cd shinyForce3D

# Install as an R package
R -e "devtools::install_local('.')"
```

## Usage

```r
# Launch the Shiny app from R
library(shinyForce3D)
shinyForce3D::run_app()
```

### From command line

```bash
Rscript -e "shinyForce3D::run_app()"
```

## Documentation

* **Vignettes**: See `vignettes/` for end‑to‑end workflows (single trial, multi‑trial averaging, export GIF).
* **API Reference**: Auto-generated help pages live in man/ (e.g. ?computeResultant, ?plotDynamic3D, ?exportGif).

## Contributing

Contributions are welcome! Please:

1. Fork the repository and create a branch for your feature or bug fix.
2. Write tests (via **testthat** and **shinytest2**) for new functionality.
3. Submit a pull request describing your changes.

## License

This software is released under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license.

## Citation

If you use this tool in your research, please cite:

> Kadlec, D., & Vial, S. (2025). *shinyForce3D: A R/Shiny App for 3D Force Vector Visualisation of Ground Reaction Force*. (https://doi.org/10.5281/zenodo.15470997)
