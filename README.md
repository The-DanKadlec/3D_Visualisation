# 3D Visualisation of Ground Reaction Forces (GRF)

This repository contains an R/Shiny application and supporting R package for visualizing ground reaction force (GRF) vector data in both dynamic interactive 3D, static 3D, and 2D formats. The tool supports single and dual dataset comparisons, customizable sphere/arrow sizing, peak detection, and export of animated GIFs.

## Features

* **Interactive 3D**: Real‑time rotation and exploration of GRF vectors as colored spheres or arrows.
* **Static 3D**: Publication‑quality snapshots from orthogonal angles, with optional butterfly lines and peak vector highlighting.
* **2D Time‑Series**: Plot of Fx, Fy, Fz, and resultant forces over stance phase, with mean ± SD ribbons, peak markers, and customizable legends.
* **GIF Export**: Generate rotating GIF animations with all annotations.
* **Adjustable Parameters**: Slider control for sphere/arrow size multiplier, filter cutoff, and color palette selection.

## Installation

```bash
# Clone the repository
git clone https://github.com/The-DanKadlec/3D_Visualisation.git
cd 3D_Visualisation

# Install as an R package
R -e "devtools::install_local('.')"
```

## Usage

```r
# Launch the Shiny app from R
library(3D_vis)
runApp(system.file('app', package = '3D_vis'))
```

### From command line

```bash
Rscript -e "grfVisual::runApp()"
```

## Example Data

Place your CSV files (Fx, Fy, Fz columns) in the `data/` folder, then use the file upload controls in the app to load and visualize.

## Documentation

* **Vignettes**: See `vignettes/` for end‑to‑end workflows (single trial, multi‑trial averaging, export GIF).
* **API Reference**: Generated documentation in `man/` for all core functions (e.g. `computeResultant()`, `plotDynamic3D()`, `exportGif()`).

## Contributing

Contributions are welcome! Please:

1. Fork the repository and create a branch for your feature or bug fix.
2. Write tests (via **testthat** and **shinytest2**) for new functionality.
3. Submit a pull request describing your changes.

## License

This software is released under the [CC0 v1.0 Universal](LICENSE) license—public domain dedication.

## Citation

If you use this tool in your research, please cite:

> Vial, J., Kadlec, J., & Smith, A. (2025). *grfVisual: An R/Shiny application for interactive ground reaction force visualization*. Journal of Open Source Software, X(X), XXX. [https://doi.org/XX.XXXXX/joss.XXXXX](https://doi.org/XX.XXXXX/joss.XXXXX)
