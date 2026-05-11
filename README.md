# Statistics of Spatial Averages and Optimal Averaging in the Presence of Missing Data

Code and data accompanying the paper:

> Seshadri, A. K. (2018). Statistics of spatial averages and optimal averaging in the presence of missing data. *Spatial Statistics*, 25, 1–21. [doi:10.1016/j.spasta.2018.04.002](https://doi.org/10.1016/j.spasta.2018.04.002)

## Overview

This paper derives approximate estimators for bias and variance of spatial averages (defined as weighted sums of point observations) when observations may be missing at random. The spatial average is a ratio of random variables, and its statistics are obtained using the delta method (truncated Taylor series expansions). The framework extends classical optimal averaging (OA) to account for missing data, and is applied to spatially averaged rainfall over India using gridded rain-gauge products from the India Meteorological Department (IMD).

## Repository Structure

```
├── README.md
├── .gitignore
├── data/
│   ├── indiadat.mat          # 1° × 1° gridded daily rainfall (1901–2011, 357 locations)
│   ├── indialatlon.mat       # Latitude/longitude coordinates for the 1° grid
│   └── meanhrraindat.mat     # 0.25° × 0.25° mean high-resolution rainfall data
└── src/
    ├── PFig1.m               # Figure 1: Verification of bias/variance models vs Monte Carlo
    ├── PFig2_3.m             # Figures 2–3: Optimal weights for min bias, variance, MSE
    ├── PFig4.m               # Figure 4: ISMR time-series and cumulative weight distributions
    ├── PFig5.m               # Figure 5: Standard error estimates (April–November)
    ├── PFig6abc.m            # Figure 6a–c: Robustness analysis with 1° data as truth
    ├── PFig6d.m              # Figure 6d: Standard error from robustness analysis
    ├── getcumdist.m          # Helper: compute cumulative distribution
    └── vect2matgeophys.m     # Helper: reshape vector data to lat/lon matrix
```

## Requirements

- **MATLAB** (tested with R2017b; should work with R2016a or later)
- MATLAB's Optimization Toolbox (for `quadprog`, used in optimal weight computation)

## Usage

1. Open MATLAB and set the working directory to the `src/` folder:
   ```matlab
   cd('/path/to/spatial-averages-missing-data/src')
   ```

2. Run any figure script directly. For example, to reproduce Figure 1:
   ```matlab
   PFig1
   ```

   Each script loads the required data from `../data/` and calls any necessary helper functions.

### Figure Descriptions

| Script | Paper Figure | Description |
|--------|-------------|-------------|
| `PFig1.m` | Fig. 1 | Compares analytical bias/variance estimators (Eqs. 22–23) with Monte Carlo simulations across different availability levels α |
| `PFig2_3.m` | Figs. 2–3 | Maps optimal weights for minimizing bias, variance, and MSE of Indian Summer Monsoon Rainfall (ISMR) |
| `PFig4.m` | Fig. 4 | ISMR time-series under optimal averaging, plus cumulative distributions of optimal weights |
| `PFig5.m` | Fig. 5 | Standard error of all-India rainfall for individual months (April–November) |
| `PFig6abc.m` | Fig. 6a–c | Robustness check: optimal weights when 1° data is treated as ground truth |
| `PFig6d.m` | Fig. 6d | Standard error from the robustness analysis |

## Data Sources

- **1° × 1° gridded rainfall**: Rajeevan, M., Bhate, J., Kale, J. D., & Lal, B. (2006). High resolution daily gridded rainfall data for the Indian region: Analysis of break and active monsoon spells. *Current Science*, 91, 296–306.
- **0.25° × 0.25° gridded rainfall**: Pai, D. S., et al. (2014). Development of a new high spatial resolution (0.25° × 0.25°) long period (1901–2010) daily gridded rainfall dataset over India and its comparison with existing datasets over the region. *Mausam*, 65, 1–18.

## Citation

If you use this code or data, please cite:

```bibtex
@article{Seshadri2018,
  title   = {Statistics of spatial averages and optimal averaging in the presence of missing data},
  author  = {Seshadri, Ashwin K.},
  journal = {Spatial Statistics},
  volume  = {25},
  pages   = {1--21},
  year    = {2018},
  doi     = {10.1016/j.spasta.2018.04.002}
}
```

## License

Please refer to the journal's data sharing policy. Code is provided as supplementary material to the paper.
