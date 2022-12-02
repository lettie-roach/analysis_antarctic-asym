# analysis_antarctic-asym

## Code for the figures and analysis in `Asymmetry in the seasonal cycle of Antarctic sea ice driven by insolation'

Paper citation: Roach, L. A., Eisenman, I., Wagner, T. J., Blanchard-Wrigglesworth, E., and Bitz, C. M.. Asymmetry in the seasonal cycle of Antarctic sea ice due to insolation (2022). Nature Geoscience. https://doi.org/10.1038/s41561-022-00913-6. [Download paper]. [Download Supplementary]

Analysis code developed by LR, with help from IE and TW. The WE15 model was developed by Wagner and Eisenman (2015, Journal of Climate) and is available in its original form at http://eisenman.ucsd.edu/code.html.

This repository contains:

- asym_funcs.py - a file containing Python functions used in multiple notebooks. This includes a Python function for the WE15 model
- processing/ - a directory containing Python notebooks that were used to process data from NSIDC, CMIP and NOAA into the formats needed for this analysis 
- processed/ - netcdf files produced by the notebooks in processing/
- we_15_runs/ - a directory containing the scripts used to run various configurations of the WE15 model as well as the netcdf output
- aasym*.ipynb - the notebooks used to produce all figures in the above paper. Anyone should be able to run these after downloading the repository, since they only used input from processed/ and we_15_runs/
- matlab/ - a directory containing Matlab scripts that were used to confirm the results of calculations above
