# analysis_antarctic-asym

## Code for the figures and analysis in `Asymmetry in the seasonal cycle of Antarctic sea ice driven by insolation', acccepted at Nature Geoscience

Paper citation: TBU

Analysis code developed by LR, with help from IE and TW. The WE15 model was developed by Wagner and Eisenman (2015, Journal of Climate) and is available in its original form at http://eisenman.ucsd.edu/code.html.

This repository contains:

- asym_funcs.py - a file containing Python functions used in multiple notebooks. This includes a Python function for the WE15 model
- processing/ - a directory containing Python notebooks that were used to process data from NSIDC, CMIP and NOAA into the formats needed for this analysis 
- processed/ - netcdf files produced by the notebooks in processing/
- we_15_runs/ - a directory containing the scripts used to run various configurations of the WE15 model as well as the netcdf output
- aasym*.ipynb - the notebooks used to produce all figures in the above paper. Anyone should be able to run these after downloaded the repository, since they only used input from processed/ and we_15_runs
