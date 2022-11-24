# Integration benchmark plots

This folder contains all files needed to generate the integration benchmark plots.

To generate the figures, run `run_plotting_Lisa.R`, all other R scripts are used in the background, together with the img folder and the metrics (metrics_scgen_added.csv) that are the output from the integration benchmark. Make sure to update file paths where necessary. The summary scores files are output from the plotting scripts.

The integration itself was done with scIB (see [the scIB GitHub page](https://github.com/theislab/scib)) and the HLCA paper methods, together with the config file in the configs folder of this repo. 