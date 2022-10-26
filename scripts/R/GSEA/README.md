# Lung cell atlas gene sets

This repository contains gene set analysis for the Human Lung Cell Atlas.

## Directory structure

* `R/` - R function files
  * `gene_sets.R` - Functions for running gene set analysis and summarising results
  * `load.R` - Functions for loading datasets
  * `plotting.R` - Functions for creating plots
  * `save.R` - Functions for saving output files
* `renv/` - **{renv}** directory
* `renv.lock` - **{renv}** lockfile
* `_targets.R` - **{targets}** pipeline file
* `_targets_packages.R` - **{targets}** packages file
* `.gitignore` - git ignore file
* `2021-LungAtlas-scRNAseq.Rproj` - RStudio project file
* `README.md` - This README
* `LICENSE.md` - License file

## Setting up

### Installing dependencies

Dependencies for this analysis are managed using **{renv}**.
You can install them by running:

```r
if (!requireNamespace("renv")) {
    install.packages("renv")
}
renv::restore()
```

### Set paths

The first three targets in `_targets.R` are used to point to the input files.
If you are running the analysis in a different environment you may need to modify these.
The first target is the path to the base input file, the second is to the sub-directory containing the input file for differential expression analysis and the third is the sub-directory containing the output files from differential expression analysis.

## Running analysis

The analysis pipeline is managed using **{targets}**.
Running the analysis should be as simple as running:

```r
targets::tar_make()
```

You can also view the analysis dependency graph using:

```r
targets::tar_visnetwork()
```

## What does the pipeline do?

To see all the steps in the analysis please refer to the `_targets.R` file but briefly what happens is:

* Load the results from the differential expression analysis
* Map gene names to ENTREZ IDs
* Download gene sets
* Run gene set analysis for each DE comparison
* Summarise gene set results
* Check for common results
* Simplify GO results
