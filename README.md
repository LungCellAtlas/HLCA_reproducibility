# HLCA reproducibility
This repository contains the code that was used for the Human Lung Cell Atlas project. 

## File organization:
You can find the main code used in the HLCA project in the [notebooks](./notebooks) of this repository. We subdivided the notebooks into three main categories:
- 1: building and annotating the atlas core  
- 2: downstream analysis of the HLCA core  
- 3: atlas extension<br>
The notebooks folder is organized based on these three categories, with every sub-folder containing the used code. Notebooks are numbered to clarify the order in which to run the notebooks.

The data and results folders as used in the code are not included in this repository, as the size of the files is too large, but will be uploaded to figshare soon.

## Figures and tables:
To help you find the notebooks you need, we here specify where to find the code to generate each figure. 

Main figures:<br>
1: manually created <br>
2a-c: notebook folder 1<br>
2d: scripts/R/integration_benchmark_plots/<br> 
3: notebook folder 1<br>
4: notebook folder 2<br>
5a-c: notebook folder 3<br>
5d: notebook folder 2<br>
5e: deconvolution analysis COMING SOON in scripts/R folder<br>
5f: notebook folder 3<br>
6: notebook folder 3<br>
7: COMING SOON disease analysis Ciro<br>

Extended Data figures:<br>
ED1: scripts/R/integration_benchmark_plots/<br>
ED2-8: notebook folder 1<br>
ED9-11: notebook folder 2<br>
ED 12-14: notebook folder 3<br>
ED 15-17: notebook folder 2<br>
ED 18-22: notebook folder 3<br>
ED 23: COMING SOON disease analysis Ciro<br>

Extended Data tables:<br>
ED table 1-2: notebook folder 3<br>
ED table 3: manually generated<br>
ED table 4-5: notebook folder 3<br>
ED table 6: notebook folder 1<br>
ED table 7-10: notebook folder 2<br>
ED table 11: deconvolution analysis COMING SOON in scripts/R folder<br>
ED table 12-13: notebook folder 3<br>
ED table disease analysis: COMING SOON disease analysis Ciro<br>
ED table X (marker gene selection parameters): notebook folder 1<br>

### Custom generation of cell type signature matrices for deconvolution
To generate custom cell type signature matrices for deconvolution of bulk samples, use the [deconvolution notebook](./notebooks/2_downstream_analysis_of_HLCA_core/7_cell_type_gene_signature_generation_for_bulk_deconvolution.ipynb)  also used for the deconvolution done in the HLCA paper.

## References:
HLCA pre-print: [Sikkema et al., bioRxiv 2022](https://www.biorxiv.org/content/10.1101/2022.03.10.483747v1.full)

## In cases of questions:
Submit an issue to this repository.
