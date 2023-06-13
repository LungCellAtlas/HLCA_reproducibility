# HLCA reproducibility
This repository contains the code that was used for the Human Lung Cell Atlas project. 

## File organization:
You can find the main code used in the HLCA project in the [notebooks](./notebooks) of this repository. We subdivided the notebooks into four main categories:
- 1: building and annotating the atlas core  
- 2: downstream analysis of the HLCA core  
- 3: atlas extension<br>
- 4: disease analysis (analysis of IPF across datasets, and cross-disease analysis)<br>

The notebooks folder is organized based on these four categories, with every sub-folder containing the used code. Notebooks are numbered to clarify the order in which to run the notebooks.

The data and results folders as used in the code are not included in this repository, as the size of the files is too large, but will be uploaded to figshare soon. Please also feel free to contact us if you are looking for a specific file.

## Figures and tables:
To help you find the notebooks you need, we here specify where to find the code to generate each figure. Further details are included in the listed folders.

__Main figures:__<br>
1: manually created <br>
2a-c: notebook folder 1<br>
2d: scripts/R/integration_benchmark_plots/<br> 
3: notebook folder 1<br>
4: notebook folder 2<br>
5a-c: notebook folder 3<br>
5d: notebook folder 2<br>
5e: scripts/R/bulk_deconvolution/process_bulk_deconvolution.R (preparation in notebook folder 2)<br>
5f: notebook folder 3<br>
6a-e: notebook folder 3<br>
6f-k: notebook folder 4<br>

__Extended Data figures:__<br>
ED1-4: notebook folder 1<br>
ED5: notebook folder 2<br>
ED 6-9: notebook folder 3<br>
ED 10: notebook folder 4<br>

__Supplementary figures:__<br>
S1: scripts/R/integration_benchmark_plots/<br>
S2-3: notebook folder 1<br>
S4-5: notebook folder 2<br>
S6: notebook folder 3<br>
S7-8: notebook folder 2<br>
S9-10: notebook folder 3<br>

__Supplementary Data tables:__<br>
S table 1-2: notebook folder 3<br>
S table 3: manually generated<br>
S table 4-5: notebook folder 3<br>
S table 6: notebook folder 1<br>
S table 7-11: notebook folder 2<br>
S table 12-13: notebook folder 3<br>
S table 14-15: notebook folder 4<br>
S table 16: notebook folder 1<br>
S table 17: notebook folder 3<br>

### Custom generation of cell type signature matrices for deconvolution
To generate custom cell type signature matrices for deconvolution of bulk samples, use the [deconvolution notebook](./notebooks/2_downstream_analysis_of_HLCA_core/11_EXAMPLE_NOTEBOOK_signature_matrix_generation_for_deconvolution.ipynb)  also used for the deconvolution done in the HLCA paper.

## References:
HLCA paper: [Sikkema et al., Nature Medicine 2023](https://www.nature.com/articles/s41591-023-02327-2)

## In cases of questions:
Submit an issue to this repository.

## Acknowledgements:
We would like to thank Ciro Ramirez Suastegui (@cramsuig), Tessa Gillett (@TessaGillet), Daniel Strobl (@danielStrobl) and Luke Zappia (@lazappi) for important contributions to the code.
