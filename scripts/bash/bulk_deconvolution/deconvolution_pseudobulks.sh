# Pseudobulk airway:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/CIBERSORTx_output_airway/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture pseudobulk_airway.txt --refsample airway_subsampled_matrix_max200cells_HUGO.txt --single_cell TRUE

# Pseudobulk nose:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/CIBERSORTx_output_nose/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture pseudobulk_nose.txt --refsample nose_subsampled_matrix_max200cells_HUGO.txt --single_cell TRUE

# Pseudobulk parenchyma:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/CIBERSORTx_output_parenchyma/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture pseudobulk_parenchyma.txt --refsample parenchyma_subsampled_matrix_max200cells_HUGO.txt --single_cell TRUE
