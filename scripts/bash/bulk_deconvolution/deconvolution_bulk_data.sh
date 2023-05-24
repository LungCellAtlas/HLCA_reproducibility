# Bulk nose:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/CIBERSORTx_output_nose/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture nose_bulk.txt --refsample nose_subsampled_matrix_max200cells_ENSG.txt --single_cell TRUE

# Bulk airway *with* NORM:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/CIBERSORTx_output_airway_withNORM/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture airway_bulk_withNORM.txt --refsample airway_subsampled_matrix_max200cells_ENSG.txt --single_cell TRUE

# Bulk parenchyma:

docker run -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/:/src/data -v /home/tessa/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/CIBERSORTx_output_parenchyma/:/src/outdir cibersortx/fractions --username t.e.gillett@umcg.nl --token  93381f107e2ffc18d2edf71e84b64897 --mixture parenchyma_bulk_noNA.txt --refsample parenchyma_subsampled_matrix_max200cells_HUGO.txt --single_cell TRUE --QN TRUE
