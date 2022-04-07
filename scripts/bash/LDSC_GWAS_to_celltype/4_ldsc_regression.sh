#! /bin/bash
# choose GWAS results of interest:
# choose from: "Shrine_lung_function", "Sakornsakolpat_COPD", "McKay_lung_cancer", "Howard_depression", "Han_asthma" 
GWAS_study="Han_asthma"
# set matching dataset name, which is often the same but in some cases not:
dataset_name=${GWAS_study}
if [ ${GWAS_study} = "Shrine_lung_function" ]; then dataset_name="Shrine_lung_function_fvc"; fi
if [ ${GWAS_study} = "McKay_lung_cancer" ]; then dataset_name="McKay_lung_cancer_adenocarcinoma"; fi
# print for logging:
echo "Working on GWAS study $GWAS_study, with dataset name $dataset_name"
# adjust as needed:
ldsc_path=$HOME/software/ldsc
# activate ldsc environment
conda activate ldsc
${ldsc_path}/ldsc.py \
--h2-cts ../../../data_GWAS_non_public/$GWAS_study/${dataset_name}.sumstats.gz \
--ref-ld-chr ../../../results/LDSC_GWAS_to_celltype/ext_files/1000G_EUR_Phase3_baseline/baseline. \
--out ../../../results/LDSC_GWAS_to_celltype/ldsc_results/${dataset_name}_1000genes \
--ref-ld-chr-cts ../../../results/LDSC_GWAS_to_celltype/cts.ldcts \
--w-ld-chr ../../../results/LDSC_GWAS_to_celltype/ext_files/weights_hm3_no_hla/weights. \
# h2-cts flag specifies that we want to do a celltype specific analysis. Put the sumstats file path here
# ref-ld-chr: location of baseline model files (chr part to specify that files are split per chromosome)
# ref-ld-chr-cts Path to file that specifies location of cell type an control annotation files 
# w-ld-chr location of LD scores used for regression weights. Using hapmap3 file is recommended
