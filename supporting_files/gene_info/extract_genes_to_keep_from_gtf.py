# written by Lisa, 2021.03.25
from gtfparse import read_gtf
# set working directory
work_dir = "./"
# genes to keep:
# based on https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#grch38_1.2.0
gene_biotypes_to_keep = [
    "protein_coding",
    "lincRNA",
    "antisense",
    "IG_LV_gene",
    "IG_V_gene",
    "IG_V_pseudogene",
    "IG_D_gene",
    "IG_J_gene",
    "IG_J_pseudogene",
    "IG_C_gene",
    "IG_C_pseudogene",
    "TR_V_gene",
    "TR_V_pseudogene",
    "TR_D_gene",
    "TR_J_gene",
    "TR_J_pseudogene",
    "TR_C_gene",
]
# returns GTF with essential columns such as "feature", "seqname", "start", "end"
# alongside the names of any optional keys which appeared in the attribute column
df = read_gtf(f"{work_dir}Homo_sapiens.GRCh38.84.gtf")
# select only genes that have one of the listed biotypes:
df_genes_to_keep = df.loc[df.gene_biotype.isin(gene_biotypes_to_keep), :].copy()
# extract gene names and ids
gene_ids_to_keep = sorted(set(df_genes_to_keep.gene_id))
gene_names_to_keep = sorted(set(df_genes_to_keep.gene_name))
# store to files:
with open(
    f"{work_dir}Homo_sapiens_GRCh38_84_gene_ids_to_keep_cellrangerbased.txt",
    "w",
) as output:
    for row in gene_ids_to_keep:
        output.write(str(row) + "\n")
with open(
    f"{workdir}Homo_sapiens_GRCh38_84_gene_names_to_keep_cellrangerbased.txt",
    "w",
) as output:
    for row in gene_names_to_keep:
        output.write(str(row) + "\n")
