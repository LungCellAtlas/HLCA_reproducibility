import numpy as np
import scanpy as sc
import anndata
import pandas as pd
import matplotlib.pyplot as plt
import collections
import glob

def get_gene_renamer_dict(gene_names):
    """This function detects genes that are misnamed. Note that it should only
    be run on adata.var.index after removing duplicate gene names! 
    It detects 1) genes with double dashes ("--") and translates them 
    to double underscores ("__"); 2) genes that occur both with a 
    "1" suffix (".1", "-1" or "_1") and without suffix, and translated 
    the with-suffix genes to the without-suffix versions, and 3)
    it detects genes that occur with corresponding but different suffixes 
    (".[suffix]" and "-[suffix]") and translates them to either 1) a version 
    without suffix, if the suffix is ".[single_digit]", or else 2) to the dash-version
    of the suffix, since this is usually the ensembl version. 
    Input: list of gene names
    Returns a dictionary with original gene names as keys, and proposed 
    translations as values."""
    # first check if there are no duplicate gene names:
    if len(gene_names) > len(set(gene_names)):
        print()
        raise ValueError("There are still duplicate gene names in your gene list, \
this function does not work properly with duplicate names. Exiting.")
    # rename genes with double dash to genes with double underscore:
    genes_with_double_dash = [gene for gene in gene_names if "--" in gene]
    genes_to_change_suffix_dict = dict()
    verbose = False
    # find genes that have a suffixed and non-suffixed version
    # (this will detect both genes that only have a non-suffixed version,
    # and genes that have variable suffixes and a non-suffixed version)
    for gene in gene_names:
        # check if it has a suffix:
        if len(gene) > 1:
            if gene[-2:] in [".1", "-1", "_1"]:
                # check if it has a non-suffix version:
                if gene[:-2] in gene_names:
                    genes_to_change_suffix_dict[gene] = gene[:-2]
    # find all genes that have both a .[suffix] and a -[suffix] version
    # (this differs from the loop above, since it will find genes that don't
    # have a short version in the data):
    # a = [x.replace("-", ".") for x in gene_names]
    a = [
        rreplace(
            original_gene_name, old="-", new=".", occurrence=1, only_if_not_followed_by="."
        )
        for original_gene_name in gene_names
    ]
    from_dash_to_dot = [
        item for item, count in collections.Counter(a).items() if count > 1
    ]
    for gene in from_dash_to_dot:
        # get gene name without suffix
        gene_minus_suffix_list = gene.split(".")
        gene_minus_suffix = ".".join(gene_minus_suffix_list[:-1])
        # map "-" version of name to new name
        dot_position = gene.rfind(".")
        gene_as_list = list(gene)
        gene_as_list_dash = gene_as_list
        gene_as_list_dash[dot_position] = "-"
        gene_dash_name = "".join(gene_as_list_dash)
        # if the suffix is a .[single_digit], change it to the name without suffix:
        if gene.split(".")[-1] in [str(digit) for digit in range(0, 10)]:
            # map "." version of name to new name
            genes_to_change_suffix_dict[gene] = gene_minus_suffix
            genes_to_change_suffix_dict[gene_dash_name] = gene_minus_suffix
            # else, change the dot version to the dash-version of the name
            # (this is usually the ensembl version):
        else:
            genes_to_change_suffix_dict[gene] = gene_dash_name
    # now create the remapping dictionary
    genes_remapper = dict()
    for gene in genes_with_double_dash:
        genes_remapper[gene] = gene.replace("--", "__")
    for gene_name_old, gene_name_new in genes_to_change_suffix_dict.items():
        genes_remapper[gene_name_old] = gene_name_new
    return genes_remapper


def subset_and_pad_adata(gene_set, adata):
    """
    This function uses a gene list provided as a Pandas dataframe with gene symbols and
    Ensembl IDs and subsets a larger Anndata object to only the genes in this list. If
    Not all genes are found in the AnnData object, then zero-padding is performed.
    """
    # Example inputs:
    # genes_filename = '/storage/groups/ml01/workspace/hlca_lisa.sikkema_malte.luecken/genes_for_mapping.csv'
    # data_filename = '/storage/groups/ml01/workspace/hlca_lisa.sikkema_malte.luecken/ready/adams.h5ad'
    # gene_set = pd.read_csv(genes_filename)
    # adata = sc.read(data_filename)

    # Prep objects
    if 'gene_symbols' in gene_set.columns:
        gene_set.index = gene_set['gene_symbols']

    else:
        raise ValueError('The input gene list was not of the expected type!\n'
                         'Gene symbols and ensembl IDs are expected in column names:\n'
                         '\t`gene_symbols` and `Unnamed: 0`')

    # Subset adata object
    common_genes = [gene for gene in gene_set['gene_symbols'].values if gene in adata.var_names]
    if len(common_genes) == 0:
        print("WARNING: YOU SHOULD PROBABLY SWITCH YOUR ADATA.VAR INDEX COLUMN TO GENE NAMES"
                  " RATHER THAN IDS! No genes were recovered.")
        return

    adata_sub = adata[:,common_genes].copy()

    # Pad object with 0 genes if needed
    if len(common_genes) < len(gene_set):
        diff = len(gene_set) - len(common_genes)
        print(f'not all genes were recovered, filling in 0 counts for {diff} missing genes...')
        
        # Genes to pad with
        genes_to_add = set(gene_set['gene_symbols'].values).difference(set(adata_sub.var_names))
        new_var = gene_set.loc[genes_to_add]
        
        if 'Unnamed: 0' in new_var.columns:
            # Assumes the unnamed column are ensembl values
            new_var['ensembl'] = new_var['Unnamed: 0']
            del new_var['Unnamed: 0']

        df_padding = pd.DataFrame(data=np.zeros((adata_sub.shape[0],len(genes_to_add))), index=adata_sub.obs_names, columns=new_var.index)
        adata_padding = sc.AnnData(df_padding, var=new_var)
        
        # Concatenate object
        adata_sub = anndata.concat([adata_sub, adata_padding], axis=1, join='outer', index_unique=None, merge='unique')

    # Ensure ensembl IDs are available
    adata_sub.var['ensembl'] = gene_set['Unnamed: 0']

    return adata_sub


def add_up_duplicate_gene_name_columns(adata, print_gene_names=True, verbose=False):
    """ This function finds duplicate gene names in adata.var (i.e. duplicates 
    in the index of adata.var). For each cell, it adds up the counts of columns 
    with the same gene name, removes the duplicate columns, and replaces the 
    counts of the remaining column with the sum of the duplicates.
    Returns anndata object."""
    print("TO DO: STORE ENSEMBL IDS OF MERGED IDS AND MATCHING COUNTS IN A SEPARATE COLUMN!!!")
    duplicate_boolean = adata.var.index.duplicated()
    duplicate_genes = adata.var.index[duplicate_boolean]
    print("Number of duplicate genes: " + str(len(duplicate_genes)))
    if print_gene_names == True:
        print(duplicate_genes)
    columns_to_replace = list()
    columns_to_remove = list()
    new_columns_array = np.empty((adata.shape[0], 0))
    for gene in duplicate_genes:
        if verbose:
            print("Calculating for gene", gene)
        # get rows in adata.var with indexname equal to gene
        # indexing zero here is to get rid of tuple output and access array
        gene_colnumbers = np.where(adata.var.index == gene)[0]
        # isolate the columns
        gene_counts = adata.X[:, gene_colnumbers].copy()
        # add up gene counts and add new column to new_columns_array
        new_columns_array = np.hstack((new_columns_array, np.sum(gene_counts, axis=1)))
        # store matching column location in real adata in list:
        columns_to_replace.append(gene_colnumbers[0])
        # store remaining column locations in columns to remove:
        columns_to_remove = columns_to_remove + gene_colnumbers[1:].tolist()
    # replace first gene column with new col:
    adata.X[:, columns_to_replace] = new_columns_array
    # remove remaining duplicate columns:
    columns_to_keep = [
        i for i in np.arange(adata.shape[1]) if i not in columns_to_remove
    ]
    adata = adata[:, columns_to_keep].copy()
    if verbose:
        print("Done!")
    return adata
    

def add_cell_annotations(adata, var_index):
    """ This function adds annotation to anndata:  
    cell level:  
    total_counts, log10_total_counts, n_genes_detected, mito_frac, ribo_frac,   
    compl(exity)  
    gene_level:  
    n_cells 

    Arguments:
        adata - anndata object, raw (unnormalized!)
        var_index < "gene_symbols", "gene_ids" > - set to which type of gene
            naming is used in adata.var.index

    Returns:
        anndata object with annotations  
    """
    # cells:
    # total transcript count per cell
    adata.obs['total_counts'] = np.sum(adata.X, axis=1)
    adata.obs['log10_total_counts'] = np.log10(adata.obs['total_counts'])
    # number of genes expressed
    # translate matrix to boolean (True if count is larger than 0):
    boolean_expression = adata.X > 0
    adata.obs['n_genes_detected'] = np.sum(boolean_expression, axis=1)
    # fraction mitochondrial transcripts per cell
    if var_index == "gene_symbols":
        mito_genes = [
            gene for gene in adata.var.index if gene.lower().startswith("mt-")
        ]
    elif var_index == "gene_ids":
        mito_genes = [
            gene_id
            for gene_id, gene_symbol in zip(adata.var.index, adata.var.gene_symbols)
            if gene_symbol.lower().startswith("mt-")
        ]
    else:
        raise ValueError(
            "var_index argument should be set to either gene_symbols or gene_ids"
        )
    # conversion to array in line below is necessary if adata.X is sparse
    adata.obs['mito_frac'] = np.array(np.sum(
        adata[:,mito_genes].X, axis=1)).flatten() / adata.obs['total_counts']
    # fraction ribosomal transcripts per cell
    if var_index == "gene_symbols":
        ribo_genes = [
            gene
            for gene in adata.var.index
            if (gene.lower().startswith("rpl") or gene.lower().startswith("rps"))
        ]
    elif var_index == "gene_ids":
        ribo_genes = [
            gene_id
            for gene_id, gene_symbol in zip(adata.var.index, adata.var.gene_symbols)
            if (
                gene_symbol.lower().startswith("rpl")
                or gene_symbol.lower().startswith("rps")
            )
        ]
    adata.obs['ribo_frac'] = np.array(np.sum(
        adata[:,ribo_genes].X, axis=1)).flatten() / adata.obs['total_counts']
    # cell complexity (i.e. number of genes detected / total transcripts)
    adata.obs['compl'] = adata.obs['n_genes_detected']\
    / adata.obs['total_counts']
    # genes
    adata.var['n_cells'] = np.sum(boolean_expression, axis=0).T
    return adata


def get_sample_annotation_table_LCA(data_dir, verbose=True):
    """reads in metadata as collected throught LCA_metadata tables.
    args:
        project_dir - path to directory that has LCA_metadata folder
    Returns:
        pandas dataframe with one row per sample, and all matching metadata."""
    # get paths to all LCA metatables available
    file_paths_long = glob.glob(f"{data_dir}/LCA_metadata*.csv")
    # store paths together with file names (some files are stored in two different
    # places, since the same lab had multiple datasets, and they filled out only
    # one table)
    path_to_name_dir = {
        file_path: file_path.split("/")[-1] for file_path in file_paths_long
    }
    files_read = list()
    meta_tables = dict()
    # store metatables for each unique file name
    for file_path, file_name in path_to_name_dir.items():
        if file_name not in files_read:
            print(file_name)
            # read csv
            meta = pd.read_csv(file_path, index_col=2)
            # remove example rows
            meta = meta.loc[[inst != "EXAMPLE INSTITUTE" for inst in meta.Institute], :]
            # remove trailing spaces from column names:
            col_renamer = {col:col.strip(" ") for col in meta.columns}
            meta.rename(col_renamer, axis=1, inplace=True)
            # store
            meta_tables[file_path] = meta
            files_read.append(file_name)
    # merge tables into one
    metadata = pd.concat(meta_tables.values())
    # remove rows with NaN as index
    print(
        "number of rows without rowname/sample name (will be removed):",
        sum(metadata.index.isnull()),
    )
    metadata = metadata.loc[metadata.index.isnull() == False, :]
    # check if sample ids abd donor ids are unique
    sample_names_unique = len(metadata.index) == len(set(metadata.index.tolist()))
    if not sample_names_unique or verbose:
        print("Sample IDs unique?", sample_names_unique)
    # print number of samples without donor ID
    n_samples_without_subject_ID = sum(metadata.subject_ID.isnull())
    if n_samples_without_subject_ID != 0 or verbose:
        print("Number of samples without donor ID:", n_samples_without_subject_ID)
    # remove spaces and commas from column names:
    col_renamer = {
        col: col.replace(",", "").replace(" ", "_").strip("_") for col in metadata.columns
    }
    metadata.rename(col_renamer, axis=1, inplace=True)
    # return result
    return metadata



def plot_QC(anndata_dict, project_name, scale=1):
    n_cols = 4
    n_rows = len(anndata_dict.keys())
    fig = plt.figure(figsize=(n_cols * 3 * scale, n_rows * 2.5 * scale))
    plt.suptitle('QC ' + project_name, fontsize=16, y=1.02)
    title_fontsize = 8
    ax_dict = dict()
    ax = 1
    for sample_ID in sorted(anndata_dict.keys()):
        adata = anndata_dict[sample_ID]
        # total counts
        ax_dict[ax] = fig.add_subplot(n_rows, n_cols, ax)
        ax_dict[ax].hist(adata.obs['log10_total_counts'], bins=50, range=(2,5))
        ax_dict[ax].set_title('log10 total counts per cell', 
                              fontsize=title_fontsize)
        ax_dict[ax].set_xlabel('log10 total counts')
        ax_dict[ax].set_ylabel(sample_ID, fontsize=16)
        ax = ax + 1
        # mitochondrial transcript percentage
        ax_dict[ax] = fig.add_subplot(n_rows, n_cols, ax)
        ax_dict[ax].hist(adata.obs['mito_frac'], bins=50, range=(0,1))
        ax_dict[ax].set_title('fraction of mitochondrial transcripts/cell', 
                      fontsize=title_fontsize)
        ax_dict[ax].set_xlabel('fraction of mito transcripts')
        ax_dict[ax].set_ylabel('frequency')
        ax = ax + 1
        # log10 number of genes detected
        ax_dict[ax] = fig.add_subplot(n_rows, n_cols, ax)
        ax_dict[ax].hist(np.log10(adata.obs['n_genes_detected']), 
                         bins=50, range=(2,5))
        ax_dict[ax].set_title('log10 number of genes detected', 
                              fontsize=title_fontsize)
        ax_dict[ax].set_xlabel('log10 number of genes')
        ax_dict[ax].set_ylabel('frequency')
        ax = ax + 1
        # cell complexity, scatter
        ax_dict[ax] = fig.add_subplot(n_rows, n_cols, ax)
        mappable_ax = str(ax) + '_m'
        ax_dict[mappable_ax] = ax_dict[ax].scatter(adata.obs['log10_total_counts'],
                    np.log10(adata.obs['n_genes_detected']), 
                    c=adata.obs['mito_frac'], s=1)
        ax_dict[ax].set_title('cell complexity \ncolored by mitochondrial fraction', 
                      fontsize=title_fontsize)
        ax_dict[ax].set_xlabel('log10 total counts')
        ax_dict[ax].set_ylabel('log10 number of genes detected')
        fig.colorbar(mappable=ax_dict[mappable_ax],ax=ax_dict[ax])
        ax = ax + 1
    plt.tight_layout()
    plt.show()
    return fig


def age_converter(age, age_range, verbose=False):
    """takes in two values (age and age range (format: [number]-[number])) 
    where only one has a value, and the other is NaN. Returns 
    either the original age or the mean of the age range."""
    if isinstance(age, float):
        if np.isnan(age):
            if verbose:
                print(age_range)
            if isinstance(age_range, float):
                if np.isnan(age_range):
                    if verbose:
                        print("no age or age range available")
                    return np.nan
            else:
                age_lower = np.float(age_range.split("-")[0])
                age_higher = np.float(age_range.split("-")[1])
                mean_age = np.mean([age_lower, age_higher])
            if verbose:
                print(mean_age)
            return mean_age
        else:
            return age
    else:
        return age


def generate_r2python_gene_mapping(r_genes, python_genes):
    """function to map genes from r (that replaces "-" with ".") to python genes.
    Returns dictionary with mapping."""
    r2python_gene_mapping = dict()
    for gene in r_genes:
        if gene in python_genes:
            r2python_gene_mapping[gene] = gene
        elif "-".join(gene.rsplit(".", 1)) in python_genes:
            r2python_gene_mapping[gene] = "-".join(gene.rsplit(".", 1))
        elif "-".join(gene.split(".", 1)) in python_genes:
            r2python_gene_mapping[gene] = "-".join(gene.split(".", 1))
        elif "-".join(gene.split(".")) in python_genes:
            r2python_gene_mapping[gene] = "-".join(gene.split("."))
        else:
            print("no equivalent found for this gene:", gene)
    return r2python_gene_mapping



def update_adata(adata, cov1, cov2, mapping):
    """
    Updates anndata obs column based on mapping from one column to the 
    other. 

    Arguments:
    adata - AnnData to update
    cov1 - column name of covariate to map from
    cov2 - column name of covariate to map to
    mapping - dictionary with mapping from covariate 1 to covariate 2

    Returns:
    adata with updated adata.obs[cov2]
    """
    cov1_in_data = adata.obs[cov1].unique()
    for cat in mapping.keys():
        if cat not in cov1_in_data:
            raise ValueError(f"{cat} not a category of {cov1} in your data.")
    cov2_updates = adata.obs[cov1].map(mapping).dropna()
    # set cov2 column to list so that we can add new categories if necessary
    adata.obs[cov2] = adata.obs[cov2].tolist()
    adata.obs.loc[cov2_updates.index, cov2] = cov2_updates.values
    # make into category again
    adata.obs[cov2] = pd.Categorical(adata.obs[cov2])
    # and return updated anndata
    return adata


# HELPER FUNCTIONS


# helper function for get_gene_renamer_dict function
def rreplace(s, old, new, occurrence, only_if_not_followed_by=None):
    """replaces occurences of character, counted from last to first, with new 
    character. 
    Arguments:
    s - string
    old - character/string to be replaced
    new - character/string to replace old with
    occurence - nth occurence up to which instances should be replaced, counting 
    from end to beginning. e.g. if set to 2, the two last instances of "old" 
    will be replaced
    only_if_not_followed_by - if this is set to a string, [old] is only replaced
    by [new] if [old] is not followed by [only_if_not_followed_by] in [s]. Useful
    for gene names like abc-1.d, in that case we do not want to replace the dash.
    Returns:
    string with replaced character
    """
    if old not in s:
        return s
    elif only_if_not_followed_by != None:
        # test if latest instance of old is followed by [only_if_not_followed_by]
        last_inst_old = s.rfind(old)
        last_inst_oinfb = s.rfind(only_if_not_followed_by)
        if last_inst_oinfb > last_inst_old:
            return s
    li = s.rsplit(old, occurrence)
    return new.join(li)

