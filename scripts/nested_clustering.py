import numpy as np 
import pandas as pd 
import scanpy as sc
import scanpy.external as sce

def create_cluster_annotation_overview(
    adata,
    n_levels,
    cluster_label,
    min_fraction_for_dominancy=0.80,
    min_fraction_annotated=0.5,
    compartment_of_interest=None,
):
    """Function to calculate for each cluster, for each annotation level, if it is 
    dominated by a cell type.
    Args:
    adata - scanpy AnnData object
    n_levels - number of annotation levels (named "ann_level_[number]" in adata.obs)
    cluster_label - column name of cluster column in adata.obs
    min_fraction_for_dominancy - minimum fraction of annotated cells to belong to 
                                one cell type, to be called "dominant". Should be 
                                higher than 0.5.
    min_fraction_annotated - minumum fraction of cells in cluster that need to be
                            annotated, before "dominancy" analysis is possible
    compartment_of_interest - ann_level_1 compartment to which to limit the cluster
                            analysis. Only clusters that belong to multiple 
                            compartments or this specific compartment are included
                            in the output df.

    Returns:
    cluster df - pandas dataframe with for each cluster information on what is the 
                dominant cluster (if there is one), and the fraction of annotated
                cells belonging to the dominant cluster
    """
    cluster_df = pd.DataFrame(
        index=adata.obs[cluster_label].cat.categories,
        columns=zip(
            ["ann{}_dom_type".format(level) for level in range(1, n_levels + 1)],
            ["ann{}_dom_fract".format(level) for level in range(1, n_levels + 1)],
        ),
    )
    for level in range(1, n_levels + 1):
        level_name = "ann_level_" + str(level)
        clust_cell_types = adata.obs.groupby([cluster_label, level_name]).agg(
            {level_name: "count"}
        )
        # count fraction of cells that is annotated at this level:
        clust_cell_types["annotated"] = [
            "no" if celltype[:2] in ["1_", "2_", "3_", "4_"] else "yes"
            for celltype in clust_cell_types.index.get_level_values(1)
        ]
        number_annotated = clust_cell_types.groupby([cluster_label, "annotated"]).agg(
            {level_name: "sum"}
        )
        fraction_annotated = number_annotated.groupby(level=0).apply(
            lambda x: x / float(x.sum())
        )
        # keep only cells that are annotated at this level:
        rows_to_keep = [
            rownumber
            for rownumber, rowname in enumerate(clust_cell_types.index.get_level_values(1))
            if not rowname[:2] in ["1_", "2_", "3_", "4_"]
        ]
        clust_cell_types = clust_cell_types.iloc[rows_to_keep, :]
        # convert to proportions
        clust_cell_types = clust_cell_types.groupby(level=0)[level_name].apply(
            lambda x: x / float(x.sum())
        )
        # add "dominant" annotation:
        dominant_types = clust_cell_types.index[
            clust_cell_types > min_fraction_for_dominancy
        ]
        dominant_fractions = clust_cell_types[clust_cell_types > min_fraction_for_dominancy]
        # copy dominant types to cluster_df:
        cluster_df.loc[
            dominant_types.get_level_values(0), "ann{}_dom_type".format(level)
        ] = dominant_types.get_level_values(1)
        # copy dominance fractions to cluster_df
        cluster_df.loc[
            dominant_fractions.index.get_level_values(0), "ann{}_dom_fract".format(level)
        ] = dominant_fractions.values
        # set underannotated entries to "underann"
        # first, make sure columns are not categorical (they would not accept new cat)
        for cat in ["ann{}_dom_type".format(level), "ann{}_dom_fract".format(level)]:
            cluster_df[cat] = cluster_df[cat].tolist()
        idx = pd.IndexSlice
        underannotated_boolean = (
            fraction_annotated.loc[idx[:, "yes"], :] < min_fraction_annotated
        )
        cluster_df.loc[
            underannotated_boolean[level_name].values,
            ["ann{}_dom_type".format(level), "ann{}_dom_fract".format(level)],
        ] = "underann"
    if compartment_of_interest != None:
        # subset epithelial and split clusters
        cluster_df = cluster_df.loc[
            [
                main_type == compartment_of_interest or split_cluster
                for main_type, split_cluster in zip(
                    cluster_df.ann1_dom_type, cluster_df.ann1_dom_type.isnull()
                )
            ],
            :,
        ]
    return cluster_df



def add_nested_clustering(
    adata,
    cluster_df,
    cluster_label_previous,
    cluster_label_new,
    cluster_res=0.2,
    min_cluster_size=100,
    verbose=True,
):
    """Function that goes through one round of clustering of already existing 
    clusters, based on the input cluster df. All clusters that don't have a 
    dominant cluster yet at all levels in the df (as long as they are 
    sufficiently annotated) will be reclustered individually.
    Returns adata with new clustering (under adata.obs[cluster_label_new].
    """
    # copy original clustering
    adata.obs[cluster_label_new] = adata.obs[cluster_label_previous].tolist()
    for cluster in cluster_df.index:
        if verbose:
            print("Cluster:", cluster)
        dom_types = cluster_df.loc[cluster, :]
        if dom_types.isnull().any():
            subadata = adata[adata.obs[cluster_label_previous] == cluster, :].copy()
            if subadata.shape[0] < min_cluster_size:
                if verbose:
                    print("cluster size smaller than", min_cluster_size, "\n")
                continue
            if verbose:
                print("reclustering...\n")
            sc.tl.pca(subadata)
            sc.tl.leiden(subadata, resolution=cluster_res, key_added=cluster_label_new)
            subadata.obs[cluster_label_new] = [
                "{}.{}".format(cluster, new_cluster)
                for new_cluster in subadata.obs[cluster_label_new]
            ]
            adata.obs.loc[subadata.obs.index, cluster_label_new] = subadata.obs[
                cluster_label_new
            ]
        else:
            if verbose:
                print("clustered to full resolution!\n")
    # order categories "numerically" (so not 1, 10, 11 but 1, 2, 3... 10, 11):
    cluster_numbers = list(sorted(set(adata.obs[cluster_label_new])))
    prefix_cluster = [float(x.split(".")[0]) for x in cluster_numbers]
    cluster_numbers_ordered = [
        cluster_numbers[idx] for idx in np.argsort(prefix_cluster)
    ]
    adata.obs[cluster_label_new] = pd.Categorical(
        adata.obs[cluster_label_new], categories=cluster_numbers_ordered
    )

    return adata

def add_nested_clustering_blind(
    adata,
    cluster_label_previous,
    cluster_label_new,
    use_rep,
    cluster_alg="leiden",
    cluster_res=0.2,
    cluster_k=30,
    min_cluster_size=50,
    redo_pca=True,
    verbose=True,
):
    """Function that goes through one round of clustering of already existing
    clusters, based on the input cluster df. All clusters will be reclustered 
    individually. ("blind" because we don't take into account annotation
    purity of clusters.) 
    Args:
        adata - anndata object to be clustered
        cluster_label_previous - parent cluster label
        cluster_label_new - label for new clustering
        use_rep - name of .obsm object to be used for neighbor graph
        cluster_alg - <"leiden","phenograph">
        cluster_res - only applicable when using "leiden" as cluster_alg
        cluster_k - only applicable when using "phenograph" as cluster_alg.
        min_cluster_size - only applicable when using "phenograph" as cluster_alg
            Make sure that cluster_k < min_cluster_size
        redo_pca - boolean. whether to re-calculate PCA for subclusters
        verbose - boolean
    Returns adata with new clustering (under adata.obs[cluster_label_new].
    """
    # copy original clustering
    clusters_previous = adata.obs[cluster_label_previous].tolist()
    adata.obs[cluster_label_new] = clusters_previous
    if not redo_pca:
        print("Not re-doing pca before nested clustering iterations!")
    for cluster in sorted(set(clusters_previous)):
        if verbose:
            print("Cluster:", cluster)
        subadata = adata[adata.obs[cluster_label_previous] == cluster, :].copy()
        if subadata.shape[0] < min_cluster_size:
            if verbose:
                print("cluster size smaller than", min_cluster_size, "\n")
            continue
        if verbose:
            print("reclustering...\n")
        if redo_pca:
            if verbose:
                print("running pca...")
            sc.tl.pca(subadata)
        if cluster_alg == "leiden":
            if verbose:
                print("calculating 30 nearest neighbors")
                print("using rep:", use_rep)
            sc.pp.neighbors(subadata, n_neighbors=30, use_rep=use_rep)
            if verbose:
                print("clustering")
            sc.tl.leiden(subadata, resolution=cluster_res, key_added=cluster_label_new)
        elif cluster_alg == "phenograph":
            subadata.obs[cluster_label_new] = pd.Categorical(
                sce.tl.phenograph(subadata.obsm[use_rep], k=cluster_k)[0]
            )
        else:
            raise ValueError("Your cluster_alg argument is incorrect.")
        subadata.obs[cluster_label_new] = [
            "{}.{}".format(cluster, new_cluster)
            for new_cluster in subadata.obs[cluster_label_new]
        ]
        adata.obs.loc[subadata.obs.index, cluster_label_new] = subadata.obs[
            cluster_label_new
        ]
    # order categories "numerically" (so not 1, 10, 11 but 1, 2, 3... 10, 11):
    # convert all cluster names to strings, instead of a mix of strings and ints:
    adata.obs[cluster_label_new] = [
        str(clust) for clust in adata.obs[cluster_label_new]
    ]
    cluster_numbers = list(sorted(set(adata.obs[cluster_label_new])))
    prefix_cluster = [float(x.split(".")[0]) for x in cluster_numbers]
    cluster_numbers_ordered = [
        cluster_numbers[idx] for idx in np.argsort(prefix_cluster)
    ]
    adata.obs[cluster_label_new] = pd.Categorical(
        adata.obs[cluster_label_new], categories=cluster_numbers_ordered
    )

    return adata

def get_cluster_markers(adata, cluster_label, marker_ref, ngenes=100, verbose=True):
    """
    Calculates markers for every cluster, using either all other cells or 
    the parent cluster as a reference (i.e. for cluster 00.00.01, it 
    uses all clusters starting with 00.00 as reference. For cluster 
    00, it uses all cells as reference).
    sc.tl.rank_genes is used for marker gene calculation.

    Arguments:
    adata - AnnData object
    cluster_label - string
        label in adata.obs that contains nested-cluster names
    marker_ref - either "all" or "sisters". Which clusters to compare with.
    ngenes - number of marker genes to get per cluster
    
    Returns:
    cluster_markers - pd.DataFrame
        dataframe with, for each cluster, 100 highest scoring genes,
        plus matching logfc and adj pvalue
    """
    # input check:
    if marker_ref == "all":
        print("Doing one versus all differential expression analysis.")
    elif marker_ref == "sisters":
        print("Doing one versus sisters differential expression analysis.")
    else:
        raise ValueError("marker_ref argument should be set to either 'all' or 'sisters'.")
    # convert clusters to strings:
    adata.obs[cluster_label] = [str(cl) for cl in adata.obs[cluster_label]]
    # store cluster set
    clusters = sorted(set(adata.obs[cluster_label]))
    colnames_nested = [
        [clust + "_gene", clust + "_logfc", clust + "_pval_adj"] for clust in clusters
    ]
    colnames = [item for sublist in colnames_nested for item in sublist]
    cluster_markers = pd.DataFrame(index=range(100), columns=colnames)
    parents_tested = list()
    for clust in clusters:
        clust_depth = len(clust.split("."))
        if clust_depth == 1:
            parent = "all"
            if parent not in parents_tested:
                if verbose:
                    print("ranking genes for parent group", parent)
                parents_tested.append(parent)
                sc.tl.rank_genes_groups(adata, groupby=cluster_label, n_genes=ngenes)
                # store results for all clusters from this parent
                # i.e. all clusters of depth 1
                for d1_cluster in [
                    clust for clust in clusters if len(clust.split(".")) == 1
                ]:
                    # create a subdf that will allow us to sort genes per cluster
                    submarker_df = pd.DataFrame(
                        index=range(ngenes),
                        columns=[
                            d1_cluster + "_gene",
                            d1_cluster + "_logfc",
                            d1_cluster + "_pval_adj",
                        ],
                    )
                    submarker_df[d1_cluster + "_gene"] = adata.uns["rank_genes_groups"][
                        "names"
                    ][d1_cluster]
                    submarker_df[d1_cluster + "_logfc"] = adata.uns[
                        "rank_genes_groups"
                    ]["logfoldchanges"][d1_cluster]
                    submarker_df[d1_cluster + "_pval_adj"] = adata.uns[
                        "rank_genes_groups"
                    ]["pvals_adj"][d1_cluster]
                    # sort values:
                    submarker_df.sort_values(
                        by=[d1_cluster + "_pval_adj", d1_cluster + "_logfc"],
                        ascending=[True, False],
                        inplace=True,
                    )
                    submarker_df = submarker_df.reset_index().drop(columns="index")
                    # and add to big dataframe
                    cluster_markers.loc[
                        submarker_df.index, submarker_df.columns
                    ] = submarker_df.values
        else:
            parent = ".".join(clust.split(".")[: clust_depth - 1])
            if parent not in parents_tested:
                # depending on reference choice, use whole adata as reference
                # or only the parent cluster.
                if marker_ref == "all":
                    subadata = adata
                elif marker_ref == "sisters":
                    subadata = adata[[cl.startswith(parent) for cl in adata.obs[cluster_label]],:].copy()
                if verbose:
                    print("ranking genes for parent group", parent)
                parents_tested.append(parent)
                siblings = [c for c in clusters if c.startswith(parent)]
                if len(siblings) < 2 and marker_ref == "sisters":
                    print("Cluster {} has only one subcluster. Skipping DEA for this parent.".format(parent))
                else:
                    sc.tl.rank_genes_groups(subadata, groupby=cluster_label, groups=siblings, n_genes=ngenes)
                    for same_depth_sibling in [
                        sib for sib in siblings if len(clust.split(".")) == clust_depth
                    ]:
                        # create a subdf that will allow us to sort genes per cluster
                        submarker_df = pd.DataFrame(
                            index=range(ngenes),
                            columns=[
                                same_depth_sibling + "_gene",
                                same_depth_sibling + "_logfc",
                                same_depth_sibling + "_pval_adj",
                            ],
                        )
                        submarker_df[same_depth_sibling + "_gene"] = subadata.uns[
                            "rank_genes_groups"
                        ]["names"][same_depth_sibling]
                        submarker_df[same_depth_sibling + "_logfc"] = subadata.uns[
                            "rank_genes_groups"
                        ]["logfoldchanges"][same_depth_sibling]
                        submarker_df[same_depth_sibling + "_pval_adj"] = subadata.uns[
                            "rank_genes_groups"
                        ]["pvals_adj"][same_depth_sibling]
                        # sort values:
                        submarker_df.sort_values(
                            by=[
                                same_depth_sibling + "_pval_adj",
                                same_depth_sibling + "_logfc",
                            ],
                            ascending=[True, False],
                            inplace=True,
                        )
                        submarker_df = submarker_df.reset_index().drop(columns="index")
                        # add to big dataframe
                        cluster_markers.loc[
                            submarker_df.index, submarker_df.columns
                        ] = submarker_df.values
    return cluster_markers



def create_cluster_mapping_overview(
    adata,
    n_levels,
    cluster_label_to_decompose,
    cluster_label_to_count_prefix,
    min_fraction_for_dominancy=0.5,
    index_name=None,
):
    """Function to calculate for a new clustering, which clusters from an old
    clustering are the dominant ones in your new clustering (or vice versa).
    Args:
    adata - scanpy AnnData object
    n_levels - number of annotation levels (named "ann_level_[number]" in adata.obs)
    cluster_label_to_decompose - column name of cluster column in adata.obs
        for which we want to know of what clusters it consists
    cluster_label_to_count_prefix - column name (excluding level number) of
        clusters by which we want to define our cluster-to-decompose
    min_fraction_for_dominancy - minimum fraction of annotated cells to belong to
                                one cell type, to be called "dominant". Should be
                                higher than 0.5.
    index_name - name to give to index column

    Returns:
    cluster df - pandas dataframe with for each cluster information on what is the
                dominant cluster (if there is one), and the fraction of annotated
                cells belonging to the dominant cluster
    """
    # set up dataframe with one row per new cluster
    cluster_df = pd.DataFrame(
        index=adata.obs[cluster_label_to_decompose].cat.categories,
        columns=zip(
            [
                f"{cluster_label_to_count_prefix}{level}_dom_type"
                for level in range(1, n_levels + 1)
            ],
            [
                f"{cluster_label_to_count_prefix}{level}_dom_fract"
                for level in range(1, n_levels + 1)
            ],
        ),
    )
    # loop through cluster-to-count levels
    for level in range(1, n_levels + 1):
        cluster_to_count_level_name = f"{cluster_label_to_count_prefix}{level}"
        clust_cell_types = adata.obs.groupby(
            [cluster_label_to_decompose, cluster_to_count_level_name]
        ).agg({cluster_to_count_level_name: "count"})
        # convert to proportions
        clust_cell_types = clust_cell_types.groupby(level=0)[
            cluster_to_count_level_name
        ].apply(lambda x: x / float(x.sum()))
        # add "dominant" annotation:
        dominant_types = clust_cell_types.index[
            clust_cell_types > min_fraction_for_dominancy
        ]
        dominant_fractions = clust_cell_types[
            clust_cell_types > min_fraction_for_dominancy
        ]
        # copy dominant types to cluster_df:
        cluster_df.loc[
            dominant_types.get_level_values(0),
            f"{cluster_to_count_level_name}_dom_type",
        ] = dominant_types.get_level_values(1)
        # copy dominance fractions to cluster_df
        cluster_df.loc[
            dominant_fractions.index.get_level_values(0),
            f"{cluster_to_count_level_name}_dom_fract",
        ] = dominant_fractions.values
    if not pd.isnull(index_name):
        cluster_df.index.name = index_name
    return cluster_df
