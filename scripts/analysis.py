import scanpy as sc
import pandas as pd
from scipy.stats import entropy
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import normalized_mutual_info_score
from sklearn.linear_model import LinearRegression
import utils # self-built module in same folder as this module


def calculate_cluster_entropy(
    adata,
    cluster_variable,
    obs_variable,
    min_frac_annotated=0.2,
    ignore_unannotated=True,
    verbose=True,
):
    """Calculates the per cluster entropy of a given obs variable.
    Arguments:
        adata - an anndata object
        cluster_variable - the name of the cluster variable in adata.obs
        obs_variable - the name of the variable of which to calculate the
            entropy. Ensure that unannotated cells are set to "None" (string)
            and not to None or NaN.
        min_frac_annotated - if less than this fraction of cells is annotated
            for a given cluster, the entropy will be set to "None"
        ignore_unannotated - (Boolean) whether to include or exclude cells with "None"
            annotation in the entropy calculation.
        verbose - (Boolean) whether to print relevant into.

    Returns:
        adata with an adata.obs variable
            "entropy_[obs_variable]_[cluster_variable]" containing the entropy
            values
    """
    # check if there are NaNs or Nones (non-string) in the observations
    # if so, throw error
    if pd.isnull(adata.obs[obs_variable]).any():
        raise ValueError(
            "There are None or NaN values in your obs_variable. Set these to 'None' (string) instead! Exiting."
        )
    cluster_obs_counts = (
        adata.obs.groupby([cluster_variable, obs_variable])
        .agg({obs_variable: "count"})
        .unstack()
    )
    # clean up column and index after unstacking:
    cluster_obs_counts.columns = cluster_obs_counts.columns.droplevel(0)
    cluster_obs_counts.columns.set_names(None, inplace=True)
    # calculate fraction of unannotated cells per cluster:
    if "None" in cluster_obs_counts.columns:
        frac_unannotated = cluster_obs_counts["None"] / cluster_obs_counts.sum(axis=1)
    else:
        frac_unannotated = pd.Series(np.zeros(cluster_obs_counts.shape[0]),index=cluster_obs_counts.index)
    # remove "None" cells (cells without a label for the obs_variable) from df:
    if ignore_unannotated:
        cluster_obs_counts = cluster_obs_counts.loc[
            :, cluster_obs_counts.columns != "None"
        ]
    # calculate row-wise entropies:
    entr = cluster_obs_counts.apply(entropy, axis=1)
    # set clusters with less than min_annotated_frac unannotated to None
    entr[frac_unannotated.loc[entr.index] > (1 - min_frac_annotated)] = None
    # map entropies to adata obs variable:
    entropy_obs_col = f"entropy_{obs_variable}_{cluster_variable}"
    if verbose:
        print(f"Storing entropy values in obs column {entropy_obs_col}")
    adata.obs[entropy_obs_col] = np.array(adata.obs[cluster_variable].map(entr).values)
    return adata


def calculate_cluster_mean(
    adata,
    cluster_variable,
    obs_variable,
    verbose=True,
):
    """Calculates the per cluster mean of a given obs variable.
    Arguments:
        adata - an anndata object
        cluster_variable - the name of the cluster variable in adata.obs
        obs_variable - the name of the variable of which to calculate the
            mean
        verbose - (Boolean) whether to print relevant into.

    Returns:
        adata - adata with an adata.obs variable
            "mean_[obs_variable]_[cluster_variable]" containing the mean
            values per cluster
        cluster_means - pandas dataframe with cluster mean values
    """
    cluster_means = adata.obs.groupby(cluster_variable).agg({obs_variable: "mean"})
    cl_to_mean_mapping = dict(zip(cluster_means.index, cluster_means[obs_variable]))
    obs_col = f"mean_{obs_variable}_{cluster_variable}"
    if verbose:
        print(f"Storing values in obs column {obs_col}")
    adata.obs[obs_col] = np.array(
        adata.obs[cluster_variable].map(cl_to_mean_mapping).values
    )
    return adata, cluster_means


def get_correlation_or_mi(cov1, cov2, obs_df, verbose=True):
    """Function to get either a linear correlation
    (sqrt(variance(predicted)/variance(observed))) or an normalized mututal
    information value quantifying the dependence of two variables.
    Arguments:
        cov1 - string of column name of covariate of interest 1
        cov1 - string of column name of covariate of interest 2
        obs_df - adata.obs dataframe with columns named as cov1 and cov2

    Returns:
        corrtype, corrvalue

        corrtype - type of correlation. Either lin_regr or nmi, or
                "single_category" if one of the variables has only one value.
        corrvalue - nmi value, or (sqrt(variance(predicted)/variance(observed))),
                or np.nan if one of the variables has only one value.

    """
    cov1_unfiltered = obs_df[cov1].values  # .copy()
    cov1_nans = np.vectorize(utils.check_if_nan)(cov1_unfiltered)
    #             var_explained.loc[comp, "overall"] = np.var(y_true_unfiltered)
    # check if covariate is categorical or numerical
    if cov1_unfiltered.dtype in ["float32", "float", "float64"]:
        cov1_type = "num"
    else:
        cov1_type = "cat"
        if len(set(cov1_unfiltered)) == 1:
            if verbose:
                print(cov1, "has only 1 category. Setting to nan.")
            corrtype = "single_category"
            corrvalue = np.nan
            return corrtype, corrvalue
        elif verbose:
            print(f"converting {cov1} to dummy variable")
    cov2_unfiltered = obs_df[cov2].values  # .copy()
    cov2_nans = np.vectorize(utils.check_if_nan)(cov2_unfiltered)
    if cov2_unfiltered.dtype in ["float32", "float", "float64"]:
        cov2_type = "num"
        # print that we are treating as numerical (only for first loop,
        # so that we don't print the same thing many times)
        if verbose:
            print(f"treating {cov2} as continuous variable")
    else:
        cov2_type = "cat"
        if len(set(cov2_unfiltered)) == 1:
            corrtype = "single_category"
            corrvalue = np.nan
            if verbose:
                print(cov2, "has only 1 category. Setting to nan.")
            return corrtype, corrvalue
        elif verbose:
            print(f"converting {cov2} to dummy variable")
    # now filter:
    cells_to_keep = [
        ~cov1_nan and ~cov2_nan for cov1_nan, cov2_nan in zip(cov1_nans, cov2_nans)
    ]
    cov1_values = cov1_unfiltered[cells_to_keep]
    cov2_values = cov2_unfiltered[cells_to_keep]
    # and covert to dummies if categoricals
    if cov1_type == "cat":
        cov1_values_undummied = cov1_values
        cov1_values = pd.get_dummies(cov1_values, drop_first=True)
    if cov2_type == "cat":
        cov2_values_undummied = cov2_values
        cov2_values = pd.get_dummies(cov2_values, drop_first=True)
    # now do the regression or calculate mutual information,
    # depending on the types of data:
    if cov1_type == "num" and cov2_type == "num":
        if verbose:
            print(
                f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, linear regr., {cov2} response var"
            )
        lrf = LinearRegression(fit_intercept=True).fit(
            cov1_values.reshape(-1, 1),
            cov2_values,
        )
        cov2_pred = lrf.predict(cov1_values.reshape(-1, 1))
        sqrt_frac_var_explained = np.sqrt(np.var(cov2_pred) / np.var(cov2_values))
        corrtype = "lin_regr"
        corrvalue = sqrt_frac_var_explained

    elif cov1_type == "cat" and cov2_type == "num":
        if verbose:
            print(
                f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, linear regr., {cov2} response var"
            )
        lrf = LinearRegression(fit_intercept=True).fit(
            cov1_values,
            cov2_values,
        )
        cov2_pred = lrf.predict(cov1_values)
        sqrt_frac_var_explained = np.sqrt(np.var(cov2_pred) / np.var(cov2_values))
        corrtype = "lin_regr"
        corrvalue = sqrt_frac_var_explained
        # same as:
    #         frac_std_explained = np.std(y_pred) / np.std(y_true)
    # same as:
    #         r, p = pearsonr(x, y_true)
    #         print(r)
    #         print(frac_std_explained, sqrt_frac_var_explained)
    elif cov1_type == "num" and cov2_type == "cat":
        if verbose:
            print(
                f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, linear regr., {cov1} response var"
            )
        lrf = LinearRegression(fit_intercept=True).fit(
            cov2_values,
            cov1_values,
        )
        cov1_pred = lrf.predict(cov2_values)
        sqrt_frac_var_explained = np.sqrt(np.var(cov1_pred) / np.var(cov1_values))
        corrtype = "lin_regr"
        corrvalue = sqrt_frac_var_explained
    elif cov1_type == "cat" and cov2_type == "cat":
        # check if one of the cats has only two categories (len(set)==1 after dummying).
        # This can be converted to numerical:
        if len(set(cov1_values)) == 1:
            if verbose:
                print(
                    f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, linear regr., {cov1} repsonse var, because it has only 2 cats"
                )
            lrf = LinearRegression(fit_intercept=True).fit(
                cov2_values,
                cov1_values,
            )
            cov1_pred = lrf.predict(cov2_values)
            sqrt_frac_var_explained = np.sqrt(np.var(cov1_pred) / np.var(cov1_values))
            corrtype = "lin_regr"
            corrvalue = sqrt_frac_var_explained[
                cov1_values.columns[0]
            ]  
        elif len(set(cov2_values)) == 1:
            if verbose:
                print(
                    f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, linear regr., {cov2} response var, because it has 2 cats"
                )
            lrf = LinearRegression(fit_intercept=True).fit(
                cov1_values,
                cov2_values,
            )
            cov2_pred = lrf.predict(cov1_values)
            sqrt_frac_var_explained = np.sqrt(np.var(cov2_pred) / np.var(cov2_values))
            corrtype = "lin_regr"
            corrvalue = sqrt_frac_var_explained[cov2_values.columns[0]]  
        else:
            # mutual info
            if verbose:
                print(f"{cov1}: {cov1_type}, {cov2}: {cov2_type}, mutual info")
            nmi = normalized_mutual_info_score(
                cov1_values_undummied, cov2_values_undummied
            )
            corrtype = "nmi"
            corrvalue = nmi
    else:
        raise ValueError("something went wrong with cat and num ass.")
    return corrtype, corrvalue
