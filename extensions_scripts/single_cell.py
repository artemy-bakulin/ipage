import scipy
import numpy as np


def get_difexp(anndata, column, cluster_A_name, cluster_B_name):
    booleans_A = anndata.obs[column] == str(cluster_A_name)
    cells_A = anndata.raw.X.toarray()[booleans_A]

    booleans_B = anndata.obs[column] == str(cluster_B_name)
    cells_B = anndata.raw.X.toarray()[booleans_B]

    genes = list(anndata.raw.var.index)

    difexp = count_difexp_with_u_test(cells_A, cells_B)
    return genes, difexp


def count_difexp_with_t_test(cells_A, cells_B):
    sign = np.sign(cells_A.mean(axis=0)-cells_B.mean(axis=0))
    p_values = scipy.stats.ttest_ind(cells_A,cells_B,  axis=0)[1]
    difexp = sign*(1-p_values)
    return difexp


def count_difexp_with_u_test(cells_A, cells_B):
    n_genes = cells_A.shape[1]
    sign = np.sign(cells_A.mean(axis=0)-cells_B.mean(axis=0))
    p_values = np.empty(n_genes)
    p_values[:] = np.nan
    for i in range(n_genes):
        if not np.array_equal(np.unique(cells_A[:,i]), np.unique(cells_B[:,i])):
            p_values[i] = scipy.stats.mannwhitneyu(cells_A[:,i],cells_B[:,i])[1]
    difexp=sign*(1-p_values)
    return difexp
