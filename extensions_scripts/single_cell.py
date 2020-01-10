import scipy
import pandas as pd
import numpy as np
from patsy import dmatrices
import statsmodels.api as sm
import statsmodels.formula.api as smf


def get_difexp(anndata, column, cluster_A_name, cluster_B_name, test='nb'):
    booleans_A = anndata.obs[column] == str(cluster_A_name)
    cells_A = anndata.raw.X.toarray()[booleans_A]

    booleans_B = anndata.obs[column] == str(cluster_B_name)
    cells_B = anndata.raw.X.toarray()[booleans_B]

    genes = list(anndata.raw.var.index)

    if test == 'nb':
        difexp = count_difexp_with_nb_regression(cells_A, cells_B)
    elif test == 't':
        difexp = count_difexp_with_t_test(cells_A, cells_B)
    elif test == 'u':
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
    difexp = sign*(1-p_values)
    return difexp


# Done it as it's explained here: https://dius.com.au/2017/08/03/using-statsmodels-glms-to-model-beverage-consumption/
def ct_response(row):
    "Calculate response observation for Cameron-Trivedi dispersion test"
    y = row['expression']
    m = row['expression_mu']
    return ((y - m)**2 - y) / m


def calculate_alpha(response, predictors, data):
    po_results = sm.GLM(response, predictors, family=sm.families.Poisson()).fit()
    ct_data = data.copy()
    ct_data['expression_mu'] = po_results.mu
    ct_data['ct_resp'] = ct_data.apply(ct_response, axis=1)
    ct_results = smf.ols('ct_resp ~ expression_mu - 1', ct_data).fit()
    # alpha_ci95 = ct_results.conf_int(0.05).loc['expression_mu']
    alpha = ct_results.params[0]
    return alpha


def calculate_p_value(data):
    formula = "expression ~ cluster"
    response, predictors = dmatrices(formula, data, return_type='dataframe')
    alpha = calculate_alpha(response, predictors, data)
    print(alpha)
    nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
    return nb_results.pvalues[1]


def count_difexp_with_nb_regression(cells_A, cells_B):
    pvalues = np.empty(cells_A.shape[1])
    for i in range(cells_A.shape[1]):
        X = cells_A.shape[0] * ['A'] + cells_B.shape[0] * ['B']
        Y = np.concatenate((cells_A[:, i], cells_B[:, i]))
        df = pd.DataFrame({'cluster': X, 'expression': Y})
        try:
            pvalues[i] = calculate_p_value(df)
        except:
            pvalues[i] = np.nan
        if (i + 1) % 1000 == 0:
            print('%.f000 genes processed' % i)

    sign = np.sign(cells_A.mean(axis=0) - cells_B.mean(axis=0))
    difexp = sign * (1 - pvalues)

    return difexp
