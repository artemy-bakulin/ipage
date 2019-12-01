import filter_db
import heatmap
import stat_ipage
import preprocess
import MI
import numpy as np
import pickle
import scipy.sparse as sparse


def preprocess_db(database_names_file, first_col_is_genes, database_index_file, filter_redundant, min_pathway_length,
                  child_unique_genes,
                  parent_unique_genes, tmp='tmp_ipage'):
    database_name = database_index_file.split('/')[-1].split('.')[0]
    db_file = '%s.ipage' % database_name

    db_names, db_profiles, db_annotations, db_genes = preprocess.get_profiles(database_index_file, first_col_is_genes,
                                                                              database_names_file)

    if filter_redundant:
        db_names, db_annotations, db_profiles = filter_db.non_redundancy_sort_pre(db_names, db_annotations, db_profiles,
                                                                                  min_pathway_length,
                                                                                  child_unique_genes,
                                                                                  parent_unique_genes)
    with open("{0}/{1}.pickle".format(tmp, db_file), "wb") as f:
        pickle.dump(db_names, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_annotations, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_genes, f, pickle.HIGHEST_PROTOCOL)
    sparse_profiles = sparse.csr_matrix(db_profiles)
    sparse.save_npz("{0}/{1}.npz".format(tmp, db_file), sparse_profiles, compressed=True)


def process_input(expression_file, database_name, input_format, output_format, expression_bins=10, abundance_bins=3,
                  sep='\t', expression_column=1, tmp='tmp_ipage'):
    genes, expression_profile = preprocess.get_expression_profile(expression_file, expression_bins,
                                                                  input_format=input_format,
                                                                  output_format=output_format, sep=sep,
                                                                  expression_column=expression_column, tmp=tmp)

    db_file = '%s.ipage' % database_name

    with open("{0}/{1}.pickle".format(tmp, db_file), 'rb') as f:
        db_names = pickle.load(f)
        db_annotations = pickle.load(f)
        db_genes = pickle.load(f)
    sparse_profiles = sparse.load_npz("{0}/{1}.npz".format(tmp, db_file))
    db_profiles = np.array(sparse_profiles.todense())
    intersected_genes = set(genes) & set(db_genes)
    db_genes_bool = [gene in intersected_genes for gene in db_genes]
    db_genes = [gene for gene in db_genes if gene in intersected_genes]
    genes_bool = [gene in intersected_genes for gene in genes]
    genes = [gene for gene in genes if gene in intersected_genes]

    db_profiles = db_profiles[:, db_genes_bool]
    expression_profile = expression_profile[genes_bool]

    indices = np.array([db_genes.index(gene) for gene in genes])
    db_profiles = db_profiles[:, indices]
    abundance_profile = db_profiles.sum(axis=0)
    abundance_profile = MI.discretize(abundance_profile, abundance_bins)
    return expression_profile, db_names, db_profiles, db_annotations, abundance_profile, genes


def count_cmi_for_profiles(expression_profile, db_profiles, abundance_profile, expression_bins, db_bins,
                           abundance_bins):
    cmis = []
    for profile in db_profiles:
        cmi = MI.cond_mut_info(expression_profile, profile, abundance_profile, expression_bins, db_bins, abundance_bins)
        cmis.append(cmi)
    cmis = np.array(cmis)
    return cmis


def statistical_testing(cmis, expression_profile, db_profiles, abundance_profile, expression_bins, db_bins,
                        abundance_bins):
    indices = np.argsort(cmis)[::-1]
    rev_indices = np.argsort(indices)
    db_profiles_ = db_profiles[indices]
    accepted_db_profiles = np.array([False] * len(db_profiles))
    z_scores = np.array([0] * len(db_profiles))
    i = 0
    false_hits = 0
    for profile in db_profiles_:
        z_score, vector_accepted = stat_ipage.test_cond_mi(expression_profile, profile, abundance_profile,
                                                           expression_bins,
                                                           db_bins, abundance_bins, shuffles=1000)
        if not vector_accepted:
            accepted_db_profiles[i] = False
            false_hits += 1
        else:
            z_scores[i] = z_score
            accepted_db_profiles[i] = True
            false_hits = 0
        if false_hits > 5:
            break
        i += 1
    return accepted_db_profiles[rev_indices], z_scores[rev_indices]


def get_rbp_expression(genes, input_format, expression_profile, accepted_db_profiles, db_annotations):
    rbp_names = [db_annotations[i] for i in range(len(db_annotations))
                 if accepted_db_profiles[i]]
    genes_symbols = preprocess.change_accessions(genes, input_format, 'gene_symbol')
    di = {el1: el2 for el1 in genes_symbols for el2 in rbp_names if el1 in el2}

    rbp_expression = dict(zip([di[gene] for gene in genes_symbols if any(gene in name for name in rbp_names)],
                              expression_profile[[any(gene in name for name in rbp_names) for gene in genes_symbols]]))
    medium = sum(rbp_expression.values()) / len(rbp_expression)  # this is due to the possible lack of information
    rbp_expression.update({el: medium for el in rbp_names if el not in rbp_expression})
    return rbp_expression


def get_rbp_expression(genes, output_format, expression_profile, accepted_db_profiles, db_annotations):
    rbp_names = [db_annotations[i] for i in range(len(db_annotations))
                 if accepted_db_profiles[i]]
    rbp_names_clear = [name.replace('_', ' ').split()[0] for name in rbp_names]
    genes_symbols = preprocess.change_accessions(genes, output_format, 'gene_symbol')

    rbp_expression = dict(zip(rbp_names,
                              [expression_profile[genes_symbols.index(name)] for name in rbp_names_clear if
                               name in genes_symbols]))

    if len(rbp_expression) != 0:
        medium = sum(rbp_expression.values()) / len(rbp_expression)  # this is due to the possible lack of information
    else:
        medium = 0
    rbp_expression.update({el: medium for el in rbp_names if el not in rbp_expression})
    return rbp_expression


def visualize_output(accepted_db_profiles, db_profiles, db_annotations, cmis, draw_bins, max_draw_output, output_name,
                     rbp_expression=None):
    p_values = {}
    for i in range(len(db_profiles)):
        if accepted_db_profiles[i]:
            p_values[db_annotations[i]] = stat_ipage.get_p_values(db_profiles[i], draw_bins)
    max_draw_output = min(max_draw_output, len(p_values))
    up_regulated_func = lambda x: sum(p_values[x][:len(p_values[x]) // 2]) <= sum(p_values[x][len(p_values[x]) // 2:])
    down_regulated_func = lambda x: sum(p_values[x][:len(p_values[x]) // 2]) >= sum(p_values[x][len(p_values[x]) // 2:])
    order_to_cmi = lambda x: cmis[db_annotations.index(x)]
    up_regulated = list(sorted(filter(up_regulated_func, p_values), key=order_to_cmi, reverse=True))[
                   :max_draw_output // 2]
    down_regulated = list(sorted(filter(down_regulated_func, p_values), key=order_to_cmi, reverse=True))[
                     :max_draw_output // 2]
    p_names = up_regulated + down_regulated

    for i in range(draw_bins):
        p_names = sorted(p_names, key=lambda x: p_values[x][i], reverse=True)
    p_values = [p_values[name] for name in p_names]
    if rbp_expression:
        rbp_expression = [rbp_expression[name] for name in p_names]

    if len(p_values) != 0:
        heatmap.draw_heatmap(p_names, p_values, output_name, rbp_expression)
