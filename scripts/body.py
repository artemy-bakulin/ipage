import filter_db
import heatmap
import stat_ipage
import preprocess
import MI
import numpy as np
import pandas as pd
import pickle
import scipy.sparse as sparse


def preprocess_db(database_names_file, first_col_is_genes, database_index_file, filter_redundant, min_pathway_length,
                  child_unique_genes, tmp):
    database_name = database_index_file.split('/')[-1].split('.')[0]

    db_names, db_profiles, db_annotations, db_genes = preprocess.get_profiles(database_index_file, first_col_is_genes,
                                                                              database_names_file)

    if filter_redundant:
        db_names, db_annotations, db_profiles = filter_db.non_redundancy_sort_pre(db_names, db_annotations, db_profiles,
                                                                                  min_pathway_length,
                                                                                  child_unique_genes)
    with open("{0}/{1}.ipage.pickle".format(tmp, database_name), "wb+") as f:
        pickle.dump(db_names, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_annotations, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_genes, f, pickle.HIGHEST_PROTOCOL)
    sparse_profiles = sparse.csr_matrix(db_profiles)
    sparse.save_npz("{0}/{1}.ipage.npz".format(tmp, database_name), sparse_profiles, compressed=True)


def process_input(expression_level, genes, database_index_file, input_format, output_format, expression_bins, abundance_bins,
                  species, tmp):

    expression_profile, genes = preprocess.get_expression_profile(expression_level, genes, expression_bins,
                                                                  input_format, output_format, species, tmp)
    database_name = database_index_file.split('/')[-1].split('.')[0]
    with open("{0}/{1}.ipage.pickle".format(tmp, database_name), 'rb') as f:
        db_names = pickle.load(f)
        db_annotations = pickle.load(f)
        db_genes = pickle.load(f)
    sparse_profiles = sparse.load_npz("{0}/{1}.ipage.npz".format(tmp, database_name))
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

    # delete genes that do not belong to any regulon
    non_zero_pos = np.where(abundance_profile != 0)[0]
    expression_profile = expression_profile[non_zero_pos]
    db_profiles = db_profiles[:, non_zero_pos]
    abundance_profile = abundance_profile[non_zero_pos]
    genes = [genes[i] for i in non_zero_pos]

    abundance_profile = MI.discretize_equal_size(abundance_profile, abundance_bins)
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


def get_rbp_expression(genes, output_format, expression_profile, accepted_db_profiles, db_names, db_annotations, species, tmp):
    rbp_names = [db_names[i].replace('_', ' ').split()[0] for i in range(len(db_names)) if accepted_db_profiles[i]]
    rbp_annotations = [db_annotations[i] for i in range(len(db_names)) if accepted_db_profiles[i]]

    genes_symbols = preprocess.change_accessions(genes, output_format, 'gene_symbol', species, tmp)
    rbp_expression = dict(zip(rbp_annotations,
                              [expression_profile[genes_symbols.index(name)] for name in rbp_names
                               if name in genes_symbols]))
    rbp_expression.update({el: np.nan for el in rbp_annotations if el not in rbp_expression})
    return rbp_expression


def produce_output(accepted_db_profiles, db_profiles, db_names, db_annotations, cmis, z_scores,
                   draw_bins, max_draw_output, output_name, rbp_expression=None):
    p_values = {}
    for i in range(len(db_profiles)):
        if accepted_db_profiles[i]:
            p_values[db_annotations[i]] = stat_ipage.get_p_values(db_profiles[i], draw_bins)
    max_draw_output = min(max_draw_output, len(p_values))
    up_regulated_func = lambda x: sum(p_values[x][:len(p_values[x]) // 2]) <= sum(p_values[x][len(p_values[x]) // 2:])
    down_regulated_func = lambda x: sum(p_values[x][:len(p_values[x]) // 2]) >= sum(p_values[x][len(p_values[x]) // 2:])
    order_to_cmi = lambda x: cmis[db_annotations.index(x)]
    up_regulated = list(sorted(filter(up_regulated_func, p_values), key=order_to_cmi, reverse=True))
    up_regulated_portion = up_regulated[:max_draw_output // 2]
    down_regulated = list(sorted(filter(down_regulated_func, p_values), key=order_to_cmi, reverse=True))
    down_regulated_portion = down_regulated[:max_draw_output // 2]
    p_names = up_regulated_portion + down_regulated_portion

    for i in range(draw_bins-2):
        p_names = sorted(p_names, key=lambda x: sum(p_values[x][i:i+3]), reverse=True)
    p_values = [p_values[name] for name in p_names]
    if rbp_expression:
        rbp_expression = [rbp_expression[name] for name in p_names]

    if len(p_values) != 0 and output_name != 'shut':
        heatmap.draw_heatmap(p_names, p_values, output_name, rbp_expression)
    output = pd.DataFrame(columns=['Group', 'CMI', 'Z-score', 'Regulation'])
    j = 0
    for name in up_regulated:
        i = db_annotations.index(name)
        output.loc[j] = [db_names[i], str(cmis[i])[:8], str(z_scores[i]), 'UP']
        j += 1
    for name in down_regulated:
        i = db_annotations.index(name)
        output.loc[j] = [db_names[i], str(cmis[i])[:8], str(z_scores[i]), 'DOWN']
        j += 1
    if output_name != 'stdout' and output_name != 'shut':
        output.to_csv(output_name + '.out', index=False, sep='\t')
    return output
