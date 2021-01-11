import os
import numpy as np
import pybiomart
import pickle
import MI
import pandas as pd
import scipy.sparse as sparse


def change_accessions(ids, input_format, output_format, species, tmp):  # refseq->ensemble->entrez;
    if input_format != output_format:
        mart_file = 'biomart_%s_%s_%s.ipage.pickle' % (species, input_format, output_format)
        mart_file = os.path.join(tmp, mart_file)
        if os.path.isfile(mart_file) and os.stat(mart_file).st_size != 0:
            with open(mart_file, 'rb') as f:
                input_to_output = pickle.load(f)

        else:
            if species == 'mouse':
                dataset = pybiomart.Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
            elif species == 'human':
                dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
            # print(*dataset.attributes.keys(), sep='\n')
            mart_attributes = {'enst': ['ensembl_transcript_id'], 'ensg': ['ensembl_gene_id'],
                               'refseq': ['refseq_mrna', 'refseq_mrna_predicted', 'refseq_ncrna',
                                          'refseq_ncrna_predicted'], 'entrez': ['entrezgene_id'],
                               'gs': ['entrezgene_accession'], 'ext': ['external_gene_name']}
            input_to_output = {}
            output_attributes = mart_attributes[output_format]
            if output_format == 'refseq':
                output_attributes = [output_attributes[0]]
            for mart in mart_attributes[input_format]:
                df1 = dataset.query(attributes=[mart] + output_attributes)
                df1 = df1[df1.iloc[:, 0].notna()]
                df1 = df1[df1.iloc[:, 1].notna()]
                if input_format == 'entrez' or output_format == 'entrez':
                    df1['NCBI gene ID'] = df1['NCBI gene ID'].apply(lambda x: '%.f' % x)
                if input_format == 'gene_symbol' or output_format == 'gene_symbol':
                    upper = lambda x: x.upper() if type(x) == str else x
                    df1['NCBI gene accession'] = df1['NCBI gene accession'].apply(upper)
                input_to_output = {**input_to_output, **dict(zip(df1.iloc[:, 0], df1.iloc[:, 1]))}
            with open(mart_file, 'wb') as f:
                pickle.dump(input_to_output, f, pickle.HIGHEST_PROTOCOL)
        new_ids = []
        for id_ in ids:
            if id_ in input_to_output.keys():
                new_ids.append(input_to_output[id_])
            else:
                new_ids.append('-')
        return new_ids
    else:
        return ids


def get_expression_profile(expression_level, genes, expression_bins, input_format, output_format,
                           species, tmp, symmetric_expression):
    df = pd.DataFrame({'genes': genes, 'expression_level': expression_level})
    df = df[df.iloc[:, 1].notna()]
    df = df.sort_values(by=df.columns[1])
    expression_level = np.array(df.iloc[:, 1])
    if symmetric_expression:
        left = MI.discretize(expression_level[expression_level < 0], expression_bins // 2)
        right = MI.discretize(expression_level[expression_level >= 0], expression_bins // 2 + expression_bins % 2)
        right += expression_bins // 2
        expression_profile = np.concatenate((left, right))
    else:
        expression_profile = MI.discretize(expression_level, expression_bins)

    genes = list(df.iloc[:, 0])
    genes = [gene.split('.')[0] for gene in genes]
    if input_format and output_format and input_format != output_format:
        genes = change_accessions(genes, input_format, output_format, species, tmp)
        gene_dict = dict(zip(genes, expression_profile))
        expression_profile = np.array([gene_dict[gene] for gene in gene_dict.keys() if gene != '-'])
        genes = [gene for gene in gene_dict.keys() if gene != '-']
    return expression_profile, genes





def get_profiles_from_table(table, sep, first_col_is_genes=True):
    df = pd.read_csv(table, sep=sep, index_col=0)
    if first_col_is_genes:
        df = df.T
    df[df != 1] = 0
    db_names, db_genes = list(df.index), list(df.columns)
    db_annotations = db_names
    db_profiles = np.array(df)
    db_genes = [el.split(',')[0] for el in db_genes]
    return db_names, db_profiles, db_annotations, db_genes


def get_profiles(db_index_file, first_col_is_genes, db_names_file=None):
    row_names = []
    column_names = set()
    for line in open(db_index_file):
        els = line.rstrip().split('\t')
        row_names.append(els[0])
        els.pop(0)
        if 'http://' in els[0]:
            els.pop(0)
        column_names |= set(els)
    column_names = list(column_names)

    db_profiles = np.zeros((len(row_names), len(column_names)), dtype=int)

    i = 0
    for line in open(db_index_file):
        els = line.rstrip().split('\t')[1:]
        if 'http://' in els[0]:
            els.pop(0)
        indices = [column_names.index(el) for el in els]
        db_profiles[i, indices] = 1
        i += 1

    if first_col_is_genes:
        db_profiles = db_profiles.T
        db_genes = row_names
        db_names = column_names
    else:
        db_genes = column_names
        db_names = row_names

    if db_names_file:
        df_annotations = pd.read_csv(db_names_file, sep='\t', header=None, index_col=0)
        df_annotations = df_annotations.reindex(db_names)
        db_annotations = list(df_annotations.iloc[:, 0])
        db_annotations = [pair[0] + '; ' + pair[1] for pair in zip(db_names, db_annotations)]
    else:
        db_annotations = db_names
    return db_names, db_profiles, db_annotations, db_genes


def dump_database(db_names, db_annotations, db_genes, db_profiles, database_name, tmp):
    pickle_file = os.path.join(tmp, "{0}.ipage.pickle".format(database_name))
    with open(pickle_file, "wb+") as f:
        pickle.dump(db_names, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_annotations, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_genes, f, pickle.HIGHEST_PROTOCOL)
    npz_file = os.path.join(tmp, "{0}.ipage.npz".format(database_name))
    sparse_profiles = sparse.csr_matrix(db_profiles)
    sparse.save_npz(npz_file, sparse_profiles, compressed=True)


def load_database(database_name, tmp):
    pickle_file = os.path.join(tmp, "{0}.ipage.pickle".format(database_name))
    with open(pickle_file, 'rb') as f:
        db_names = pickle.load(f)
        db_annotations = pickle.load(f)
        db_genes = pickle.load(f)
    npz_file = os.path.join(tmp, "{0}.ipage.npz".format(database_name))
    sparse_profiles = sparse.load_npz(npz_file)
    db_profiles = np.array(sparse_profiles.todense())
    return db_names, db_annotations, db_genes, db_profiles


def sort_genes(genes, db_genes, expression_profile, db_profiles, delete_genes_not_in_expression=True,
               delete_genes_not_in_db=False):

    expression_profile = np.atleast_2d(expression_profile)
    db_profiles = np.atleast_2d(db_profiles)

    genes, unique_e_inds = np.unique(genes, return_index=True)
    genes = genes.tolist()
    expression_profile = expression_profile[:, unique_e_inds]

    db_genes, unique_db_inds = np.unique(db_genes, return_index=True)
    db_genes = db_genes.tolist()
    db_profiles = db_profiles[:, unique_db_inds]

    genes_not_in_db_genes = set(genes) - set(db_genes)
    genes_not_in_genes = set(db_genes) - set(genes)
    genes += list(genes_not_in_genes)
    db_genes += list(genes_not_in_db_genes)

    expression_profile_supl = np.zeros((expression_profile.shape[0], len(genes_not_in_genes)))
    db_profiles_supl = np.zeros((db_profiles.shape[0], len(genes_not_in_db_genes)))
    expression_profile = np.concatenate((expression_profile, expression_profile_supl), axis=1)
    db_profiles = np.concatenate((db_profiles, db_profiles_supl), axis=1)

    db_indices = np.argsort(db_genes)
    indices = np.argsort(genes)
    genes = sorted(genes)
    expression_profile = expression_profile[:, indices]
    db_profiles = db_profiles[:, db_indices]

    if delete_genes_not_in_expression:
        indices = np.where(~np.isin(genes, list(genes_not_in_genes)))[0]
        expression_profile = expression_profile[:, indices]
        db_profiles = db_profiles[:, indices]
        genes = [genes[i] for i in indices]

    if delete_genes_not_in_db:
        db_indices = np.where(~np.isin(genes, genes_not_in_db_genes))[0]
        expression_profile = expression_profile[:, db_indices]
        db_profiles = db_profiles[:, db_indices]
        genes = [genes[i] for i in db_indices]

    if expression_profile.shape[0] == 1:
        expression_profile = expression_profile.flatten()
        ranking_indices = np.argsort(expression_profile)
        expression_profile = expression_profile[ranking_indices]
        db_profiles = db_profiles[:, ranking_indices]
        genes = [genes[i] for i in ranking_indices]

    return genes, expression_profile, db_profiles

