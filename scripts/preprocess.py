import os
import numpy as np
import pybiomart
import pickle
import MI
import pandas as pd
import scipy.sparse as sparse


def change_accessions(ids, input_format, output_format, species, tmp):  # refseq->ensemble->entrez;
    if input_format != output_format:
        mart_file = '%s/biomart_%s%s_%s.ipage.pickle' % (tmp, species, input_format, output_format)
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


def get_expression_profile(expression_level, genes, expression_bins, input_format, output_format, species, tmp):
    df = pd.DataFrame({'genes': genes, 'expression_level': expression_level})
    df = df[df.iloc[:, 1].notna()]
    df = df.sort_values(by=df.columns[1])
    expression_level = np.array(df.iloc[:, 1])
    if np.count_nonzero(expression_level < 0) != 0 and np.count_nonzero(expression_level >= 0) != 0:
        left = MI.discretize(expression_level[expression_level < 0], expression_bins // 2)
        right = MI.discretize(expression_level[expression_level >= 0], expression_bins // 2)
        right += expression_bins // 2
        expression_profile = np.concatenate((left, right))
    else:
        expression_profile = MI.discretize(expression_level[expression_level >= 0], expression_bins // 2)

    genes = list(df.iloc[:, 0])
    genes = [gene.split('.')[0] for gene in genes]
    if input_format and output_format and input_format != output_format:
        genes = change_accessions(genes, input_format, output_format, species, tmp)
        gene_dict = dict(zip(genes, expression_profile))
        expression_profile = np.array([gene_dict[gene] for gene in gene_dict.keys() if gene != '-'])
        genes = [gene for gene in gene_dict.keys() if gene != '-']
    return expression_profile, genes


def get_profiles(db_index_file, first_col_is_genes, db_names_file=None):
    with open(db_index_file) as f:
        lines = filter(None, (line.rstrip() for line in f))
        id_a = set()
        id_b = set()
        for line in lines:
            id_a |= set([line.split('\t')[0]])
            id_b |= set(line.split('\t')[1:])
        id_a.add('')
        id_a.remove('')
        id_b.add('')
        id_b.remove('')
        id_a = list(id_a)
        id_b = list(filter(lambda x: 'http://' not in x, id_b))

    profiles = np.zeros((len(id_b), len(id_a)), dtype='float64')
    id_a_number = dict(zip(id_a, range(len(id_a))))
    id_b_number = dict(zip(id_b, range(len(id_b))))
    with open(db_index_file) as f:
        lines = filter(None, (line.rstrip() for line in f))
        for line in lines:
            ids = line.split()
            a = ids[0]
            for b in ids[1:]:
                if 'http://' not in b:
                    profiles[id_b_number[b]][id_a_number[a]] = 1

    if first_col_is_genes:
        db_genes = id_a
        db_names = id_b
    else:
        profiles = np.transpose(profiles)
        db_genes = id_b
        db_names = id_a
    if db_names_file:
        db_annotations = db_names.copy()
        with open(db_names_file) as f:
            lines = filter(None, (line.rstrip() for line in f))
            for line in lines:
                name = line.split('\t')[0]
                if name in db_names:
                    if len(line.split('\t')) > 1:
                        annotation = '; '.join(line.split('\t')[:2])
                    else:
                        annotation = name
                    db_annotations[db_names.index(name)] = annotation
    else:
        db_annotations = db_names

    genes_bool = np.sum(profiles, axis=0) != 0
    profiles = profiles[:, genes_bool]
    db_genes = [db_genes[i] for i in range(len(db_genes)) if genes_bool[i]]
    return db_names, profiles, db_annotations, db_genes


def dump_database(db_names, db_annotations, db_genes, db_profiles, database_name, tmp):
    with open("{0}/{1}.ipage.pickle".format(tmp, database_name), "wb+") as f:
        pickle.dump(db_names, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_annotations, f, pickle.HIGHEST_PROTOCOL)
        pickle.dump(db_genes, f, pickle.HIGHEST_PROTOCOL)
    sparse_profiles = sparse.csr_matrix(db_profiles)
    sparse.save_npz("{0}/{1}.ipage.npz".format(tmp, database_name), sparse_profiles, compressed=True)


def load_database(database_name, tmp):
    with open("{0}/{1}.ipage.pickle".format(tmp, database_name), 'rb') as f:
        db_names = pickle.load(f)
        db_annotations = pickle.load(f)
        db_genes = pickle.load(f)
    sparse_profiles = sparse.load_npz("{0}/{1}.ipage.npz".format(tmp, database_name))
    db_profiles = np.array(sparse_profiles.todense())
    return db_names, db_annotations, db_genes, db_profiles

# to be deleted
'''def sort_genes(genes, db_genes, expression_profile, db_profiles, delete_zero_genes=True):
    intersected_genes = list(set(genes) & set(db_genes))
    db_genes_bool = np.isin(db_genes, intersected_genes)
    db_genes = np.array(db_genes)[db_genes_bool].tolist()
    genes_bool = np.isin(genes, intersected_genes)
    genes = np.array(genes)[genes_bool].tolist()

    # slightly artificial maintenance of function polymorphism
    if len(expression_profile.shape) == 2:
        expression_profile = expression_profile.T
    if len(db_profiles.shape) == 2:
        db_profiles = db_profiles.T

    db_profiles = db_profiles[db_genes_bool]
    expression_profile = expression_profile[genes_bool]

    indices = np.array([db_genes.index(gene) for gene in genes])
    db_profiles = db_profiles[indices]

    if len(expression_profile.shape) == 2:
        expression_profile = expression_profile.T
    if len(db_profiles.shape) == 2:
        db_profiles = db_profiles.T

    if delete_zero_genes:
        non_zero_pos = np.where(db_profiles.sum(0) != 0)[0]
        expression_profile = expression_profile[non_zero_pos]
        db_profiles = db_profiles[:, non_zero_pos]
        genes = [genes[i] for i in non_zero_pos]

    return genes, expression_profile, db_profiles'''


def sort_genes(genes, db_genes, expression_profile, db_profiles, delete_zero_genes=True):
    genes_not_in_db_genes = set(genes) - set(db_genes)
    genes_not_in_genes = set(db_genes) - set(genes)
    genes += list(genes_not_in_genes)
    db_genes += list(genes_not_in_db_genes)

    expression_profile = np.atleast_2d(expression_profile)
    db_profiles = np.atleast_2d(db_profiles)

    expression_profile_supl = np.zeros((expression_profile.shape[0], len(genes_not_in_genes)))
    expression_profile_supl[:] = np.nan
    db_profiles_supl = np.zeros((db_profiles.shape[0], len(genes_not_in_db_genes)))
    db_profiles_supl[:] = np.nan
    expression_profile = np.concatenate((expression_profile, expression_profile_supl), axis=1)
    db_profiles = np.concatenate((db_profiles, db_profiles_supl), axis=1)

    indices = np.array([db_genes.index(gene) for gene in genes])
    db_profiles = db_profiles[:, indices]

    if expression_profile.shape[0] == 1:
        expression_profile = expression_profile.flatten()

    if delete_zero_genes:
        non_zero_pos = np.where(np.isnan(db_profiles).sum(0) == 0)[0]
        expression_profile = expression_profile[non_zero_pos]
        db_profiles = db_profiles[:, non_zero_pos]
        genes = [genes[i] for i in non_zero_pos]

        non_zero_pos = np.where(~np.isnan(expression_profile))[0]
        expression_profile = expression_profile[non_zero_pos]
        db_profiles = db_profiles[:, non_zero_pos]
        genes = [genes[i] for i in non_zero_pos]

    return genes, expression_profile, db_profiles



