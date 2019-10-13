import os
import numpy as np
import pybiomart
import pickle
import MI
import pandas as pd


def change_accessions(ids, input_format, output_format):  # refseq->ensemble->entrez;
    if input_format != output_format:
        mart_file = 'biomart_%s_%s.ipage.pickle' % (input_format, output_format)
        if os.path.isfile(mart_file) and os.stat(mart_file).st_size != 0:
            with open(mart_file, 'rb') as f:
                input_to_output = pickle.load(f)

        else:
            dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
            mart_attributes = {'enst': ['ensembl_transcript_id'], 'ensg': ['ensembl_gene_id'],
                               'refseq': ['refseq_mrna', 'refseq_mrna_predicted', 'refseq_ncrna',
                                          'refseq_ncrna_predicted'], 'entrez': ['entrezgene_id'],
                               'gene_symbol': ['entrezgene_accession']}
            input_to_output = {}
            for mart in mart_attributes[input_format]:
                df1 = dataset.query(attributes=[mart] + mart_attributes[output_format])
                df1 = df1[df1.iloc[:, 0].notna()]
                df1 = df1[df1.iloc[:, 1].notna()]
                if input_format == 'entrez' or output_format == 'entrez':
                    df1['NCBI gene ID'] = df1['NCBI gene ID'].apply(lambda x: '%.f' % x)
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


def get_expression_profile(expression_file, nbins=10, sep='\t', id_column=0, expression_column=1, input_format=None,
                           output_format=None):
    df = pd.read_csv(expression_file, sep=sep)
    df = df.sort_values(by=df.columns[expression_column])
    expression_level = np.array(df.iloc[:, expression_column])
    expression_profile = MI.discretize(expression_level, nbins)
    genes = list(df.iloc[:, id_column])
    if input_format and output_format and input_format != output_format:
        genes = change_accessions(genes, input_format, output_format)
        gene_dict = dict(zip(genes, expression_profile))
        expression_profile = np.array([gene_dict[gene] for gene in gene_dict.keys() if gene != '-'])
        genes = [gene for gene in gene_dict.keys() if gene != '-']
    return genes, expression_profile


def get_profiles(db_index_file, first_col_is_genes, db_names_file=None):
    with open(db_index_file) as f:
        lines = filter(None, (line.rstrip() for line in f))
        id_a = set()
        id_b = set()
        for line in lines:
            id_a |= set(line.split('\t')[0])
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
        db_annotations = [''] * len(db_names)
        with open(db_names_file) as f:
            lines = filter(None, (line.rstrip() for line in f))
            for line in lines:
                if len(line.split('\t')) == 3:
                    annotation = line.split('\t')[1]
                else:
                    annotation = line.split('\t')[0]
                db_annotations[db_names.index(line.split('\t')[0])] = annotation
    else:
        db_annotations = db_names
    genes_bool = np.sum(profiles, axis=0) != 0
    profiles = profiles[:, genes_bool]
    db_genes = [db_genes[i] for i in range(len(db_genes)) if genes_bool[i]]
    return db_names, profiles, db_annotations, db_genes
