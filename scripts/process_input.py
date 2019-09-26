import os
import numpy as np
import pybiomart
import pickle
def change_accessions(ids,input_format,output_format): #refseq->ensemble->entrez;
    if input_format!=output_format:
        mart_file='biomart_%s_%s.pickle' % (input_format,output_format)
        if os.path.isfile(mart_file) and os.stat(mart_file).st_size != 0:
            with open(mart_file, 'rb') as f:
                input_to_output = pickle.load(f)

        else:
            dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl',host='http://www.ensembl.org')
            mart_attributes={'enst':['ensembl_transcript_id'],'ensg':['ensembl_gene_id'],'refseq':['refseq_mrna','refseq_mrna_predicted','refseq_ncrna','refseq_ncrna_predicted'],'entrez':['entrezgene_id'],'gene_symbol':['entrezgene_accession']}
            input_to_output={}
            for mart in mart_attributes[input_format]:
                df1=dataset.query(attributes=[mart]+mart_attributes[output_format])
                df1=df1[df1.iloc[:,0].notna()]
                df1=df1[df1.iloc[:,1].notna()]
                if input_format=='entrez' or output_format=='entrez':
                    df1['NCBI gene ID']=df1['NCBI gene ID'].apply(lambda x:'%.f'% x)
                input_to_output = {**input_to_output, **dict(zip(df1.iloc[:,0],df1.iloc[:,1]))}
            with open(mart_file, 'wb') as f:
                pickle.dump(input_to_output, f, pickle.HIGHEST_PROTOCOL)
        new_ids=[]
        for id_ in ids:
            try:
                new_ids.append(input_to_output[id_])
            except:
                new_ids.append('-')
        return new_ids
    else:
        return ids
def get_expression_level(expression_file): #get_expression_level
    expression_level={}
    with open(expression_file) as f:
        next(f)
        lines = filter(None, (line.rstrip() for line in f))
        for line in lines:
            expression_level[line.split()[0]]=float(line.split()[1])
        return expression_level

#first variant which makes nearly equal bins
def get_discrete_profile(v,nbins):
    bin_size=len(v)//nbins
    unique=set(v)
    unique.discard(9999)
    unique=sorted(unique)
    bins=[]
    cat=0
    value_score={}
    filling=0
    for i in unique:
        filling+=v.count(i)
        value_score[i]=str(cat)
        if filling>=bin_size and cat<nbins-1:
            filling=0
            cat+=1
    value_score[9999]=str(9999)  
    return np.array([int(i) for i in [value_score[i] for i in v]])

def get_expression_profile(expression_level,nbins,input_format, output_format):
    genes=sorted(expression_level,key=lambda x : expression_level[x])
    genes=change_accessions(genes,input_format,output_format)
    expression=sorted(expression_level.values())
    profile=get_discrete_profile(expression,nbins)
    return genes, profile

def get_names(db_names_file):
    db_names={}
    with open(db_names_file) as f:
        for line in f:
            name=line.split()[0]
            db_names[name]=' '.join(line.split()[1:-1])+'; '+name
    return db_names

def get_indexes_1f(db_file):
    db_indexes={}
    with open(db_file) as f:
        lines = filter(None, (line.rstrip() for line in f))
        for line in lines:
            db_indexes[line.split('\t')[0]]=list(filter(lambda x:'http://' not in x, line.split('\t')[2:]))
    return db_indexes

def get_profiles_1f(indexes, genes):
    profiles={}
    for name in indexes:
        #print(indexes)
        profiles[name]=np.array([1 if gene in indexes[name] else 0 for gene in genes])
    return profiles

def get_profiles_2f(db_index_file,db_names, genes):
    x=np.zeros((len(db_names), len( genes)), dtype='i1')
    gene_number=dict(zip(genes,range(len(genes))))
    ann_number=dict(zip(db_names.keys(),range(len(db_names))))
    with open(db_index_file) as f:
        lines = filter(None, (line.rstrip() for line in f))
        for line in lines:
            ids=line.split()
            transcript=ids[0]
            if transcript in genes:
                for id_ in ids[1:]:
                    if id_ in db_names.keys():
                        x[ann_number[id_]][gene_number[transcript]]=1
    x=dict(zip(db_names.keys(),x))
    return x