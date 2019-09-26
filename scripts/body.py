from filter_db import *
from heatmap import *
from stat_ipage import *
from process_input import *
import MI
import numpy as np
def process_input(expression_file=None, database_names_file=None, database_index_file=None, input_format=None, output_format=None, expression_bins=0,filter_redundant=False):
    expression_level=get_expression_level(expression_file)
    genes, expression_profile=get_expression_profile(expression_level,expression_bins,input_format,output_format) 
    if database_names_file and database_index_file:
        db_names=get_names(database_names_file)
        db_profiles=get_profiles_2f(database_index_file,db_names, genes)
    elif database_index_file:
        db_indexes=get_indexes_1f(database_index_file)
        db_profiles=get_profiles_1f(db_indexes, genes)
        db_names={item:item for item in db_profiles.keys()}
    if filter_redundant:
        db_profiles=non_redundancy_sort_pre(db_profiles)
    abundance_profile=list(np.sum([db_profiles[name] for name in db_profiles],axis=0))
    abundance_bins=3
    abundance_profile=get_discrete_profile(abundance_profile,abundance_bins)
    return expression_profile, db_profiles, abundance_profile, db_names
def count_cmi_for_profiles(expression_profile,db_profiles,abundance_profile,expression_bins,db_bins,abundance_bins):
    cmi_dic={}
    for name in db_profiles:
        cmi = MI.cond_mut_info(expression_profile,db_profiles[name],abundance_profile, expression_bins, db_bins, abundance_bins)
        cmi_dic[name] = cmi
    return cmi_dic
def statistical_testing(cmi_dic,expression_profile, db_profiles, abundance_profile,expression_bins, db_bins, abundance_bins,filter_redundant=False):
    tested=0
    accepted=0
    cmi_output={}
    cmi_sorted=sorted(cmi_dic,key=lambda x: cmi_dic[x], reverse=True)
    for name in cmi_sorted:
        vector_accepted=test_cond_mi(expression_profile,db_profiles[name],abundance_profile, expression_bins, db_bins, abundance_bins)
        if vector_accepted:
            cmi_output[name]=cmi_dic[name]
            accepted+=1
        if accepted>50 or tested>200:
            break        
    if filter_redundant==True:
        cmi_output=non_redundancy_sort_post(cmi_output,db_profiles)
    return cmi_output

def visualize_output(cmi_output,db_profiles, db_names, draw_bins, max_draw_output):
    number=max_draw_output if len(cmi_output)>=max_draw_output else len(cmi_output)
    top_cmi=sorted(cmi_output, key=lambda x: cmi_output[x], reverse=True)
    p_values={}
    for name in top_cmi[:number]:
        p_values[db_names[name]]=get_p_values(db_profiles[name],draw_bins)
    draw_heatmap(p_values)