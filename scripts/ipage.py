import sys
import argparse
from body import *
if __name__=='__main__':
    input_=sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-e_file','--expression_file', type=str,
                        help='The name of the expression file')
    parser.add_argument('-db_names','--database_names_file', type=str,
                        help='The name of the database\'s file with names')
    parser.add_argument('-db_index','--database_index_file', type=str,
                        help='The name of the database\'s file with indexes')
    parser.add_argument('-e_bins','--expression_bins', type=int, default=10,
                        help='The number of the expression bins')
    parser.add_argument('-d_bins','--draw_bins', type=int, default=15,
                        help='The number of the bins in the heatmap')
    parser.add_argument('-db_format','--database_format', type=str,
                        help='The database format')
    parser.add_argument('-e_format','--expression_file_format', type=str,
                        help='The expression file format')
    parser.add_argument('-max_draw','--max_draw_output', type=int, default=50,
                        help='The number of results to be presented on the heatmap')
    parser.add_argument('-f','--filter_redundant', action="store_true",
                        help='Filter redundant pathways')
    args = parser.parse_args(input_)

    db_bins=2
    abundance_bins=3
    expression_bins=args.expression_bins 
    draw_bins=args.draw_bins 
    filter_redundant=args.filter_redundant
    max_draw_output=args.max_draw_output

    input_args={
          'expression_file':args.expression_file,
          'database_names_file':args.database_names_file,
          'database_index_file':args.database_index_file,
          'expression_bins':expression_bins,
          'input_format':args.expression_file_format,
          'output_format':args.database_format,
          'filter_redundant':filter_redundant
          }

    expression_profile, db_profiles, abundance_profile, db_names=process_input(**input_args)
    cmi_dic=count_cmi_for_profiles(expression_profile,db_profiles,abundance_profile,expression_bins,db_bins,abundance_bins)
    cmi_output=statistical_testing(cmi_dic, expression_profile, db_profiles, abundance_profile,expression_bins,db_bins,abundance_bins,filter_redundant=filter_redundant)
    visualize_output(cmi_output,db_profiles, db_names, draw_bins, max_draw_output)
    with open('task_output','w+') as f:
    	for el in sorted(cmi_output,key=lambda x:cmi_dic[x],reverse=True):
    		f.write(el+'\t'+str(cmi_output[el])+'\n')


