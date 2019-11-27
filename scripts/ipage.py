import sys
import argparse
import body
import os
import pandas as pd

if __name__ == '__main__':
    input_ = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--expression_file', type=str,
                        help='The name of the expression file')

    parser.add_argument('-n', '--database_names_file', type=str,
                        help='The name of the database\'s file with names')
    parser.add_argument('-i', '--database_index_file', type=str,
                        help='The name of the database\'s file with indexes')

    parser.add_argument('-e_sep', '--separator_in_expression_file', type=str, default='\t',
                        help='The type of the separator in expression file')
    parser.add_argument('-e_col', '--column_with_stability', type=str, default='2',
                        help='Column in expression file of required stability values')

    parser.add_argument('-db_format', '--database_format', type=str,
                        help='The database format')
    parser.add_argument('-e_format', '--expression_file_format', type=str,
                        help='The expression file format')

    parser.add_argument('-o', '--output_name', type=str,
                        help='The name of the output file')

    parser.add_argument('-e_bins', '--expression_bins', type=int, default=10,
                        help='The number of the expression bins')
    parser.add_argument('-s_bins', '--sum_bins', type=int, default=3,
                        help='The number of bins in the sum profile')
    parser.add_argument('-h_bins', '--draw_bins', type=int, default=15,
                        help='The number of the bins in the heatmap')

    parser.add_argument('-max_draw', '--max_draw_output', type=int, default=50,
                        help='The number of results to be presented on the heatmap')

    parser.add_argument('-f', '--filter_redundant', action="store_true",
                        help='Filter redundant pathways')
    parser.add_argument('-g', '--first_col_is_genes', action="store_true",
                        help='First column is genes')
    parser.add_argument('-preprocess', '--preprocess', action="store_true",
                        help='Preprocess database')

    parser.add_argument('-clip ', '--eclip', action="store_true",
                        help='Adds a column with rbp\'s expression')

    args = parser.parse_args(input_)

    database_name = args.database_index_file.split('/')[-1].split('.')[0]

    db_bins = 2
    expression_file = args.expression_file
    expression_bins = args.expression_bins
    database_names_file = args.database_names_file
    database_index_file = args.database_index_file
    abundance_bins = args.sum_bins
    draw_bins = args.draw_bins
    child_unique_genes = 0.2
    parent_unique_genes = 0.4
    min_pathway_length = 6
    filter_redundant = args.filter_redundant
    max_draw_output = args.max_draw_output
    first_col_is_genes = args.first_col_is_genes
    sep = args.separator_in_expression_file
    input_format = args.expression_file_format
    output_format = args.database_format
    expression_columns = args.column_with_stability
    eclip = args.eclip
    tmp = 'tmp_ipage'

    if not os.path.isdir(tmp):
        os.mkdir(tmp)

    if args.preprocess:
        body.preprocess_db(database_names_file, first_col_is_genes, database_index_file, filter_redundant,
                           min_pathway_length, child_unique_genes, parent_unique_genes, tmp)
    if expression_file and database_index_file:

        if args.output_name:
            output_dir = args.output_name
            output_name = args.output_name.split('/')[-1]
        else:
            output_dir = 'output_ipage/'
            output_name = args.expression_file.split('/')[-1].split('.')[0]

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        if args.column_with_stability == 'all':
            expression_columns = range(1, pd.read_csv(expression_file, sep='\t').shape[1])
        else:
            expression_columns = (int(el) - 1 for el in expression_columns.split(','))

        for expression_column in expression_columns:

            output_name += '.' + str(expression_column) if len(list(expression_columns)) > 1 else ''
            output_name = '/'.join([output_dir, output_name])

            expression_profile, db_names, db_profiles, db_annotations, abundance_profile, genes = body.process_input(
                expression_file,
                database_name,
                input_format,
                output_format,
                expression_bins,
                abundance_bins,
                sep,
                expression_column,
                tmp)
            cmis = body.count_cmi_for_profiles(expression_profile, db_profiles, abundance_profile, expression_bins,
                                               db_bins, abundance_bins)
            accepted_db_profiles, z_scores = body.statistical_testing(cmis, expression_profile, db_profiles,
                                                                      abundance_profile,
                                                                      expression_bins, db_bins, abundance_bins)
            if eclip:
                rbp_expression = body.get_rbp_expression(genes, output_format, expression_profile, accepted_db_profiles,
                                                         db_annotations)
            else:
                rbp_expression = None

            body.visualize_output(accepted_db_profiles, db_profiles, db_annotations, cmis, draw_bins, max_draw_output,
                                  output_name, rbp_expression)

            with open(output_name + '.out', 'w+') as f:
                for i in range(len(accepted_db_profiles)):
                    if accepted_db_profiles[i]:
                        f.write(db_names[i] + '\t' + str(cmis[i]) + '\t' + str(z_scores[i]) + '\n')
