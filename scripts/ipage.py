import sys
import argparse
import body

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
    parser.add_argument('-e_col', '--column_with_stability', type=int, default=2,
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
    expression_column = args.column_with_stability

    if args.preprocess:
        body.preprocess_db(database_names_file, first_col_is_genes, database_index_file, filter_redundant,
                           min_pathway_length, child_unique_genes, parent_unique_genes)
    if expression_file and database_index_file:
        output_name = args.output_name if args.output_name else args.expression_file.split('/')[-1].split('.')[0]
        expression_profile, db_names, db_profiles, db_annotations, abundance_profile = body.process_input(
            expression_file,
            database_name,
            input_format,
            output_format,
            expression_bins,
            abundance_bins,
            sep,
            expression_column)
        cmis = body.count_cmi_for_profiles(expression_profile, db_profiles, abundance_profile, expression_bins, db_bins,
                                           abundance_bins)
        accepted_db_profiles, z_scores = body.statistical_testing(cmis, expression_profile, db_profiles,
                                                                  abundance_profile,
                                                                  expression_bins, db_bins, abundance_bins)
        body.visualize_output(accepted_db_profiles, db_profiles, db_annotations, cmis, draw_bins, max_draw_output,
                              output_name)

        with open(output_name + '.out', 'w+') as f:
            for i in range(len(accepted_db_profiles)):
                if accepted_db_profiles[i]:
                    f.write(db_names[i] + '\t' + str(cmis[i]) + '\t' + str(z_scores[i]) + '\n')
