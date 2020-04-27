import sys
import argparse
import body
import os
import pandas as pd


def preprocess(database_index_file, database_names_file=None, first_col_is_genes=True, filter_redundant=False,
               tmp='tmp_ipage', min_pathway_length=6):
    child_unique_genes = 0.2
    if not os.path.isdir(tmp):
        os.mkdir(tmp)
    body.preprocess_db(database_names_file, first_col_is_genes, database_index_file, filter_redundant,
                       min_pathway_length, child_unique_genes, tmp)


def read_expression_file(expression_file, sep='\t', column=1):
    id_column = 0
    expression_column = column
    df = pd.read_csv(expression_file, sep=sep, skiprows=1, header=None)
    genes = df.iloc[:, id_column]
    expression_level = df.iloc[:, expression_column]
    return genes, expression_level


def ipage(genes, expression_level, database_name, output_name='stdout', e_ft=None, db_ft=None,
          e_bins=10, freq_bins=3, species='human', draw_bins=15, max_draw_output=50, regulator=False,
          tmp='tmp_ipage', function='cmi', p_value=0.01):

    db_bins = 2

    if output_name != 'stdout' and output_name != 'shut':
        if len(output_name.split('/')) != 1:
            output_dir = '/'.join(output_name.split('/')[:-1])
        else:
            output_dir = output_name
            output_name += '/' + output_name

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)


    expression_profile, db_names, db_profiles, db_annotations, abundance_profile, genes = body.process_input(
        expression_level, genes, database_name, e_ft, db_ft, e_bins, freq_bins,
        species, tmp)

    cmis = body.count_cmi_for_profiles(expression_profile, db_profiles, abundance_profile, e_bins,
                                       db_bins, freq_bins, function)
    accepted_db_profiles, z_scores = body.statistical_testing(cmis, expression_profile, db_profiles,
                                                              abundance_profile,
                                                              e_bins, db_bins, freq_bins, function, p_value)
    if regulator:
        regulator_expression = body.get_rbp_expression(genes, db_ft, expression_profile,
                                                       accepted_db_profiles, db_names, db_annotations, species,
                                                       tmp)
    else:
        regulator_expression = None

    output = body.produce_output(accepted_db_profiles, db_profiles, db_names, db_annotations, cmis, z_scores,
                        draw_bins, max_draw_output, output_name, regulator_expression)
    # if output_name == 'stdout' or output_name == 'shut':
    return output


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

    parser.add_argument('-db_ft', '--database_format', type=str, default='ensg',
                        help='The database format')
    parser.add_argument('-e_ft', '--expression_file_format', type=str, default='ensg',
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

    parser.add_argument('-sp', '--species', type=str, default='human',
                        help='Adds a column with rbp\'s expression')

    parser.add_argument('-reg', '--regulator', action="store_true",
                        help='Adds a column with regulator\'s expression')

    parser.add_argument('-func', '--function', type=str,
                        help='Add information function to use')

    args = parser.parse_args(input_)

    database_name = args.database_index_file.split('/')[-1].split('.')[0]

    expression_file = args.expression_file
    expression_bins = args.expression_bins
    database_names_file = args.database_names_file
    database_index_file = args.database_index_file
    abundance_bins = args.sum_bins
    draw_bins = args.draw_bins
    filter_redundant = args.filter_redundant
    max_draw_output = args.max_draw_output
    first_col_is_genes = args.first_col_is_genes
    sep = args.separator_in_expression_file
    input_format = args.expression_file_format
    output_format = args.database_format
    expression_columns = args.column_with_stability
    species = args.species
    regulator = args.regulator

    tmp = 'tmp_ipage'
    if not os.path.isdir(tmp):
        os.mkdir(tmp)

    if args.preprocess:
        preprocess(database_names_file, database_index_file, first_col_is_genes, filter_redundant, tmp)
    if expression_file and database_index_file:

        if args.column_with_stability == 'all':
            expression_columns = range(1, pd.read_csv(expression_file, sep=sep).shape[1])
        else:
            expression_columns = (int(el) for el in expression_columns.split(','))

        for expression_column in expression_columns:
            if args.output_name:
                output_name = args.output_name
            else:
                output_name = 'output_ipage/' + args.expression_file.split('/')[-1].split('.')[0]
            output_name += '.' + str(expression_column) if len(list(expression_columns)) > 1 else ''

            genes, expression_level = read_expression_file(expression_file, sep, expression_column)

            ipage(expression_level, genes, database_index_file, output_name, input_format, output_format, expression_bins,
                  abundance_bins, species, draw_bins, max_draw_output, regulator, tmp)

