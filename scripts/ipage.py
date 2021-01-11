import sys
import argparse
import body
import os
import pandas as pd
import preprocess
import filter_db


def process_annotation(annotation_index_file=None, istable=None, sep='\t',
                  annotation_names_file=None, first_col_is_genes=True,
                  filter_redundant=False, child_unique_genes=0.2, min_pathway_length=20,
                  annotation_dir='annotation_dir', annotation_name=None):

    if not os.path.isdir(annotation_dir):
        os.mkdir(annotation_dir)

    if istable:
        annotation_table = annotation_index_file
        db_names, db_profiles, db_annotations, db_genes = \
            preprocess.get_profiles_from_table(annotation_table, sep=sep, first_col_is_genes=first_col_is_genes)

    else:
        db_names, db_profiles, db_annotations, db_genes = \
            preprocess.get_profiles(annotation_index_file, first_col_is_genes, annotation_names_file)

    if filter_redundant:
        db_names, db_annotations, db_profiles = filter_db.non_redundancy_sort_pre(db_names, db_annotations, db_profiles,
                                                                                  min_pathway_length,
                                                                                  child_unique_genes)
    if not annotation_name:
        annotation_name = annotation_index_file.split('/')[-1].split('.')[0]
    preprocess.dump_database(db_names, db_annotations, db_genes, db_profiles, annotation_name, annotation_dir)


def read_expression_file(expression_file, sep='\t', id_column=0, de_column=1):
    expression_column = de_column
    df = pd.read_csv(expression_file, sep=sep, skiprows=1, header=None)
    genes = df.iloc[:, id_column]
    expression_level = df.iloc[:, expression_column]
    return genes, expression_level


def ipage(de_genes, de_profile, annotation_name, annotation_dir='annotation_dir', output_name='stdout',
          de_ft=None, ann_ft=None, species='human', symmetric_expression=True,
          de_bins=10, a_bins=3, function='cmi', alpha=0.01, holm_bonferroni=False,
          max_heatmap_rows=20, export_heatmap=False, heatmap_bins=15,
          regulator=False, cmap_main='RdBu_r', cmap_reg='RdBu_r'):

    db_bins = 2
    if output_name != 'stdout' and output_name != 'shut':
        if len(output_name.split('/')) != 1:
            output_dir = os.path.join(*output_name.split('/')[:-1])
        else:
            output_dir = output_name
            output_name += '/' + output_name

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

    de_profile_discr, ann_names, ann_profiles, ann_annotations, abundance_profile, genes = body.process_input(
        de_profile, de_genes, annotation_name, de_ft, ann_ft, de_bins, a_bins,
        species, annotation_dir, symmetric_expression)

    cmis = body.count_cmi_for_profiles(de_profile_discr, ann_profiles, abundance_profile, de_bins,
                                       db_bins, a_bins, function)
    accepted_ann_profiles, z_scores = body.statistical_testing(cmis, de_profile_discr, ann_profiles,
                                                               abundance_profile, de_bins, db_bins, a_bins,
                                                               function, alpha, holm_bonferroni)
    if regulator:
        regulator_expression = body.get_rbp_expression(genes, ann_ft, de_profile,
                                                       accepted_ann_profiles, ann_names, ann_annotations, species,
                                                       annotation_dir)
    else:
        regulator_expression = None

    output = body.produce_output(accepted_ann_profiles, ann_profiles, ann_names, ann_annotations, cmis, z_scores,
                                 heatmap_bins, max_heatmap_rows, output_name, cmap_main, cmap_reg, regulator_expression,
                                 export_heatmap)
    if output_name == 'stdout' or output_name == 'shut':
        return output


if __name__ == '__main__':
    input_ = sys.argv[1:]

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-annotation_dir', '--annotation_dir', type=str, default='./',
                        help='The folder where annotations are stored')
                            
    parser.add_argument('-preprocess', '--preprocess', action="store_true",
                        help='Preprocess annotation')
    parser.add_argument('-annotation', '--annotation', type=str,
                        help='''If preprocess is called then annotation should specify the name of the file
                        with gene annotations, otherwise simply provide the name of the preprocessed annotation''')
    parser.add_argument('-aliases', '--annotation_aliases', type=str,
                        help='Aliases for gene set names from gene annotation, provide only if preprocess is called')
    parser.add_argument('-istable', '--annotation_is_table', action="store_true",
                        help='Specifies if annotation is in table format')
    parser.add_argument('-f', '--filter_redundant', action="store_true",
                        help='Filter redundant gene sets')
    parser.add_argument('-g', '--first_element_is_gene', action="store_true",
                        help='The first element in each row is a gene name')
    parser.add_argument('-length', '--min_pathway_length', type=int, default=20,
                        help='The minimal pathway length')
    parser.add_argument('-fs', '--filter_stringency', type=float, default=0.2,
                        help='Flitering stringency')
    parser.add_argument('-asep', '--annotation_separator', type=str, default='\t',
                        help='The separator in annotation file')
    #
    #
    parser.add_argument('-e', '--expression_file', type=str,
                        help='The name of the expression file')
    parser.add_argument('-e_sep', '--separator_in_expression_file', type=str, default='\t',
                        help='The separator in expression file')
    parser.add_argument('-e_col', '--expression_column', type=str, default='1',
                        help='''Number of the column in the expression file with values, can be a number, 
                        a comma delimited list of numbers. If \'all\' specified then all columns are used''')

    parser.add_argument('-ann_ft', '--annotation_format', type=str, default='ensg',
                        help='The annotation gene format (ensg, enst, gs, refseq, entrez, ext)')
    parser.add_argument('-e_ft', '--expression_format', type=str, default='ensg',
                        help='The expression file gene format (ensg, enst, gs, refseq, entrez, ext)')    
    parser.add_argument('-sp', '--species', type=str, default='human',
                        help='Adds a column with rbp\'s expression (human, mouse)')
    
    parser.add_argument('-e_bins', '--expression_bins', type=int, default=10,
                        help='The number of the expression bins')
    parser.add_argument('-asymmetric', '--symmetric_expression', action="store_false",
                        help='Provide if the expression is not zero centered')      
    parser.add_argument('-a_bins', '--abundance_bins', type=int, default=3,
                        help='The number of bins in the sum profile')
    parser.add_argument('-heatmap_bins', '--heatmap_bins', type=int, default=15,
                        help='The number of the bins in the heatmap')
    
    parser.add_argument('-max_draw', '--max_draw_output', type=int, default=20,
                        help='The number of rows to be presented on the heatmap')
    parser.add_argument('-reg', '--regulator', action="store_true",
                        help='Adds a column with regulator\'s expression to the heatmap')

    parser.add_argument('-func', '--function', type=str, default='cmi',
                        help='Add information function to use (cmi, mi)')
    parser.add_argument('-alpha', '--alpha', type=float, default=0.01,
                        help='Alpha threshold for hypothesis testing')

    parser.add_argument('-o', '--output_name', type=str,
                        help='The name of the output file')

    args = parser.parse_args(input_)

    annotation_name = args.annotation.split('/')[-1].split('.')[0]
    if args.preprocess:
        process_annotation(args.annotation, args.annotation_is_table, args.annotation_separator,
                           args.annotation_aliases, args.first_element_is_gene,
                           args.filter_redundant, args.filter_stringency, args.min_pathway_length,
                           args.annotation_dir, annotation_name)

    if args.expression_file:

        if args.expression_column == 'all':
            expression_columns = list(range(1, pd.read_csv(args.expression_file, sep=sep).shape[1]))
        else:
            expression_columns = list((int(el) for el in args.expression_column.split(',')))

        if args.output_name:
            output_name_template = args.output_name
        else:
            output_name_template = os.path.join('output_ipage', args.expression_file.split('/')[-1].split('.')[0])

        for expression_column in expression_columns:
            output_name = output_name_template + ('.' + str(expression_column) if len(expression_columns) > 1 else '')

            de_genes, de_profile = read_expression_file(args.expression_file, sep=args.separator_in_expression_file,
                                                        id_column=0, de_column=expression_column)

            ipage(de_genes, de_profile, annotation_name, args.annotation_dir, output_name,
                  args.expression_format, args.annotation_format, args.species, args.symmetric_expression,
                  args.expression_bins, args.abundance_bins, args.function, args.alpha, False,
                  args.max_draw_output, False, args.heatmap_bins, args.regulator)

