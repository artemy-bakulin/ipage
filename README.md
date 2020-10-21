# iPAGE-2

iPAGE-2 is a tool for the pathway-level analysis of differential gene expression, which implements the ideas of
information theory.

Introduction to the pathway-level analysis requires the clarification of the thesaurus used in the field.

So, first, the term **pathway** – by this word we mean practically any set of genes: a metabolic chain,
genes regulated by a transcription factor, targets of an RNA binding protein. There are various sources of **pathways'
annotations**, files that contain the description of which genes to which pathways belong to. We provide a few preprocessed
pathways' annotations with our package but one should feel free to use annotation files of their own.

The predictions of iPAGE are based on the **differential expression**: the assessment of the difference of genes' expression
distributions between two conditions. Differential expression should be precomputed before using iPAGE (e.g. by Deseq2).
Virtually any metric can be used with iPAGE as it implements non-parametric testing therefore its output should be robust
to mathematical transformation of input data but we recommend that differential expression is computed as a p-value of the significance of
transcript's abundance change subtracted from one taken as a positive value if the abundance increases and negative if
vice versa ('sign(fold_change) * (1-p)'). Further we refer to the vector of differential expression values as the
**differential expression profile**.

## Download and installation

Currently iPAGE is accessible from this github page. The package can be cloned to local environment using a link:
"https://github.com/artemy-bakulin/ipage.git".

Installation of dependency packages may also be required "pybiomart", "statsmodels", "scipy", "argparse", "numpy", "pandas"
(use pip or conda for this purpose).

## The typical session

So having cloned iPAGE-2, installed all dependencies and computed differential expression, one could immediately go ahead
exploring pathways' regulation.

iPAGE-2 can be used both as a python package and as a command line tool. Here we explore the typical python session (also consult 'test/test_session.ipynb').

Start with inserting iPAGE to your path.

```
import sys
sys.path.insert(1,'direction/to/ipage/scripts/')
import ipage
```

Then download differential expression data in any convenient manner, for example:

```
 de_genes, de_profile = ipage.read_expression_file('expression_file.csv', sep=',')
```
'de_genes', 'de_profile' are two lists, numpy arrays or pandas series.

And finally iPAGE command can be executed.

```
annotation_name = 'msig_db'
de_gene_format = 'ensg'
annotation_gene_format = 'gs'
annotation_dir = 'direction/to/tmp_ipage'
ipage.ipage(de_genes, de_profile, annotation_name, output_name='test', de_ft=de_gene_format, ann_ft=annotation_gene_format, annotation_dir=annotation_dir)
```
The program creates two output files: a tab-delimited table with all significant hits ('test/test.csv') and a heatmap
with top hits ('test/test.jpg').

<img src="test/test.jpg" width="400">

## A comment on the analysis of RBPs

iPAGE can also be used for the analysis of RBPs' deregulation. Though in this case not differential expression but differential stability should be used.
We advise the use of REMBRANDTS (https://github.com/csglab/REMBRANDTS) package to calculate genes' stability.
Then for each gene t test (z test) should be used to calculate the p values. Differential stability should be computed as
(1 - p.value) * sign_of_stability_change.
This package provides RBP regulons annotation: 'hybrid_CLIP'.



## Functions' specification

### ipage.ipage

```
ipage.ipage(de_genes, de_profile, annotation_name, annotation_dir='annotation_dir', output_name='stdout',
          de_ft=None, ann_ft=None, species='human', symmetric_expression=True,
          de_bins=10, a_bins=3, function='cmi', alpha=0.01, holm_bonferroni=False,
          max_draw_output=20, export_heatmap=False, heatmap_bins=15,
          regulator=False, cmap_main='RdBu_r', cmap_reg='YlOrBr')
```

This function takes gene expression data and discovers significantly deregulated pathways.
It produces a table and a heatmap.

Parameters:

	de_genes : list, pd.Series, np.array
				List of genes in expression profile.

	de_profile : list, pd.Series, np.array
				List of DE values

	annotation_name : str
				Name of the pathways' database
				(Should be preprocessed before launch with ipage.preprocess_db)

	output_name : str
				Name of the output.
				If provided "stdout", the output is directed to stdout

	de_ft : str
				Gene accession type in expression profile.
				Available accessions: ensg, enst, refseq, entrez, gs (gene symbol)

	ann_ft : str
				Gene accession type in annotation
				Available accessions: ensg, enst, refseq, entrez, gs (gene symbol)

	de_bins: int
				The number of bins by which expression is discretized

	a_bins: int
				The number of bins by which gene abundance is discretized

	species: str
				Species ('human'/'mouse') which genes are analyzed.
				Should be specified if annotation's and differential expression's file use different accessions.

	heatmap_bins: int
				Number of columns on a heatmap.

	max_heatmap_rows: int
				Number of rows on a heatmap.

	regulator: bool
				Can be added if regulons are named with gene symbols after another gene.
				Adds an additional column with expression of the regulator to heatmap.

	annotation_dir: str
				The folder where gene annotations and other intersession files are stored.

	function: str
				The function which iPAGE uses ('mi' for mutual information, 'cmi for conditional mutual information').
				Conditional mutual information is recommended for biased set annotations.

	alpha: float
				The threshold for the p-value.

	holm_bonferroni: bool
				Specifies if Holm-Bonferroni correction should be used or not.

	cmap_main: str
				Colormap from matplotlib collection to be used on heatmap.

	cmap_reg: str
				Colormap from matplotlib collection to be used on heatmap in the expression column.

	export_heatmap: bool
				If specified exports heatmap data in pickle format.


	symmetric_expression: bool
				If expression is restricted and symmetrical around zero.

### ipage.process_annotation

```
ipage.process_annotation(annotation_table=None, sep='\t',
                annotation_index_file=None, annotation_names_file=None, first_col_is_genes=True,
                filter_redundant=False, child_unique_genes=0.2, min_pathway_length=20,
                annotation_dir='annotation_dir', annotation_name=None)
```

This function processes annotation files to the compatible with iPAGE form.
It is common that annotation files are formatted either as a binary table (where rows are genes and columns are pathways,
at in each cell of this table there is either 1 or 0 if gene belongs to this pathway or not) or a file where on each line
a pathway and its corresponding genes are listed separated by delimeter ("index format"). "ipage.process_annotation" works with both formats.

Parameters:

    annotation_table: str
                Name of the annotation file in binary table format.
                Either annotation_table or annotation_index_file should be specified.

    sep: str
                Separator in the annotation_table file.

    annotation_index_file: str
                Name of the annotation file in the index format.
                Either annotation_table or annotation_index_file should be specified.

    annotation_names_file: str
                Sometimes there is a file with description of pathways names in the annotation_index_file, this parameter
                sets the name of this file. Should be specified only if annotation_index_file is specified but it also may be omitted.

    first_col_is_genes: bool
                annotation_index_file may be formatted so that the first element in each row is a gene or that it is a pathway name,
                the parameter specifies this format issue.

    filter_redundant: bool
                It is not unusual of pathways annotation to contain multiple almost identical entries,
                this parameter specifies whether they should be sorted.

    min_pathway_length: int
                The parameter specifies the minimal pathway length.

    annotation_dir: str
                Specifies the directory where the processed annotations are stored.

    annotation_name: str
                Specifies the name of the processed annotations.

### ipage.read_expression_file
```
ipage.read_expression_file(expression_file, sep='\t', id_column=0, de_column=1)
```
This function reads differential expression file.

Parameters:

    expression_file: str
               Specifies the name of the expression file.

    sep: str
                Specifies the file's separator.str

    id_column: int
                Specifies gene id column's number in file.

    de_column: int
                Specifies differential expression column's number in file.