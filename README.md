# ipage
Conditional mutual information based version of ipage

To test the programme run the following command:

python3.6 scripts/ipage.py -e_file test/bladder.exp -db_names test/human_ensembl_names.txt -db_index test/human_ensembl_index.txt -db_format ensg -e_format refseq -e_bins 10 -d_bins 15 -f -max_draw 50

-e_file expression file

-db_names database names file

-db_index database index file

-e_bins expression_bins

-d_bins bins in heatmap

-db_format database format (ensg, enst, refseq,entrez,gene_symbol)

-e_format expression file format

-max_draw maximum number of results presented in graphical output

-f filter_redundant
