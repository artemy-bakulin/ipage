import pandas as pd
import os


# You can find tx2gene here: https://zenodo.org/record/1324497#.Xef4NKdeOqA
def deseq(files_a, files_b, conditions=('A', 'B'), name='output.deseq.tsv', detailed=False, keep_file=False,
          tx2gene='gencode.v28.annotation.tx2gene.csv'):
    files = ','.join(files_a+files_b)
    conditions = ','.join([conditions[0]]*len(files_a)+[conditions[1]]*len(files_b))
    cmd = 'Rscript deseq.R %s %s %s %s %s' % (name, files, conditions, bool(detailed), tx2gene)
    os.system(cmd)
    df = pd.read_csv(name, sep='\t')
    if not keep_file:
        os.remove(name)
    return df

