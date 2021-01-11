library("tximport")
library("readr")
library("DESeq2")
args = commandArgs(trailingOnly=TRUE)
name=args[1]
files = unlist(strsplit(args[2], ','))
conditions = unlist(strsplit(args[3], ','))
detailed = as.logical(args[4])
tx2gene_file = args[5]
tx2gene <- read_csv(file.path(tx2gene_file))
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
colData=data.frame('SampleName'= files, condition=conditions)
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)
res<-results(dds)
if (detailed) {
  write.table(res, file=name, sep="\t",row.names = FALSE)
} else {
  df <- cbind(rownames(res), unlist(lapply(res$log2FoldChange,sign))*(1-res$pvalue))
  colnames(df)<-c('genes', 's*(1-p)')
  write.table(df, file=name, sep="\t",row.names = FALSE)
}
