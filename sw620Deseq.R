library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(rlang)
library(readxl)
library(apeglm)

#DESeq2 analysis of SW620 cells

##contents

#data, metadata upload, filter, transform
#sample heatmap
#gene heatmap
#PCA


exp_counts <- read.table(file="clipboard", sep = "\t",
                         header = T, row.names=1)
#file = "clipboard" after copying the entire file from excel sheet.
sample_data <- read.table(file="clipboard", sep = "\t", header = T,
                          row.names=1)
#here the metadata was made manually looking at the data
#Creating a DESeqDataSet from the Input Data

data_deseq <- DESeqDataSetFromMatrix(countData = exp_counts,
                                     colData = sample_data, design = ~ 1)
#while making metadata copy paste row names from col names of count data, 
#so that there's no small differences in clm names.

nrow(data_deseq) # = 58830
data_deseq = data_deseq[rowSums(counts(data_deseq))>1]
nrow(data_deseq) # = 33379
head(counts(data_deseq))

#Transforming the dataset to take care of different abundance of genes

rld = rlog(data_deseq,blind = FALSE)

#sample heatmap
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
rownames(sampleDistMatrix) <- paste( rld$cell_type, rld$Sample,
                                     sep="-" )

colnames(sampleDistMatrix) <- paste( rld$cell_type,rld$Sample,
                                     sep="-" )
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists)



#PCA



plotPCA(rld,intgroup = c("Sample"))
colnames(sample_data)

#gene heatmap

geneVars = rowVars(assay(rld))
geneVarsOrdered = order(geneVars, decreasing = TRUE)
topVarGenes = head(geneVarsOrdered,200)

mat = assay(rld)[topVarGenes,]
mat = mat - rowMeans(mat)
head(mat,10)
df <- as.data.frame(colData(rld)[,c("Sample")])
clear_col_names <- paste( rld$cell_type, rld$Sample, sep="."
)
topGenesHeatmap <- pheatmap(mat,fontsize = 1)
head(assay(rld)[topVarGenes,])


#design cell_type

data_deseq_group <- DESeqDataSetFromMatrix(countData = exp_counts,
                                           colData = sample_data, design = ~ Cell_type)
nrow(data_deseq_group)
data_deseq_group = data_deseq_group[rowSums(counts(data_deseq_group))>1]
#filtering the dataset for genes with 0 counts
nrow(data_deseq_group)

#from video
dds = DESeq(data_deseq_group)
res = results(dds)
BiocManager::install("apeglm")



reslfc = lfcShrink(dds, coef = 2)

(resOrdered <- res[order(res$padj),])
summary(res)
#results: out of 33379 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 2949, 8.8%
#LFC < 0 (down)     : 2680, 8%
#outliers [1]       : 6, 0.018%
#low counts [2]     : 14883, 45%
#(mean count < 7)

flt_vs_gc = as.data.frame(res$log2FoldChange)
head(flt_vs_gc)

plotMA(reslfc, ylim = c(-1,1))

write.csv(as.data.frame(resOrdered), file = "deseq_de_paired_from_video.csv")

#data transformation

dds <- estimateSizeFactors(dds)
#se == std error
#se <- SummarizedExperiment(log2(counts(dds, normalize = TRUE)+ 1),colData = colData(dds))
se2 = rlog(dds, blind = FALSE)
se = se2
plotPCA(DESeqTransform(se),intgroup = "Sample")

#sample heatmap
sampleDists3 <- dist(t(assay(se)))
sampleDists3
sampleDistMatrix3 <- as.matrix(sampleDists3)
head(sampleDistMatrix3)
rownames(sampleDistMatrix3) <- paste( se$cell_type, se$Sample,
                                      sep="-" )

colnames(sampleDistMatrix3) <- paste( se$cell_type,se$Sample,
                                      sep="-" )
pheatmap(sampleDistMatrix3, clustering_distance_rows=sampleDists3,
         clustering_distance_cols=sampleDists3)




library(AnnotationDbi)

library(org.Hs.eg.db)

foldchanges <- as.data.frame(res$log2FoldChange, row.names = row.names(res))

head(foldchanges)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "SYMBOL",
                     keytype = "GENENAME",
                     multiVals = "first")


res$entrez <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = "ENTREZID",
                     keytype = "GENENAME",
                     multiVals = "first")

head(res)


#gene heatmap
geneVars3 = rowVars(assay(se))
head(assay(se))
?rowVars
geneVarsOrdered3 = order(geneVars3, decreasing = TRUE)
topVarGenes3 = head(geneVarsOrdered3,50)

mat3 = assay(se)[topVarGenes3,]
mat3 = mat3 - rowMeans(mat3)
head(mat3,10)
df3 <- as.data.frame(colData(se)[,c("Sample")])
clear_col_names3 <- paste( se$cell_type, se$Sample, sep="."
)
topGenesHeatmap3 <- pheatmap(mat3,fontsize = 4)














#plotPCA(DESeqTransform(se),intgroup = "Sample")

#Transforming the dataset to take care of different abundance of genes
#uses regularized log transformation
rld2 = rlog(data_deseq_group,blind = FALSE)

#pca
plotPCA(rld2,intgroup = c("Sample"))
colnames(sample_data)

#sample heatmap
sampleDists2 <- dist(t(assay(rld2)))
sampleDists2
sampleDistMatrix2 <- as.matrix(sampleDists2)
head(sampleDistMatrix2)
rownames(sampleDistMatrix2) <- paste( rld2$cell_type, rld2$Sample,
                                      sep="-" )

colnames(sampleDistMatrix2) <- paste( rld2$cell_type,rld2$Sample,
                                      sep="-" )
pheatmap(sampleDistMatrix2, clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2)

#gene heatmap
geneVars2 = rowVars(assay(rld2))
geneVarsOrdered2 = order(geneVars2, decreasing = TRUE)
topVarGenes2 = head(geneVarsOrdered2,50)

mat2 = assay(rld2)[topVarGenes2,]
mat2 = mat2 - rowMeans(mat2)
head(mat2,10)
df2 <- as.data.frame(colData(rld2)[,c("Sample")])
clear_col_names2 <- paste( rld2$cell_type, rld2$Sample, sep="."
)
topGenesHeatmap2 <- pheatmap(mat2,fontsize = 4)
