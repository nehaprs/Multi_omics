#==================
#rechecking RNAseq
#==================
'
rnaseq in multiomics: ADAM9 FDR suspiciosly low. Redoing the part to see whats up
'
setwd("~/BINF/o2pls/for manuscript/revision")

library(dplyr)
library(edgeR)
library(readxl)

mrna <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
counts = mrna[c(2:13)]

group = factor(c(rep("Control", 6), rep("Treatment", 6)))

dge = DGEList(counts = counts, group = group)
dge = calcNormFactors(dge)
dge = estimateDisp(dge)
fit = glmFit(dge)
lrt = glmLRT(fit, coef=2)
results = lrt$table
results$PValue <- lrt$table$PValue     
results$FDR <- p.adjust(results$PValue, method = "BH") 
mrna = cbind(mrna, results)
writexl::write_xlsx(mrna,"mrna_rechecked.xlsx")

#######################
#try glmQLFit
########################


mrna <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
counts = mrna[c(2:13)]

group = factor(c(rep("Control", 6), rep("Treatment", 6)))

dge = DGEList(counts = counts, group = group)
dge = calcNormFactors(dge)
dge = estimateDisp(dge)
#fit = glmFit(dge)
#lrt = glmLRT(fit, coef=2)
fit = glmQLFit(dge)
lrt = glmQLFTest(fit, coef=2)
results = lrt$table
results$PValue <- lrt$table$PValue     
results$FDR <- p.adjust(results$PValue, method = "BH") 
mrna = cbind(mrna, results)
writexl::write_xlsx(mrna,"mrna_glmqlfit.xlsx")

