#===================
#MCP revision round 2
#12.18.2025
#=====================

'
 comment 1: Some proteins listed in Supplemental Table S2 were identified based on a single unique peptide. 
 The corresponding proteins should be removed from Table S2.

'

library(dplyr)
library(readxl)
library(writexl)

s2 <- read_excel("~/BINF/o2pls/for manuscript/revision/second-revision/Table S2.xlsx")

s2.2 = s2[s2$...12 != 1,]
write_xlsx(s2.2,"tableS2_v2.xlsx")

'
 comment 2: I would suggest using as a background the set of genes that are expressed in HCT116
 WT cells under their experimental conditions.

'

#create this background set
#what is 'expressed'?
# in the list of 12k genes, at least one control > 0

#df from binfold/rnaseq/counting
#counts ≥ 1 in ≥ 1  (how many columns have a value ≥ 10)
df <- read_excel("~/BINF/o2pls/for manuscript/revision/second-revision/controlsv1.xlsx")


df = as.data.frame(df)


num_cols = sapply(df, is.numeric )
df1 = df[rowSums(df[,num_cols] >10) > 0,]
write_xlsx(df1, "controls_thresh2.xlsx")
#thresh: v1: 0, v2 =10

