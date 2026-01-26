#===================
#MCP revision round 2
#12.18.2025
#1.6.2026: gsea
#=====================

'
 comment 1: Some proteins listed in Supplemental Table S2 were identified based on a single unique peptide. 
 The corresponding proteins should be removed from Table S2.

'
install.packages("msigdbr")
library(dplyr)
library(tibble)
library(readxl)
library(writexl)
library(fgsea)
library(msigdbr)


s2 <- read_excel("~/BINF/o2pls/for manuscript/revision/second-revision/tableS2_v2.xlsx")

s2.2 = s2[ s2$...12 != 0, ]
write_xlsx(s2.2,"tableS2_v3.xlsx")

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
df2 = df[rowSums(df[,num_cols] >5) > 0,]
write_xlsx(df2, "controls_thresh3.xlsx")
#thresh: v1: 0, v2 =10, v3 = 5

##############
#gsea for second revision: incomplete
##############
setwd("~/BINF/o2pls/for manuscript/revision/second-revision")

bg_genes <- read_excel("controls/controls_thresh2.xlsx", 
                          sheet = "Sheet2")
bg_genes <- as.character(bg_genes[[1]])

de_table <- read_excel("degs.xlsx")
de_table$gene = de_table$GeneName
str(de_table$PValue)
de_table$PValue <- as.numeric(de_table$PValue)
de_table$stat <- sign(de_table$logFC) * (-log10(de_table$PValue))

dt <- de_table %>%
  filter(!is.na(gene)) %>%
  distinct(gene, .keep_all = TRUE)

dt <- dt %>% filter(gene %in% bg_genes)

ranks <- dt$stat
names(ranks) <- dt$gene

m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)

pathways <- split(m_df$gene_symbol, m_df$gs_name)

universe <- names(ranks)
pathways <- lapply(pathways, intersect, universe)

fg <- fgseaMultilevel(
  pathways = pathways,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500
)

fg <- fg %>% arrange(padj)
head(fg, 10)
#msigdb: ony 2 pathways, 

library(clusterProfiler)
library(org.Hs.eg.db)

m <- bitr(names(ranks),
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = org.Hs.eg.db)

# keep only mapped genes, ensure 1-to-1
ranks2 <- ranks[m$SYMBOL]
names(ranks2) <- m$ENTREZID

# drop duplicates (keep the strongest stat per Entrez)
ranks2 <- tapply(ranks2, names(ranks2), function(x) x[which.max(abs(x))])
ranks2 <- sort(ranks2, decreasing = TRUE)

class(ranks2)
is.numeric(ranks2)
is.vector(ranks2)


ranks2 <- as.numeric(ranks2)
names(ranks2) <- names(tapply(ranks2, names(ranks2), identity)) 

