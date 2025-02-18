###surface annotate pttpossig

library(dplyr)
library(OmicsPLS)
library(readr)
library(readxl)
library(writexl)
library(openxlsx)
library(readxl)
library(pracma)
library(tidyr)

pttpos_idmap <- read_excel("analysis after integration/pttpos_idmap.xlsx")
pttpos_idmap = as.data.frame(pttpos_idmap)

pttpos_idmap$`Subcellular location [CC]` <- gsub(
  "\\{.*?\\}"," ",pttpos_idmap$`Subcellular location [CC]`
)

pttpos_idmap$`Subcellular location [CC]` <- 
  gsub("Note.*?\\.", "", pttpos_idmap$`Subcellular location [CC]`)

pttpos_idmap$`Subcellular location [CC]` <- 
  gsub("SUBCELLULAR LOCATION:", "", pttpos_idmap$`Subcellular location [CC]`)

#make subset of 'membrane' and 'secreted' annotated
uniprot_suf_secr <- pttpos_idmap %>% 
  filter(grepl("membrane", pttpos_idmap$`Subcellular location [CC]`,
               ignore.case = TRUE) |
           grepl("secreted", pttpos_idmap$`Subcellular location [CC]`,
                 ignore.case = TRUE, ))

write_xlsx(uniprot_suf_secr,"surface annotation/uniprot_suf_secr.xlsx")

#FC cut off
setwd("~/BINF/o2pls/tr-sec/output/analysis after integration/fc1.5")
ttpos1.5 = ttpossigBH[abs(ttpossigBH$log2FC) > log2(1.5), ]
write_xlsx(ttpos1.5, "ttpos1.5.xlsx")

ttneg1.5 = ttnegsigBH[abs(ttnegsigBH$log2FC) > log2(1.5), ]
write_xlsx(ttneg1.5, "ttneg1.5.xlsx")

##comparing correlations between rwa data and joint components
#finding correlations based on raw data



pttpos1.5 = pttpossigBH[abs(pttpossigBH$log2FC) > log2(1.5), ]
write_xlsx(pttpos1.5, "pttpos1.5.xlsx")

pttneg1.5 = pttnegsigBH[abs(pttnegsigBH$log2FC) > log2(1.5), ]
write_xlsx(pttneg1.5, "pttneg1.5.xlsx")


