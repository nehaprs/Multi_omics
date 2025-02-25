
library(dplyr)
library(readxl)
library(writexl)

#list all substrates with FC and FDR

setwd("~/BINF/o2pls/compare datasets/substrates with FDR and FC")

#list all the references
rna <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx")
cp <- read_excel("~/BINF/BINF_old/cellular proteomics/results_volcano_cleaned.xlsx")
sec <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/secretomics/perseus output/secretome_results.xlsx")

cp$pval = 10^(-cp$`-Log(P-value)`)
cp$FDR = p.adjust(cp$pval, method = "fdr")

#select ref columns
rna2 = rna[,c(1,14,17)]
cp2 = cp[,c(1,3,5,16)]
sec2 = sec[,c(1,3,2,17)]
colnames(rna2)[2] = "log2FC"
#colnames(cp2) = c("ID", "GeneName", "log2FC", "FDR")

#label datatype to FC and FDR columns
colnames(rna2) = c("GeneName", "log2FC.tr", "FDR.tr")
colnames(cp2) = c("ID", "GeneName", "log2FC.cp", "FDR.cp")
colnames(sec2) = c("ID", "GeneName", "log2FC.sec", "FDR.sec")

rm(rna, cp, sec)


#load substrate datasets

subs.tr_cp <- read_excel("~/BINF/o2pls/for manuscript/tr-cp substrates.xlsx", 
                             sheet = "all 188 membrane proteins")

subs.tr_sec <- read_excel("~/BINF/o2pls/tr-sec/output/analysis after integration/surface annotation/uniprot_suf_secr.xlsx", 
                                   sheet = "Sheet2")

subs.cp_sec <- read_excel("~/BINF/o2pls/cp-sec/output/analysis after integration/surface annotation/uniprot_suf_secr.xlsx")

#clean substrate datasets
subs.tr_cp2 = subs.tr_cp[,-c(2,3,5)]
colnames(subs.tr_cp2) = c("ID", "GeneName", "Subcellular location [CC]")

subs.tr_sec2 = subs.tr_sec[,c(1,4,5)]
colnames(subs.tr_sec2) = c("ID", "GeneName", "Subcellular location [CC]")

subs.cp_sec2 = subs.cp_sec[,c(1,3,5)]
colnames(subs.cp_sec2) = c("ID", "GeneName", "Subcellular location [CC]")

#add columns FDR, FC

FCFDR.subs.tr_cp = inner_join(subs.tr_cp2, rna2, by = "GeneName")
FCFDR.subs.tr_cp = inner_join(FCFDR.subs.tr_cp, cp2, by = c("ID", "GeneName"))
write_xlsx(FCFDR.subs.tr_cp,"FCFDR.subs.tr_cp.xlsx")

FCFDR.subs.tr_sec = inner_join(subs.tr_sec2, rna2, by = "GeneName")
FCFDR.subs.tr_sec = inner_join(FCFDR.subs.tr_sec, sec2, by = c( "GeneName"))
FCFDR.subs.tr_sec = FCFDR.subs.tr_sec[, - 1]
colnames(FCFDR.subs.tr_sec)[5] = "ID"
write_xlsx(FCFDR.subs.tr_sec, "FCFDR.subs.tr_sec.xlsx")

FCFDR.subs.cp_sec = inner_join(subs.cp_sec2, cp2, by = c("ID", "GeneName"))
FCFDR.subs.cp_sec = inner_join(FCFDR.subs.cp_sec, sec2, by = c( "ID","GeneName"))
write_xlsx(FCFDR.subs.cp_sec, "FCFDR.subs.cp_sec.xlsx")


#FDRFC of common proteins
setwd("~/BINF/o2pls/compare datasets/common proteins")

'''
comm_all3 <- read_excel("~/BINF/o2pls/compare datasets/comm_all3.xlsx")
FCFDR.comm_all3 = inner_join(comm_all3, rna2, by = "GeneName")
FCFDR.comm_all3 = inner_join(FCFDR.comm_all3, cp2, by = c( "ID","GeneName"))
FCFDR.comm_all3 = inner_join(FCFDR.comm_all3, sec2, by = c( "ID","GeneName"))
write_xlsx(FCFDR.comm_all3, "FCFDR.comm_all3.xlsx")
'''

xlsx_files = list.files(pattern = "xlsx$")

for(file in xlsx_files){
  df = read_excel(file)
  
  result <- df %>%
    inner_join(rna2, by = "GeneName") %>%
    inner_join(cp2, by = c("ID", "GeneName")) %>%
    inner_join(sec2, by = c("ID", "GeneName"))
  
  output_filename = paste0("FCFDR.", tools::file_path_sans_ext(basename(file)), ".xlsx")
  write_xlsx(result,output_filename )
  
}