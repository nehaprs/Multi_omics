
library(dplyr)
library(readxl)
library(writexl)

setwd("~/BINF/o2pls/compare datasets")
#compare substrates from tr-cp and tr-sec

subs_tr_sec <- read_excel("uniprot_suf_secr.xlsx", 
                             sheet = "Sheet2")

subs_tr_cp <- read_excel("~/BINF/o2pls/for manuscript/tr-cp substrates.xlsx", 
                           sheet = "all 188 membrane proteins", 
                           col_types = c("text", "text", "skip", 
                                             +         "text", "text", "text"))

comm_trcp_trsec = inner_join(subs_tr_cp, subs_tr_sec, by = "Entry")

write_xlsx(comm_trcp_trsec, "comm_trcp_trsec.xlsx")

#compare substrates cp-sec with from tr-cp and tr-sec

subs_cp_sec <- read_excel("~/BINF/o2pls/cp-sec/output/analysis after integration/surface annotation/uniprot_suf_secr.xlsx")
comm_trcp_cpsec = inner_join(subs_tr_cp, subs_cp_sec, by = "Entry")
write_xlsx(comm_trcp_cpsec, "comm_trcp_cpsec.xlsx")

comm_trsec_cpsec = inner_join(subs_tr_sec, subs_cp_sec, by = "Entry")
write_xlsx(comm_trsec_cpsec, "comm_trsec_cpsec.xlsx")

comm_all3 = inner_join(subs_tr_cp, comm_trsec_cpsec,by = "Entry" )
write_xlsx(comm_all3, "comm_all3.xlsx")
