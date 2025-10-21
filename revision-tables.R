################
#resolve discrepency in numbers for ProteinGroup.txt
##################
proteinGroups <- read_excel("~/BINF/o2pls/for manuscript/revision/round3/proteinGroups.xlsx")


proteinGroups$genes = sub(";.*","", proteinGroups$`Gene names`)
head(proteinGroups$genes)

cp <- read_excel("cp_all.xlsx")

df_subset <- df[df$genes %in% df2$`Gene Names`, ]

proteinGroups2 = proteinGroups[proteinGroups$genes %in% cp$`Gene Names`, ]

##################
proteinGroups <- read_excel("~/BINF/o2pls/for manuscript/revision/round3/proteinGroups.xlsx")


proteinGroups$ID = sub(";.*","", proteinGroups$`Majority protein IDs`)
head(proteinGroups$ID)

cp <- read_excel("cp_all.xlsx")

proteinGroups2 = proteinGroups[proteinGroups$ID %in% cp$ID, ]
writexl::write_xlsx(proteinGroups2,"proteinGroups_filtered.xlsx")
