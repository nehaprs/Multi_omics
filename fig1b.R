#all common proteins
alltargets_bh1xy <- read_excel("~/BINF/o2pls/tr-cp redo for foxo3/alltargets_bh1xy.xlsx")
#get logFC
cp_All_results <- read_excel("~/BINF/o2pls/tr-cp redo for foxo3/cp_All_results.xlsx")
cpfdr = cp_All_results[cp_All_results$FDR < 0.05,]
cpfdr = cpfdr[abs(cpfdr$Log2FC) > log2(1.2), ]
rm(cp_All_results)
cpfdr = cpfdr[,c(1,3,15,16)]
alltargetsDEP = inner_join(alltargets_bh1xy, cpfdr, by = "ID")
alltargetsDEP =alltargetsDEP[,-10]
alltargetsDEP = alltargetsDEP %>% rename(FDR.cp = FDR)
alltargetsDEP = rename(alltargetsDEP, log2FC.cp = Log2FC)

#mrnalog2fc

mRNA_All_Results <- read_excel("~/BINF/BINF_old/multi-omics/ADAM9 KD HCT116 proteomics/Protein and mRNA/DEGs and DEPs/DEGs/mRNA All Results.xlsx", 
                               col_types = c("text", "skip", "skip", 
                                             "skip", "skip", "skip", "skip", "skip", 
                                             "skip", "skip", "skip", "skip", "skip", 
                                             "numeric", "numeric", "numeric", 
                                             "numeric"))

alltargetsDEP2 = inner_join(alltargetsDEP, mRNA_All_Results, by = "GeneName")
rm(alltargetsDEP)
alltargetsDEP = alltargetsDEP2
rm(alltargetsDEP2)
alltargetsDEP = alltargetsDEP[,-c(13,14)]
alltargetsDEP = alltargetsDEP %>% rename(FDR.tr = FDR, log2FC.tr = logFC)


####

#plot logFC scatterplots

#id transcriptional and ptt

tt = alltargetsDEP %>% filter(Corr_xyJnt > 0, bh1.xy < 0.05)
#323 tts

ptt = alltargetsDEP %>% filter(Corr_xyJnt < 0 | (Corr_xyJnt > 0 & bh1.xy > 0.05))
#1963 ptt
genes_to_label <- c("ADAM9","ALCAM","IGSF8","IL10RB","FOXO3")

################################
#make plot
ggplot(alltargetsDEP, aes(x = log2FC.tr, y = log2FC.cp) )+
  geom_hline(yintercept = 0, color = "black")+
  geom_vline(xintercept = 0, color = "black")+
  # base layer: all points in grey, no legend
  geom_point(color = "grey", size = 2) +
  
  # colored points for transcriptional targets and post-transcriptional targets
  
  geom_point(data = tt, aes(color = "Candidate Transcriptional Targets"), size = 2)+
  geom_point(data = ptt, aes(color = "Candidate Post-Transcriptional Targets"), size = 2)+
  
  # manual scale: map the two labels to colors
  scale_color_manual(name   = "Target type",
                     values = c(
                       "Candidate Transcriptional Targets"     = "blue",
                       "Candidate Post-Transcriptional Targets" = "grey"
                     )) +
  
  # clean theme & axis labels
  #theme_minimal(base_size = 14) +
  labs(x = expression(log[2]*"FC mRNA"),
       y = expression(log[2]*"FC protein"))+
  theme(
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25)
  )



'
+
  geom_text_repel(
    data = alltargetsDEP %>% filter(GeneName %in% genes_to_label),
    aes(x = log2FC.tr, y = log2FC.cp, label = GeneName),
    size = 4,
    min.segment.length = 0,
    segment.color = "black",
    #fontface = "bold",
    box.padding = 1,
    point.padding = 0.5,
    max.overlaps = Inf,
    inherit.aes = FALSE,
    na.rm = TRUE
  )

'
