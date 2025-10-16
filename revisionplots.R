library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)
library(ggbreak)
library(tidyr)
library(ggsignif)
setwd("~/BINF/o2pls/for manuscript/revision")

#####################
#scatterplot logfc vc logfc
#last edit: 10.6.2025
#####################

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
                       "Candidate Post-Transcriptional Targets" = "red"
                     )) +
  
  # clean theme & axis labels
  #theme_minimal(base_size = 14) +
  labs(
    x = "log2FC of proteome",
    y = "log2FC of transcriptome"
  ) +
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

###########

#gene name vs joint correlation
ggplot(alltargetsDEP, aes(x = GeneName , y = Corr_xyJnt) )+
  
  # base layer: all points in grey, no legend
  geom_point(color = "grey", size = 2) +
  geom_text_repel(
    data = alltargetsDEP %>% filter(GeneName %in% genes_to_label),
    aes(label = GeneName),
    size = 3
  ) +
  
  geom_point(data = tt, aes(color = "Candidate Transcriptional Targets"), size = 2)+
  geom_point(data = ptt, aes(color = "Candidate Post-Transcriptional Targets"), size = 2)+
  
  # manual scale: map the two labels to colors
  scale_color_manual(name   = "Target type",
                     values = c(
                       "Candidate Transcriptional Targets"     = "blue",
                       "Candidate Post-Transcriptional Targets" = "red"
                     )) 
  
#heatmap
df2 = alltargetsDEP
Gene1 = df2$GeneName
Gene2 = df2$GeneName
# 1. Pivot to a wide matrix form
mat_df <- df2 %>%
  pivot_wider(
    names_from  = Gene2,
    values_from = Corrxy_jnt,
    values_fill = NA
  )

# 2. Convert back to long form for ggplot (if needed)
plot_df <- mat_df %>%
  pivot_longer(
    -Gene1,
    names_to  = "Gene2",
    values_to = "Corrxy_jnt"
  ) %>%
  rename(GeneX = Gene1, GeneY = Gene2)

# 3. Draw the heatmap
ggplot(plot_df, aes(x = GeneX, y = GeneY, fill = Corrxy_jnt)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low    = "blue",
    mid    = "white",
    high   = "red",
    midpoint = 0,
    na.value = "grey95",
    name   = "Corrxy_jnt"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title  = element_blank()
  )

#heatmap: HUGE matrix with the same info as the scatterplot with only the diagonals colored
#cant do jnt comp cp vd joint comp tr: joint comp are multi-dimensional.

##################################
#bar graphs adam9 expression
##################################

library(ggplot2)


#adam9 data mrna

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5", "CT6",
              "KD1","KD2","KD3","KD4","KD5","KD6")),
  Condition = rep(c("Control", "Knockdown"), each = 6),
  Expression = c(17153,	16111,	12290,	15711,	14619,	16318,	#control
                 2839,	3551,	3747,	5170,	3878,	3673) #kd
)

#adam9 data proteome

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5",
              "KD1","KD2","KD3","KD4","KD5")),
  Condition = rep(c("Control", "Knockdown"), each = 5),
  Expression = c(17.79033553,	17.83873419,	17.5847391,	17.82758708,	17.81835974,
                 15.40862074,	15.61941198,	15.70166099,	15.626183,	15.44307587
) 
)

#################################

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5", "CT6",
              "KD1","KD2","KD3","KD4","KD5","KD6")),
  Condition = rep(c("Control", "Knockdown"), each = 6),
  Expression = c(17153,	16111,	12290,	15711,	14619,	16318,	#control
                 2839,	3551,	3747,	5170,	3878,	3673) #kd
)
#adam9 data proteome

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5",
              "KD1","KD2","KD3","KD4","KD5")),
  Condition = rep(c("Control", "Knockdown"), each = 5),
  Expression = c(17.79033553,	17.83873419,	17.5847391,	17.82758708,	17.81835974,
                 15.40862074,	15.61941198,	15.70166099,	15.626183,	15.44307587
  ) 
)

ctrl_mean <- mean(df$Expression[df$Condition == "Control"], na.rm = TRUE)
kd_mean   <- mean(df$Expression[df$Condition == "Knockdown"], na.rm = TRUE)

df$Expression_norm <- (df$Expression - kd_mean) / (ctrl_mean - kd_mean)

ggplot(df, aes(x = Condition, y = Expression_norm)) +
  geom_jitter(aes(shape = Condition), width = 0.12, height = 0,
              size = 2.5, alpha = 0.8) +
  scale_shape_manual(values = c("Control" = 16, "Knockdown" = 17)) +
  
  # horizontal min and max lines
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar",
               width = 0.4, size = 0.8, color = "black") +
  
  # horizontal mean line (thicker)
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 0,
               color = "black", size = 1.5) +
  
  # significance bar
  geom_signif(comparisons = list(c("Control", "Knockdown")),
              map_signif_level = TRUE, test = "t.test",
              annotations = "p-value < 0.01") +
  
  labs(x = NULL, y = "Normalized ADAM9 Expression in Proteome") +
  
  #theme(panel.border = element_rect(color = "grey", fill = NA))+
  annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf,
           color = "grey")+
  theme_minimal()+
  NoLegend()+
  scale_y_continuous(
    #limits = c(0, 1.1),
    #expand = c(0, 0),
    breaks = seq(0, 1.2, by = 0.2)
  ) +
  theme(
    #panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 0.3),
    axis.ticks.length.y = unit(2, "mm"),
    axis.line = element_line(color = "grey"),
    legend.position = "none"
  )
# run test and extract p-value
p_val <- t.test(Expression_norm ~ Condition, data = df)$p.value
#1.60975e-06 for mrna
#2.74156e-09 for proteome


###############################

df$Expression_norm <- (df$Expression - min(df$Expression)) /
  (max(df$Expression) - min(df$Expression))


ggplot(df, aes(x = Condition, y = Expression_norm)) +
  geom_jitter(aes(shape = Condition), width = 0.12, height = 0,
              size = 2.5, alpha = 0.8) +
  scale_shape_manual(values = c("Control" = 16, "Knockdown" = 17)) +
  
  # horizontal min and max lines
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar",
               width = 0.4, size = 0.8, color = "black") +
  
  # horizontal mean line (thicker)
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 0,
               color = "black", size = 1.5) +
  
  # significance bar
  geom_signif(comparisons = list(c("Control", "Knockdown")),
              map_signif_level = TRUE, test = "t.test") +
  
  labs(x = NULL, y = "Normalized ADAM9 Expression in transcriptome") +
  
  #theme(panel.border = element_rect(color = "grey", fill = NA))+
 annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf,
             color = "grey")+
  theme_minimal()+
  NoLegend()+
  
theme(panel.grid = element_blank(), axis.line = element_line(color = "grey"))


######################
##latest version of scatterplot
#####################
#adam9 data mrna

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5", "CT6",
              "KD1","KD2","KD3","KD4","KD5","KD6")),
  Condition = rep(c("Control", "Knockdown"), each = 6),
  Expression = c(17153,	16111,	12290,	15711,	14619,	16318,	#control
                 2839,	3551,	3747,	5170,	3878,	3673) #kd
)

df$rawExp = df$Expression
#adam9 data proteome

df = data.frame(
  Sample = (c("CT1", "CT2", "CT3", "CT4", "CT5",
              "KD1","KD2","KD3","KD4","KD5")),
  Condition = rep(c("Control", "Knockdown"), each = 5),
  Expression = c(17.79033553,	17.83873419,	17.5847391,	17.82758708,	17.81835974,
                 15.40862074,	15.61941198,	15.70166099,	15.626183,	15.44307587
  ) 
)

#expressions here are log2 transformed. make them raw counts
df$rawExp = 2^(df$Expression)
#########################
control_mean <- mean(df$rawExp[df$Condition == "Control"])
df$Expression_norm <- df$rawExp / control_mean


ggplot(df, aes(x = Condition, y = Expression_norm)) +
  geom_jitter(aes(shape = Condition), width = 0.12, height = 0,
              size = 2.5, alpha = 0.8) +
  scale_shape_manual(values = c("Control" = 16, "Knockdown" = 17)) +
  
  # horizontal min and max lines
  stat_summary(fun.min = min, fun.max = max, geom = "errorbar",
               width = 0.4, size = 0.8, color = "black") +
  
  # horizontal mean line (thicker)
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, fatten = 0,
               color = "black", size = 1.5) +
  
  # significance bar
  geom_signif(comparisons = list(c("Control", "Knockdown")),
              map_signif_level = TRUE, test = "t.test",
              annotations = "") +
  
  labs(x = NULL, y = "Normalized ADAM9 Expression in Transcriptome") +
  
  #theme(panel.border = element_rect(color = "grey", fill = NA))+
  #annotate("segment", x = 0, xend = 0, y = -Inf, yend = Inf,
  #color = "grey")+
  theme_minimal()+
  NoLegend()+
  scale_y_continuous(
    limits = c(0, 1.2),
    expand = c(0, 0),
    breaks = seq(0, 1.2, by = 0.2)) +
  
  theme(
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(color = "black", size = 0.3),
    axis.ticks.length.y = unit(2, "mm"),
    axis.line = element_line(color = "grey"),
    legend.position = "none",
    panel.grid = element_blank()
  ) 
# run test and extract p-value
p_val <- t.test(Expression_norm ~ Condition, data = df)$p.value
#1.60975e-06 for mrna
#5.388286e-06 for proteome

############################
##scatterplot for flowchart
############################

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
                       "Candidate Post-Transcriptional Targets" = "red"
                     )) +
  
  # clean theme & axis labels
  #theme_minimal(base_size = 14) +
  labs(
    x = "log2FC of transcriptome",
    y = "log2FC of proteome"
  ) +
  theme(
    axis.title.x = element_text(size = 23),
    axis.title.y = element_text(size = 23)
  )