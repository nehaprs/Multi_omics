library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)
library(ggbreak)
setwd("~/BINF/o2pls/for manuscript/revision/volcano")

#volcanoplot for cp

res <- read_excel("cp_res.xlsx")
# cut‐offs
fc.cut   <- log2(1.2)
fdr.cut1 <- 0.05
fdr.cut2 <- 0.01

# Prep data
res2 <- res %>%
  mutate(
    negLogFDR = -log10(FDR),
    group = case_when(
      log2FC >  0 & FDR < fdr.cut1 ~ "Upregulated",
      log2FC < 0 & FDR < fdr.cut1 ~ "Downregulated",
      TRUE                                ~ "Not significant"
    ),
    group = factor(group, levels = c("Downregulated","Not significant","Upregulated"))
  )

# Compute y‐positions for the cut‐off lines
y1 <- -log10(fdr.cut1)
y2 <- -log10(fdr.cut2)

p = ggplot(res2, aes(x = log2FC, y = negLogFDR, color = group)) +
  
  # Points
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Downregulated"   = "blue",
    "Not significant" = "grey70",
    "Upregulated"     = "red"
  )) +
  xlim(-4,4) +
  # Vertical FC cutoffs
  geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dotted") +
  
  # Horizontal FDR cutoffs
  geom_hline(yintercept = y1, linetype = "dashed") +
  geom_hline(yintercept = y2, linetype = "dashed") +
  
  #label axes
  labs(
    x = expression(log[2](KD/WT)),
    y = expression(-log[10](FDR)),
    color = NULL
  )+

  # Label the horizontal lines
  annotate("text", x = -3, y = y1 + 0.1,
           label = "FDR = 0.05", hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = -3, y = y2 + 0.1,
           label = "FDR = 0.01", hjust = 0, vjust = 0, size = 3) +
  
  
  # Label the vertical lines
  annotate("text", x = -fc.cut -1, y = 0 - 0.1,
           label = "FC = -1.2", hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = fc.cut + 0.1, y = 0 - 0.1,
           label = "FC = 1.2", hjust = 0, vjust = 0, size = 3) +
  
  # Highlight & label ADAM9, nudged upward
  geom_text_repel(
    data = filter(res2, gene == "ADAM9"),
    aes(label = gene),
    nudge_y = 0.3,
    direction = "y",
    segment.size = 0.2,
    size = 4,
    color = "black"
  ) 

###############

#volcano plot for rnaseq

res <- read_excel("mrnares.xlsx")
# cut‐offs
fc.cut   <- log2(1.5)
fdr.cut1 <- 0.05
fdr.cut2 <- 0.01

# Prep data
res2 <- res %>%
  mutate(
    negLogFDR = -log10(FDR),
    group = case_when(
      log2FC >  0 & FDR < fdr.cut1 ~ "Upregulated",
      log2FC < 0 & FDR < fdr.cut1 ~ "Downregulated",
      TRUE                                ~ "Not significant"
    ),
    group = factor(group, levels = c("Downregulated","Not significant","Upregulated"))
  )

# Compute y‐positions for the cut‐off lines
y1 <- -log10(fdr.cut1)
y2 <- -log10(fdr.cut2)

 ggplot(res2, aes(x = log2FC, y = negLogFDR, color = group)) +
  
  # Points
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Downregulated"   = "blue",
    "Not significant" = "grey70",
    "Upregulated"     = "red"
  )) +
  xlim(-3,3) +
  # Vertical FC cutoffs
  geom_vline(xintercept = c(-fc.cut, fc.cut), linetype = "dotted") +
  
  # Horizontal FDR cutoffs
  geom_hline(yintercept = y1, linetype = "dashed") +
  #geom_hline(yintercept = y2, linetype = "dashed") +
  
  #label axes
  labs(
    x = expression(log[2](KD/WT)),
    y = expression(-log[10](FDR)),
    color = NULL
  )+
  
  # Label the horizontal lines
  annotate("text", x = -3, y = y1 + 1,
           label = "FDR = 0.05", hjust = 0, vjust = 0, size = 3) +
  #annotate("text", x = -3, y = y2 - 0.1,
   #        label = "FDR = 0.01", hjust = 0, vjust = 0, size = 3) +
  
  
  # Label the vertical lines
  annotate("text", x = -fc.cut -1, y = 0 - 1,
           label = "FC = -1.5", hjust = 0, vjust = 0, size = 3) +
  annotate("text", x = fc.cut + 0.1, y = 0 - 1,
           label = "FC = 1.5", hjust = 0, vjust = 0, size = 3) +
  
  #break on y-axis
    scale_y_break(c(50,130)) +
  
  # Highlight & label ADAM9, nudged upward
  geom_text(
    data = filter(res2, gene == "ADAM9", negLogFDR > 130),
    aes(label = gene),
    nudge_y = 3, direction = "y", segment.size = 0.2,
    color = "black"
  ) 
