library(dplyr)
library(ggplot2)


rm(list = ls())
options(stringsAsFactors = F)

Base_edit_table <- read.table("/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/A107_plot_v2.txt",
                              sep = "\t",header = F)
names(Base_edit_table) <- c("chr","position","ref","rate","num","mutate","type","Sample")

level <- c("C>A","C>T","C>G","G>A","G>T","G>C","A>T","A>C","A>G","T>A","T>C","T>G",
           "C->A","C->T","C->G","G->A","G->T","G->C","A->T","A->C","A->G","T->A","T->C","T->G")
Base_edit_table$mutate <- factor(Base_edit_table$mutate, level = level)

Scatter <- ggplot(data=Base_edit_table, aes(x=mutate, y=rate,color=mutate)) +
  geom_jitter(position=position_jitter(0.1)) +
  geom_boxplot(outlier.colour=NA, width=.3,position= position_nudge(x=-0)) +
  theme_classic() +
  theme(legend.position = "none") + 
  ylab("Ratio") +
  xlab("") +
  geom_hline(aes(yintercept=0.002),linetype="dashed") +
  scale_color_manual(values=c("A->G" = "blue", "T->C" = "blue",
                              "A>G" = "red", "T>C" = "red")) + 
  theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=12),
    axis.line = element_line(linewidth =0),
    axis.text.x = element_text(size=12, color = "black",angle = 90,vjust = 0.5),
    axis.text.y = element_text(size=12, color = "black"),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA)) + 
  geom_vline(xintercept=c(12.5),size=.4,) 

Scatter
ggsave(Scatter,filename="/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/Per_base_mutate_ratio.pdf",
       width = 8,height = 3.5)
