# Order by PM_P
library(scales)
library(patchwork)

rm(list = ls())
options(stringsAsFactors = F)

Base_edit_function <- function(edit_table, base_type) {
  edit_table <- edit_table[which(edit_table$mutate %in% base_type),]
  averages <- edit_table %>%
    group_by(Sample) %>%
    mutate(every_four = ceiling((row_number() ) / 4)) %>%  # 计算每四行的新分组
    group_by(Sample, every_four) %>%
    summarise(avg_value = mean(rate, na.rm = TRUE),
              every_four_num = mean(Pos, na.rm = TRUE)) %>% 
    ungroup() %>%
    select(Sample, every_four_num, avg_value)  # 选择需要展示的列
  
  ggplot(data=edit_table, aes(x=Pos, y=rate,fill=Sample)) +
    ylab(paste0("Edit rate")) +
    xlab(paste0("Position")) +
    ggtitle(base_type) +
    geom_line(data=averages[which(averages$Sample %in% "WXM027Vector_Prime"),], 
              aes(x=every_four_num, y=avg_value),size=0.5,colour="#0000ff") + 
    geom_line(data=averages[which(averages$Sample %in% "A107WXM027EST96h_LA"),], 
              aes(x=every_four_num, y=avg_value),size=0.5,colour="#ffa500",alpha=1) + 
    geom_point(data=edit_table[which(edit_table$Sample %in% "WXM027Vector_Prime"),],
               aes(x=Pos, y=rate),size=0.5,colour="#0000ff") +
    geom_point(data=edit_table[which(edit_table$Sample %in% "A107WXM027EST96h_LA"),],
               aes(x=Pos, y=rate),size=0.5,colour="#ffa500",alpha=1) +   
    scale_x_continuous(limits = c(0,220),breaks = seq(0, 220, 50),expand = c(0,0)) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5,color="black", size=14,face="bold"),
      axis.title.x = element_text(size=12),
      axis.title.y = element_text(size=12),
      axis.text = element_text(color = "black", size = 10),
      axis.line = element_line(linewidth =0),
      legend.position = "none",
      panel.border = element_rect(color = "black", size = 0.8, fill = NA))
}

# All base edit rate
edit_table <- read.table("/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/Merge_position_mutate.results",
                         sep = "\t",header = T)
edit_table <- edit_table[-c(225),]

edit_table$Pos <- as.numeric(edit_table$Pos)
edit_table$rate <- as.numeric(edit_table$rate)

averages <- edit_table %>%
  group_by(Sample) %>%
  mutate(every_four = ceiling((row_number() ) / 4)) %>%  # 计算每四行的新分组
  group_by(Sample, every_four) %>%
  summarise(avg_value = mean(rate, na.rm = TRUE),
            every_four_num = mean(Pos, na.rm = TRUE)) %>% 
  ungroup() %>%
  select(Sample, every_four_num, avg_value)  # 选择需要展示的列

Scatter <- ggplot(data=edit_table, aes(x=Pos, y=rate,fill=Sample)) +
  # geom_point(stat="identity",size=1) +
  ylab(paste0("Edit rate")) +
  xlab(paste0("Position")) +
  geom_line(data=averages[which(averages$Sample %in% "WXM027Vector_Prime"),], 
            aes(x=every_four_num, y=avg_value),size=0.5,colour="#0000ff") + 
  geom_line(data=averages[which(averages$Sample %in% "A107WXM027EST96h_LA"),], 
            aes(x=every_four_num, y=avg_value),size=0.5,colour="#ffa500",alpha=1) + 
  geom_point(data=edit_table[which(edit_table$Sample %in% "WXM027Vector_Prime"),],
             aes(x=Pos, y=rate),size=0.5,colour="#0000ff") +
  geom_point(data=edit_table[which(edit_table$Sample %in% "A107WXM027EST96h_LA"),],
             aes(x=Pos, y=rate),size=0.5,colour="#ffa500",alpha=1) +   
  scale_x_continuous(limits = c(0,230),breaks = seq(0, 220, 50),expand = c(0,0)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_classic() + 
  theme(
    axis.title.x = element_text(size=12),
    legend.position = "none",
    axis.title.y = element_text(size=12),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
Scatter
ggsave(Scatter,filename="/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/All_edit_rate_line.pdf",
       width = 7,height = 4)

# Per base edit rate
edit_table <- read.table("/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/Merge_position.results",
                         sep = "\t",header = F)
names(edit_table) <- c("Chr","Pos","Ref","rate","num","mutate","type","Sample")
# edit_table[which(edit_table$Pos > 160 & edit_table$Pos < 230), 4] <- 0

# C
C_to_T_table <- c("C->T")
C_to_A_table <- c("C->A")
C_to_G_table <- c("C->G")
C_to_T_plot <- Base_edit_function(edit_table,C_to_T_table)
C_to_A_plot <- Base_edit_function(edit_table,C_to_A_table)
C_to_G_plot <- Base_edit_function(edit_table,C_to_G_table)

# A
A_to_T_table <- c("A->T")
A_to_C_table <- c("A->C")
A_to_G_table <- c("A->G")
A_to_T_plot <- Base_edit_function(edit_table,A_to_T_table)
A_to_C_plot <- Base_edit_function(edit_table,A_to_C_table)
A_to_G_plot <- Base_edit_function(edit_table,A_to_G_table)

# G
G_to_T_table <- c("G->T")
G_to_C_table <- c("G->C")
G_to_A_table <- c("G->A")
G_to_T_plot <- Base_edit_function(edit_table,G_to_T_table)
G_to_C_plot <- Base_edit_function(edit_table,G_to_C_table)
G_to_A_plot <- Base_edit_function(edit_table,G_to_A_table)

# T 
T_to_G_table <- c("T->G")
T_to_C_table <- c("T->C")
T_to_A_table <- c("T->A")
T_to_G_plot <- Base_edit_function(edit_table,T_to_G_table)
T_to_C_plot <- Base_edit_function(edit_table,T_to_C_table)
T_to_A_plot <- Base_edit_function(edit_table,T_to_A_table)

Merge_plot <- C_to_A_plot + C_to_T_plot + C_to_G_plot +
  G_to_A_plot + G_to_T_plot + G_to_C_plot +
  T_to_A_plot + T_to_C_plot + T_to_G_plot +
  A_to_C_plot + A_to_T_plot + A_to_G_plot +
  plot_layout(nrow = 4)
Merge_plot

ggsave(Merge_plot,filename="/home/qzh/Analysis/Shengwang/YY583_base_edit/02-result/plot/Per_base_edit_rate_line.pdf",
       width = 12,height = 12)
