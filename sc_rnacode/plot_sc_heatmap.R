#scale data
setwd('../rawdata/')
load('sc_percent.Rdata')
library(tidyverse)
up_df_scale <- up_df_percent %>% 
  column_to_rownames('Sub_Cluster') %>% 
  scale()

down_df_scale <- down_df_percent %>% 
  column_to_rownames('Sub_Cluster') %>% 
  scale()
df_scale <- df_percent %>% 
  column_to_rownames('Sub_Cluster') %>% 
  scale()

#plot
library(colorspace)
library(ComplexHeatmap)
library(circlize)
#hcl_palettes(plot = TRUE)
col_fun <- colorRamp2(
  c(-4,-2,0,2,4), 
  c('#4A6FE3','#9DA8E2','#F9F9F9','#E495A5','#D33F6A')
)
pdf(file = '../results/up_stagegenes_heatmap.pdf',width = 7,height = 9)
Heatmap(up_df_scale,cluster_rows = T,cluster_columns = T,show_column_names = F,
        show_row_dend = F,show_column_dend = F,col = col_fun)
dev.off()


pdf(file = '../results/down_stagegenes_heatmap.pdf',width = 7,height = 9)
Heatmap(down_df_scale,cluster_rows = T,cluster_columns = T,show_column_names = F,
        show_row_dend = F,show_column_dend = F,col = col_fun)
dev.off()

