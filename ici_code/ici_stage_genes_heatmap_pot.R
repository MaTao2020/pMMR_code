library(tidyverse)
library(Seurat)
load('../../2.scrna/rawdata/scs.Rdata')
load('../../2.scrna/rawdata/sc_rna.Rdata')

sc <- merge(scrna,scs)
metadata <- sc@meta.data
metadata <- metadata %>% 
  filter(!Sub_Cluster %in% c('hC01_P0413','hC02_P0825','hC03_P0411','hC04_P1212'))

  genes <- read_csv('genes.csv')
gene_name <- genes$x  

df <- FetchData(object=sc,vars=gene_name)
df <- df[rownames(metadata),]

df[df>0] <- 1

df <- metadata %>% 
  select(Sub_Cluster) %>% 
  cbind(df) %>% 
  select(Sub_Cluster, everything() )

df_percent <- df %>% 
  group_by(Sub_Cluster) %>% 
  summarise(n = n()
            ,across(where(is.numeric),sum)/n*100) %>% 
  select(-n) 


#plot
df_scale <-df_percent %>% 
  column_to_rownames('Sub_Cluster') %>% 
  scale()


library(colorspace)
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(
  c(-4,-2,0,2,4), 
  c('#4A6FE3','#9DA8E2','#F9F9F9','#E495A5','#D33F6A')
)

Heatmap(df_scale,cluster_rows = T,cluster_columns = T,show_column_names = F,
        show_row_dend = F,show_column_dend = F,col = col_fun)






















