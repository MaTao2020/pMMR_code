#percent heatmap

library(tidyverse)
library(Seurat)
getwd()
setwd('../rawdata/')
load('scs.Rdata')#smart
load('sc_rna.Rdata')#10x

sc <- merge(scrna,scs)
metadata <- sc@meta.data
metadata <- metadata %>% 
  filter(!Sub_Cluster %in% c('hC01_P0413','hC02_P0825','hC03_P0411','hC04_P1212'))

genes <- read_csv('up_down_genes_vst_age50_Pvalue0.05.csv')
gene_name <- genes$gene_name

df <- FetchData(object=sc,vars=gene_name)#701个基因没被匹配到，901个匹配
df <- df[rownames(metadata),]
###### %-based heatmaps and tables
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

write.csv(df_percent,file = '../results/heatmap_up_down_percent.csv')



up_gene <- genes %>% 
  filter(cluster == '1') %>% 
  select(gene_name) %>% 
  pull() %>% 
  intersect(colnames(df_percent))

up_df_percent <- df_percent %>% 
 select(Sub_Cluster,up_gene)

down_gene <- setdiff(colnames(df_percent),colnames(up_df_percent))

down_df_percent <- df_percent %>% 
  select(Sub_Cluster,all_of(down_gene))

write.csv(up_df_percent,file = '../results/heatmap_up_percent.csv')
write.csv(down_df_percent,file = '../results/heatmap_down_percent.csv')

save(df_percent,up_df_percent,down_df_percent,file = 'sc_percent.Rdata')

# 
# metadata <- scs@meta.data
# 
# MSI <- c('P0825','P0413')
# 
# metadata <- metadata %>% 
#   filter(!Sub_Cluster %in% c('hC01_P0413','hC02_P0825','hC03_P0411','hC04_P1212'))
# 
# meta_msi <- filter(metadata,Sample %in% MSI)
# meta_mss <- filter(metadata,!Sample %in% MSI)
# 
# msi_meta <- meta_msi %>% 
#   group_by(Sub_Cluster) %>% 
#   summarise(count=n(),percent_msi = (count/3418)*100) %>% 
#   mutate(group_msi = 'msi') %>% 
#   select(type = Sub_Cluster,percent_msi,group_msi)
# 
# 
# mss_meta <- meta_mss %>% 
#   group_by(Sub_Cluster) %>% 
#   summarise(count=n(),percent_mss = (count/6781)*100) %>% 
#   mutate(group_mss = 'mss') %>% 
#   select(type = Sub_Cluster,percent_mss,group_mss)
# 
# 
# meta <- full_join(msi_meta,mss_meta,by = 'type')
# 
# meta <- replace_na(meta,list(percent_msi = 0,group_msi = 'msi'))
# 
# meta %>% 
#   pivot_longer(cols = starts_with('percent'),names_to = 'group',values_to = 'percent') %>% 
#   ggplot(aes(x=type,y=percent,fill = group))+geom_col(position = position_dodge())+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1,size = 10),axis.text.y = element_text(size = 10))
# 
# 
# meta <- mutate(meta,percent = percent_msi/percent_mss)
# 
# meta %>%
#   ggplot(aes(x=type,y=percent))+geom_col(position = position_dodge())+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90,vjust = 1,hjust = 1,size = 10),axis.text.y = element_text(size = 10))
# 
# 
# 
# 










