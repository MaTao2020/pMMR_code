library(tidyverse)
setwd('../rawdata/')

load('crc_rawdata.Rdata')
load('crc_tpm.Rdata')


df <- counts_df %>% 
  column_to_rownames('gene_id')

#counts_df <- df
#过滤在所有样本都是零表达的基因 
counts_df <- df[apply(df,1, function(x) sum(x>1) > 1),]

tpm <- filter(tpm,gene_id %in% rownames(counts_df))
rm(df)


#load clinical data
library(readxl)
clin <-  read_excel('../result/crc_clin_clean.xlsx')

clin <- distinct(clin,sample_id,.keep_all = T) #全部是pMMR

clin <- clin %>% 
  replace_na(list(days_to_last_follow_up = 0,days_to_death = 0)) %>% 
  mutate(time = ifelse(days_to_last_follow_up >= days_to_death , days_to_last_follow_up,
                       days_to_death)) %>% 
  select(sample_id,gender,age, status ,time,pathologic_M,pathologic_N,pathologic_T,stage,
         tumor_type,tumor_subtype = disease_type) %>% 
  mutate(group = ifelse(str_sub(sample_id,14,15) == '11','Normal','Tumor')) %>% 
  filter(group == 'Tumor')

clin2 <- inner_join(clin,metadata,by = c('sample_id' = 'TCGA_id'))

clin <- filter(clin,sample_id %in% clin2$sample_id)

colnames(counts_df) <- str_sub(colnames(counts_df),1,16)
colnames(tpm) <- str_sub(colnames(tpm),1,16)

counts_tu <- select(counts_df,clin2$sample_id)
tpm_tu <- select(tpm,gene_id,clin2$sample_id)

#共计398例
save(clin,counts_tu,tpm_tu,file = 'crc_counts_tpm_clin.Rdata')












