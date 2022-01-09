library(tidyverse)
setwd('../../rawdata/')

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
clin <-  read_excel('../result/crc_clin_clean_dmmr.xlsx')

clin <- distinct(clin,sample_id,.keep_all = T) #全部是dMMR

clin <- clin %>% 
  replace_na(list(days_to_last_follow_up = 0,days_to_death = 0)) %>% 
  mutate(time = ifelse(days_to_last_follow_up >= days_to_death , days_to_last_follow_up,
                       days_to_death)) %>% 
  select(sample_id,gender,age, status ,time,pathologic_M,pathologic_N,pathologic_T,stage) %>% 
  mutate(group = ifelse(str_sub(sample_id,14,15) == '11','Normal','Tumor')) %>% 
  filter(group == 'Tumor')

clin2 <- inner_join(clin,metadata,by = c('sample_id' = 'TCGA_id'))

clin <- filter(clin,sample_id %in% clin2$sample_id)

colnames(counts_df) <- str_sub(colnames(counts_df),1,16)
colnames(tpm) <- str_sub(colnames(tpm),1,16)

counts_tu <- select(counts_df,clin2$sample_id)
tpm_tu <- select(tpm,gene_id,clin2$sample_id)

# ID CONVERSATION
gtf_v22 <- read_tsv("gencode.gene.info.v22.tsv")

id2symbol <- gtf_v22 %>%
  dplyr::select(1, 2)

tpm_tu <- tpm_tu %>%
  inner_join(id2symbol, by = "gene_id") %>%
  relocate("gene_name", .before = "gene_id") %>%
  select(-gene_id)

tpm_tu <- tpm_tu %>%
  group_by(gene_name) %>%
  summarise_all(mean)
tpm_tu <- rename(tpm_tu, "gene symbols" = "gene_name")

write.table(tpm_tu, file = "../result/tpm_tu_dmmr.txt", sep = "\t", row.names = F, col.names = T)
write.csv(tpm_tu, file = "../result/tpm_tu_dmmr.csv", row.names = F, col.names = T)
#clin

clin$age <- case_when(clin$age %in% c(30:50) ~ '30_50',  
                      clin$age %in% c(51:90) ~ '51_90',
                      TRUE ~ as.character(clin$age)
)

sum(is.na(clin))

clin$age <- factor(clin$age)
clin$gender <- factor(clin$gender,levels = c('female','male'),labels = c('1','0'))
clin$status <- factor(clin$status,levels = c('Dead','Alive'),labels = c('1','0'))
clin$pathologic_M <- as.factor(clin$pathologic_M)
clin$pathologic_N <- as.factor(clin$pathologic_N)
clin$pathologic_T <- as.factor(clin$pathologic_T)
clin$stage <- as.factor(clin$stage)
write.csv(clin,'../code/dmmr/clin_dmmr.csv')










