library(tidyverse)
library(readxl)
setwd('../../rawdata/')
clin <-  read_excel('../result/crc_clin_clean_dmmr.xlsx')
clin <- distinct(clin,sample_id,.keep_all = T) #全部是dMMR
clin <- clin %>% 
  replace_na(list(days_to_last_follow_up = 0,days_to_death = 0)) %>% 
  mutate(time = ifelse(days_to_last_follow_up >= days_to_death , days_to_last_follow_up,
                       days_to_death)) %>% 
  select(sample_id,gender,age, status ,time,pathologic_M,pathologic_N,pathologic_T,stage) %>% 
  mutate(group = ifelse(str_sub(sample_id,14,15) == '11','Normal','Tumor')) %>% 
  filter(group == 'Tumor')
library(compareGroups)
clin <- select(clin,-group)
tab <- descrTable(~.,data = clin)
export2word(tab, file='../result/table1——dmmr.docx')
