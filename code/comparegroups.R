#绘制临床基线表
library(tidyverse)
library(compareGroups)
setwd('../rawdata/')
load('crc_counts_tpm_clin.Rdata')

clin$age <- case_when(clin$age %in% c(30:50) ~ '30_50',  
                      clin$age %in% c(51:90) ~ '51_90',
                      TRUE ~ as.character(clin$age)
)

sum(is.na(clin))#在pathologic_M有两个NA  复制为MX

clin$pathologic_M <- replace_na(clin$pathologic_M,'MX')

clin <- select(clin,-group)

tab <- descrTable(~.,data = clin)
export2word(tab, file='../result/table1.docx')






















