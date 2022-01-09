library(tidyverse)
library(DESeq2)
setwd('../rawdata/')
load('crc_counts_tpm_clin.Rdata')

clin$age <- case_when(clin$age %in% c(30:50) ~ '30_50',  
                      clin$age %in% c(51:90) ~ '51_90',
                      TRUE ~ as.character(clin$age)
)

sum(is.na(clin))#在pathologic_M有两个NA  复制为MX

clin$pathologic_M <- replace_na(clin$pathologic_M,'MX')
clin$age <- factor(clin$age)
clin$gender <- factor(clin$gender,levels = c('female','male'),labels = c('1','0'))
clin$status <- factor(clin$status,levels = c('Dead','Alive'),labels = c('1','0'))
clin$pathologic_M <- as.factor(clin$pathologic_M)
clin$pathologic_N <- as.factor(clin$pathologic_N)
clin$pathologic_T <- as.factor(clin$pathologic_T)
clin$stage <- as.factor(clin$stage)
clin$tumor_subtype <- factor(clin$tumor_subtype,levels = c('Adenomas and Adenocarcinomas','Complex Epithelial Neoplasms',
                                                           'Cystic, Mucinous and Serous Neoplasms','Epithelial Neoplasms, NOS'),
                             labels = c('1','2','3','4'))


clin <- column_to_rownames(clin,'sample_id')

identical(rownames(clin),colnames(counts_tu))

dds <- DESeqDataSetFromMatrix(countData=counts_tu, colData=clin, 
                              design=~gender+age+pathologic_M+pathologic_N+pathologic_T+stage) 

dds <- DESeq(dds,test="LRT", 
             reduced =~gender+age+pathologic_M+pathologic_N+pathologic_T,
             parallel = T)

resultsNames(dds)

save(clin,counts_tu,tpm_tu,dds,file = 'DESQ2_dds_counts_tpm_age50.Rdata')




#dds <- readRDS('desq_dds.Rds')




