library(tidyverse)
setwd('../rawdata/')
#load clin data
coad_clin <- read_tsv('TCGA-COAD.GDC_phenotype.tsv/TCGA-COAD.GDC_phenotype.tsv')
read_clin <- read_tsv('TCGA-READ.GDC_phenotype.tsv/TCGA-READ.GDC_phenotype.tsv')

#merge
setdiff(colnames(coad_clin),colnames(read_clin))

coad_clin <- select(coad_clin,-c(dbgap_registration_code,program))

crc_clin <- rbind(coad_clin,read_clin)

write.csv(crc_clin,file = '../result/crc_clin.csv',row.names = F,col.names = T,na = ' ')
