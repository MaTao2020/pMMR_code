###
library(tidyverse)
library(GEOquery)
library(data.table)
GSE179351_rawdata <- getGEO('GSE179351',
                           
                           GSEMatrix=TRUE,
                           
                           getGPL = F)

###提取临床数据
GSE179351_pdata=pData(GSE179351_rawdata[[1]])
head(GSE179351_pdata) ##查看临床数据,发现这个是经过RMA标准化的


##提取表达数据
GSE179351_exprSet=exprs(GSE179351_rawdata[[1]])#无数据

GSE179351 <- fread('../rawdata/GSE179351_RawCountsForAllSamples.txt/GSE179351_RawCountsForAllSamples.txt')
GSE179351 <- select(GSE179351,gene.symbol,gene.type,starts_with('X'))
GSE179351 <- GSE179351 %>% 
  filter(GSE179351$gene.type == 'protein_coding') %>% 
  select(gene.symbol,starts_with('X'))
 #其中只需要4号病人
gse_clin <- filter(GSE179351_pdata,GSE179351_pdata$title %in% c('4-Pre-Tx','4-Pre-xRT'))

gse_expr <- select(GSE179351,gene.symbol,starts_with('X4.Pre'))

save(gse_clin,gse_expr,file = '../rawdata/gse179351_counts_clin.Rdata')

