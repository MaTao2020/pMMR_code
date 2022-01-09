library(tidyverse)
library(GEOquery)
library(data.table)
library(DESeq2)
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

pre_r <- c('50-Pre-Tx','33-Pre-Tx','41-Pre-Tx')
pre_nr <- c('52-Pre-Tx','63-Pre-Tx','3-Pre-Tx')

expr <- GSE179351 %>% 
  rename_with(~ str_sub(.,2,str_length(colnames(GSE179351))),everything()) %>% 
  rename_with(~str_replace_all(.,'\\.','-'),everything()) %>% 
  select(`ene-symbol`,pre_r,pre_nr) %>% 
  group_by(`ene-symbol`) %>% 
  summarise_all(median) %>% 
  column_to_rownames('ene-symbol') %>% 
  round(digits = 0)

df <- data.frame(sample_name=c('50-Pre-Tx','33-Pre-Tx','41-Pre-Tx',
                               '52-Pre-Tx','63-Pre-Tx','3-Pre-Tx'),
                 group = c('r','r','r','nr','nr','nr'))

df <- column_to_rownames(df,'sample_name')

identical(colnames(expr),rownames(df))

dds <- DESeqDataSetFromMatrix(countData = expr,
                              colData = df,
                              design = ~group)


dds <- DESeq(dds, parallel = T)

resultsNames(dds)

res <- results(dds)
des_gene <- res %>% 
  as.data.frame(.) %>% 
  filter(padj < 0.05)

write.csv(des_gene,file = '../results/pdac_diff_genes.csv')
