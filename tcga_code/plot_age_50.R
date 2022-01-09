library(DESeq2)
library(tidyverse)
library(DEGreport)
library(RColorBrewer)
library(colorspace)    
setwd('../rawdata/')
load('DESQ2_dds_counts_tpm_age50.Rdata')


vsd = vst(dds,blind=FALSE)
vsd_df <- assay(vsd) 
 #plotPCA(vsd,"gender")
 #plotPCA(vsd,'stage')

#batch
vsd_df <- limma::removeBatchEffect(vsd_df, coveriate =model.matrix(~gender+age+pathologic_M+pathologic_N+pathologic_T,data = clin),
                                   design = model.matrix(~stage,data=clin))


#identical(colnames(vsd),rownames(clin))

res_LRT <- results(dds)
res_LRT

res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigLRT_genes <- res_LRT_tb %>% 
  filter(pvalue < 0.05)

cluster_rlog <- vsd_df[sigLRT_genes$gene, ]
clusters <- degPatterns(cluster_rlog, metadata = clin, time = "stage", col=NULL,minc = 15)

cluster_groups <- clusters$normalized
up_gene <- filter(cluster_groups,cluster_groups$cluster == '1')
down_gene <- filter(cluster_groups,cluster_groups$cluster =='2')
up_down_gene <- filter(cluster_groups,cluster_groups$cluster %in% c(1,2))

#plot
#scale_color_discrete_qualitative(palette = "Dark3")
library(ggsci)
degPlotCluster(up_down_gene,time="stage",color="stage")+ 
  geom_line(aes_string(group="genes"),alpha = 0.05)+
  scale_color_nejm()+
  theme_bw()

ggsave(filename = '../result/stage_updown.pdf',width = 8,height = 4)



#ID CONVERSATION
gtf_v22 <- read_tsv("gencode.gene.info.v22.tsv")

id2symbol <- gtf_v22 %>% 
  dplyr::select(1, 2)


up_down_gene <- distinct(up_down_gene,genes,.keep_all = T)
id2symbol <- dplyr::rename(id2symbol,genes=gene_id)
up_down_gene <- up_down_gene %>% 
  inner_join(id2symbol,by='genes')
up_gene <- filter(up_down_gene,cluster == '1')
down_gene <- filter(up_down_gene,cluster == '2')


write.csv(up_down_gene,file = '../result//up_down_genes_vst_age50_Pvalue0.05.csv',row.names = F)
write.csv(up_gene,file = '../result/up_genes_vst_age50_Pvalue0.05.csv',row.names = F)
write.csv(down_gene,file = '../result/down_genes_vst_age50_Pvalue0.05.csv',row.names = F)

###









