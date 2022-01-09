library(tidyverse)

stage_genes <- read_csv('../rawdata/up_down_genes_vst_age50_Pvalue0.05.csv')
ici_genes <- read_csv('../rawdata/diff_genes.csv')

genes <- intersect(stage_genes$gene_name,ici_genes$X1)
genes1 <- inner_join(stage_genes,ici_genes,by = c('gene_name'='X1'))

#write.csv(genes1,file = '../results/ici_updown_genes.csv')

#plot venn
library(ggvenn)
library(ggsci)
a <- list(`ICI genes` = ici_genes$X1,
          `stage genes` = stage_genes$gene_name)

ggvenn(a,c("ICI genes","stage genes"),show_percentage = F,
       stroke_color = "white",stroke_size = 1,
       fill_color = c("#0173b5","#bb3c2a"),
       set_name_color =c("#0173b5","#bb3c2a"),
       text_size = 5)
ggsave(filename = '../results/ici_stage_venn.pdf',height = 4,width = 6)
