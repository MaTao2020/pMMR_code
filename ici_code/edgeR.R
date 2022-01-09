# 
library(tidyverse)
library(edgeR)
#edgeR包可以做无重复的差异分析，不过需要认为指定一个dispersion值(设置BCV值)，这样得到的结果比较主观，
#不同的人就可以有不同的结果。通常如果是实验控制的好的人类数据，那么选择BCV=0.4，比较好的模式生物选择BCV=0.1。
getwd()
setwd('../rawdata/')

load('gse179351_counts_clin.Rdata')


gse_expr <- gse_expr %>% 
  group_by(gene.symbol) %>% 
  summarise_all(median) %>% 
  column_to_rownames('gene.symbol')


group_list <- factor(c(rep("Contral",1),rep("Treat",1)))

#数据预处理
#（1）构建 DGEList 对象
dgelist <- DGEList(counts = gse_expr, group = group_list)
#（2）过滤 low count 数据，例如 CPM 标准化（推荐）
keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')

bcv = 0.4  #设置BCV值
et <- exactTest(dgelist, dispersion=bcv^2)

genes <- topTags(et,n = nrow(dgelist$counts))
genes <- genes$table
diff_genes <- genes %>% 
  filter(FDR < 0.05)

write.csv(diff_genes,file = '../results/diff_genes.csv')

# plot Volcano


DEG_all <- genes %>% 
  rownames_to_column("symbol") %>% 
  dplyr::select(symbol, logFC, Pvalue = FDR) %>% #选择三列
  mutate(direction = factor(ifelse(Pvalue < 0.05 & abs(logFC) > 0,#添加direction一列
                                   ifelse(logFC > 0, "Up", "Down"),"NS"),
                            levels=c('Up','Down','NS')))


library(ggrepel)
library(ggsci)
ggplot(data = DEG_all, aes(x = logFC, y = -log10(Pvalue), colour = direction)) + #数据映射
  geom_point(alpha = 0.6) +#散点图，alpha就是点的透明度
  scale_color_manual(values = c("#bc3d2a", "#0173b5", "#808080")) + #手动调颜色
  geom_text_repel(data = DEG_all %>% filter(Pvalue < 0.01, abs(logFC) > 1),#加注释，筛选出差异非常显著的基因
                  aes(label = symbol),#便签就是基因名 geom_text_repel：adds text directly to the plot.
                  size = 3,
                  segment.color = "black", #连接线的颜色，就是名字和点之间的线
                  show.legend = FALSE) + 
  theme_bw() +#设定主题
  theme(legend.title = element_blank()) +#不想要direction这个title
  ylab(expression(-log[10]("Adjusted P Value"))) +#expression的作用就是让log10的10下标
  xlab(expression(log[2]("Fold Change"))) +
  xlim(-7, 7) +
  ylim(0, 7) +
  geom_vline(xintercept = c(-1, 1), #加垂直线，在-1和+1之间划线
             lty = 2,
             col = "black",
             lwd = 0.6) +
  geom_hline(yintercept = -log10(0.05),
             lty = 2,
             col = "black",
             lwd = 0.6)  

ggsave(filename = '../results/ICI_diffgenes_vocalno.pdf',width = 5,height = 4)













