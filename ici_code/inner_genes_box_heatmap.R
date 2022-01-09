
library(tidyverse)
library(GEOquery)
library(data.table)
GSE179351_rawdata <- getGEO("GSE179351",
  GSEMatrix = TRUE,
  getGPL = F
)

### 提取临床数据
GSE179351_pdata <- pData(GSE179351_rawdata[[1]])
head(GSE179351_pdata) ## 查看临床数据,发现这个是经过RMA标准化的


## 提取表达数据
GSE179351_exprSet <- exprs(GSE179351_rawdata[[1]]) # 无数据

GSE179351 <- fread("../rawdata/GSE179351_RawCountsForAllSamples.txt/GSE179351_RawCountsForAllSamples.txt")
GSE179351 <- select(GSE179351, gene.symbol, gene.type, starts_with("X"))
GSE179351 <- GSE179351 %>%
  filter(GSE179351$gene.type == "protein_coding") %>%
  select(gene.symbol, starts_with("X"))

# 4\74\73\68\18\30\ Pre-Tx  Pre-xRT

gse_expr <- select(
  GSE179351, gene.symbol, X4.Pre.Tx, X4.Pre.xRT, X73.Pre.Tx, X73.Pre.xRT,
  X74.Pre.Tx, X74.Pre.xRT, X68.Pre.Tx, X68.Pre.xRT, X18.Pre.Tx,
  X18.Pre.xRT, X30.Pre.Tx, X30.Pre.xRT
)


sample_id <- c(
  "X4.Pre.Tx", "X4.Pre.xRT", "X73.Pre.Tx", "X73.Pre.xRT",
  "X74.Pre.Tx", "X74.Pre.xRT", "X68.Pre.Tx", "X68.Pre.xRT", "X18.Pre.Tx",
  "X18.Pre.xRT", "X30.Pre.Tx", "X30.Pre.xRT"
)

sample_id <- sample_id %>%
  str_replace_all("\\.", "-") %>%
  str_sub(2, str_length(.))


gse_clin <- filter(GSE179351_pdata, GSE179351_pdata$title %in% sample_id)

genes <- read_csv("../results/ici_updown_genes.csv")

genes <- genes %>%
  select(gene_name) %>%
  pull()

sample_id1 <- c(
  "73-Pre-Tx", "74-Pre-Tx", "68-Pre-Tx", "18-Pre-Tx",
  "30-Pre-Tx"
)

sample_id2 <- c("73-Pre-xRT", "74-Pre-xRT", "68-Pre-xRT", "18-Pre-xRT", "30-Pre-xRT")

clin <- gse_clin %>%
  select(title) %>%
  mutate(group = case_when(
    title == "4-Pre-Tx" ~ "R-pre",
    title == "4-Pre-xRT" ~ "R-post",
    title %in% sample_id1 ~ "NR-pre",
    title %in% sample_id2 ~ "NR-post"
  ))


## tpm
gtf_v22 <- read_tsv("gencode.gene.info.v22.tsv")

tpm <- gse_expr %>%
  inner_join(gtf_v22, by = c("gene.symbol" = "gene_name")) %>%
  mutate(across(
    .cols = starts_with("X"),
    .fns = ~ (((. * 10^9) / (sum(.[gene_type == "protein_coding"]) * exon_length)) * 10^6 /
      (sum((. * 10^9) / (sum(.[gene_type == "protein_coding"]) * exon_length))))
  )) %>%
  select(colnames(gse_expr))

genes <- read_csv("genes.csv")

genes <- genes %>%
  pull(x)

tpm <- tpm %>%
  filter(gene.symbol %in% genes) %>%
  rename_with(~ str_replace_all(., "\\.", "-"), everything()) %>%
  column_to_rownames("gene-symbol")

tpm <- log2(tpm + 1)

tpm <- tpm %>%
  t(.) %>%
  as.data.frame() %>%
  scale()

rownames(tpm) <- str_sub(rownames(tpm), 2, str_length(rownames(tpm)))
tpm <- t(tpm)

tpm1 <- tpm %>%
  as.data.frame(.) %>%
  dplyr::select(contains("Tx"))



library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(
  c(-4, -2, 0, 2, 4),
  c("#4A6FE3", "#9DA8E2", "#F9F9F9", "#E495A5", "#D33F6A")
)
row_anno <- HeatmapAnnotation(gp = gpar(col = c('#b6503d','#3f7bb1','#de975f','#68946f')),
  group = c(
    "R-pre", "R-post", "NR-pre", "NR-post",
    "NR-pre", "NR-post", "NR-pre", "NR-post",
    "NR-pre", "NR-post", "NR-pre", "NR-post"
  ), which = "colum",
  show_annotation_name = T
)
pdf(file = '27genes_R_NR_heatmap.pdf',width =8,height = 8 )
Heatmap(tpm,
  cluster_columns = F, show_row_dend = F, show_column_names = F,
  show_row_names = T,top_annotation = row_anno)
dev.off()
