setwd("../rawdata/")
library(tidyverse)
load("crc_counts_tpm_clin.Rdata")

# ID CONVERSATION
gtf_v22 <- read_tsv("gencode.gene.info.v22.tsv")

id2symbol <- gtf_v22 %>%
  dplyr::select(1, 2)

tpm_tu <- tpm_tu %>%
  inner_join(id2symbol, by = "gene_id") %>%
  relocate("gene_name", .before = "gene_id") %>%
  select(-gene_id)

tpm_tu <- tpm_tu %>%
  group_by(gene_name) %>%
  summarise_all(mean)
tpm_tu <- rename(tpm_tu, "gene symbols" = "gene_name")

write.table(tpm_tu, file = "../result/tpm_tu.txt", sep = "\t", row.names = F, col.names = T)
write.csv(tpm_tu, file = "../result/tpm_tu.csv", row.names = F, col.names = T)

# save(clin,tpm_tu,file = 'tpm_clin.Rdata')
