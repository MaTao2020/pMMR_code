library(tidyverse)
setwd('../rawdata/')

json_file <- 'metadata.cart.2021-09-14.json'
metadata <- jsonlite::read_json(path = json_file, simplifyVector = T)
metadata <- tibble::tibble(
  file_name = metadata$file_name,
  md5sum = metadata$md5sum,
  TCGA_id_full = bind_rows(metadata$associated_entities)$entity_submitter_id,
  TCGA_id = stringr::str_sub(TCGA_id_full, 1, 16),
  patient_id = stringr::str_sub(TCGA_id, 1, 12),
  tissue_type_id = stringr::str_sub(TCGA_id, 14, 15),
  tissue_type = sapply(tissue_type_id, function(x) {
    switch(x,
           "01" = "Primary Solid Tumor",
           "02" = "Recurrent Solid Tumor",
           "03" = "Primary Blood Derived Cancer - Peripheral Blood",
           "05" = "Additional - New Primary",
           "06" = "Metastatic",
           "07" = "Additional Metastatic",
           "11" = "Solid Tissue Normal")}),   
  group = ifelse(tissue_type_id == "11", "Normal", "Tumor"))

#expression data
file_path <- 'gdc_download_20210914_071640.000645'
counts_files <- list.files(path = file_path, 
                           pattern = "counts.gz", 
                           full.names = T, 
                           recursive = T)
head(counts_files)
library(tools)
all(tools::md5sum(counts_files) %in% metadata$md5sum)

test <- readr::read_tsv(counts_files[1], col_names = F)
head(test)
tail(test)#delete tail 5

# merged counts data 
counts_df <- counts_files %>% 
  lapply(function(x) {
    tmp <- read_tsv(x, col_names = F) %>% 
      purrr::set_names("gene_id", basename(x))
    cat(which(counts_files == x), "of", length(counts_files), "\n")
    return(tmp)
  }) %>%
  reduce(function(x, y) full_join(x, y, by = "gene_id")) %>% 
  dplyr::select(gene_id, metadata$file_name) %>% 
  set_names("gene_id", metadata$TCGA_id_full) %>% 
  dplyr::slice(1:(nrow(.)-5))


head(counts_df)[1:5]



##tpm
gtf_v22 <- read_tsv("gencode.gene.info.v22.tsv")

head(counts_df)[,1:6]

tpm <- counts_df %>% 
  inner_join(gtf_v22, by = "gene_id") %>% 
  mutate(across(.cols = starts_with("TCGA"),
                .fns = ~(((.*10^9)/(sum(.[gene_type == "protein_coding"])*exon_length))*10^6/
                           (sum((.*10^9)/(sum(.[gene_type == "protein_coding"])*exon_length))))
  )
  ) %>% 
  select(colnames(counts_df))

save(counts_df,metadata,file = 'crc_rawdata.Rdata')

save(tpm,file = 'crc_tpm.Rdata')



#ID conversation
# 
# id2symbol <- gtf_v22 %>%
#   dplyr::select(1, 2)
# 
# counts_df <- id2symbol %>%
#   inner_join(counts_df, by = "gene_id")
# 
# head(counts_df) [,1:6] 
# 
# 
# tpm <- id2symbol %>% 
#   inner_join(tpm, by = "gene_id")
# 
# head(tpm) [,1:6]

