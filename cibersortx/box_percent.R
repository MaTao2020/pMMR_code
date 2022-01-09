library(tidyverse)
library(ggsci)
library(colorspace)
setwd("../rawdata/")
clin <- read_csv("clin_ciber.csv")

clin <- column_to_rownames(clin, "X1")

df <- read_csv("CIBERSORTx_Job22_Results.csv")

df <- df %>%
  dplyr::select("Mixture":"hT18_CD8-LAYN") %>%
  column_to_rownames("Mixture")

percent <- colSums(df != 0) / 398
percent <- as.data.frame(percent)

percent %>%
  rownames_to_column("Celltype") %>%
  mutate(data = percent * 100) %>%
  mutate(Celltype = fct_reorder(Celltype, data)) %>%
  ggplot(aes(Celltype, data, fill = Celltype)) +
  geom_col() +
  theme_bw() +
  theme(legend.position = "none") +
  coord_flip() +
  geom_hline(yintercept = 50, linetype = 2) +
  scale_fill_discrete_qualitative(palette = "Warm")+
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.1)))
  

ggsave(filename = '../results/box_percent_plot.pdf',width = 6,height = 8)
