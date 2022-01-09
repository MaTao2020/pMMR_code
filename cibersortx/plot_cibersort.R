#plot cibersortx OLR
library(tidyverse)
setwd('../results/')
df <- read_csv('cibersort_OLR_satge.csv')
df <- df[,-1]

df <- df %>% 
  mutate(logP =-log10(p_value))  %>% 
  mutate(group = case_when(value > 0 & p_value < 0.05 ~ 'up',
                           value < 0 & p_value < 0.05 ~ 'down',
                           TRUE ~ 'other')) %>% 
  arrange(p_value) %>% 
  mutate(celltype=fct_reorder(celltype,p_value,.desc = T))


library(ggsci)
ggplot(df,aes(x=value,y=celltype))+geom_vline(xintercept=0,linetype=2,color="gray44")+
  geom_errorbar(aes(xmin=lowerCI,xmax=upperCI),width=0.2)+geom_point(aes(color=group,size=logP))+
  theme_bw()+
  scale_size_continuous(name = "-log pvalue",range=c(1,5)) +
  scale_color_manual(values=c("#0072b5","black","#bc3c29"))
  

df %>% mutate(celltype= fct_reorder(celltype,p_value))
 
