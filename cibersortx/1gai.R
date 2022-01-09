#
library(tidyverse)

setwd('../rawdata/')
clin <- read_csv('clin_ciber.csv')

clin <- column_to_rownames(clin,'...1')

df <- read_csv('CIBERSORTx_Job22_Results.csv')

df <- df %>% 
  dplyr::select('Mixture':'hT18_CD8-LAYN') %>% 
  column_to_rownames('Mixture')

df <- df[,((colSums(df !=0)/398) >= 0.5)]

colnames(df) <- str_replace(colnames(df),'-','_')
data <- cbind(df,clin)
#
library(MASS)
df1 <- data %>%
  dplyr::mutate(
    across(colnames(df) | gender | age |  pathologic_M | pathologic_N| pathologic_T, as.factor),
    across(colnames(df), ~ fct_inseq(., ordered = TRUE))
  )

mydf<-matrix(nrow=(ncol(df1)-8),ncol=5)

for (i in 1:(ncol(df1)-8)){
  c<-colnames(df1)[i]
  d <- paste('gender','age','pathologic_M', 'pathologic_N','pathologic_T','stage',sep = '+')
  fml <- as.formula(paste(c ,'~',d))
  
  mod_mass=polr(fml,data = df1,Hess = TRUE,
                method = c("logistic"))
  ctable <- coef(summary(mod_mass))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  ctable <- cbind(ctable, "p value" = p)
  mydf[i,1]<-c
  mydf[i,2:5]<-ctable["stage",]
}

mydf <- mydf %>% as.data.frame %>% 
  set_names('celltype','value','std_error','t_value','p_value')

mydf <- mydf %>% 
  mutate(across(colnames(mydf)&!celltype,as.numeric)) %>% 
  mutate(lowerCI = value - 1.96*std_error,
         upperCI = value + 1.96*std_error) %>% 
  dplyr::select(-t_value)

#write.csv(mydf,file = '../results/cibersort_OLR_satge.csv')

df <- mydf %>% 
  mutate(logP =-log10(p_value))  %>% 
  mutate(group = case_when(value > 0 & p_value < 0.05 ~ 'up',
                           value < 0 & p_value < 0.05 ~ 'down',
                           TRUE ~ 'other')) %>% 
  arrange(p_value) %>% 
  mutate(celltype=fct_reorder(celltype,p_value,.desc = T))

library(ggsci)

ggplot(df,aes(x=value,y=celltype))+geom_vline(xintercept=0,linetype=2,color="gray44")+
  geom_errorbar(aes(xmin=lowerCI,xmax=upperCI),width=0.5)+geom_point(aes(color=group,size=logP))+
  theme_bw()+
  scale_color_manual(values=c("#00AFBB","black","#FC4E07")) +
  scale_size_continuous(name = "-log pvalue",range=c(1,5))+
  scale_color_jco()

df %>% mutate(celltype= fct_reorder(celltype,p_value))



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

ggsave('../results/cibersotx_2.pdf',width = 8,height = 8)




