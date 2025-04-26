library(GEOquery)
library(tidyverse)
library(ggpubr)

# get GSE47460 data
eSet <- getGEO('GSE47460', destdir=".",
               AnnotGPL = F,
               getGPL = F)
b <- eSet[[1]]
raw_exprSet <- as.data.frame(exprs(b))
pd <- pData(b)
colnames(pd)[50] <- "d_type"

# select COPD samples only
pd <- filter(pd,d_type == 'Chronic Obstructive Lung Disease')
raw_exprSet <- raw_exprSet[,c(as.character(pd$geo_accession))]

# expression matrix annotation to gene symbol
gpl <- getGEO('GPL14550', destdir=".")
id2symbol <- Table(gpl)[,c(1,7)]
raw_exprSet$ID <- rownames(raw_exprSet)
raw_exprSet <- raw_exprSet %>%
  inner_join(id2symbol,by = "ID") %>%
  dplyr::select(-ID) %>% 
  dplyr::select(GENE_SYMBOL,1:(length(raw_exprSet)-1)) %>% 
  mutate(rowMean =rowMeans(.[grep("GSM", names(.))])) %>%
  arrange(desc(rowMean)) %>% 
  distinct(GENE_SYMBOL,.keep_all = T) %>% 
  dplyr::select(-rowMean) 
rownames(raw_exprSet) <- raw_exprSet[,1]
exprSet <- raw_exprSet[,-1]

plot_dat <- exprSet[c('CD207','CD8A'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('geo_accession') %>% 
  inner_join(pd,by = 'geo_accession') %>% 
  select(2,3,45:50) %>% 
  rename(CD8 = 2,emph = 3, DLCO = 4, FEV1_post = 5,FEV1_pre = 6, FVC_post = 7, FVC_pre = 8) %>% 
  mutate(across(c(3:8),as.numeric))

# correlation plot
plot_dat %>% 
  ggscatter(x = 'CD8', y = 'emph', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "%Emphysema")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))

plot_dat %>% 
  ggscatter(x = 'CD8', y = 'DLCO', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "DLCO%pred")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))


plot_dat %>% 
  ggscatter(x = 'CD8', y = 'FEV1_post', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "FEV1%pred (post-bd)")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))

plot_dat %>% 
  ggscatter(x = 'CD8', y = 'FEV1_pre', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "FEV1%pred (pre-bd)")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))


plot_dat %>% 
  ggscatter(x = 'CD8', y = 'FVC_post', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "FVC%pred (post-bd)")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))


plot_dat %>% 
  ggscatter(x = 'CD8', y = 'FVC_pre', 
            add = "reg.line", conf.int = TRUE, cor.coef = TRUE,
            size = 0.5,
            add.params = list(color = "#3d6cb2",
                              fill = "#737373",size = 1),
            cor.coeff.args = list(method = "spearman", 
                                  color = "#3d6cb2",
                                  size =3,
                                  cor.coef.name = 'r'),
            xlab = "CD8A expression", ylab = "FVC%pred (pre-bd)")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))