library(tidyverse)
library(ggpubr)
library(rstatix)

# GSE18385 ----
rm(list = ls())
load('C:\\Files\\geo_traning\\Xuan\\GSE18385\\raw_exp_data.rda')
pd <- read.csv('C:\\Files\\geo_traning\\Xuan\\GSE18385\\clinical_data.csv')
pd <- pd[,c(2,10:12)]
names(pd) <- c('sample','group','smoke_pack_years','airway')
pd <- pd %>% mutate(airway = case_when(
  airway == 'large airways' ~ 'Large aiway',
  TRUE ~ 'Small airway'
))

dat_tmp <- as.data.frame(t(raw_exprSet['CD207',]))
all(rownames(dat_tmp) == pd$sample)
dat_tmp <- cbind(pd,dat_tmp)


# large airway CD207
boxplot_data <- filter(dat_tmp,airway == 'Large aiway') %>% 
  select(2,5) %>% 
  rename(value = 2) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 11)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_fill_manual(values = c("#b3c5e3","#F5DAB7")) +
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CD207')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE18385_CD207_large_aiway.pdf',width = 0.8,height = 1.2)




# small airway CD207
boxplot_data <- filter(dat_tmp,airway == 'Small airway') %>% 
  select(2,5) %>% 
  rename(value = 2) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 11)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CD207')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE18385_CD207_small_aiway.pdf',width = 0.8,height = 1.2)





rm(list = ls())
load('C:\\Files\\geo_traning\\Xuan\\GSE18385\\raw_exp_data.rda')
pd <- read.csv('C:\\Files\\geo_traning\\Xuan\\GSE18385\\clinical_data.csv')
pd <- pd[,c(2,10:12)]
names(pd) <- c('sample','group','smoke_pack_years','airway')
pd <- pd %>% mutate(airway = case_when(
  airway == 'large airways' ~ 'Large aiway',
  TRUE ~ 'Small airway'
))

dat_tmp <- as.data.frame(t(raw_exprSet['CSF2',]))
all(rownames(dat_tmp) == pd$sample)
dat_tmp <- cbind(pd,dat_tmp)


# large airway CD207
boxplot_data <- filter(dat_tmp,airway == 'Large aiway') %>% 
  select(2,5) %>% 
  rename(value = 2) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 10)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_fill_manual(values = c("#b3c5e3","#F5DAB7")) +
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CSF2')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE18385_CSF2_large_aiway.pdf',width = 0.8,height = 1.2)




# small airway CD207
boxplot_data <- filter(dat_tmp,airway == 'Small airway') %>% 
  select(2,5) %>% 
  rename(value = 2) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 10)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CSF2')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE18385_CSF2_small_aiway.pdf',width = 0.8,height = 1.2)




# GSE5058 ----
rm(list = ls())
load('../../GSE5058/GSE5058_eSet.rda')
b=eSet[[1]]
raw_exprSet=as.data.frame(exprs(b))
pd <- read.csv('../../GSE5058/clinical_info.csv') %>%
  rename(Age = 3, Sex =4,Ethnic = 5) %>% 
  mutate(
    Age = as.numeric(str_extract(Age, "\\d+")), 
    Sex = str_remove(Sex, ".*:\\s"),    
    Ethnic = str_replace(Ethnic, "(?i)^ethnic group\\s*:\\s*", "")
  ) %>% 
  mutate(Group = if_else(str_detect(characteristics_ch1.3, "early-COPD"), "early-COPD", Group)) %>% 
  dplyr::select(geo_accession,Age,Sex,Ethnic,Group) %>% 
  filter(Group != 'early-COPD') %>%
  mutate(dataset = 'GSE5058')

load('../../GSE5058/GSE8545_eSet.rda')
b2=eSet2[[1]]
raw_exprSet2=as.data.frame(exprs(b2))
pd2 <- read.csv('../../GSE5058/clinical_info_GSE8545.csv') %>% 
  rename(Age = 3, Sex =4,Ethnic = 5) %>% 
  mutate(
    Age = as.numeric(str_extract(Age, "\\d+")), 
    Sex = str_remove(Sex, ".*:\\s"),
    Ethnic = str_replace(Ethnic, "(?i)^ethnic group\\s*:\\s*", "")
  ) %>% 
  mutate(Group = if_else(str_detect(characteristics_ch1.3, "early-COPD"), "early-COPD", Group)) %>% 
  dplyr::select(geo_accession,Age,Sex,Ethnic,Group) %>% 
  mutate(dataset = 'GSE8545')

common_cols <- intersect(colnames(raw_exprSet),colnames(raw_exprSet2))
raw_exprSet2 <- raw_exprSet2[,-which(colnames(raw_exprSet2) %in% common_cols)]
exp = cbind(raw_exprSet,raw_exprSet2)

pd2 <- pd2 %>% 
  filter(!geo_accession %in% common_cols)
dat <- pd %>% bind_rows(pd2)

load('C:/Files/geo_traning/Xuan/GSE5058/exprSet.Rda')
tmp <-  exprSet[c('CD207','CSF2'),] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('geo_accession') %>% 
  inner_join(dat %>% select(geo_accession,Group), by = 'geo_accession') %>% 
  filter(Group != 'COPD')


# CD207
boxplot_data <- tmp %>% 
  select(group = Group,value = CD207) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 9)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CD207')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE5058_CD207.pdf',width = 0.8,height = 1.2)




# CSF2
boxplot_data <- tmp %>% 
  select(group = Group,value = CSF2) %>% 
  mutate(group = str_to_sentence(group)) %>% 
  mutate(group = fct_relevel(group,c('Non-smoker','Smoker')))

sem <- function(x, na.rm = FALSE) {
  out <-sd(x, na.rm = na.rm)/sqrt(length(!is.na(x)))
  return(out)
}
boxplot_dat <- boxplot_data %>%  
  group_by(group) %>% 
  summarise(value_mean=mean(value,na.rm = T),sd=sem(value,na.rm = T))


stat_dat <- boxplot_data %>% 
  wilcox_test(value~group) %>%  
  add_significance('p') %>% 
  add_xy_position() %>% 
  mutate(y.position = 8)


set.seed(123)
p <- ggplot(boxplot_data,aes(group,value,color = group))+
  geom_boxplot(lwd = 0.2,outlier.colour = NA)+ 
  geom_jitter(aes(fill=group,color = group), width = 0.3,size=0.2,shape = 16)+
  scale_color_manual(values = c("#576fa0","#FFA448"))+
  xlab(NULL)+ylab('Expression levels')+
  labs(title = 'CSF2')+
  theme_classic() + 
  theme(legend.position = "none",
        panel.grid = element_line(size = 0.1, linetype = 2),
        plot.title = element_text(size=6,hjust=0.5,margin = margin(b=1,t=2,unit = "pt")),
        # legend.margin=margin(l=-10,unit = "pt"),
        plot.margin = margin(t=0,r=0,b=0,l=0,unit = "pt"), # 默认都是5.5，即base_size/2，默认base_size是11
        axis.line = element_line(linewidth = 0.2, color = "black"),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(margin = margin(r=2,unit = "pt")),
        axis.title.x = element_text(margin = margin(t=2,unit = "pt")),
        axis.ticks = element_line(linewidth = 0.2, color = "black"),
        axis.ticks.length = unit(.06, "cm"),
        axis.text.y = element_text(size=6, colour = 'black',margin = margin(r=1,unit = "pt")),
        axis.text.x = element_text(size=5, colour = 'black',margin = margin(t=0,unit = "pt"),angle = 20, hjust = 0.8))+
  # guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  stat_pvalue_manual(
    stat_dat,label = 'p.signif',
    size = 1.8,
    tip.length = 0,
    bracket.size = 0.2,
    bracket.shorten = 0.1,
    bracket.nudge.y = 0
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.09)))
p
ggsave(p,filename = 'tmp/GSE5058_CSF2.pdf',width = 0.8,height = 1.2)

