### code to read in previously published chemoreceptor expression data 
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tidytext)

# Slotman 2020 - Anopheles quadriannulatus and Anopheles coluzzii
coluzzii_quadriannulatus_M <- read_excel(path='../raw data/13071_2020_4085_MOESM4_ESM.xlsx')
coluzzii_M_F <- read_excel(path='../raw data/13071_2020_4085_MOESM6_ESM.xlsx')
quadriannulatus_M_F <- read_excel(path='../raw data/13071_2020_4085_MOESM7_ESM.xlsx')


# Matthews 2018 - Aedes aegypti neurotranscriptome with AaegL5

download.file(url='https://github.com/VosshallLab/AGWG-AaegL5/raw/master/AGWG%20AaegL5%20Chemoreceptor%20TPM.xlsx',
              destfile = 'Aedes_ntx_raw.xlsx')

aedes_ntx_L5_an <- read_excel(path='Aedes_ntx_raw.xlsx') %>% 
  select(Gene_name,FeAnBF,FeAnO,FeAnSF,MaAn) %>%
  arrange(desc(FeAnSF)) %>% 
  mutate(Gene_name=factor(Gene_name, levels=Gene_name))

aedes_ntx_L5_pa <- read_excel(path='Aedes_ntx_raw.xlsx') %>% 
  select(Gene_name,FePaSF,MaRs,FeRsBF,FeRsSF) %>%
  arrange(desc(FePaSF)) %>% 
  mutate(Gene_name=factor(Gene_name, levels=Gene_name))

aedes_ntx_L5.tidy.an <- gather(aedes_ntx_L5_an,tissue,value,-Gene_name) 
aedes_ntx_L5.tidy.an$family <- str_sub(aedes_ntx_L5.tidy.an$Gene_name,1,2)

aedes_ntx_L5.tidy.pa <- gather(aedes_ntx_L5_pa,tissue,value,-Gene_name) 
aedes_ntx_L5.tidy.pa$family <- str_sub(aedes_ntx_L5.tidy.pa$Gene_name,1,2)

anPlot <- aedes_ntx_L5.tidy.an %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=61))

paPlot <- aedes_ntx_L5.tidy.pa %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=5))

anPlot_noGR <- aedes_ntx_L5.tidy.an %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  filter(family %in% c("Or","Ir")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=61))

paPlot_noGR <- aedes_ntx_L5.tidy.pa %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  filter(family %in% c("Or","Ir")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(aes(yintercept=0)) + geom_vline(aes(xintercept=5))


## summarise counts
aedes_ntx_L5.tidy.an %>% filter(value > 1) %>% count(tissue)
aedes_ntx_L5.tidy.pa %>% filter(value > 1) %>% count(tissue)

aedes_ntx_L5.tidy.an %>%filter(family %in% c("Or","Ir")) %>% 
  filter(value > 1) %>% count(tissue)
aedes_ntx_L5.tidy.pa %>%filter(family %in% c("Or","Ir")) %>%
  filter(value > 1) %>% count(tissue)



