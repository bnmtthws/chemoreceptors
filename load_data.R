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

aedes_ntx_L5 <- read_excel(path='Aedes_ntx_raw.xlsx') %>% 
  select(Gene_name,FeAnBF,FeAnO,FeAnSF,FePaSF,MaAn) %>%
  arrange(desc(FeAnSF)) %>% 
  mutate(Gene_name=factor(Gene_name, levels=Gene_name))

aedes_ntx_L5.tidy <- gather(aedes_ntx_L5,tissue,value,-Gene_name) 
aedes_ntx_L5.tidy 
aedes_ntx_L5.tidy$family <- str_sub(aedes_ntx_L5.tidy$Gene_name,1,2)

aedes_ntx_L5.tidy %>% mutate(Gene_name = reorder_within(Gene_name, value, tissue)) %>% 
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_hline(aes(yintercept=0)) + geom_hline(aes(yintercept=))
