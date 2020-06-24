### code to read in previously published chemoreceptor expression data 
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tidytext)
library(cowplot)

# Slotman 2020 - Anopheles quadriannulatus and Anopheles coluzzii
coluzzii_quadriannulatus_M <- read_excel(path='../raw data/13071_2020_4085_MOESM4_ESM.xlsx')
coluzzii_M_F <- read_excel(path='../raw data/13071_2020_4085_MOESM6_ESM.xlsx')
quadriannulatus_M_F <- read_excel(path='../raw data/13071_2020_4085_MOESM7_ESM.xlsx')

Agam.LUT <- read_excel('Agam-LUT.xlsx')
slotman.LUT <- read_excel('slotman-annotations.xlsx')
merge(Agam.LUT,slotman.LUT,by='Accession',all=TRUE) -> all.LUT
write.csv(all.LUT,'Agam.LUT.csv')

Agam.LUT <- read_excel('Agam.LUT.all.xlsx')

my.agam.summary <- read.table('~/Documents/bioinfo/anopheles_chemo/counts/summarize_defaults_multi',header=TRUE)

colnames(my.agam.summary) <- c('Geneid','Chr','Start','End','Strand','Length','quad.palp.2','quad.palp.1','quad.ant.2','quad.ant.1','col.palp.2','col.palp.1','col.ant.2','col.ant.1')

# convert counts to TPM

my.agam.summary <- my.agam.summary %>% mutate(col.palp.2.rpk = 1000 * col.palp.2 / Length)
my.agam.summary <- my.agam.summary %>% mutate(col.palp.1.rpk = 1000 * col.palp.1 / Length)
my.agam.summary <- my.agam.summary %>% mutate(col.ant.2.rpk = 1000 * col.ant.2 / Length)
my.agam.summary <- my.agam.summary %>% mutate(col.ant.1.rpk = 1000 * col.ant.1 / Length) %>%
  select(Geneid,Length,col.palp.2.rpk,col.palp.1.rpk,col.ant.2.rpk,col.ant.1.rpk)

my.agam.summary <- my.agam.summary %>% mutate(col.palp.2.tpm = col.palp.2.rpk / (sum(my.agam.summary$col.palp.2.rpk) / 1000000))
my.agam.summary <- my.agam.summary %>% mutate(col.palp.1.tpm = col.palp.1.rpk / (sum(my.agam.summary$col.palp.1.rpk) / 1000000))
my.agam.summary <- my.agam.summary %>% mutate(col.ant.2.tpm = col.ant.2.rpk / (sum(my.agam.summary$col.ant.2.rpk) / 1000000))
my.agam.summary <- my.agam.summary %>% mutate(col.ant.1.tpm = col.ant.1.rpk / (sum(my.agam.summary$col.ant.1.rpk) / 1000000))

## pull out just chemoreceptors

coluzzii_chemo <- coluzzii_M_F[coluzzii_M_F$Geneid %in% Agam.LUT$Accession,]
coluzzii_chemo <- merge(coluzzii_chemo, Agam.LUT,by=1)

my.coluzzii <- my.agam.summary[my.agam.summary$Geneid %in% Agam.LUT$Accession,]
my.coluzzii <- merge(my.coluzzii,Agam.LUT,by=1)

my.coluzzii.ant <- my.coluzzii %>% select(Geneid,col.ant.1.tpm,col.ant.2.tpm,Gene,Family)
my.coluzzii.palp <- my.coluzzii %>% select(Geneid,col.palp.1.tpm,col.palp.2.tpm,Gene,Family)
#my.quad.ant <- my.coluzzii %>% select(V1,V9,V10,Gene,Family)
#my.quad.palp <- my.coluzzii %>% select(V1, V7, V8,Gene,Family)

#my.coluzzii.ant.mult <- as_tibble(my.coluzzii.mult) %>% select(V1,V13,V14,Gene,Family) %>% rename(Accession=V1,Rep1=V13,Rep2=V14)
#my.coluzzii.palp.mult <- my.coluzzii.mult %>% select(V1,V11,V12,Gene,Family) %>% rename(Accession=V1,Rep1=V11,Rep2=V12)
#my.quad.ant.mult <- my.coluzzii.mult %>% select(V1,V9,V10,Gene,Family) %>% rename(Accession=V1,Rep1=V9,Rep2=V10)
#my.quad.palp.mult <- my.coluzzii.mult %>% select(V1, V7, V8,Gene,Family) %>% rename(Accession=V1,Rep1=V7,Rep2=V8)

agam.plot.an <- my.coluzzii.ant %>% mutate(RepMean = (as.numeric(col.ant.1.tpm) + as.numeric(col.ant.2.tpm)) / 2) %>% mutate(Gene = reorder(Gene,-RepMean)) %>%
  ggplot(data=,aes(x=Gene,y=log10((RepMean)+0.01),colour=Family)) + geom_point() + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') + ggtitle('Anopheles antenna') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
 geom_vline(aes(xintercept=67))#+geom_hline(aes(yintercept=0))

agam.plot.pa <- my.coluzzii.palp %>% mutate(RepMean = (as.numeric(col.palp.1.tpm) + as.numeric(col.palp.2.tpm)) / 2) %>% mutate(Gene = reorder(Gene,-RepMean)) %>%
  ggplot(data=,aes(x=Gene,y=log10((RepMean)+0.01),colour=Family)) + geom_point() + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') + ggtitle('Anopheles palp') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
 geom_vline(aes(xintercept=5))#+geom_hline(aes(yintercept=0))

agam.test.plot.an <- coluzzii_chemo %>% mutate(Geneid = reorder(Geneid, -CFA)) %>%
ggplot(data=,aes(x=Geneid,y=log10(CFA+0.01),colour=Family)) + geom_point() + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)')+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
 geom_vline(aes(xintercept=61))#+geom_hline(aes(yintercept=0))

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

anPlot <- aedes_ntx_L5.tidy.an %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% filter(tissue == 'FeAnSF') %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') + ggtitle('Aedes antenna') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
  geom_vline(aes(xintercept=61))#+geom_hline(aes(yintercept=0))

paPlot <- aedes_ntx_L5.tidy.pa %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% filter(tissue == 'FePaSF') %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') + ggtitle('Aedes palp') + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ ylim(-2,4)+
   geom_vline(aes(xintercept=5))#+geom_hline(aes(yintercept=0))

anPlot_noGR <- aedes_ntx_L5.tidy.an %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  filter(family %in% c("Or","Ir")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
  geom_vline(aes(xintercept=61))#+geom_hline(aes(yintercept=0))

paPlot_noGR <- aedes_ntx_L5.tidy.pa %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  filter(family %in% c("Or","Ir")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
   geom_vline(aes(xintercept=6)) #+geom_hline(aes(yintercept=0))

## summarise counts
aedes_ntx_L5.tidy.an %>% filter(value > 1) %>% count(tissue)
aedes_ntx_L5.tidy.pa %>% filter(value > 1) %>% count(tissue)

aedes_ntx_L5.tidy.an %>%filter(family %in% c("Or","Ir")) %>% 
  filter(value > 1) %>% count(tissue)
aedes_ntx_L5.tidy.pa %>%filter(family %in% c("Or","Ir")) %>%
  filter(value > 1) %>% count(tissue)


## Drosophila
dmel <- read_excel(path='../raw data/Menuz-2014-summary.xlsx') %>% 
  select(Gene_name,ato,cs) %>%
  arrange(desc(cs)) %>% 
  mutate(Gene_name=factor(Gene_name, levels=Gene_name))

dmel.tidy <- gather(dmel,tissue,value,-Gene_name)
dmel.tidy$family <- str_sub(dmel.tidy$Gene_name,1,2)

dmel_plot <- dmel.tidy %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% filter(tissue=='cs') %>%
  filter(family %in% c("Or","Ir","Gr")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (RPKM + 0.01)') + ggtitle('Drosophila antenna') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4)+
  geom_hline(aes(yintercept=0.8)) + geom_vline(aes(xintercept=54))

dmel_plot_noGr <- dmel.tidy %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% 
  filter(family %in% c("Or","Ir")) %>%
  ggplot(data=,aes(x=Gene_name,y=log10(value+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (RPKM + 0.01)') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ ylim(-2,4) + 
  geom_hline(aes(yintercept=0.8)) + geom_vline(aes(xintercept=43))

dmel_plot.tpm <- dmel.tidy %>% mutate(Gene_name = reorder_within(Gene_name, -value, tissue)) %>% filter(tissue=='cs') %>%
  filter(family %in% c("Or","Ir","Gr")) %>% mutate(tpm = 10^6 * (value/704923.6998)) %>% 
  ggplot(data=,aes(x=Gene_name,y=log10(tpm+0.01),colour=family)) + geom_point() + 
  facet_grid(.~tissue, scales="free_x") + 
  scale_x_reordered() + xlab("Chemoreceptor genes") + 
  ylab('log10 (TPM + 0.01)') + ggtitle('Drosophila antenna') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(-2,4) +
  geom_vline(aes(xintercept=51)) + 
  geom_vline(aes(xintercept=58)) +  geom_hline(aes(yintercept=log10(0.617087867364482 + 0.01))) +
  geom_hline(aes(yintercept=log10(8.09378021047491 + 0.01)))
## summarise counts
dmel.tidy %>% filter(value > 10^0.8) %>% count(tissue)
dmel.tidy %>%filter(family %in% c("Or","Ir")) %>%
  filter(value > 10^0.8) %>% count(tissue)




## final plots

aedes_grid <- plot_grid(anPlot,paPlot,nrow=1)
ano_grid <- plot_grid(agam.plot.an,agam.plot.pa,nrow=1)

an_only_plot <- plot_grid(anPlot,agam.plot.an,dmel_plot,nrow=3)
pa_only_plot <- plot_grid(paPlot,agam.plot.pa,nrow=2)

all_together <- plot_grid(aedes_grid,ano_grid,plot_grid(dmel_plot.tpm,dmel_plot.tpm,nrow=1),nrow=3,rel_widths = c(2, 2, 1))
