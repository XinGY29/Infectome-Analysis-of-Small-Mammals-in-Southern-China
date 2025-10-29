####figure 3
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
##install packages
install.packages('ggthemes')
install.packages('tidyr')
install.packages('ggplot2')
install.packages('pairwiseAdonis')
install.packages('ggpubr')
install.packages('ggvenn')
install.packages('patchwork')
install.packages('RColorBrewer')
install.packages('dplyr')
install.packages('vegan')
install.packages('grid')
install.packages('gridExtra')
###import library
library(ggthemes)
library(tidyr)
library(ggplot2)
library(pairwiseAdonis)
library(ggpubr)
library(ggvenn)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(cowplot)
library(grid)
library(gridExtra) 
##Work Path
setwd(data_path)

###figure 3a
distribution.num <- read.csv(paste0(data_path,'/fig3_pathogens_distribution.csv',sep=''),check.names = F,row.names = 1)
fig3a <- ggplot(distribution.num, aes(x=Tissue,y = Num,fill=Type)) +
  geom_bar(stat = "identity",
           width = 0.8, size = 0.25,
           alpha = 1)+ 
  labs(x="",y="Num of pathogens")+theme_few()+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))+
  scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
                                      c("RNA Virus","DNA Virus","Bacteria","Eukaryote")))
###figure 3b
tissue.table.long.taxo <- read.csv(paste0(data_path,'/fig3_abundance.csv',sep=''),check.names = F,row.names = 1)
Lung.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung')$Pathogens)
Spleen.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen')$Pathogens)
Gut.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut')$Pathogens)
Tissue <- list(Lung=Lung.Pathogens,Spleen=Spleen.Pathogens,Gut=Gut.Pathogens)
fig3b <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),set_name_size = 3,text_size = 2,show_percentage=FALSE)

###figure 3c
richness.data <- read.csv(paste0(data_path,'/fig3_richness_violin.csv',sep=''),check.names = F,row.names = 1)
richness.data$Type <- factor(richness.data$Type,levels = c("RNA virus","DNA virus","Bacteria","Eukaryote"))
fig3c1<-ggplot(richness.data,aes(Tissue,Richness,color=Tissue,fill=Tissue)) + 
  geom_violin(trim=FALSE)+
  stat_compare_means(comparisons=list(c("Lung","Spleen"),c("Lung","Gut"),
                                      c("Spleen","Gut")),size=2,
                     label = "p.signif",method = 'wilcox.test')+
  facet_grid(. ~ Type) +
  ylim(0,10)+
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),strip.text = element_text(size= 6),
        axis.text = element_text(size = 5))+
  scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
  scale_fill_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))

tissue.table.long.taxo <- subset(tissue.table.long.taxo,RPM > 0)
tissue.table.long.taxo$Type <- factor(tissue.table.long.taxo$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
fig3c2<-ggplot(tissue.table.long.taxo,aes(Tissue,log10(RPM+1),color=Tissue,fill=Tissue)) + 
  geom_violin(trim=FALSE)+
  stat_compare_means(comparisons=list(c("Lung","Spleen"),c("Lung","Gut"),
                                      c("Spleen","Gut")),size=2,
                     label = "p.signif",method = 'wilcox.test')+
  facet_grid(. ~ Type) +
  ylim(0,8)+
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),strip.text = element_text(size= 6),
        axis.text = element_text(size = 5))+
  scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
  scale_fill_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))

f3c <- fig3c1/fig3c2

###figure 3d
Tissue <- read.csv(paste0(data_path,'/fig3_Tissue_Pathogens_Composition.csv',sep=''),check.names = F,row.names = 1)

fig3d <- ggplot(data = Tissue,aes(x=Family,y=RPM,fill=Tissue))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
  theme_few()+
  geom_point(aes(x=Family,y=0,color=Type),size=1,shape=15)+
  scale_color_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),c("RNA Virus","DNA Virus","Bac","Euk")))+
  labs(x="",y='Relative RPM percentage(%)')+
  scale_x_discrete(limits=unique(factor(Tissue$Family)))+
  theme(legend.key.size = unit(3, 'mm'),
        legend.position = 'top',
        text=element_text(size= 6),
        axis.text.x=element_text(angle=90,hjust = 1))+coord_flip()

###figure 3e
tissue.table.long.indi <- read.csv(paste0(data_path,'/fig3_12viruses.csv',sep=''),check.names = F,row.names = 1)
selected.pathogens <- c('Rat hepatitis E virus','Wenzhou virus',
                        'Betacoronavirus 1',
                        'Orthohantavirus seoulense',
                        'Norwalk virus','Bartonella kosoyi',
                        'Klebsiella variicola',
                        'Angiostrongylus cantonensis',
                        'Cryptosporidium ubiquitum',
                        'Nippostrongylus brasiliensis',
                        'Pneumocystis carinii','Trypanosoma lewis')
###selected_Table
selected.table <- subset(tissue.table.long.indi,Pathogens%in%selected.pathogens)
selected.table$Pathogens <- factor(selected.table$Pathogens,levels = selected.pathogens)

fig3e <- ggplot(subset(selected.table,RPM>0),aes(Tissue,log10(RPM+1),color=Tissue))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
  theme_few()+
  scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
  scale_x_discrete(limits = c('Gut','Lung','Spleen'))+
  facet_wrap(vars(Pathogens),nrow =2,scales = "free")+
  theme(legend.position = 'none',axis.title=element_text(size= 6),
        strip.text = element_text(size= 6),
        plot.title = element_text(hjust = 0.5 ),
        text=element_text(size= 5),
        axis.text = element_text(size = 5),
        axis.text.x=element_text(angle=60,hjust = 1))
fig3e

###save figs

ggsave(filename = paste0(out_path,'/fig 3a.pdf',sep=''), plot = fig3a, width = 80, height = 60,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3b.pdf',sep=''), plot = fig3b, width = 60, height = 60,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3c.pdf',sep=''), plot = f3c, width = 140, height = 85,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3d.pdf',sep=''), plot = fig3d, width = 60, height = 120,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3e.pdf',sep=''), plot = fig3e, width = 180, height = 80,units = 'mm')

