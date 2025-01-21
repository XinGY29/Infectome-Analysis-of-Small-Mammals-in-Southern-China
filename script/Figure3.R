####figure3
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','raw_data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")
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
###Tissue
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
##read data
Library_Tissue_RPM <- read.csv('Library_Pathogen_Mean.csv',check.names = F)[,-1]
Library_Tissue_RPM$Type <- paste(Library_Tissue_RPM$Collected_Region,Library_Tissue_RPM$Collected_Time,
                                 Library_Tissue_RPM$Tissue,sep = '-')
Ann_Type <- unique(Library_Tissue_RPM[,c("Type","Tissue")])
Library_Tissue_RPM_long <- gather(subset(Library_Tissue_RPM,select = -c(
  Library,Hos_Taxo,Tissue,Collected_Time,Collected_Region)),key = "Taxo",value = RPM, -"Type")
Library_RPM_Tissue <- aggregate(Library_Tissue_RPM_long$RPM,
                                list(Library_Tissue_RPM_long$Type,Library_Tissue_RPM_long$Taxo),mean)
colnames(Library_RPM_Tissue) <- c("Type","Taxo","RPM")
#Anno
Ann_Table <- read.csv('Virus_Name_Anno.csv')
Ann_row <- unique(Ann_Table[,c("classify",'Level')])
Family <- unique(Ann_row$classify)

###colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
            "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
color20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

Family_color <- color37[1:length(Family)]
filll_color=setNames(Family_color,Family)

### bar graph
Library_RPM_Tissue_1 <- merge(Library_RPM_Tissue,Ann_Type)
Ann_Table <- read.csv('F:/Mouse_Result_24/RPM_Information/rmCrossLib/Virus_Name_Anno.csv')
Ann_Table <- subset(Ann_Table,Taxo != 'Desulfovibrio piger')
Library_RPM_Tissue_1 <- merge(Library_RPM_Tissue_1,Ann_Table[,c('Taxo','classify','Level')])
Tissue <- aggregate(Library_RPM_Tissue_1$RPM,
                    list(Library_RPM_Tissue_1$classify,Library_RPM_Tissue_1$Tissue,
                         Library_RPM_Tissue_1$Level),mean)
colnames(Tissue) <- c("Family","Tissue","Type","RPM")
Tissue$Type <- factor(Tissue$Type,levels = c("RNA Virus","DNA Virus","Bac","Euk"))
Tissue <- Tissue %>%
  filter(!(Family %in% c('Picobirnaviridae','Anelloviridae')))

Tissue.sum <-  aggregate(Tissue$RPM,list(Tissue$Family),sum)
colnames(Tissue.sum) <- c('Family','RPM.sum')
Tissue <- merge(Tissue,Tissue.sum,all.x = TRUE)
Tissue$rate <- Tissue$RPM/Tissue$RPM.sum
Tissue.wide <- spread(Tissue[,c('Family','Tissue','rate')],key = 'Tissue',value = 'rate',fill = 0)
Tissue.wide <- merge(Tissue.wide,unique(Tissue[,c('Family','Type')]))
Tissue.wide <- Tissue.wide[order(Tissue.wide$Lung,decreasing = TRUE),]
Tissue.wide <- Tissue.wide[order(Tissue.wide$Spleen,decreasing = TRUE),]
Tissue.wide <- Tissue.wide[order(Tissue.wide$Gut,decreasing = TRUE),]
Tissue.wide <- Tissue.wide[order(Tissue.wide$Type),]
Tissue$Family <- factor(Tissue$Family,levels = Tissue.wide$Family)
Tissue <- Tissue[order(Tissue$Family),]

Relative.plot <- ggplot(data = Tissue,aes(x=Family,y=RPM,fill=Tissue))+
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
Relative.plot
# ###病原在组织之间分布的多样性
# Library_RPM_Tissue_1 <- Library_RPM_Tissue_1 %>%
#   filter(!(classify %in% c('Picobirnaviridae','Anelloviridae')))
# Tissue_RPM <- aggregate(Library_RPM_Tissue_1$RPM,
#                         list(Library_RPM_Tissue_1$Tissue,Library_RPM_Tissue_1$Taxo,
#                              Library_RPM_Tissue_1$Level),sum)
# colnames(Tissue_RPM) <- c("Tissue","Taxo","Level","RPM")
# Tissue_RPM_wide <- spread(Tissue_RPM[,c("Tissue","Taxo","RPM")],key=Tissue,value=RPM,fill=0)
# 
# rownames(Tissue_RPM_wide) <- Tissue_RPM_wide$Taxo
# Tissue_RPM_wide <- Tissue_RPM_wide[,-1]
# Ann_Table$Level <- factor(Ann_Table$Level,levels = c("RNA Virus","DNA Virus","Bac","Euk"))
# Ann_Table <- Ann_Table[order(Ann_Table$Level),]
# rownames(Ann_Table) <- Ann_Table$Taxo
# Tissue_RPM_wide <- Tissue_RPM_wide[intersect(Ann_Table$Taxo,rownames(Tissue_RPM_wide)),]
# # ##PCA
# df.dist = vegdist(Tissue_RPM_wide,method='bray')
# 
# df.dist[is.na(df.dist)] <- 0
# pca <- prcomp(log(Tissue_RPM_wide+1))
# pca$x
# plot_data <- data.frame(pca$x)[1:2]
# plot_data$Taxo <- rownames(plot_data)
# names(plot_data)[1:2] <- c('PCA1', 'PCA2')
# plot_data <- merge(plot_data, Ann_Table, by = 'Taxo', all.x = TRUE)
# focus=c('Trypanosoma lewis','Angiostrongylus cantonensis','Bartonella kosoyi',
#         'Porcine bocavirus','Rat minute virus 2a','Guangdong rodent torque teno virus 1',
#         'Guangdong rodent dependoparvovirus 4','Guangdong rodent dependoparvovirus 1',
#         'Guangdong rodent anellovirus 3','Rodent pegivirus','Wenzhou Apodemus agrarius hepacivirus 1',
#         'Wnezhou virus','Beilong virus','Guangdong Rattus norvegicus arterivirus virus',
#         'Guangdong Rattus norvegicus arterivirus virus 2','Guangdong rodent arterivirus 1',
#         'Guangdong rodent arterivirus 2','Norway rat pestivirus','Orthohantavirus seoulense',
#         'Parabovirus A','Rat hepatitis E virus','Rodent astrovirus 3')
# labels <- subset(plot_data,Taxo %in% focus)
# 
# PCA_Tissue <- ggplot() +
#   geom_point(data = plot_data, aes(x=PCA1, y=PCA2,color=Level),alpha=.7, size=1) +
#   # labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#   #      y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
#   scale_colour_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
#                                         c("RNA Virus","DNA Virus","Bac","Euk")))+
#   # stat_ellipse(aes(fill = Level,color=Level), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
#   geom_text(data=labels,aes(x=PCA1,y=PCA2,label=Taxo),size=5,size.unit='pt',vjust=1) +
#   theme_few()+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'))+
#   coord_fixed(1.5)
# PCA_Tissue

###病毒相关的丰度信息
tissue.table <- read.csv('Individual_Mammal-related_Pathogen_Library.csv',check.names = F)
tissue.table <- subset(tissue.table,select = -c(`Desulfovibrio piger`))
colnames(tissue.table)[1] <- 'Library'
Library_infor <- read.csv('Library.csv')
Pathogen.info <- read.csv('Pathongen_Taxon_Table.txt',sep='\t')
Mammal.related <- subset(Pathogen.info,Host_Type == 'Mammal-related')
Mammal.related <- unique(Mammal.related[,c('species','family','Type')])
Mammal.related$Type <- ifelse(Mammal.related$Type=='Bac','Bacteria',
                              ifelse(Mammal.related$Type=='Euk','Eukaryote',
                                     ifelse(Mammal.related$Type=='RNA_Virus','RNA Virus',
                                            ifelse(Mammal.related$Type=='DNA_Virus','DNA Virus',Mammal.related$Type))))
Mammal.related <- Mammal.related %>%
  filter(!(family %in% c('Picobirnaviridae','Anelloviridae')))

tissue.table.long <- gather(tissue.table,key = "Pathogens",value = RPM, -"Library")
tissue.table.long.taxo <- merge(tissue.table.long,Library_infor[,c('Library','Tissue')])
Library.anno <- unique(tissue.table.long.taxo[,c('Library','Tissue')])
colnames(Mammal.related) <- c('Pathogens','Family','Type')
tissue.table.long.taxo <- merge(tissue.table.long.taxo,Mammal.related)
RNA.Virus.long <- subset(tissue.table.long.taxo,Type == 'RNA Virus')
DNA.Virus.long <- subset(tissue.table.long.taxo,Type == 'DNA Virus')
Bac.long <- subset(tissue.table.long.taxo,Type == 'Bacteria')
Euk.long <- subset(tissue.table.long.taxo,Type == 'Eukaryote')

RNA.virus.wide <- spread(RNA.Virus.long[,c("Library","Pathogens","RPM")],key = Pathogens,value = RPM,fill=0)
rownames(RNA.virus.wide) <- RNA.virus.wide$Library
RNA.virus.wide <- RNA.virus.wide[,-1]
Shannon <- diversity(t(RNA.virus.wide), index = "shannon", MARGIN = 2, base = exp(1))
Richness <- specnumber(t(RNA.virus.wide), MARGIN = 2)#spe.rich =sobs
RNAvirus.richness <- as.data.frame(cbind(Shannon, Richness))
RNAvirus.richness$Library <- rownames(RNAvirus.richness)
RNAvirus.richness$Type <- 'RNA virus'

DNA.virus.wide <- spread(DNA.Virus.long[,c("Library","Pathogens","RPM")],key = Pathogens,value = RPM,fill=0)
rownames(DNA.virus.wide) <- DNA.virus.wide$Library
DNA.virus.wide <- DNA.virus.wide[,-1]
Shannon <- diversity(t(DNA.virus.wide), index = "shannon", MARGIN = 2, base = exp(1))
Richness <- specnumber(t(DNA.virus.wide), MARGIN = 2)#spe.rich =sobs
DNAvirus.richness <- as.data.frame(cbind(Shannon, Richness))
DNAvirus.richness$Library <- rownames(DNAvirus.richness)
DNAvirus.richness$Type <- 'DNA virus'

Euk.wide <- spread(Euk.long[,c("Library","Pathogens","RPM")],key = Pathogens,value = RPM,fill=0)
rownames(Euk.wide) <- Euk.wide$Library
Euk.wide <- Euk.wide[,-1]
Shannon <- diversity(t(Euk.wide), index = "shannon", MARGIN = 2, base = exp(1))
Richness <- specnumber(t(Euk.wide), MARGIN = 2)#spe.rich =sobs
Euk.richness <- as.data.frame(cbind(Shannon, Richness))
Euk.richness$Library <- rownames(Euk.richness)
Euk.richness$Type <- 'Eukaryote'

Bac.wide <- spread(Bac.long[,c("Library","Pathogens","RPM")],key = Pathogens,value = RPM,fill=0)
rownames(Bac.wide) <- Bac.wide$Library
Bac.wide <- Bac.wide[,-1]
Shannon <- diversity(t(Bac.wide), index = "shannon", MARGIN = 2, base = exp(1))
Richness <- specnumber(t(Bac.wide), MARGIN = 2)#spe.rich =sobs
Bac.richness <- as.data.frame(cbind(Shannon, Richness))
Bac.richness$Library <- rownames(Bac.richness)
Bac.richness$Type <- 'Bacteria'

richness.data <- rbind(Bac.richness,RNAvirus.richness,DNAvirus.richness,Euk.richness)
richness.data <- merge(richness.data,Library.anno)
richness.data$Type <- factor(richness.data$Type,levels = c("RNA virus","DNA virus","Bacteria","Eukaryote"))
  
Richness.plot <- ggplot(richness.data,aes(Tissue,Richness,color=Tissue))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness')+ geom_jitter(size=0.5)+
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
Richness.plot
richness.data$Type <- factor(richness.data$Type,levels = )
Shannon.plot <- ggplot(richness.data,aes(Tissue,Shannon,color=Tissue))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='Shannon')+ geom_jitter(size=0.5)+
  stat_compare_means(comparisons=list(c("Lung","Spleen"),c("Lung","Gut"),
                                      c("Spleen","Gut")),size=2.5,
                     label = "p.signif",method = 'wilcox.test')+
  facet_grid(. ~ Type) +
  ylim(0,2)+
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5))+
  scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
  scale_fill_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))
Shannon.plot
# ###Abundance Plot sum
# tissue.table.long.taxo$Type <- factor(tissue.table.long.taxo$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
# Abundance.plot.sum <- ggplot(subset(tissue.table.long.taxo,RPM > 0),aes(Tissue,log10(RPM+1),color=Tissue))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='Abundance')+ geom_jitter(size=0.5)+
#   stat_compare_means(comparisons=list(c("Lung","Spleen"),c("Lung","Gut"),
#                                       c("Spleen","Gut")),size=2.5,
#                      label = "p.signif",method = 'wilcox.test')+
#   facet_grid(. ~ Type) +
#   ylim(0,7.5)+
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),strip.text = element_text(size= 6),
#         axis.text = element_text(size = 5))+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   scale_fill_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))
# Abundance.plot.sum
# 
# 
# Diversity <- Shannon.plot/Richness.plot
# Diversity
# Ann.Pathogens <- unique(tissue.table.long.taxo[,c('Pathogens','Type','Family')])
# Ann.Pathogens <- Ann.Pathogens[order(Ann.Pathogens$Family),]
# Ann.Pathogens <- Ann.Pathogens[order(Ann.Pathogens$Type),]
# ##Abundance
# Abundance.plot <- ggplot(subset(tissue.table.long.taxo,RPM > 0),aes(Pathogens,log10(RPM+1),color=Type))+
#   geom_boxplot()+labs(x="",y='Abundance of Tissues')+ geom_jitter(size=0.5)+
#   facet_grid(Tissue~.)+theme_few()+
#   scale_x_discrete(limits=unique(factor(Ann.Pathogens$Pathogens)))+
#   theme(legend.position = 'bottom',
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# Abundance.plot

##Venn
Lung.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung')$Pathogens)
Spleen.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen')$Pathogens)
Gut.Pathogens <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut')$Pathogens)
Tissue <- list(Lung=Lung.Pathogens,Spleen=Spleen.Pathogens,Gut=Gut.Pathogens)
Venn.plot <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),set_name_size = 3,text_size = 2,show_percentage=FALSE)
Venn.plot
# ###Pathogens Clade Venn
# Lung.Pathogens.dna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung' & Type == 'DNA Virus')$Pathogens)
# Spleen.Pathogens.dna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen' & Type == 'DNA Virus')$Pathogens)
# Gut.Pathogens.dna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut' & Type == 'DNA Virus')$Pathogens)
# Tissue <- list(Lung=Lung.Pathogens.dna,Spleen=Spleen.Pathogens.dna,Gut=Gut.Pathogens.dna)
# Venn.plot.dna <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),
#                         set_name_size = 3,text_size = 2,stroke_size=0.5,show_percentage=FALSE)
# Venn.plot.dna
# 
# Lung.Pathogens.rna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung' & Type == 'RNA Virus')$Pathogens)
# Spleen.Pathogens.rna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen' & Type == 'RNA Virus')$Pathogens)
# Gut.Pathogens.rna <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut' & Type == 'RNA Virus')$Pathogens)
# Tissue <- list(Lung=Lung.Pathogens.rna,Spleen=Spleen.Pathogens.rna,Gut=Gut.Pathogens.rna)
# Venn.plot.rna <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),
#                         set_name_size = 3,text_size = 2,stroke_size=0.5,show_percentage=FALSE)
# Venn.plot.rna
# 
# Lung.Pathogens.BAC <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung' & Type == 'Bacteria')$Pathogens)
# Spleen.Pathogens.BAC <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen' & Type == 'Bacteria')$Pathogens)
# Gut.Pathogens.BAC <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut' & Type == 'Bacteria')$Pathogens)
# Tissue <- list(Lung=Lung.Pathogens.BAC,Spleen=Spleen.Pathogens.BAC,Gut=Gut.Pathogens.BAC)
# Venn.plot.BAC <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),
#                         set_name_size = 3,text_size = 2,stroke_size=0.5,show_percentage=FALSE)
# Venn.plot.BAC
# 
# Lung.Pathogens.EUK <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung' & Type == 'Eukaryote')$Pathogens)
# Spleen.Pathogens.EUK <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen' & Type == 'Eukaryote')$Pathogens)
# Gut.Pathogens.EUK <- unique(subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut' & Type == 'Eukaryote')$Pathogens)
# Tissue <- list(Lung=Lung.Pathogens.EUK,Spleen=Spleen.Pathogens.EUK,Gut=Gut.Pathogens.EUK)
# Venn.plot.EUK <- ggvenn(Tissue,fill_color=c('#E1C855','#51B1B7','#E07B54'),
#                         set_name_size = 3,text_size = 2,stroke_size=0.5,show_percentage=FALSE)
# Venn.plot.EUK

##distribution
pathogens.distribution <- unique(subset(tissue.table.long.taxo,RPM >0)[,c('Pathogens','Tissue','Type')])
distribution.num <- aggregate(pathogens.distribution$Pathogens,
                              list(pathogens.distribution$Tissue,pathogens.distribution$Type),length)
colnames(distribution.num)  <- c('Tissue','Type','Num')
distribution.plot <- ggplot(distribution.num, aes(x=Tissue,y = Num,fill=Type)) +
  geom_bar(stat = "identity",
           width = 0.8, size = 0.25,
           alpha = 1)+ 
  labs(x="",y="Num of pathogens")+theme_few()+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))+
  scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
                                        c("RNA Virus","DNA Virus","Bacteria","Eukaryote")))

###linkage of virus
# tissue.table.long.indi <- merge(tissue.table.long.taxo,Library_infor[,c('Library','Individual')])
# ##wenzhou virus
# wzv.RPM <- subset(tissue.table.long.indi,Pathogens=='Wnezhou virus' & RPM >0)
# wzv.RPM.wide <- spread(wzv.RPM[c('Individual','RPM','Tissue')],key=Tissue,value = RPM,fill = 0)
# wzvall.plot <- ggplot(wzv.RPM,aes(Tissue,log10(RPM+1),color=Tissue))+
#   geom_boxplot()+labs(x="",y='Abundance of Tissues/log10(RPM+1)')+ geom_jitter(size=0.5)+
#   theme_few()+
#   ggtitle("Wenzhou virus")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# 
# wzvGS.plot <- ggplot(wzv.RPM.wide,aes(x=log10(Gut+1),y=log10(Spleen+1)))+
#   geom_point()+geom_vline(xintercept= 0,color='gray',linetype='dashed')+
#   geom_hline(yintercept = 0,color='gray',linetype='dashed')+
#   stat_smooth(method="lm",se=TRUE)+
#   stat_cor(data=wzv.RPM.wide, method = "pearson")+
#   theme_few()+
#   ggtitle("Wenzhou virus")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# 
# wzvLS.plot <- ggplot(wzv.RPM.wide,aes(x=log10(Lung+1),y=log10(Spleen+1)))+
#   geom_point()+geom_vline(xintercept= 0,color='gray',linetype='dashed')+
#   geom_hline(yintercept = 0,color='gray',linetype='dashed')+
#   stat_smooth(method="lm",se=TRUE)+
#   stat_cor(data=wzv.RPM.wide, method = "pearson")+
#   theme_few()+
#   ggtitle("Wenzhou virus")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# wzc.plot <- wzvall.plot+wzvGS.plot+wzvLS.plot

##Guangdong rodent  arterivirus 1
# Gra.RPM <- subset(tissue.table.long.indi,Pathogens=='Guangdong rodent arterivirus 1' & RPM >0)
# Gra.RPM.wide <- spread(Gra.RPM[c('Individual','RPM','Tissue')],key=Tissue,value = RPM,fill = 0)
# Gral.plot <- ggplot(Gra.RPM,aes(Tissue,log10(RPM+1),color=Tissue))+
#   geom_boxplot()+labs(x="",y='Abundance of Tissues/log10(RPM+1)')+ geom_jitter(size=0.5)+
#   theme_few()+
#   ggtitle("Guangdong rodent arterivirus 1")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# Gral.plot
# 
# GralGS.plot <- ggplot(Gra.RPM.wide,aes(x=log10(Gut+1),y=log10(Spleen+1)))+
#   geom_point()+geom_vline(xintercept= 0,color='gray',linetype='dashed')+
#   geom_hline(yintercept = 0,color='gray',linetype='dashed')+
#   stat_smooth(method="lm",se=TRUE)+
#   stat_cor(data=Gra.RPM.wide, method = "pearson")+
#   theme_few()+
#   ggtitle("Guangdong rodent arterivirus 1")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# GralGS.plot
# 
# GralLS.plot <- ggplot(Gra.RPM.wide,aes(x=log10(Lung+1),y=log10(Spleen+1)))+
#   geom_point()+geom_vline(xintercept= 0,color='gray',linetype='dashed')+
#   geom_hline(yintercept = 0,color='gray',linetype='dashed')+
#   stat_smooth(method="lm",se=TRUE)+
#   stat_cor(data=Gra.RPM.wide, method = "pearson")+
#   theme_few()+
#   ggtitle("Guangdong rodent arterivirus 1")+
#   scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# GralLS.plot
# 
# Gra.plot <- Gral.plot+GralGS.plot+GralLS.plot
# 
# example.polt <- wzc.plot/Gra.plot

###All Pathogen Polt
# tissue.table.long.indi$Pathogens <- ifelse(tissue.table.long.indi$Pathogens=='Wnezhou virus','Wenzhou virus',tissue.table.long.indi$Pathogens)
# rna.virus.data <- subset(tissue.table.long.indi,Type == 'RNA Virus' & RPM >0)
# rna.virus.plot <- ggplot(rna.virus.data,aes(Tissue,log10(RPM+1)))+geom_boxplot()+
#   geom_point()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, color="white") +
#   geom_line(aes(group=Individual), linetype="dashed", col="skyblue")+
#   facet_wrap(vars(Pathogens))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# rna.virus.plot
# 
# dna.virus.data <- subset(tissue.table.long.indi,Type == 'DNA Virus' & RPM >0)
# dna.virus.plot <- ggplot(dna.virus.data,aes(Tissue,log10(RPM+1)))+geom_boxplot()+
#   geom_point()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, color="white") +
#   geom_line(aes(group=Individual), linetype="dashed", col="skyblue")+
#   facet_wrap(vars(Pathogens))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# dna.virus.plot
# 
# bacteria.data <- subset(tissue.table.long.indi,Type == 'Bacteria' & RPM >0 )
# bacteria.plot <- ggplot(bacteria.data,aes(Tissue,log10(RPM+1)))+geom_boxplot()+
#   geom_point()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, color="white") +
#   geom_line(aes(group=Individual), linetype="dashed", col="skyblue")+
#   facet_wrap(vars(Pathogens))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# bacteria.plot
# 
# Eukaryote.data <- subset(tissue.table.long.indi,Type == 'Eukaryote' & RPM >0)
# Eukaryote.plot <- ggplot(Eukaryote.data,aes(Tissue,log10(RPM+1)))+geom_boxplot()+
#   geom_point()+
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.7, color="white") +
#   geom_line(aes(group=Individual), linetype="dashed", col="skyblue")+
#   facet_wrap(vars(Pathogens))+
#   theme(legend.position = 'none',axis.title=element_text(size= 6),
#         text=element_text(size= 5),
#         axis.text = element_text(size = 5),
#         axis.text.x=element_text(angle=60,hjust = 1))
# Eukaryote.plot

###Library View of Richness
# library.tissue.rpm <- gather(Library_Tissue_RPM,key = "Pathogens",value = RPM, -"Library",-'Hos_Taxo',
#                                     -'Tissue',-'Collected_Time',-'Collected_Region')
# library.tissue.rpm <- merge(library.tissue.rpm,Ann.Pathogens)
# library.tissue.rpm <- transform(library.tissue.rpm,RPM = as.numeric(RPM))
# library.tissue.rpm1 <- aggregate(library.tissue.rpm$RPM,
#           list(library.tissue.rpm$Library,library.tissue.rpm$Family,library.tissue.rpm$Tissue),mean)
# colnames(library.tissue.rpm1) <- c("Library",'family','Tissue','RPM')
# 
# 
# library.tissue.rpm1 <- library.tissue.rpm1[order(library.tissue.rpm1$RPM),]
# library.tissue.rpm1.sum <- aggregate(library.tissue.rpm1$RPM,
#                                 list(library.tissue.rpm1$Library),sum)
# library.tissue.rpm1.sum <- library.tissue.rpm1.sum[order(library.tissue.rpm1.sum$x,decreasing = T),]
# library.tissue.rpm1$Library <- factor(library.tissue.rpm1$Library,levels = library.tissue.rpm1.sum$Group.1)
# hhh <- Ann_row
# colnames(hhh) <- c('family','Type')
# library.tissue.rpm1 <- merge(library.tissue.rpm1,hhh,all.x = TRUE)
# library.tissue.rpm1 <- aggregate(library.tissue.rpm1$RPM,
#                                  list(library.tissue.rpm1$Type,library.tissue.rpm1$Library,library.tissue.rpm1$Tissue),sum)
# colnames(library.tissue.rpm1) <- c('Type','Library','Tissue','RPM')
# gut <- subset(library.tissue.rpm1,Tissue=='Gut' & RPM>0)
# gut.wide <- spread(gut[,c('Library','Type','RPM')],key = 'Type',value = 'RPM')
# gut.wide <- gut.wide[order(gut.wide$Bac,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$Euk,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`DNA Virus`,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`RNA Virus`,decreasing = TRUE),]
# gut$Library <- factor(gut$Library,levels = unique(gut.wide$Library))
# gut.plot <- ggplot(data = gut,aes(x=Library,y=log10(RPM+1),fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill")+
#   scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
#                                       c("RNA Virus","DNA Virus","Bac","Euk")))+theme_few()+
#   labs(x="")+ggtitle('Gut')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# gut.plot
# 
# gut <- subset(library.tissue.rpm1,Tissue=='Spleen' & RPM >0)
# gut.wide <- spread(gut[,c('Library','Type','RPM')],key = 'Type',value = 'RPM')
# gut.wide <- gut.wide[order(gut.wide$Bac,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$Euk,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`DNA Virus`,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`RNA Virus`,decreasing = TRUE),]
# gut$Library <- factor(gut$Library,levels = unique(gut.wide$Library))
# Spleen.plot <- ggplot(data = gut,aes(x=Library,y=log10(RPM+1),fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill")+
#   scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
#                                       c("RNA Virus","DNA Virus","Bac","Euk")))+theme_few()+
#   labs(x="")+ggtitle('Spleen')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# Spleen.plot
# 
# gut <- subset(library.tissue.rpm1,Tissue=='Lung'& RPM >0)
# gut.wide <- spread(gut[,c('Library','Type','RPM')],key = 'Type',value = 'RPM')
# gut.wide <- gut.wide[order(gut.wide$Bac,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$Euk,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`DNA Virus`,decreasing = TRUE),]
# gut.wide <- gut.wide[order(gut.wide$`RNA Virus`,decreasing = TRUE),]
# gut$Library <- factor(gut$Library,levels = unique(gut.wide$Library))
# Lung.plot <- ggplot(data = gut,aes(x=Library,y=log10(RPM+1),fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill")+
#   scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
#                                       c("RNA Virus","DNA Virus","Bac","Euk")))+theme_few()+
#   labs(x="")+ggtitle('Lung')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# Lung.plot

# log.library.plot <- ggplot(data = library.tissue.rpm1,aes(x=Library,y=log10(RPM+1),fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill")+
#   scale_fill_manual(values = setNames(c('#DB432C','#438870','#838AAF','#C4B797'),
#                                       c("RNA Virus","DNA Virus","Bac","Euk")))+theme_few()+
#   facet_grid(. ~ Tissue)+ labs(x="")+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# 
# Euk.richness <- merge(Euk.richness,Library.anno)
# Euk.richness <- Euk.richness[order(Euk.richness$Richness,decreasing = T),]
# Euk.richness <- Euk.richness[order(Euk.richness$Tissue,decreasing = T),]
# Euk.richness$Library <- factor(Euk.richness$Library,levels = Euk.richness$Library)
# Euk.library.plot <- ggplot(data = Euk.richness,aes(x=Library,y=Richness,fill=Tissue))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25)+
#   scale_fill_manual(values = c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen'))+theme_few()+
#   labs(x="",y='Eukaryote Richness')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# 
# Bac.richness <- merge(Bac.richness,Library.anno)
# Bac.richness <- Bac.richness[order(Bac.richness$Richness,decreasing = T),]
# Bac.richness <- Bac.richness[order(Bac.richness$Tissue,decreasing = T),]
# Bac.richness$Library <- factor(Bac.richness$Library,levels = Bac.richness$Library)
# Bac.library.plot <- ggplot(data = Bac.richness,aes(x=Library,y=Richness,fill=Tissue))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25)+
#   scale_fill_manual(values = c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen'))+theme_few()+
#   labs(x="",y='Bacteria Richness')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# 
# RNAvirus.richness <- merge(RNAvirus.richness,Library.anno)
# RNAvirus.richness <- RNAvirus.richness[order(RNAvirus.richness$Richness,decreasing = T),]
# RNAvirus.richness <- RNAvirus.richness[order(RNAvirus.richness$Tissue,decreasing = T),]
# RNAvirus.richness$Library <- factor(RNAvirus.richness$Library,levels = RNAvirus.richness$Library)
# RNAvirus.library.plot <- ggplot(data = RNAvirus.richness,aes(x=Library,y=Richness,fill=Tissue))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25)+
#   scale_fill_manual(values = c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen'))+theme_few()+
#   labs(x="",y='RNA virus Richness')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# 
# DNAvirus.richness <- merge(DNAvirus.richness,Library.anno)
# DNAvirus.richness <- DNAvirus.richness[order(DNAvirus.richness$Richness,decreasing = T),]
# DNAvirus.richness <- DNAvirus.richness[order(DNAvirus.richness$Tissue,decreasing = T),]
# DNAvirus.richness$Library <- factor(DNAvirus.richness$Library,levels = DNAvirus.richness$Library)
# DNAvirus.library.plot <- ggplot(data = DNAvirus.richness,aes(x=Library,y=Richness,fill=Tissue))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25)+
#   scale_fill_manual(values = c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen'))+theme_few()+
#   labs(x="",y='DNA virus Richness')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_blank())
# Library.plot <- DNAvirus.library.plot/RNAvirus.library.plot/Bac.library.plot/Euk.library.plot

###Polt 12 filgures
selected.pathogens <- c('Rat hepatitis E virus','Wenzhou virus',
                        'Betacoronavirus 1',
                        'Orthohantavirus seoulense',
                        'Norwalk virus','Bartonella kosoyi',
                        'Klebsiella variicola',
                        'Angiostrongylus cantonensis',
                        'Cryptosporidium ubiquitum',
                        'Nippostrongylus brasiliensis',
                        'Pneumocystis carinii','Trypanosoma lewis')

function_name <- function(R){
  hhh <- subset(tissue.table.long.indi,Pathogens==R)
  ggplot(subset(hhh,RPM>0),aes(Tissue,log10(RPM+1),color=Tissue))+
    geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
    theme_few()+
    ggtitle(R)+
    ylab('')+
    scale_color_manual(values = setNames(c('#E07B54','#E1C855','#51B1B7'),c('Gut','Lung','Spleen')))+
    scale_x_discrete(limits = c('Gut','Lung','Spleen'))+
    theme(legend.position = 'none',axis.title=element_text(size= 6),
          plot.title = element_text(hjust = 0.5 ),
          text=element_text(size= 5),
          axis.text = element_text(size = 5),
          axis.text.x=element_text(angle=60,hjust = 1))}

sumplot <- lapply(selected.pathogens,function_name)

plot_sum <- plot_grid(plotlist = sumplot, align = "h", 
          nrow = 2)
plot_sum
###selected_Table
selected.table <- subset(tissue.table.long.indi,Pathogens%in%selected.pathogens)
selected.table$Pathogens <- factor(selected.table$Pathogens,levels = selected.pathogens)

selected.pa.plot <- ggplot(subset(selected.table,RPM>0),aes(Tissue,log10(RPM+1),color=Tissue))+
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
selected.pa.plot


###coninfection
# Tissue.Pathogens.lung <- subset(tissue.table.long.taxo,RPM >0 & Tissue=='Lung')
# lung1 <- aggregate(Tissue.Pathogens.lung$Pathogens,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung1) <- c('Library','Pathogen_Num')
# lung2 <- aggregate(Tissue.Pathogens.lung$Type,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung2) <- c('Library','Type_Num')
# lung <- merge(lung1,lung2)
# lung$Type <- ifelse(lung$Pathogen_Num ==1 ,'Single Pathogen',
#                     ifelse(lung$Type_Num ==1 ,'Single Pathogen Type','Multi Pathogen Type'))
# color_set <- setNames(c('#D3D3D3','#DEBF80','#DCCD5B','#6179A7'),
#                       c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# lung.sum <- aggregate(lung$Library,
#                       list(lung$Type),length)
# lung.sum[nrow(lung.sum) + 1,] = c('No Pathogen', 693-sum(lung.sum$x))
# colnames(lung.sum) <- c('Type','Num')
# lung.sum <- lung.sum[order(lung.sum$Num,decreasing = TRUE),]
# lung.sum$Type <- factor(lung.sum$Type,levels = c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# lung.sum <- transform(lung.sum,Num = as.numeric(Num))
# lung.sum <- lung.sum[order(lung.sum$Type),]
# lung.pie <- ggplot(lung.sum, aes(x="", y=Num, fill=Type))+
#   geom_bar(width = 1, stat = "identity")+labs(x='',y='')+
#   scale_fill_manual(values = color_set)+
#   coord_polar("y", start=0,direction = -1)+theme_minimal()+
#   geom_text(aes(x=1.2,y=sum(lung.sum$Num)-cumsum(lung.sum$Num)+lung.sum$Num/2 ,
#                 label=as.character(lung.sum$Num)),size=3)
# lung.pie
# 
# Tissue.Pathogens.lung <- subset(tissue.table.long.taxo,RPM >0 & Tissue=='Spleen')
# lung1 <- aggregate(Tissue.Pathogens.lung$Pathogens,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung1) <- c('Library','Pathogen_Num')
# lung2 <- aggregate(Tissue.Pathogens.lung$Type,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung2) <- c('Library','Type_Num')
# lung <- merge(lung1,lung2)
# lung$Type <- ifelse(lung$Pathogen_Num ==1 ,'Single Pathogen',
#                     ifelse(lung$Type_Num ==1 ,'Single Pathogen Type','Multi Pathogen Type'))
# color_set <- setNames(c('#D3D3D3','#DEBF80','#DCCD5B','#6179A7'),
#                       c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# s.sum <- aggregate(lung$Library,
#                       list(lung$Type),length)
# s.sum[nrow(s.sum) + 1,] = c('No Pathogen', 693-sum(s.sum$x))
# colnames(s.sum) <- c('Type','Num')
# s.sum <- s.sum[order(s.sum$Num,decreasing = TRUE),]
# s.sum$Type <- factor(s.sum$Type,levels = c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# s.sum$pro <- (s.sum$Num/693) * 100
# s.sum <- transform(s.sum,Num = as.numeric(Num))
# s.sum <- s.sum[order(s.sum$Type),]
# spleen.pie <- ggplot(s.sum, aes(x="", y=Num, fill=Type))+
#   geom_bar(width = 1, stat = "identity")+labs(x='',y='')+
#   scale_fill_manual(values = color_set)+
#   coord_polar("y", start=0,direction = -1)+theme_minimal()+
#   geom_text(aes(x=1.2,y=sum(s.sum$Num)-cumsum(s.sum$Num)+s.sum$Num/2 ,
#                 label=as.character(s.sum$Num)),size=3)
# spleen.pie
# 
# Tissue.Pathogens.lung <- subset(tissue.table.long.taxo,RPM >0 & Tissue=='Gut')
# lung1 <- aggregate(Tissue.Pathogens.lung$Pathogens,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung1) <- c('Library','Pathogen_Num')
# lung2 <- aggregate(Tissue.Pathogens.lung$Type,
#                    list(Tissue.Pathogens.lung$Library),length)
# colnames(lung2) <- c('Library','Type_Num')
# lung <- merge(lung1,lung2)
# lung$Type <- ifelse(lung$Pathogen_Num ==1 ,'Single Pathogen',
#                     ifelse(lung$Type_Num ==1 ,'Single Pathogen Type','Multi Pathogen Type'))
# color_set <- setNames(c('#D3D3D3','#DEBF80','#DCCD5B','#6179A7'),
#                       c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# g.sum <- aggregate(lung$Library,
#                       list(lung$Type),length)
# g.sum[nrow(g.sum) + 1,] = c('No Pathogen', 693-sum(g.sum$x))
# colnames(g.sum) <- c('Type','Num')
# g.sum <- g.sum[order(g.sum$Num,decreasing = TRUE),]
# g.sum$Type <- factor(g.sum$Type,levels = c('No Pathogen','Single Pathogen','Single Pathogen Type','Multi Pathogen Type'))
# g.sum <- transform(g.sum,Num = as.numeric(Num))
# g.sum <- g.sum[order(g.sum$Type),]
# gut.pie <- ggplot(g.sum, aes(x="", y=Num, fill=Type))+
#   geom_bar(width = 1, stat = "identity")+labs(x='',y='')+
#   scale_fill_manual(values = color_set)+
#   coord_polar("y", start=0,direction = -1)+theme_minimal()+
#   geom_text(aes(x=1.2,y=sum(g.sum$Num)-cumsum(g.sum$Num)+g.sum$Num/2 ,
#                 label=as.character(g.sum$Num)),size=3)
# gut.pie
# plot <- gut.pie+lung.pie+spleen.pie+plot_layout(guides = 'collect')
# 
# ggsave(
#   filename = paste('F:/Mouse_Result_24/Figure/TMP1/','Figure s3.6 coinfection.pdf',sep=''),
#   plot,
#   width = 180,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )


##richness & Abudance

###save figures
ggsave(filename = paste0(out_path,'/fig 3a.pdf',sep=''), distribution.plot, width = 60, height = 40,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3b.pdf',sep=''), Venn.plot, width = 40, height = 40,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3c.pdf',sep=''),Relative.plot,  width = 60, height = 120,units = 'mm')
f3d <- plot_grid(
               Richness.plot+theme(legend.position = "none"),
               Abundance.plot.sum+theme(legend.position = "none"),align = "v", axis = "l",
               ncol = 1,rel_heights = c(0.7,0.7)
               )
ggsave(filename = paste0(out_path,'/fig 3d.pdf',sep=''),f3d,  width = 120, height = 80,units = 'mm')
ggsave(filename = paste0(out_path,'/fig 3e.pdf',sep=''),selected.pa.plot,  width = 150, height = 70,units = 'mm')
P1 <- plot_grid(distribution.plot+theme(legend.position = "none"),
                Venn.plot,ncol = 2,rel_widths = c(1,1))+draw_plot_label(
                  c("a","b"),
                  c(0,0.5),
                  c(1,1),
                  size = 8
                )
P2 <- plot_grid(P1,
                Richness.plot+theme(legend.position = "none"),
                Abundance.plot.sum+theme(legend.position = "none"),align = "v", axis = "l",
                ncol = 1,rel_heights = c(1,0.7,0.7))+draw_plot_label(
                  c("d"),
                  c(0),
                  c(0.66),
                  size = 8
                )
P2
P3 <- plot_grid(Relative.plot+theme(legend.position = "none"),NULL,
                ncol = 1,rel_heights = c(15,1))+
  draw_plot_label(
    c("c"),
    c(0),
    c(1),
    size = 8
  )
P4 <- plot_grid(P2,P3,
                ncol = 2,rel_widths = c(2,1.5))
P4

P5 <- plot_grid(P4,selected.pa.plot+theme(legend.position = "none"),
                ncol = 1,rel_heights = c(3,2))+
  draw_plot_label(
    c("e"),
    c(0),
    c(0.4),
    size = 8
  )
P5

ggsave(
  filename = paste0(out_path,'/fig 3.pdf',sep=''),
  P5,
  width = 150,             # 宽
  height = 180,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

# SUB.VENN <- plot_grid(Venn.plot.rna,Venn.plot.dna,Venn.plot.BAC,Venn.plot.EUK,
#                       nrow = 1,rel_widths = c(1,1,1,1))+draw_plot_label(
#   c("RNA Virus",'DNA Virus','Bacteria','Eukaryote'),
#   c(0.05,0.05+0.25,0.05+0.5,0.05+0.75),
#   c(0.1,0.1,0.1,0.1),
#   size = 8
# )
# SUB.VENN
# ggsave(
#   filename = paste('F:/Mouse_Result_24/Figure/TMP1/','Figure s3.4 venn.pdf',sep=''),
#   SUB.VENN,
#   width = 120,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )

ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP1/','Figure3.test.pdf',sep=''),
  P2,
  width = 150*2/3.5,             # 宽
  height = 120,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

