####figure4
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','raw_data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")
##install package
# install.packages(ggthemes)
# install.packages(tidyr)
# install.packages(ggplot2)
# install.packages(pairwiseAdonis)
# install.packages(cowplot)
# install.packages(dplyr)
# install.packages(RColorBrewer)
# install.packages(scales)
# install.packages(ade4)   # 用于计算PcoA
# install.packages(vegan)  # 用于计算距离
# install.packages(ggpubr)
# install.packages(ComplexUpset)
# install.packages(ggsci)
# install.packages(tidygraph)
# install.packages(ggrepel)
# install.packages(showtext)
# install.packages(circlize)
# install.packages(igraph)
# install.packages(ggraph)

library(ggthemes)
library(tidyr)
library(ggplot2)
library(pairwiseAdonis)
library(cowplot)
library(dplyr)
library(RColorBrewer)
library(scales)
library(ade4)   # 用于计算PcoA
library(vegan)  # 用于计算距离
library(ggpubr)
library(ComplexUpset)
library(ggsci)
library(tidygraph)
library(ggrepel)
library(showtext)
library(circlize)
library(igraph)
library(ggraph)
##Read Pathogen file
RMP_Table <- read.csv('Pathogens_Mammal-related_Wide_RPM_Table.csv',check.names = F)
RMP_Table <- subset(RMP_Table,select = -c(`Desulfovibrio piger`,`Oryza sativa`))
colnames(RMP_Table)[1] <- 'Individual'
Ann_Table <- read.csv('Pathongen_Taxon_Table_2.txt',sep='\t')
Ann_Table <- unique(Ann_Table)
Mammal.releted <- subset(Ann_Table,Host_Type == 'Mammal-related')
Individual.infor <- read.csv('Full_Tissue_Individual.csv')
Ann_Table$Type <- ifelse(Ann_Table$Type == 'RNA_Virus','RNA Virus',
                         ifelse(Ann_Table$Type == 'Bac','Bacteria',
                                ifelse(Ann_Table$Type == 'Euk','Eukaryote',
                                       ifelse(Ann_Table$Type == 'DNA_Virus','DNA Virus',Ann_Table$Type))))
Individual.infor$Hos_Taxo <- ifelse(Individual.infor$Hos_Taxo == 'Berylmys sp.','Berylmys bowersi',Individual.infor$Hos_Taxo)
Individual.infor$Collected_Region <- ifelse(Individual.infor$Collected_Region == 'AP','Anpu',
                                            ifelse(Individual.infor$Collected_Region == 'ZJ','Zhanjiang',
                                                   ifelse(Individual.infor$Collected_Region == 'DB','Maoming',
                                                          ifelse(Individual.infor$Collected_Region == 'HY','Heyuan',
                                                                 ifelse(Individual.infor$Collected_Region == 'JY','Jieyang',
                                                                        ifelse(Individual.infor$Collected_Region == 'SG','Shaoguan',
                                                                               ifelse(Individual.infor$Collected_Region == 'SS','Foshan',
                                                                                      ifelse(Individual.infor$Collected_Region == 'SZ','Shenzhen',
                                                                                             ifelse(Individual.infor$Collected_Region == 'YF','Yunfu',''
                                                                                             )))))))))
library(dplyr)
Mammal.releted <- Mammal.releted %>%
  filter(!(family %in% c('Picobirnaviridae','Anelloviridae')))

Pathogens <- intersect(colnames(RMP_Table),Mammal.releted$species)
RMP_Table <- RMP_Table[,c('Individual',Pathogens)]

###colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
            "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
color20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")

Host <- unique(Individual.infor$Hos_Taxo)
Host_color <- colors[1:length(Host)]
Region <- unique(Individual.infor$Collected_Region)
Region_color <- color20[1:length(Region)]
Time <- unique(Individual.infor$Collected_Time)
Time_color <- color20[(length(Region)+1):(length(Region)+length(Time))]
Type <- unique(Ann_Table$Type)
Type_color <- c("#91D1C2FF","#DC0000FF","#8491B4FF","#B09C85FF")
Family <- unique(Mammal.releted$family)
Family_color <- color37[1:length(Family)]
anno_color <- list(Host=setNames(Host_color,Host),
                   Time=setNames(Time_color,Time),
                   Region=setNames(Region_color,Region),
                   Level=setNames(Type_color,Type),
                   classify=setNames(Family_color,Family)
)
Individual.infor$Collected_Time <- factor(Individual.infor$Collected_Time,level=c("2021-Winter","2022-Spring","2022-Summer","2022-Autumn","2022-Winter"))
Ann_Table$Type <- factor(Ann_Table$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))

###bar graph
RPM.table.long <- gather(RMP_Table,key = "species",value = RPM, -"Individual")
RPM.table.long <- merge(RPM.table.long,Ann_Table)
RPM.table.long <- merge(RPM.table.long,Individual.infor)
lwd_pt <- .pt*72.27/96

RPM.table.long.sum <- aggregate(RPM.table.long$RPM,
                                list(RPM.table.long$Individual),sum)
colnames(RPM.table.long.sum) <- c("Individual","RPM")
RPM.table.long.sum <- merge(RPM.table.long.sum,Individual.infor)
RPM.table.long.sum.individual.sum <- aggregate(RPM.table.long.sum$Individual,
                                               list(RPM.table.long.sum$Hos_Taxo),length)
colnames(RPM.table.long.sum.individual.sum) <- c('Taxo','Num')
RPM.table.long.sum.individual.sum <- RPM.table.long.sum.individual.sum[order(RPM.table.long.sum.individual.sum$Num,decreasing = TRUE),]
RPM.table.long.sum$Hos_Taxo <- factor(RPM.table.long.sum$Hos_Taxo,levels = RPM.table.long.sum.individual.sum$Taxo)
RPM.table.long.sum <- RPM.table.long.sum[order(RPM.table.long.sum$RPM,decreasing = T),]
RPM.table.long.sum <- RPM.table.long.sum[order(RPM.table.long.sum$Hos_Taxo),]

RPM.table.long.sum$Individual <- factor(RPM.table.long.sum$Individual,levels = RPM.table.long.sum$Individual)
# Library.plot <- ggplot(data = RPM.table.long.sum,aes(x=Individual,y=log10(RPM+1)))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,fill='#9FC9DF')+
#   #scale_fill_manual(values = setNames(Family_color,Family))+
#   theme_few()+
#   labs(x="")+
#   geom_point(aes(x=Individual,y=-0.1,color=Hos_Taxo),shape=15,size=0.4)+
#   scale_color_manual(values = setNames(Host_color,Host))+
#   theme(legend.key.size = unit(3,'mm'),
#         text=element_text(size= 6),
#         axis.text.x=element_blank())
# Library.plot
###Abundance plot
# Abundance.plot <- ggplot(subset(RPM.table.long,RPM >0),aes(Hos_Taxo,log10(RPM+1),color=Hos_Taxo))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
#   facet_grid(. ~ Type) + scale_color_manual(values=setNames(Host_color,Host))+
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0))
# 
# ggsave(
#   filename = paste('F:/Mouse_Result_24/Figure/TMP1/','Figure 4.2 pathogens Abundance.pdf',sep=''),
#   Abundance.plot,
#   width = 120,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )

##richness
individual.rna.wide <- spread(subset(RPM.table.long,Type=='RNA Virus')[,c('Individual','species','RPM')],
                              key = 'species',value = 'RPM')
rownames(individual.rna.wide) <- individual.rna.wide$Individual
individual.rna.wide <- individual.rna.wide[,-1]
rna.richness <- specnumber(t(individual.rna.wide), MARGIN = 2)#spe.rich =sobs

individual.dna.wide <- spread(subset(RPM.table.long,Type=='DNA Virus')[,c('Individual','species','RPM')],
                              key = 'species',value = 'RPM')
rownames(individual.dna.wide) <- individual.dna.wide$Individual
individual.dna.wide <- individual.dna.wide[,-1]
dna.richness <- specnumber(t(individual.dna.wide), MARGIN = 2)#spe.rich =sobs

individual.bac.wide <- spread(subset(RPM.table.long,Type=='Bacteria')[,c('Individual','species','RPM')],
                              key = 'species',value = 'RPM')
rownames(individual.bac.wide) <- individual.bac.wide$Individual
individual.bac.wide <- individual.bac.wide[,-1]
bac.richness <- specnumber(t(individual.bac.wide), MARGIN = 2)#spe.rich =sobs

individual.euk.wide <- spread(subset(RPM.table.long,Type=='Eukaryote')[,c('Individual','species','RPM')],
                              key = 'species',value = 'RPM')
rownames(individual.euk.wide) <- individual.euk.wide$Individual
individual.euk.wide <- individual.euk.wide[,-1]
euk.richness <- specnumber(t(individual.euk.wide), MARGIN = 2)#spe.rich =sobs
richness <- as.data.frame(cbind(rna.richness,dna.richness,euk.richness,bac.richness))
richness$Individual <- rownames(richness)
richness.long <- gather(richness,key = 'Type',value = 'Num',-'Individual')
richness.long$Type <- ifelse(richness.long$Type == 'rna.richness','RNA virus',
                             ifelse(richness.long$Type == 'dna.richness','DNA virus',
                                    ifelse(richness.long$Type == 'bac.richness','Bacteria',
                                           ifelse(richness.long$Type == 'euk.richness','Eukaryote',''
                                           ))))
richness.long <- merge(richness.long,Individual.infor)
richness.long$Type <- factor(richness.long$Type,levels = c("RNA virus","DNA virus","Bacteria","Eukaryote"))
richness.long$Hos_Taxo <- factor(richness.long$Hos_Taxo,levels = RPM.table.long.sum.individual.sum$Taxo)
richness.long.sum <- aggregate(richness.long$Num,
                               list(richness.long$Individual),sum)
richness.long.sum <- richness.long.sum[order(richness.long.sum$x,decreasing = T),]
richness.long$Individual <- factor(richness.long$Individual,levels = richness.long.sum$Group.1)
richness.long <- richness.long[order(richness.long$Type),]
richness.long <- richness.long[order(richness.long$Individual),]
richness.long <- richness.long[order(richness.long$Hos_Taxo),]
richness.long$Individual <- factor(richness.long$Individual,levels = unique(richness.long$Individual))
richness.long1 <- richness.long
colnames(richness.long1) <- c("Individual","Type","Num","Collected_Time","Collected_Region","Host")
richness.plot.sum <- ggplot(data = richness.long1,aes(x=Individual,y=Num))+
  geom_bar(stat = "identity",width = 1,size = 0.25,fill='#9FC9DF')+
  # scale_fill_manual(values = setNames(c("#DC0000FF","#91D1C2FF","#8491B4FF","#B09C85FF"),
  #                                     c("RNA virus","DNA virus","Bacteria","Eukaryote")))+
  theme_few()+
  labs(x="")+
  geom_point(aes(x=Individual,y=-0.1,color=Host),shape=15,size=0.4)+
  scale_color_manual(values = setNames(Host_color,Host))+
  theme(legend.key.size = unit(2,'mm'),legend.position = 'bottom',
        text=element_text(size= 6),
        axis.text.x=element_blank())
richness.plot.sum

ggsave(
  filename = paste0(out_path,'/fig 4a.pdf',sep=''),
  richness.plot.sum,
  width = 120,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

# write.csv(richness.long1,'F:/Mouse_Result_24/Manuscripts/richness.long1.csv')

##物种柱状图  
# Species <- aggregate(RPM.table.long$RPM,
#                      list(RPM.table.long$family,RPM.table.long$Hos_Taxo,
#                           RPM.table.long$Type),mean)
# colnames(Species) <- c("Family","Host","Type","RPM")
# Species$Type <- factor(Species$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
# Species_plot <- ggplot(data = Species,aes(x=Host,y=log10(RPM+1),fill=Family))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
#   )+scale_fill_manual(values = setNames(Family_color,Family))+theme_few()+
#   facet_grid(. ~ Type)+ labs(x="")+
#   theme(legend.key.size = unit(3,'mm'),
#         text=element_text(size= 6),
#         axis.text.x=element_text(angle=30,hjust = 1))
# Species_plot

###PCA/PCoA/NSDM
# library(ade4)   # 用于计算PcoA
# library(vegan)  # 用于计算距离
# library(pairwiseAdonis)
# library(ggpubr)
# host.pathogen.table <- aggregate(RPM.table.long$RPM,
#                                  list(RPM.table.long$species,RPM.table.long$Hos_Taxo,
#                                       RPM.table.long$Type,RPM.table.long$Collected_Time,
#                                       RPM.table.long$Collected_Region),mean)
# colnames(host.pathogen.table) <- c('Pathogens','Host','Type','Time','Site','RPM')
# host.pathogen.table$Region <- paste(host.pathogen.table$Host,host.pathogen.table$Time,host.pathogen.table$Site,
#                                     sep='-')
# host.pathogen.table.wide <- spread(subset(host.pathogen.table,RPM >0)[,c('Region','Pathogens','RPM')],key = Pathogens,
#                                    value = RPM,fill = 0)
# rownames(host.pathogen.table.wide) <- host.pathogen.table.wide$Region
# host.pathogen.table.wide <- subset(host.pathogen.table.wide,select = -c(Region))
# df.dist = vegdist(host.pathogen.table.wide,method='bray')
# pcoa <- cmdscale(df.dist, k = (nrow(host.pathogen.table.wide) - 1), eig = TRUE)
# plot_data <- data.frame({pcoa$point})[1:2]
# plot_data$Region <- rownames(plot_data)
# names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
# eig = pcoa$eig
# plot_data <- merge(plot_data, unique(host.pathogen.table[,c('Region','Host')]), by = 'Region', all.x = TRUE)
# Host <- unique(host.pathogen.table$Host)
# Host_color <- colors[1:length(Host)]
# PCoA_Host <- ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2,color=Host)) +
#   geom_point(alpha=.7, size=1) +
#   labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
#        y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
#   scale_colour_manual(values = setNames(Host_color,Host))+
#   stat_ellipse(aes(fill = Host,color=Host), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
#   theme_few()+
#   theme(legend.key.size = unit(3,'mm'),
#         text=element_text(size= 6))+coord_fixed(1)
# PCoA_Host


##Abundance Polt
# Species.abun.plot <- ggplot(data = subset(Species,RPM>0),aes(x=Family,y=RPM,fill=Host))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
#   )+scale_fill_manual(values = setNames(Host_color,Host))+theme_few()+
#   labs(x="")+
#   theme(text=element_text(size= 6),legend.key.size = unit(3,'mm'),
#         axis.text.x=element_text(angle=30,hjust = 1))
# Species.abun.plot
# ggsave(
#   filename = paste('F:/Mouse_Result_24/Figure/TMP1/','Figure S4.2_Abundance_of_Host.pdf',sep=''),
#   Species.abun.plot,
#   width = 150,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )

###Upset
Upset.long <- RPM.table.long
Upset_sum <- aggregate(Upset.long$RPM, by = list(Upset.long$Hos_Taxo,Upset.long$species,Upset.long$Type), 
                       FUN = sum)
colnames(Upset_sum) <- c("Host","Pathogen","Type","RPM")
Upset_sum$RPM <- ifelse(Upset_sum$RPM > 0 ,1,0)
Upset_sum_wide <- spread(Upset_sum,key=Host,value = RPM,fill = 0)
Host <- unique(colnames(Upset_sum_wide[,-c(1,2)]))
Upset_Plot <- upset(Upset_sum_wide, Host,width_ratio=0.15,height_ratio=1,
                    matrix=(
                      intersection_matrix(
                        geom=geom_point(
                          size=1
                        ))),
                    base_annotations = list(
                      'Intersection size'=intersection_size(mapping = aes(fill=Type),
                                                            text=list(size=4/lwd_pt))+
                        
                        scale_fill_manual(values=setNames(Type_color,Type))+
                        ylab('Num of Shared Pathogens')+
                        theme(axis.title = element_text(size=5),
                              axis.text = element_text(size=5),
                              text = element_text(size=5),
                              legend.key.size = unit(3,'mm'),
                              legend.position = "left")),
                    set_sizes = (
                      upset_set_size()+geom_text(aes(label=..count..),size=4/lwd_pt, hjust=1.1, stat='count')+
                        ylab('Num of Pathogens')+
                        theme(axis.title.x = element_text(size=5),
                              axis.text = element_text(size=5),
                              text = element_text(size=5))
                    ),
                    queries=list(upset_query(set='Suncus murinus', fill='red')),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(text=element_text(size=5)),
                        'overall_sizes'=theme(axis.text.x=element_text(angle=90))
                      ))
)
Upset_Plot
ggsave(
  filename = paste(out_path,'fig 4d.pdf',sep='/'),Upset_Plot, # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 120,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
###Network
# library(tidygraph)
# library(ggraph)
noda.names <- c("Betacoronavirus 1","Beilong virus","Guangdong rodent arterivirus 1",
                "Guangdong rodent arterivirus 2","Norwalk virus",'Porcine bocavirus','Rat minute virus 2a',
                "Norway rat pestivirus","Orthohantavirus seoulense","Rat hepatitis E virus",
                "Rodent astrovirus 1","Rodent astrovirus 3","Rodent coronavirus","Wnezhou virus",
                "Guangdong murine mastadenovirus 2","Guangdong rodent dependoparvovirus 1",
                "Bartonella kosoyi","Angiostrongylus cantonensis","Cryptosporidium ubiquitum",
                "Nippostrongylus brasiliensis","Tritrichomonas sp","Chlamydia muridarum",
                "Bandicota indica","Berylmys bowersi","Mus caroli","Niviventer lotipes",
                "Rattus andamanensis","Rattus losea","Rattus norvegicus",
                "Rattus tanezumi","Suncus murinus")
graph_long <- gather(Upset_sum_wide,key = "Node",value = 'RPM', -"Pathogen",-'Type')
graph_long <- subset(graph_long,RPM > 0)
# selected_Path <- aggregate(graph_long$Node,list(graph_long$Pathogen),length)
# colnames(selected_Path) <- c("Pathogen","Num")
# Pathogens <- subset(selected_Path,Num >1)$Pathogen
# graph_long <- subset(graph_long,Pathogen %in% Pathogens)
nodes <- data.frame(name=unique(c(graph_long$Pathogen,graph_long$Pathogen)))
relations <- data.frame(from=graph_long$Node,to=graph_long$Pathogen,Host=graph_long$Node)
g <- as_tbl_graph(relations)
g <- g %>% mutate(Type = ifelse(name %in% graph_long$Node, 'Host',
                                ifelse(name %in% subset(Ann_Table,Type=='Bacteria')$species, 'Bacteria',
                                       ifelse(name %in% subset(Ann_Table,Type=='Eukaryote')$species,'Eukaryote',
                                              ifelse(name %in% subset(Ann_Table,Type=='RNA Virus')$species,'RNA Virus',
                                                     ifelse(name %in% subset(Ann_Table,Type=='DNA Virus')$species,'DNA Virus','Others')))))) %>%
  mutate(Name = ifelse(name %in% graph_long$Node, name,
                       ifelse(name %in% subset(Ann_Table,Type=='Bacteria')$species, 'Bacteria',
                              ifelse(name %in% subset(Ann_Table,Type=='Eukaryote')$species,'Eukaryote',
                                     ifelse(name %in% subset(Ann_Table,Type=='RNA Virus')$species,'RNA Virus',
                                            ifelse(name %in% subset(Ann_Table,Type=='DNA Virus')$species,'DNA Virus',"Pathogen"))))))
g <- g %>% mutate(showd_name = ifelse(name %in% noda.names, name, ''
))


Name <- unique(c(graph_long$Node,"Pathogen"))
Name_color <- colors[1:length(Name)]
Name_color <- c(Name_color,c('#9BC985','#F7D58B','#797BB7','#B595BF'))
Name <- c(Name,c('Bacteria','Eukaryote','RNA Virus','DNA Virus'))
# library(ggrepel)
# library(showtext)
font_add('Arial','/Library/Fonts/Arial.ttf') 
showtext_auto()
Network <- ggraph(g,layout = 'stress') +
  geom_edge_fan(aes(color = Host),alpha=0.5,show.legend = FALSE) +
  geom_node_point(aes(shape=Type,fill=Name,color='black'),size=3,position="jitter") + 
  scale_shape_manual(values = c('Host'=23,'Bacteria'=21,'Eukaryote'=22,
                                'RNA Virus'=24,'DNA Virus'=25))+
  scale_fill_manual(values = setNames(Name_color,Name))+
  scale_edge_color_manual(values = setNames(Name_color,Name))+
  scale_color_manual(values = setNames(Name_color,Name))+
  #geom_node_text(aes(label=showd_name),size=2.5)+
  theme_graph()+
  theme(text=element_text(size= 5),legend.position = 'right')+coord_fixed()

Network

###绘制不同病毒Richness和Abundant的多样性之间的关系
##All Abundance & Richness
# Individual.wide.rpm <- spread(RPM.table.long[,c('Individual','species','RPM')],key = species,value = RPM,fill = 0)
# rownames(Individual.wide.rpm) <- Individual.wide.rpm$Individual
# Individual.wide.rpm <- Individual.wide.rpm[,-1]
# Shannon <- vegan::diversity(t(Individual.wide.rpm), index = "shannon", MARGIN = 2, base = exp(1))
# Richness <- specnumber(t(Individual.wide.rpm), MARGIN = 2)#spe.rich =sobs
# all.diversity <- as.data.frame(cbind(Shannon, Richness))
# all.diversity$Individual <- rownames(all.diversity)
# all.diversity <- merge(all.diversity,Individual.infor)
# all.diversity <- subset()
# host.compare <- list(c("Suncus murinus","Rattus tanezumi"),c("Suncus murinus","Rattus norvegicus"),
#                      c("Suncus murinus","Bandicota indica"),c("Suncus murinus","Rattus losea"),
#                      c("Suncus murinus","Rattus andamanensis"),c("Suncus murinus","Berylmys bowersi"),
#                      c("Suncus murinus","Mus caroli"),c("Suncus murinus","Niviventer lotipes"),
#                      c("Rattus tanezumi","Rattus norvegicus"),c("Rattus tanezumi","Bandicota indica"),
#                      c("Rattus tanezumi","Rattus losea"),c("Rattus tanezumi","Rattus andamanensis"),
#                      c("Rattus tanezumi","Berylmys bowersi"),c("Rattus tanezumi","Mus caroli"),
#                      c("Rattus tanezumi","Niviventer lotipes"),c("Rattus norvegicus","Bandicota indica"),
#                      c("Rattus norvegicus","Rattus losea"),c("Rattus norvegicus","Rattus andamanensis"),
#                      c("Rattus norvegicus","Berylmys bowersi"),c("Rattus norvegicus","Mus caroli"),
#                      c("Rattus norvegicus","Niviventer lotipes"),c("Bandicota indica","Rattus losea"),
#                      c("Bandicota indica","Rattus andamanensis"),c("Bandicota indica","Berylmys bowersi"),
#                      c("Bandicota indica","Mus caroli"),c("Bandicota indica","Niviventer lotipes"),
#                      c("Rattus losea","Rattus andamanensis"),c("Rattus losea","Berylmys bowersi"),
#                      c("Rattus losea","Mus caroli"),c("Rattus losea","Niviventer lotipes"),
#                      c("Rattus andamanensis","Berylmys bowersi"),c("Rattus andamanensis","Mus caroli"),
#                      c("Rattus andamanensis","Niviventer lotipes"),c("Berylmys bowersi","Mus caroli"),
#                      c("Berylmys bowersi","Niviventer lotipes"),c("Mus caroli","Niviventer lotipes"))
# 
# Richness.plot <- ggplot(all.diversity,aes(Hos_Taxo,Richness,color=Hos_Taxo))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness')+ geom_jitter(size=0.5)+
#   stat_compare_means(comparisons=list(c("Suncus murinus","Rattus norvegicus"),
#                                       c("Suncus murinus","Bandicota indica"),
#                                       c("Rattus tanezumi","Bandicota indica"),
#                                       c("Rattus tanezumi","Rattus losea"),
#                                       c("Rattus norvegicus","Rattus losea"),
#                                       c("Bandicota indica","Rattus losea"),
#                                       c("Bandicota indica","Rattus andamanensis"),
#                                       c("Bandicota indica","Mus caroli")
#   ),size=2,
#   label = "p.signif",method = 'wilcox.test')+
#   scale_color_manual(values = setNames(Host_color,Host))+
#   theme_few()+coord_flip()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
# Richness.plot
# Shannon.plot <- ggplot(all.diversity,aes(Hos_Taxo,Shannon,color=Hos_Taxo))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='Shannon')+ geom_jitter(size=0.5)+
#   stat_compare_means(comparisons=list(
#     c("Suncus murinus","Rattus tanezumi"),
#     c("Suncus murinus","Rattus norvegicus"),
#     c("Suncus murinus","Bandicota indica"),
#     c("Rattus tanezumi","Bandicota indica"),
#     c("Rattus tanezumi","Rattus losea"),
#     c("Rattus norvegicus","Bandicota indica"),
#     c("Rattus norvegicus","Rattus losea"),
#     c("Rattus norvegicus","Mus caroli"),
#     c("Bandicota indica","Rattus losea"),
#     c("Bandicota indica","Rattus andamanensis"),
#     c("Bandicota indica","Mus caroli")
#   ),size=2,
#   label = "p.signif",method = 'wilcox.test')+
#   scale_color_manual(values = setNames(Host_color,Host))+
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
# Shannon.plot

###不同宿主之间跨物种病毒之间的比例率的情况科+大类的水平
##Dae level
# cross.host <- subset(Upset_sum,RPM > 0)
# cross.num <- aggregate(cross.host$Pathogen,
#                        list(cross.host$Pathogen),length)
# colnames(cross.num) <- c('species','Num')
# cross.num <- merge(Mammal.releted,cross.num,all.x = TRUE)
# cross.num[is.na(cross.num)] <- 1
# cross.num$Cross <- ifelse(cross.num$Num>1,'Detected in Multiple Host','Detected in Single Host')
# Dae.level <- aggregate(cross.num$Cross,list(cross.num$family,cross.num$Cross),length)
# colnames(Dae.level) <- c('family','Type','Num')
# Dae.level <- Dae.level[order(Dae.level$Num,decreasing = TRUE),]
# Dae.level$family <- factor(Dae.level$family,levels = unique(Dae.level$family))
# dae.level.plot <- ggplot(data = Dae.level,aes(x=family,y=Num,fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25
#   )+theme_few()+
#   labs(x="",y='Num of pathogens')+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
#         axis.text.x=element_text(angle=60,hjust = 1))
# ggsave(
#   filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figures4.5 V2.pdf',sep=''),
#   dae.level.plot,
#   width = 180,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )
##Clade level
# clade.level <- aggregate(cross.num$Cross,list(cross.num$Type,cross.num$Cross),length)
# colnames(clade.level) <- c('family','Type','Num')
# clade.level <- clade.level[order(clade.level$Num,decreasing = TRUE),]
# clade.level$family <- factor(clade.level$family,levels = unique(clade.level$family))
# clade.level.plot <- ggplot(data = clade.level,aes(x=family,y=Num,fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25
#   )+theme_few()+
#   labs(x="",y='Num of pathogens')+coord_flip()+
#   theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),legend.position = 'top',
#         axis.text.x=element_text(angle=60,hjust = 1))
# clade.level.plot
# ggsave(
#   filename = paste0(out_path,'/fig 4c.pdf',sep=''),
#   clade.level.plot,
#   width = 60,             # 宽
#   height = 60,            # 高
#   units = "mm",          # 单位
#   dpi = 300,              # 分辨率DPI
#   limitsize = FALSE
# )
###不同宿主不同类型的病毒的多样性的情况Richness/Shannon
# Individual.wide.rpm
# rna.virus.name <- unique(subset(Mammal.releted,Type=='RNA_Virus')$species)
# dna.virus.name <- unique(subset(Mammal.releted,Type=='DNA_Virus')$species)
# bacteria.name <- unique(subset(Mammal.releted,Type=='Bac')$species)
# eukaryote.name <- unique(subset(Mammal.releted,Type=='Euk')$species)
# library(vegan)
# RNA.virus.wide <- Individual.wide.rpm[,rna.virus.name]
# Shannon <- diversity(t(RNA.virus.wide), index = "shannon", MARGIN = 2, base = exp(1))
# Richness <- specnumber(t(RNA.virus.wide), MARGIN = 2)#spe.rich =sobs
# RNAvirus.richness <- as.data.frame(cbind(Shannon, Richness))
# RNAvirus.richness$Individual <- rownames(RNAvirus.richness)
# RNAvirus.richness$Type <- 'RNA virus'
# DNA.virus.wide <- Individual.wide.rpm[,intersect(colnames(Individual.wide.rpm),dna.virus.name)]
# Shannon <- diversity(t(DNA.virus.wide), index = "shannon", MARGIN = 2, base = exp(1))
# Richness <- specnumber(t(DNA.virus.wide), MARGIN = 2)#spe.rich =sobs
# DNAvirus.richness <- as.data.frame(cbind(Shannon, Richness))
# DNAvirus.richness$Individual <- rownames(DNAvirus.richness)
# DNAvirus.richness$Type <- 'DNA virus'
# BAC.wide <- Individual.wide.rpm[,intersect(colnames(Individual.wide.rpm),bacteria.name)]
# Shannon <- diversity(t(BAC.wide), index = "shannon", MARGIN = 2, base = exp(1))
# richness <- specnumber(t(BAC.wide), MARGIN = 2)#spe.rich =sobs
# bac.richness <- as.data.frame(cbind(Shannon, Richness))
# bac.richness$Individual <- rownames(bac.richness)
# bac.richness$Type <- 'Bacteria'
# EUK.wide <- Individual.wide.rpm[,intersect(colnames(Individual.wide.rpm),eukaryote.name)]
# Shannon <- diversity(t(EUK.wide), index = "shannon", MARGIN = 2, base = exp(1))
# richness <- specnumber(t(EUK.wide), MARGIN = 2)#spe.rich =sobs
# euk.richness <- as.data.frame(cbind(Shannon, Richness))
# euk.richness$Individual <- rownames(euk.richness)
# euk.richness$Type <- 'Eukaryote'
# 
# richness.data <- rbind(bac.richness,RNAvirus.richness,DNAvirus.richness,euk.richness)
# richness.data <- merge(richness.data,Individual.infor)
# 
# richness.plot <- ggplot(richness.data,aes(Hos_Taxo,Richness,color=Hos_Taxo))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness of Tissues')+ geom_jitter(size=0.5)+
#   facet_grid(. ~ Type) +
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90))
# shannon.plot <- ggplot(richness.data,aes(Hos_Taxo,Shannon,color=Hos_Taxo))+
#   geom_boxplot(outlier.shape = NA)+labs(x="",y='Shannon of Tissues')+ geom_jitter(size=0.5)+
#   facet_grid(. ~ Type) +
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90))

###绘制不同物种里面新的病原体的数量
Individual.wide.rpm1 <- Individual.wide.rpm
Individual.wide.rpm1$Individual <- rownames(Individual.wide.rpm1)
Individual.wide.rpm.long <- gather(Individual.wide.rpm1,key = 'species',value = RPM,-'Individual')
Individual.wide.rpm.long <- subset(Individual.wide.rpm.long,RPM > 0)
Individual.wide.rpm.long <- merge(Individual.wide.rpm.long,Individual.infor)
Individual.wide.rpm.long <- merge(Individual.wide.rpm.long,Ann_Table)
pathogens.long <- unique(Individual.wide.rpm.long[,c('species','Hos_Taxo','isExit','Type')])
pathogens.long.sum <- aggregate(pathogens.long$species,
                                list(pathogens.long$Hos_Taxo,pathogens.long$isExit,pathogens.long$Type),length)
colnames(pathogens.long.sum) <- c('Host','Type','Pathogen_Type','Num')
pathogens.long.sum$Type <- ifelse(pathogens.long.sum$Type == 'Exit','exited','new')
# host.new.virus <- ggplot(data = pathogens.long.sum,aes(x=Host,y=Num,fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25
#   )+theme_few()+facet_grid(. ~ Pathogen_Type)+scale_fill_manual(values=setNames(c('#A32A31','#407BD0'),c('new','exited')))+
#   theme(legend.key.size = unit(3,'mm'),
#         text=element_text(size= 6),
#         axis.text.x=element_text(angle=90))
# host.new.virus.sum <- ggplot(data = pathogens.long.sum,aes(x=Host,y=Num,fill=Type))+
#   geom_bar(stat = "identity",width = 0.8,size = 0.25
#   )+xlab('')+theme_few()+scale_fill_manual(values=setNames(c('#A32A31','#407BD0'),c('new','exited')))+
#   theme(legend.key.size = unit(3,'mm'),
#         text=element_text(size= 6),
#         axis.text.x=element_text(angle=90,vjust = 0,hjust = 1))
###对病原丰度情况进行拟合
all.diversity1 <- subset(all.diversity,Hos_Taxo != 'Mus caroli' & Hos_Taxo != 'Niviventer lotipes' & 
                           Hos_Taxo != 'Berylmys bowersi')
all.diversity1$Hos_Taxo <- factor(all.diversity1$Hos_Taxo,levels = c('Bandicota indica','Rattus norvegicus','Rattus andamanensis','Rattus tanezumi','Suncus murinus','Rattus losea'))
Richness.plot1 <- ggplot(all.diversity1,aes(Hos_Taxo,Richness,color=Hos_Taxo))+
  geom_boxplot()+labs(x="",y='Richness')+
  scale_color_manual(values = setNames(Host_color,Host))+
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 60,hjust = 1))
##构建函数进行拟合
all.diversity1$Hos_Taxo <- factor(all.diversity1$Hos_Taxo,levels = unique(all.diversity1$Hos_Taxo))
nb_model <- MASS::glm.nb(Richness ~ Hos_Taxo,  data = all.diversity1)
# summary(nb_model)
x2.nb_model <- sum(residuals(nb_model, type = "pearson")^2)
# x2.nb_model
# qchisq(0.975, df.residual(nb_model))

species_fitted = predict(nb_model, terms = "Hos_Taxo",type = "terms", se=T)
pred_df = tibble(Individual = all.diversity1$Individual, 
                 species=all.diversity1$Hos_Taxo, 
                 sp_effect = species_fitted$fit[,1],
                 se = species_fitted$se.fit[,1]) %>%
  group_by(species)
pred_df_1 <- aggregate(list(pred_df$sp_effect,pred_df$se), by = list(pred_df$species), 
                       FUN = mean)
colnames(pred_df_1) <- c("species","sp_effect","se")
pred_df_1$species = factor(pred_df_1$species, levels=arrange(pred_df_1, desc(pred_df_1$sp_effect))$species)

Estimated.richness.plot <- ggplot(aes(x=species, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  ylab("Estimated effect size\n(number of pathogen per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =0,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))+
  coord_flip()
Estimated.richness.plot

ggsave(
  filename = paste0(out_path,'/fig 4b.pdf',sep=''),
  Estimated.richness.plot,
  width = 70,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

###判断不同类型的病原的跨物种传播的概率
# cross.host <- RPM.table.long[,c('species','RPM','family','Type','Hos_Taxo')]
# cross.host <- subset(cross.host,RPM>0)
# cross.host <- unique(cross.host[,c('species','family','Type','Hos_Taxo')])
# cross.host.sum <- aggregate(cross.host$Hos_Taxo,
#                             list(cross.host$species,cross.host$family,cross.host$Type),length)
# colnames(cross.host.sum) <- c('Species','Family','Type','Num')
# cross.host.sum$Type <- factor(cross.host.sum$Type,levels = c('Bacteria','DNA Virus', 'RNA Virus','Eukaryote'))
# bargraph.plot1 <- ggplot(cross.host.sum,aes(Type,Num,color=Type))+
#   geom_boxplot()+labs(x="",y='Num of host shared')+geom_jitter()+
#   scale_color_manual(values=setNames(Type_color,Type))+
#   theme_few()+
#   theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
#         axis.text = element_text(size = 5),axis.text.x = element_text(angle = 60,hjust = 1))
###跨类别的情况
Host.taxo.table = data.frame(
  Hos_Taxo=c("Suncus murinus","Rattus tanezumi","Rattus norvegicus","Bandicota indica","Rattus losea","Rattus andamanensis","Berylmys bowersi","Mus caroli","Niviventer lotipes"),
  Host_Family=c('Soricidae','Muridae','Muridae','Muridae','Muridae','Muridae','Muridae','Muridae','Muridae'),
  Host_Genus=c('Suncus','Rattus','Rattus','Bandicota','Rattus','Rattus','Berylmys','Mus','Niviventer'),
  Host_order=c('Eulipoyphla','Rodentia','Rodentia','Rodentia','Rodentia','Rodentia','Rodentia','Rodentia','Rodentia')
)
cross.host.taxo <- merge(cross.host,Host.taxo.table)
cross.host.taxo <- unique(cross.host.taxo[,c('Hos_Taxo','species','family','Type','Host_Family',
                                             'Host_Genus','Host_order')])
cross.species <- aggregate(cross.host.taxo$Hos_Taxo,
                           list(cross.host.taxo$species,cross.host.taxo$family,cross.host.taxo$Type),length)
colnames(cross.species) <- c('species','family','Type','Species')
cross.host.taxo <- unique(cross.host.taxo[,c('species','family','Type','Host_Family',
                                             'Host_Genus','Host_order')])
cross.genus <- aggregate(cross.host.taxo$Host_Genus,
                         list(cross.host.taxo$species,cross.host.taxo$family,cross.host.taxo$Type),length)
colnames(cross.genus) <- c('species','family','Type','genus')
cross.host.taxo <- unique(cross.host.taxo[,c('species','family','Type','Host_Family','Host_order')])
cross.family <- aggregate(cross.host.taxo$Host_Family,
                          list(cross.host.taxo$species,cross.host.taxo$family,cross.host.taxo$Type),length)
colnames(cross.family) <- c('species','family','Type','Family')
cross.host.taxo <- unique(cross.host.taxo[,c('species','family','Type','Host_order')])
cross.order <- aggregate(cross.host.taxo$Host_order,
                         list(cross.host.taxo$species,cross.host.taxo$family,cross.host.taxo$Type),length)
colnames(cross.order) <- c('species','family','Type','Order')
cross.sum <- merge(cross.genus,cross.family)
cross.sum <- merge(cross.sum,cross.order)
cross.sum <- merge(cross.sum,cross.species)
cross.sum$crossType <- ifelse(cross.sum$Order > 1,'Cross-order transmission',
                              ifelse(cross.sum$Family > 1,'Cross-family transmission',
                                     ifelse(cross.sum$genus > 1,'Cross-genus transmission',
                                            ifelse(cross.sum$Species > 1,'Cross-species transmission','Detected in Single Species'))))
cross.sum.length <- aggregate(cross.sum$species,
                              list(cross.sum$Type,cross.sum$crossType),length)
colnames(cross.sum.length) <- c('Pathogens','Type','Num')
cross.sum.length$Type <- factor(cross.sum.length$Type,levels = c('Detected in Single Species',
                                                                 'Cross-species transmission',
                                                                 'Cross-genus transmission',
                                                                 'Cross-family transmission',
                                                                 'Cross-order transmission'
))
cross.type <- c('Detected in Single Species',
                'Cross-species transmission',
                'Cross-genus transmission',
                'Cross-family transmission',
                'Cross-order transmission'
)
cross.color <- c('#9281DD','#B7617D','#E4B112','#B7DB29','#159FD7')

cross.sum.length$Pathogens <- factor(cross.sum.length$Pathogens,levels = c('Bacteria','DNA Virus','RNA Virus','Eukaryote'))
cross.plot.2 <- ggplot(data = cross.sum.length,aes(x=Pathogens,y=Num,fill=Type))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25)+
  scale_fill_manual(values = setNames(cross.color,cross.type))+
  theme_few()+coord_flip()+
  labs(x="",y='Num of Pathogens')+
  theme(legend.key.size = unit(3,'mm'),
        text=element_text(size= 6),
        axis.text.x=element_text(angle = 90,hjust = 1))
cross.plot.2

ggsave(
  filename = paste0(out_path,'/fig 4c.pdf',sep=''),
  cross.plot.2,
  width = 100,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

##estime
# nb_model <- MASS::glm.nb(Num ~ Type,  data = cross.host.sum)
# summary(nb_model)
# x2.nb_model <- sum(residuals(nb_model, type = "pearson")^2)
# x2.nb_model
# qchisq(0.975, df.residual(nb_model))
# 
# species_fitted = predict(nb_model, terms = "Type",type = "terms", se=T)
# pred_df = tibble(Individual = cross.host.sum$Species, 
#                  species=cross.host.sum$Type, 
#                  sp_effect = species_fitted$fit[,1],
#                  se = species_fitted$se.fit[,1]) %>%
#   group_by(species)
# pred_df_1 <- aggregate(list(pred_df$sp_effect,pred_df$se), by = list(pred_df$species), 
#                        FUN = mean)
# colnames(pred_df_1) <- c("species","sp_effect","se")
# pred_df_1$species = factor(pred_df_1$species, levels=arrange(pred_df_1, desc(pred_df_1$sp_effect))$species)
# 
# Infected.polt <- ggplot(aes(x=species, y=sp_effect), data=pred_df_1) +
#   geom_point(size=1.5) +
#   geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
#   geom_hline(yintercept = 0, linetype=2) +
#   theme_bw() +
#   ylab("Estimated effect size\n(number of host pathogens infected)")+
#   xlab("")+
#   theme( axis.text.x = element_text(angle =60, hjust=1,size=5),
#          axis.text.y = element_text(size=5),
#          axis.title.y=element_text(size=5))
# Infected.polt

####Positivate Rate information
positivate.data <- RPM.table.long
positivate.data$Value <- ifelse(positivate.data$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Hos_Taxo + Value + species, data=positivate.data)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)

hhhh_1 <- aggregate(hhhh$Freq, by = list(hhhh$Hos_Taxo,hhhh$species), 
                    FUN = sum)
colnames(hhhh_1) <- c('Hos_Taxo','species','Freq.Sum')

hhhh_1 <- subset(merge(hhhh,hhhh_1,all.x = TRUE),Value == "Positive")
hhhh_1 <- subset(hhhh_1,Freq >0 )
hhhh_1$Positiveta_Rate <- hhhh_1$Freq/hhhh_1$Freq.Sum
hhhh_1 <- hhhh_1[,c('Hos_Taxo','species','Positiveta_Rate')]
colnames(hhhh_1) <- c('Node','Pathogen','Positivate.Rate')
colnames(graph_long)
hhhh_1 <- merge(graph_long,hhhh_1)

relations <- hhhh_1[,c('Node','Pathogen','Positivate.Rate')]
colnames(relations) <- c('from','to','Positivate.Rate')

actors <- data.frame(
  name = unique(c(hhhh_1$Node,hhhh_1$Pathogen
  )))

actors <- actors %>% 
  mutate(Type = ifelse(name %in% graph_long$Node, 'Host',
                       ifelse(name %in% subset(Ann_Table,Type=='Bacteria')$species, 'Bacteria',
                              ifelse(name %in% subset(Ann_Table,Type=='Eukaryote')$species,'Eukaryote',
                                     ifelse(name %in% subset(Ann_Table,Type=='RNA Virus')$species,'RNA Virus',
                                            ifelse(name %in% subset(Ann_Table,Type=='DNA Virus')$species,'DNA Virus','Others')))))) %>%
  mutate(Name = ifelse(name %in% graph_long$Node, name,
                       ifelse(name %in% subset(Ann_Table,Type=='Bacteria')$species, 'Bacteria',
                              ifelse(name %in% subset(Ann_Table,Type=='Eukaryote')$species,'Eukaryote',
                                     ifelse(name %in% subset(Ann_Table,Type=='RNA Virus')$species,'RNA Virus',
                                            ifelse(name %in% subset(Ann_Table,Type=='DNA Virus')$species,'DNA Virus',"Pathogen")))))) %>%
  mutate(PathogensType = ifelse(name %in% graph_long$Node,name,
                                ifelse(name %in% subset(cross.sum,crossType == 'Cross-order transmission')$species,'Cross-order transmission',
                                       ifelse(name %in% subset(cross.sum,crossType == 'Cross-genus transmission')$species,'Cross-genus transmission',
                                              ifelse(name %in% subset(cross.sum,crossType == 'Cross-species transmission')$species,'Cross-species transmission',
                                                     ifelse(name%in% subset(cross.sum,crossType == 'Cross-family transmission')$species,'Cross-family transmission','Detected in single species')))))) %>%
  mutate(shownames = ifelse(name %in% subset(cross.sum,crossType == 'Cross-order transmission')$species,name, 
                            ifelse(Type == 'Host',name,'')))
Species.names <- unique(c(unique(graph_long$Node),c('Niviventer lotipes')))
Species.names.color <- c('#DE7833','#912C2C','#F2BB6B','#C2ABC8','#329845','#AED185','#B43970',
                         '#A695BD','#43978F')
Type <- c("Cross-order transmission","Cross-genus transmission","Cross-species transmission","Detected in single species")
color2 <- c('#ff4d4d','#EDAE92','#A4C8D9','#1D75B5')
color.panel <- setNames(c(color2,Species.names.color),c(Type,Species.names))
order <- c("Rodentia","Eulipoyphla")
order.color <- c('#F89FA8','#F9E9A4')
order.color.type <- setNames(order.color,order)
library(igraph)
g <- graph_from_data_frame(relations, directed = FALSE, vertices = actors)
V(g)$color <- color.panel[actors$PathogensType[match(V(g)$name, actors$name)]]
V(g)$degree <- degree(g)/2
V(g)$size <- sqrt(V(g)$degree)*3
Rodentia.list <- c(subset(cross.host.taxo,Host_order == 'Rodentia')$species,subset(Host.taxo.table,Host_order == 'Rodentia')$Hos_Taxo )
Eulipoyphla.list <- c(subset(cross.host.taxo,Host_order == 'Eulipoyphla')$species,subset(Host.taxo.table,Host_order == 'Eulipoyphla')$Hos_Taxo )
Rodentia.list1 <- setdiff(Rodentia.list,Eulipoyphla.list)
Eulipoyphla.list1 <- setdiff(Eulipoyphla.list,Rodentia.list)
# pdf('F:/Mouse_Result_24/Figure/TMP/Figure_4_net.pdf', width = 7, height = 5)
# plot(g,layout=layout_nicely,vertex.color = V(g)$color,vertex.size = V(g)$size,vertex.label = NA,
#      #vertex.label = V(g)$shownames,vertex.label.cex = 1,
#      edge.width = E(g)$Positivate.Rate,
#      mark.groups = list(V(g)$name[V(g)$name %in% Rodentia.list1],
#                         V(g)$name[V(g)$name %in% Eulipoyphla.list1]),
#      mark.col = rainbow(2, alpha = 0.1),
#      mark.border = rainbow(2, alpha = 1))
# legend("topright", legend = Species.names,  col = Species.names.color,   pch = 16,  title = "Host Species")
# legend("topleft", legend = Type,  col = color2,   pch = 16,  title = "Pathogens Type")
# legend("bottomleft", legend = order,  col = rainbow(2, alpha = 1),   pch = 16,  title = "Type")
# dev.off()

##创建ggrapg对象
ig_tidy <- tidygraph::as_tbl_graph(g)

gr1_layout <- ggraph::create_layout(ig_tidy,layout = 'stress')
gr1_layout$x <- 8-gr1_layout$x
gr1_layout$x <- ifelse(gr1_layout$name %in%c(subset(relations,from == 'Suncus murinus')$to,c('Suncus murinus')),
                       gr1_layout$x-2,gr1_layout$x)
gr1_layout$x <- ifelse(gr1_layout$name %in%c(subset(relations,from == 'Niviventer lotipes')$to,c('Niviventer lotipes')),
                       gr1_layout$x-2,gr1_layout$x)
gr1_layout$x <- ifelse(gr1_layout$name %in% subset(relations,from == 'Rattus tanezumi')$to,
                       ifelse(gr1_layout$PathogensType == 'Detected in single species',
                              gr1_layout$x+0.5,gr1_layout$x),gr1_layout$x)
gr1_layout$point.type <- ifelse(gr1_layout$name %in% Eulipoyphla.list1,'Eulipoyphla',
                                ifelse(gr1_layout$name %in% Rodentia.list1,'Rodentia','Others'))
gr1_layout$x <- ifelse(gr1_layout$name == 'Hepatozoon sp',gr1_layout$x+0.5,gr1_layout$x)
gr1_layout$y <- ifelse(gr1_layout$name == 'Hepatozoon sp',gr1_layout$y+0.5,gr1_layout$y)

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

test <- subset(gr1_layout,name %in% unique(subset(relations,from == 'Suncus murinus')$to))
test <- subset(test,PathogensType=='Cross-order transmission')
test <- subset(test,name != 'Brachylaima sp')
# ggplot(data = test, aes(x = x, y = y)) +
#   geom_smooth(method = "lm", se=TRUE, 
#               color="black", formula = y ~ x) +
#   geom_point()+
#   theme_bw()+geom_text(x = 1, y = 2, 
#                        label = lm_eqn(test), parse = TRUE)

test <- subset(gr1_layout,name %in% unique(subset(relations,from == 'Suncus murinus')$to))
test <- subset(test,PathogensType=='Cross-order transmission')
test <- test[order(test$y),]
#test$x <- seq(from=0.8,by=0.2,length.out=10)
test$name
test$x <- c(1.0, 2.4, 2.6, 1.6, 1.8, 1.4, 1.2,2.0 , 2.2,0.8)
test <- test %>%  rowwise() %>% mutate(y = 3.4-x)
rownames(test) <- test$name

# test['Bartonella kosoyi',]$x

gr1_layout$x <- ifelse(gr1_layout$name %in% test$name,test[gr1_layout$name,]$x,gr1_layout$x)
gr1_layout$y <- ifelse(gr1_layout$name %in% test$name,test[gr1_layout$name,]$y,gr1_layout$y)
gr1_layout$x <- ifelse(gr1_layout$name %in%c(subset(relations,from == 'Suncus murinus')$to,c('Suncus murinus')),
                       gr1_layout$x+1,gr1_layout$x)
gr1_layout$y <- ifelse(gr1_layout$name %in%c(subset(relations,from == 'Niviventer lotipes')$to,c('Niviventer lotipes')),
                       gr1_layout$y+5,gr1_layout$y)

###转换成ggplot2绘制
point_sit <- read.csv('Site.csv')
##绘制提取宿主点的信息
host.table <- data.frame(subset(gr1_layout,Type == 'Host'))
host.table <- merge(host.table,point_sit)
##提取病原体的信息
pathogens.table <- data.frame(subset(gr1_layout,Type != 'Host'))
pathogens.table <- merge(pathogens.table,point_sit)
##提取边的信息
host.table1 <- host.table[,c('new_x','new_y','name')]
colnames(host.table1) <- c('xstart','ystart','from')

pathogens.table1 <- pathogens.table[,c('new_x','new_y','name')]
colnames(pathogens.table1) <- c('xend','yend','to')
linkage.table <- merge(relations,host.table1,all.x = TRUE)
linkage.table <- merge(linkage.table,pathogens.table1,all.x = TRUE)

gr1_layout1 <- merge(gr1_layout,point_sit)

cros <- ggplot(gr1_layout1,aes(x=new_x,y=new_y))+
  ggforce::geom_mark_hull(aes(fill = point.type, filter = point.type != 'Others'),concavity = 3)+
  geom_segment(aes(x = xstart, y = ystart, xend = xend, yend = yend,linewidth=Positivate.Rate), 
               color='gray',alpha=0.5,show.legend=TRUE,
               data = linkage.table)+
  geom_point(aes(x = new_x, y = new_y,shape=Type,size=2,fill=PathogensType,color=PathogensType),show.legend=TRUE,
             data = host.table)+
  geom_point(aes(x = new_x, y = new_y, shape=Type,size=1,fill=PathogensType,color=PathogensType),show.legend=TRUE,
             data = pathogens.table)+
  scale_shape_manual(values = c('Host'=23,'Bacteria'=21,'Eukaryote'=22,
                                'RNA Virus'=24,'DNA Virus'=25))+
  geom_text(aes(x = new_x, y = new_y,label=shownames),size.unit = 'pt',hjust=-0.5,
            data = pathogens.table)+
  geom_text(aes(x = new_x, y = new_y,label=shownames),size.unit = 'pt',hjust=-0.5,
            data = host.table)+
  scale_color_manual(values = color.panel)+
  scale_fill_manual(values = color.panel)+
  theme_void()+
  theme(text=element_text(size= 5),legend.position = 'right',legend.key.size = unit(3,'mm'))

ggsave(filename = paste0(out_path,'/fig 4e.pdf',sep=''),
       cros,
       width = 180,             # 宽
       height = 120,            # 高
       units = "mm",          # 单位
       dpi = 300,              # 分辨率DPI
       limitsize = FALSE)

###figure
P1 <- plot_grid(richness.plot.sum+theme(legend.position = 'bottom'), 
                Richness.plot1+theme(legend.position = "none"),
                Estimated.richness.plot,
                align = "h", axis = "bt",nrow = 1, rel_widths = c(2,1,1))+
  draw_plot_label(
    c("a", "b"),
    c(0, 0.5),
    c(1, 1),
    size = 8
  )
P1

P2 <- plot_grid(cross.plot.2+theme(legend.position = 'bottom'),Upset_Plot+theme(legend.position = "none"),
                axis = "bt",nrow = 1, rel_widths = c(1,2))+
  draw_plot_label(
    c("e", "f"),
    c(0, 1/3),
    c(1, 1),
    size = 8
  )
P2
P3 <- plot_grid(P1,P2,nrow = 2,align = "h", axis = "bt")
P3


P3 <- plot_grid(P3,cros,nrow = 2,align = "h", axis = "bt")
P3

ggsave(
  filename = paste0(out_path,'/fig 4.pdf',sep=''),
  P3,
  width = 180,             # 宽
  height = 240,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
