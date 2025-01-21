####figure2
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','raw_data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")

setwd(data_path)
install.packages('ggthemes')
install.packages('tidyr')
install.packages('ggplot2')
install.packages('pairwiseAdonis')
install.packages('ggpubr')
install.packages('ggpubr')
install.packages('RColorBrewer')
install.packages('pheatmap')
install.packages('ggplot2')
install.packages("circlize")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ggthemes)
library(tidyr)
library(ggplot2)
library(pairwiseAdonis)
library(ggpubr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(dplyr)

##Read Pathogen file
Ann_Table <- read.csv('Virus_Name_Anno.csv')
Ann_row <- unique(Ann_Table[,c("classify",'Level')])
Ann_Table <- Ann_Table[,c("Taxo",'classify')]
Ann_col <- read.csv('Time_Region_Host_Pathogen_Mean_Annotation.csv')
Ann_col <- Ann_col[,-1]
rownames(Ann_col) <- Ann_col[,1]
rownames(Ann_row) <- Ann_row[,1]

###colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
            "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
color20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d",
            "#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c",
            "#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d",
            "#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3",
            "#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d",
            "#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b",
            "#bd5975","#312520")

Host <- unique(Ann_col$Host)
Host_color <- colors[1:length(Host)]
Region <- unique(Ann_col$Region)
Region_color <- color20[1:length(Region)]
Time <- unique(Ann_col$Time)
Time_color <- color20[(length(Region)+1):(length(Region)+length(Time))]
Type <- unique(Ann_row$Level)
Type_color <- c("#91D1C2FF","#DC0000FF","#8491B4FF","#B09C85FF")
Family <- unique(Ann_row$classify)
Family_color <- color37[1:length(Family)]
anno_color <- list(Host=setNames(Host_color,Host),
                   Time=setNames(Time_color,Time),
                   Region=setNames(Region_color,Region),
                   Level=setNames(Type_color,Type),
                   classify=setNames(Family_color,Family)
)

##绘制阳性率的图片
individual <- read.csv('Full_Tissue_Individual.csv')
#length(unique(individual$Individual))
Positive_rate <- read.csv('Pathogen_Positive_Rate.csv')
Ann_Table <- read.csv('Virus_Name_Anno.csv')

Positive_rate_selected <- subset(Positive_rate,Taxo %in% Ann_Table$Taxo)
Type <- c("Bac","Euk","RNA Virus","DNA Virus")
Type_color <- c("#91D1C2FF","#DC0000FF","#8491B4FF","#B09C85FF")

Ann_Table$Level <- factor(Ann_Table$Level,levels = c("RNA Virus","DNA Virus","Bac","Euk"))
Ann_Table <- Ann_Table[order(Ann_Table$classify),]
Ann_Table <- Ann_Table[order(Ann_Table$Level),]

####统计在不同组织内部的病原的丰度情况
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

tissue.table.long <- gather(tissue.table,key = "Pathogens",value = RPM, -"Library")
tissue.table.long.taxo <- merge(tissue.table.long,Library_infor[,c('Library','Tissue')])
tissue.table.long.taxo$Pathogens <- ifelse(tissue.table.long.taxo$Pathogens %in% subset(
  Mammal.related,family=='Picobirnaviridae')$species,'Picobirnaviridae',tissue.table.long.taxo$Pathogens)

tissue.max <- aggregate(tissue.table.long.taxo$RPM,list(tissue.table.long.taxo$Tissue,
                                                        tissue.table.long.taxo$Pathogens),max)
colnames(tissue.max) <- c("Tissue","Pathogens","RPM")
tissue.mean <- aggregate(tissue.table.long.taxo$RPM,list(tissue.table.long.taxo$Tissue,
                                                        tissue.table.long.taxo$Pathogens),mean)
colnames(tissue.mean) <- c("Tissue","Pathogens","RPM")
tissue.sum <- merge(tissue.max,tissue.mean)
tissue.table.long.taxo <- merge(tissue.table.long,Library_infor[,c('Library','Individual','Tissue')])
tissue.table.long.taxo$Pathogens <- ifelse(tissue.table.long.taxo$Pathogens %in% subset(
  Mammal.related,family=='Picobirnaviridae')$species,'Picobirnaviridae',tissue.table.long.taxo$Pathogens)

tissue.table.long.indi <- aggregate(tissue.table.long.taxo$RPM,list(tissue.table.long.taxo$Individual,
                                                                    tissue.table.long.taxo$Pathogens),sum)
colnames(tissue.table.long.indi) <- c('Individual','Pathogens','RPM')
##Region
Time.table <- read.csv('Full_Tissue_Individual.csv')
tissue.table.long.Region <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Region')]))
tissue.table.long.Region <- merge(tissue.table.long.Region,Time.table)
tissue.table.long.Region$Collected_Region <- ifelse(tissue.table.long.Region$Collected_Region == 'AP','Anpu',
                        ifelse(tissue.table.long.Region$Collected_Region == 'ZJ','Zhanjiang',
                               ifelse(tissue.table.long.Region$Collected_Region == 'DB','Maoming',
                                      ifelse(tissue.table.long.Region$Collected_Region == 'HY','Heyuan',
                                             ifelse(tissue.table.long.Region$Collected_Region == 'JY','Jieyang',
                                                    ifelse(tissue.table.long.Region$Collected_Region == 'SG','Shaoguan',
                                                           ifelse(tissue.table.long.Region$Collected_Region == 'SS','Foshan',
                                                                  ifelse(tissue.table.long.Region$Collected_Region == 'SZ','Shenzhen',
                                                                         ifelse(tissue.table.long.Region$Collected_Region == 'YF','Yunfu',''
                                                                         )))))))))
region.max <- aggregate(tissue.table.long.Region$RPM,list(tissue.table.long.Region$Collected_Region,
                                                          tissue.table.long.Region$Pathogens),max)
colnames(region.max) <- c('Region','Pathogens','RPM')
###Time
tissue.table.long.Time <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Time')]))
tissue.table.long.Time <- merge(tissue.table.long.Time,Time.table)
time.max <- aggregate(tissue.table.long.Time$RPM,list(tissue.table.long.Time$Collected_Time,
                                                      tissue.table.long.Time$Pathogens),max)
colnames(time.max) <- c('Time','Pathogens','RPM')

##Host
tissue.table.long.Host <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Hos_Taxo')]))
host.max <- aggregate(tissue.table.long.Host$RPM,list(tissue.table.long.Host$Hos_Taxo,
                                                      tissue.table.long.Host$Pathogens),max)
colnames(host.max) <- c('Host','Pathogens','RPM')

#Max Table
tissue.wide <- spread(tissue.max,key = 'Tissue',value = RPM,fill = 0)
rownames(tissue.wide) <- tissue.wide$Pathogens
tissue.wide <- tissue.wide[,-1]
time.wide <- spread(time.max,key = 'Time',value = RPM,fill = 0)
rownames(time.wide) <- time.wide$Pathogens
time.wide <- time.wide[,-1]
region.wide <- spread(region.max,key = 'Region',value = RPM,fill = 0)
rownames(region.wide) <- region.wide$Pathogens
region.wide <- region.wide[,-1]
host.wide <- spread(host.max,key = 'Host',value = RPM,fill = 0)
rownames(host.wide) <- host.wide$Pathogens
host.wide <- host.wide[,-1]
max_table <- cbind(host.wide,tissue.wide,region.wide,time.wide)
##Positivate Rate Sum
tissue.table.long.indi.selected <- subset(tissue.table.long.indi,RPM > 0)

positivate.num <- aggregate(tissue.table.long.indi.selected$Individual,
                            list(tissue.table.long.indi.selected$Pathogens),length)
colnames(positivate.num) <- c('Taxo','Freq')
positivate.num$Samples <- length(unique(tissue.table.long.indi.selected$Individual))

positivate.num$Positive_Rate <- positivate.num$Freq/positivate.num$Samples*100
colnames(positivate.num)[1] <- 'species'
rownames(positivate.num) <- positivate.num$species
###Mamm-related Pathogens
Mammal.related$Type <- factor(Mammal.related$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
Mammal.related <- Mammal.related[order(Mammal.related$family),]
Mammal.related <- Mammal.related[order(Mammal.related$Type),]
Mammal.related[nrow(Mammal.related) + 1,] = c('Picobirnaviridae', 'Picobirnaviridae', 'RNA Virus')
positivate.num <- merge(positivate.num,Mammal.related)
positivate.num$Type <- factor(positivate.num$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
positivate.num <- positivate.num[order(positivate.num$Positive_Rate,decreasing = T),]
positivate.num$family <- factor(positivate.num$family,
                                levels = unique(positivate.num$family))
positivate.num <- positivate.num[order(positivate.num$family),]
positivate.num <- positivate.num[order(positivate.num$Type),]
positivate.num$family <- factor(positivate.num$family,
                                levels = unique(positivate.num$family))
positivate.num <- unique(positivate.num)
rownames(positivate.num)  <- positivate.num$species
max_table <- max_table[positivate.num$species,]

positivate.num <- positivate.num[intersect(positivate.num$species,rownames(max_table)),]

# RPM in echo sample
rpm.table <- spread(tissue.table.long.indi,key = Individual,value = RPM,fill = 0)
rownames(rpm.table) <- rpm.table$Pathogens
table.max <- data.matrix(rpm.table[,-1])
table.max <- table.max[intersect(Mammal.related$species,rownames(max_table)),]
table.max <- table.max[unique(rownames(positivate.num)),]

##Pisitivate Table
##Tissue
head(tissue.table.long.taxo)
tissue.table.full <- aggregate(tissue.table.long.taxo$RPM,
                               list(tissue.table.long.taxo$Library,
                                    tissue.table.long.taxo$Pathogens),sum)
colnames(tissue.table.full) <- c('Library','Pathogens','RPM')
tissue.table.wide <- spread(tissue.table.full,key = 'Pathogens',
                             value = 'RPM',fill = 0)
tissue.table.full <- merge(Library_infor,tissue.table.wide,all.x = T)
tissue.table.full <- tissue.table.full[-1,]
tissue.table.full[is.na(tissue.table.full)] <- 0
tissue.table.full.long <- gather(subset(tissue.table.full,select = -c(Hos_Taxo,RPM,Individual,Region,Time)),
                                 key = "Pathogens",value = RPM, -"Library",-"Tissue")
tissue.table.full.long$Value <- ifelse(tissue.table.full.long$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Tissue + Pathogens + Value, data=tissue.table.full.long)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh1 <- aggregate(hhhh$Freq,list(hhhh$Tissue,hhhh$Pathogens),sum)
colnames(hhhh1) <- c('Tissue','Pathogens','Sum')
hhhh <- subset(merge(hhhh,hhhh1),Value=='Positive')
hhhh$Positiveta_Rate <- hhhh$Freq/hhhh$Sum
##Tissue
tissue.positivate.wide <- spread(hhhh[,c('Tissue','Pathogens','Positiveta_Rate')],key = 'Tissue',value = Positiveta_Rate,fill = 0)
rownames(tissue.positivate.wide) <- tissue.positivate.wide$Pathogens
tissue.wide <- tissue.positivate.wide[,-1]

###Time
tissue.table.long.Time <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Time')]))
tissue.table.long.Time <- merge(tissue.table.long.Time,Time.table)

tissue.table.long.Time$Value <- ifelse(tissue.table.long.Time$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Collected_Time + Pathogens + Value, data=tissue.table.long.Time)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh1 <- aggregate(hhhh$Freq,list(hhhh$Collected_Time,hhhh$Pathogens),sum)
colnames(hhhh1) <- c('Collected_Time','Pathogens','Sum')
hhhh <- subset(merge(hhhh,hhhh1),Value=='Positive')
hhhh$Positiveta_Rate <- hhhh$Freq/hhhh$Sum

time.positivate.wide <- spread(hhhh[,c('Collected_Time','Pathogens','Positiveta_Rate')],
                               key = 'Collected_Time',value = Positiveta_Rate,fill = 0)
rownames(time.positivate.wide) <- time.positivate.wide$Pathogens
time.wide <- time.positivate.wide[,-1]

##Region
tissue.table.long.Region <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Region')]))
tissue.table.long.Region <- merge(tissue.table.long.Region,Time.table)
tissue.table.long.Region$Collected_Region <- ifelse(tissue.table.long.Region$Collected_Region == 'AP','Anpu',
                                                    ifelse(tissue.table.long.Region$Collected_Region == 'ZJ','Zhanjiang',
                                                           ifelse(tissue.table.long.Region$Collected_Region == 'DB','Maoming',
                                                                  ifelse(tissue.table.long.Region$Collected_Region == 'HY','Heyuan',
                                                                         ifelse(tissue.table.long.Region$Collected_Region == 'JY','Jieyang',
                                                                                ifelse(tissue.table.long.Region$Collected_Region == 'SG','Shaoguan',
                                                                                       ifelse(tissue.table.long.Region$Collected_Region == 'SS','Foshan',
                                                                                              ifelse(tissue.table.long.Region$Collected_Region == 'SZ','Shenzhen',
                                                                                                     ifelse(tissue.table.long.Region$Collected_Region == 'YF','Yunfu',''
                                                                                                     )))))))))
tissue.table.long.Region$Value <- ifelse(tissue.table.long.Region$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Collected_Region + Pathogens + Value, data=tissue.table.long.Region)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh1 <- aggregate(hhhh$Freq,list(hhhh$Collected_Region,hhhh$Pathogens),sum)
colnames(hhhh1) <- c('Collected_Region','Pathogens','Sum')
hhhh <- subset(merge(hhhh,hhhh1),Value=='Positive')
hhhh$Positiveta_Rate <- hhhh$Freq/hhhh$Sum

region.positivate.wide <- spread(hhhh[,c('Collected_Region','Pathogens','Positiveta_Rate')],
                               key = 'Collected_Region',value = Positiveta_Rate,fill = 0)
rownames(region.positivate.wide) <- region.positivate.wide$Pathogens
region.wide <- region.positivate.wide[,-1]


##Host
tissue.table.long.Host <- merge(tissue.table.long.indi,unique(Library_infor[c('Individual','Hos_Taxo')]))
tissue.table.long.Host$Value <- ifelse(tissue.table.long.Host$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Hos_Taxo + Pathogens + Value, data=tissue.table.long.Host)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh1 <- aggregate(hhhh$Freq,list(hhhh$Hos_Taxo,hhhh$Pathogens),sum)
colnames(hhhh1) <- c('Hos_Taxo','Pathogens','Sum')
hhhh <- subset(merge(hhhh,hhhh1),Value=='Positive')
hhhh$Positiveta_Rate <- hhhh$Freq/hhhh$Sum

host.positivate.wide <- spread(hhhh[,c('Hos_Taxo','Pathogens','Positiveta_Rate')],
                                 key = 'Hos_Taxo',value = Positiveta_Rate,fill = 0)
rownames(host.positivate.wide) <- host.positivate.wide$Pathogens
host.wide <- host.positivate.wide[,-1]

##merge
PosiR_table <- cbind(host.wide,tissue.wide[rownames(host.wide),],region.wide,time.wide)
PosiR_table <- PosiR_table*100

PosiR_table <- PosiR_table[intersect(Mammal.related$species,rownames(PosiR_table)),]
PosiR_table <- PosiR_table[unique(rownames(positivate.num)),]
PosiR_table$Pathgen <- rownames(PosiR_table)


PosiR_table <- PosiR_table %>%
  filter(!Pathgen %in% c('Picobirnaviridae','Guangdong rodent aleptorquevirus',
                         'Guangdong rodent anellovirus 3',
                         'Guangdong rodent torque teno virus 1',
                         'Guangdong rodent anellovirus 2',
                         'Guangdong rodent anellovirus 1',
                         'Rodent torque teno virus 3',
                         'Guangdong rodent torque teno virus 3'))

PosiR_table <- subset(PosiR_table,select = -c(`Pathgen`))

library(ComplexHeatmap)
col_Ann <- HeatmapAnnotation(
  Type = rep(c('Host','Tissue','Region','Time'),c(9,3,9,5)),
  col = list(Type=c('Host'='#55B7E6','Tissue'='#193E8F',
                    'Region'='#E53528','Time'='#F09739')),
  annotation_name_gp= gpar(fontsize = 5),simple_anno_size = unit(0.2, "cm"))
positivate.num <- subset(positivate.num,family != 'Anelloviridae')
positivate.num <- subset(positivate.num,family != 'Picobirnaviridae')

names <- setdiff(rownames(table.max),c('Picobirnaviridae','Guangdong rodent aleptorquevirus',
                                      'Guangdong rodent anellovirus 3',
                                      'Guangdong rodent torque teno virus 1',
                                      'Guangdong rodent anellovirus 2',
                                      'Guangdong rodent anellovirus 1',
                                      'Rodent torque teno virus 3',
                                      'Guangdong rodent torque teno virus 3'))
table.max1 <- table.max[names,]

Family <- unique(positivate.num$family)
Family_color <- color37[1:length(Family)]
Type <- unique(positivate.num$Type)
Type_color <- c('#DB432C','#438870','#838AAF','#C4B797')
row_Ann <- rowAnnotation(
  'Pathogens Type' = positivate.num$Type,
  Family = positivate.num$family,
  col = list(Family=setNames(Family_color,Family),'Pathogens Type'=setNames(Type_color,Type)),
  'Positive rate' = anno_barplot(positivate.num$Positive_Rate,add_numbers = FALSE,
                                 column_names_gp= gpar(fontsize = 5)),
  'log10(RPM+1)' = anno_boxplot(log10(table.max1+1)),
  annotation_name_gp= gpar(fontsize = 5),simple_anno_size = unit(0.2, "cm")
)

col_fun = colorRamp2(c(0, 6), c("white", "red"))
ht_opt(heatmap_column_names_gp = gpar(fontsize = 5), 
       heatmap_column_title_gp = gpar(fontsize = 5),
       heatmap_row_title_gp = gpar(fontsize = 5),
       heatmap_column_title_gp = gpar(fontsize = 5),
       legend_title_gp = gpar(fontsize = 5),
       legend_labels_gp = gpar(fontsize = 5),
       legend_grid_width = unit(3, "mm"),
       legend_grid_height =unit(3, "mm"),
       legend_gap=unit(2,'mm'),
       merge_legends=TRUE
)
PosiR_table <- as.matrix(PosiR_table)
hhh<-Heatmap(log2(PosiR_table +1), name = "log2(Positive rate % +1)",col = col_fun, cluster_rows=F,
             cluster_columns = F,
             row_split = positivate.num[,6],column_split = rep(c('Host','Tissue','Region','Time'),c(9,3,9,5)),
             row_title = NULL,column_title = NULL,rect_gp = gpar(col = 'gray', lwd = 0.1),
             width  = unit(60, "mm"), height = unit(150, "mm"),
             row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
             top_annotation = col_Ann,right_annotation = row_Ann
)
hhh
# lgd_list = list(
#   Legend(col_fun = col_fun, title = "log2(Positive rate % +1)",labels_gp = gpar(fontsize = 5),
#          direction = "horizontal",title_gp = gpar(fontsize = 6)),
#   Legend(labels = c('Host','Tissue','Region','Time'),
#          grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
#          legend_gp = gpar(fill = c('Host'='#55B7E6','Tissue'='#193E8F',
#                                    'Region'='#E53528','Time'='#F09739')), title = "Type",
#          labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 6)),
#   Legend(labels = Family,
#          legend_gp = gpar(fill = setNames(Family_color,Family)),
#          grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
#          title = "Family",
#          labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 6)),
#   Legend(labels = Type,
#          legend_gp = gpar(fill = setNames(Type_color,Type)),
#          grid_height = unit(2.5, "mm"), grid_width = unit(2.5, "mm"),
#          title = "Pathogens Type",
#          labels_gp = gpar(fontsize = 5),title_gp = gpar(fontsize = 6))
# )
pdf(paste0(out_path,'/fig 2.pdf',sep=''),height=7,width=10)
draw(hhh,merge_legend = FALSE)
# draw(hhh,merge_legend = FALSE,annotation_legend_list = lgd_list)
dev.off()

