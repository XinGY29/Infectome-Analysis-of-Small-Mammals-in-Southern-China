###figure 2
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
RPM_Table <- paste(script_dir, '..','data', "fig2_RPM_Table.csv", sep = "/")
positive_rate <- paste(script_dir, '..','data','fig2_Positive_num.csv', sep = "/")
Positive_num <- paste(script_dir, '..','data','fig2_Positive_Table.csv', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")

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
library(circlize)
library(ComplexHeatmap)
library(dplyr)

color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d",
            "#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c",
            "#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d",
            "#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3",
            "#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d",
            "#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b",
            "#bd5975","#312520")
positivate.num <- read.csv(positive_rate,check.names = F,row.names = 1)
positivate.num <- subset(positivate.num, family != '')
positivate.num$Type <- factor(positivate.num$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
positivate.num$family <- factor(positivate.num$family, levels =c(
  'Arteriviridae','Hepeviridae','Arenaviridae','Astroviridae','Coronaviridae','Flaviviridae','Picornaviridae',
  'Hantaviridae','Paramyxoviridae','Sedoreoviridae','Caliciviridae','Parvoviridae','Adenoviridae',
  'Bartonellaceae','Chlamydiaceae','Enterobacteriaceae','Coriobacteriaceae','Pneumocystaceae','Ancylostomatidae',
  'Tritrichomonadidae','Trypanosomatidae','Strongyloididae','Hexamitidae','Cryptosporidiidae','Eimeriidae',
  'Brachylaimidae','Hepatozoidae','Oxyuridae','Hypotrichomonadidae','Retortamonadidae','Oxymonadida_Unclassified',
  'Davaineidae','Babesiidae','Schistosomatidae'))

PosiR_table <- read.csv(Positive_num,check.names = F,row.names = 1)
table.max1 <- read.csv(RPM_Table,check.names = F,row.names = 1)
table.max1 <- as.matrix(table.max1)

Family <- unique(positivate.num$family)
Family_color <- color37[1:length(Family)]
Type <- unique(positivate.num$Type)
Type_color <- c('#DB432C','#438870','#838AAF','#C4B797')

###Order
positivate.num <- positivate.num[order(positivate.num$meanrpm,decreasing = TRUE),]
positivate.num <- positivate.num[order(positivate.num$family),]
positivate.num <- positivate.num[order(positivate.num$Type),]
PosiR_table <- PosiR_table[positivate.num$species,]

row_Ann <- rowAnnotation(
  'Pathogens Type' = positivate.num$Type,
  Family = positivate.num$family,
  col = list(Family=setNames(Family_color,Family),'Pathogens Type'=setNames(Type_color,Type)),
  'Positive rate' = anno_barplot(positivate.num$Positive_Rate,add_numbers = FALSE,
                                 column_names_gp= gpar(fontsize = 5)),
  'log10(RPM+1)' = anno_boxplot(log10(table.max1[positivate.num$species,]+1)),
  annotation_name_gp= gpar(fontsize = 5),simple_anno_size = unit(0.2, "cm")
)

col_fun = colorRamp2(c(0, 90), c("white", "#B22222"))
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
pdf(paste0(out_path,'/fig 2.pdf',sep=''), width = 10, height = 10)

hhh<-Heatmap(PosiR_table, name = "Positive rate",col = col_fun, cluster_rows=F,
             cluster_columns = F,
             row_split = positivate.num[,6],column_split = rep(c('Host','Tissue','Region','Time'),c(9,3,9,5)),
             row_title = NULL,column_title = NULL,rect_gp = gpar(col = 'gray', lwd = 0.1),
             width  = unit(60, "mm"), height = unit(150, "mm"),
             row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
             # top_annotation = col_Ann,
             right_annotation = row_Ann,
             show_heatmap_legend = TRUE
)
# lgd_list = list(
#   Legend(col_fun = col_fun, title = "Positivity rate",labels_gp = gpar(fontsize = 5),
#          direction = "vertical",title_gp = gpar(fontsize = 6)),
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
draw(hhh,merge_legend = TRUE)
# draw(hhh,merge_legend = FALSE,annotation_legend_list = lgd_list)
dev.off()

