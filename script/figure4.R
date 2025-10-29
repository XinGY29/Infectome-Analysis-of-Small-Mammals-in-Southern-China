##figure 4a
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','data', sep = "/")
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

###figure 4a
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
            "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
richness.long1 <- read.csv(paste0(data_path,'/fig4.richness.long1.csv',sep=''),check.names = F,row.names = 1)
# richness.long1 <- richness.long1[order(richness.long1$Type),]
# richness.long1 <- richness.long1[order(richness.long1$Individual),]
# richness.long1 <- richness.long1[order(richness.long1$Hos_Taxo),]
host.num <- aggregate(richness.long1$Individual,by = list(richness.long1$Host),FUN = length)
colnames(host.num) <- c('host','num')
host.num <- host.num[order(host.num$num,decreasing = TRUE),]

richness.num <- aggregate(richness.long1$Num,by = list(richness.long1$Individual,richness.long1$Host),FUN = sum)
colnames(richness.num) <- c('Individual','host','num')
richness.num$host <- factor(richness.num$host,levels = host.num$host)
richness.num <- richness.num[order(richness.num$num,decreasing = TRUE),]
richness.num <- richness.num[order(richness.num$host),]
richness.long1$Individual <- factor(richness.long1$Individual,levels = richness.num$Individual)
Host <- c("Suncus murinus", "Rattus tanezumi", "Rattus norvegicus", "Bandicota indica", "Rattus losea", "Rattus andamanensis", "Berylmys bowersi", "Mus caroli", "Niviventer lotipes")
Host_color <- colors[1:length(Host)]

fig4a <- ggplot(data = richness.long1,aes(x=Individual,y=Num))+
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
fig4a
ggsave(filename = paste0(out_path,'/fig 4a.pdf',sep=''), plot = fig4a, width = 120, height = 60,units = 'mm')

###figure 4d
cross.sum.length <- read.csv(paste0(data_path,'/fig4.d.csv',sep=''),check.names = F,row.names = 1)
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
ggsave(filename = paste0(out_path,'/fig 4d.pdf',sep=''), plot = cross.plot.2, width = 100, height = 60,units = 'mm')
###figure 4e
df.long <- read.csv(paste0(data_path,'/fig4.e2.csv',sep=''),check.names = F,row.names = 1)
heatmap <- ggplot(df.long,aes(x=x,y=key,fill = log2(value)))+
  geom_point(shape=21,size=8)+scale_fill_gradient2(low = "red",mid = "white", high = "blue", midpoint = log2(0.05))+
  theme_bw()+
  labs(x="",y='')+
  theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
        axis.text.x=element_text(angle=90,hjust = 1))

# Order.positivate.wide$Type <- factor(Order.positivate.wide$Type,levels = rev(unique(df.long$x)))
Order.positivate.wide <- read.csv(paste0(data_path,'/fig4.e.csv',sep=''),check.names = F,row.names = 1)
bar.plot <- ggplot(Order.positivate.wide,aes(Type))+geom_bar(aes(weight = rate),fill="#3C5488FF")+theme_bw()+
  xlab('')+ylab('Cross Order Rate')+
  theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
        axis.text.x=element_text(angle=90,hjust = 1))+coord_flip()+guides(fill = F)+scale_y_reverse()
fig4e <- bar.plot+heatmap
ggsave(filename = paste0(out_path,'/fig 4e.pdf',sep=''), plot = fig4e, width = 140, height = 60,units = 'mm')

###figure 4b
plot_data <- read.csv(paste0(data_path,'/fig4.b.dev_prop_shared_virus.csv',sep=''),check.names = F,row.names = 1)
color = c("Host evolutionary distance" = "#DB432C",
          "Climate" = "#838AAF",
          'Spatial Distance'='#C4B797',
          "Unexplained" = "gray")
plot_data <- plot_data[order(plot_data$dev_exp,decreasing = TRUE),]
legend_labels <- plot_data$var
names(legend_labels) <- plot_data$var
Richness.pie <- ggplot(plot_data) +
  geom_bar(aes(x="", y=-dev_exp, fill=var), color="black", stat="identity", width=1) +
  coord_polar("y", start= pi*0/180) +
  scale_fill_manual(values =color ,name="Variable", 
                    labels=legend_labels) + 
  #ggtitle('Explain for Richness')+
  theme(legend.position = "right", 
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        legend.key.size = unit(2, 'mm'),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=7,hjust=0.5))
ggsave(
  filename = paste0(out_path,'/fig 4b.pdf',sep=''),
  Richness.pie,
  width = 90,             # 宽
  height = 60,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)


###figure 4c
position <- read.csv(paste0(data_path,'/Site.csv',sep=''),check.names = F,row.names = 1)
networkdata <- read.csv(paste0(data_path,'/fig4c.csv',sep=''),check.names = F)

gr1_layout1 <- read.csv(paste0(data_path,'/fig4c1.gr1_layout1.csv',sep=''),check.names = F,row.names = 1)
linkage.table <- read.csv(paste0(data_path,'/fig4c1.linkage.table.csv',sep=''),check.names = F,row.names = 1)
host.table <- read.csv(paste0(data_path,'/fig4c1.host.table.csv',sep=''),check.names = F,row.names = 1)
pathogens.table <- read.csv(paste0(data_path,'/fig4c1.pathogens.table.csv',sep=''),check.names = F,row.names = 1)

Species.names <- unique(c(unique(host.table$name),c('Niviventer lotipes')))
Species.names.color <- c('#DE7833','#912C2C','#F2BB6B','#C2ABC8','#329845','#AED185','#B43970',
                         '#A695BD','#43978F')
Type <- c("Cross-order transmission","Cross-genus transmission","Cross-species transmission","Detected in single species")
color2 <- c('#ff4d4d','#EDAE92','#A4C8D9','#1D75B5')
color.panel <- setNames(c(color2,Species.names.color),c(Type,Species.names))

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

ggsave(
  filename = paste0(out_path,'/fig 4c.pdf',sep=''),
  cros,
  width = 180,             # 宽
  height = 120,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
