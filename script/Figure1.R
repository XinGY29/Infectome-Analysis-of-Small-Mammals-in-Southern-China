####figure_1
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
Geo_path <- paste(script_dir, '..','Geo_data', "中国.geojson", sep = "/")
data_path <- paste(script_dir, '..','raw_data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")
##maps
install.packages("tidyverse")
library(tidyverse)
install.packages("sf")
library(sf)
install.packages("ggspatial")
library(ggspatial)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("sjmisc")
library(sjmisc)
install.packages("tidyr")
library(tidyr)
install.packages("scatterpie")
library(scatterpie)
china_prov =  read_sf(Geo_path)
Guangdong = subset(china_prov,province=='广东省')

Host <- read.csv(paste0(data_path,'/Individual_Frequence.csv',sep=''),check.names = F)
Host_long <- gather(Host,key = 'Taxo',value = 'Num',-'Region')
Host_Taxo <- read.csv(paste0(data_path,'/Region_LongLat.csv',sep=''))
Host_Taxo$Collected_Region <- ifelse(Host_Taxo$Collected_Region == 'AP',"Anpu",
                                    ifelse(Host_Taxo$Collected_Region == 'DB',"Maoming",
                                           ifelse(Host_Taxo$Collected_Region == 'HY',"Heyuan",
                                                  ifelse(Host_Taxo$Collected_Region == 'JY',"Jieyang",
                                                         ifelse(Host_Taxo$Collected_Region == 'SG',"Shaoguan",
                                                                ifelse(Host_Taxo$Collected_Region == 'SS',"Foshan",
                                                                       ifelse(Host_Taxo$Collected_Region == 'SZ',"Shenzhen",
                                                                              ifelse(Host_Taxo$Collected_Region == 'YF',"Yunfu",
                                                                                     ifelse(Host_Taxo$Collected_Region == 'ZJ',"Zhanjiang",'')))))))))

Region <- unique(Host_Taxo[,c("Region","Collected_Region")])
longlat <- aggregate(list(Host_Taxo$Longtitude,Host_Taxo$Latitude),list(Region = Host_Taxo$Collected_Region),mean)
colnames(longlat) <- c("Collected_Region",'long',"lat")

selected <- merge(Host_long,Region)
selected <- merge(selected,longlat)
selected_long <- aggregate(selected$Num,
                           list(selected$Collected_Region,selected$Taxo,selected$long,selected$lat),sum)
colnames(selected_long) <-c("Collected_Region","Taxo","long",'lat','Num')
selected_wide <- spread(selected_long,key='Taxo',value = 'Num',fill = 0)
selected_wide$total <- rowSums(selected_wide[, 5:length(colnames(selected_wide))])
setwd(out_path)
names<-unique(selected_long$Taxo)
color = c("Bandicota indica"= "#2C6344",
          "Berylmys bowersi"= "#A4C97C",
          "Rattus tanezumi"= "#308192",
          "Niviventer lotipes"= "#AED2E2",
          "Rattus andamanensis"= "#C74D26",
          "Rattus losea"= "#E38D26",
          "Rattus norvegicus"= "#F1CC74",
          "Suncus murinus"= "#61496D",
          "Mus caroli"= "#CAC1D4")
###绘制广东地区地图
lwd_pt <- .pt*72.27/96

###绘制湛江地区的地图
Zhanjiang <- subset(china_prov,city=='湛江市')
longlat <- aggregate(list(Host_Taxo$Longtitude,Host_Taxo$Latitude),list(
  Host_Taxo$Collected_Region,Host_Taxo$Collected_Time),mean)
colnames(longlat) <- c("Collected_Region","Collected_Time",'long',"lat")
longlat3 <- subset(selected_wide,Collected_Region == 'Anpu' | Collected_Region == 'Zhanjiang')
Time <- merge(Host_long,Region)
Time <- merge(Time,unique(Host_Taxo[,c("Region","Collected_Time")]))
Time <- subset(Time,Collected_Region=='Anpu' | Collected_Region=='Zhanjiang')
Time_long <- aggregate(Time$Num,list(Time$Collected_Region,Time$Collected_Time,Time$Taxo),sum)
colnames(Time_long) <- c("Collected_Region","Collected_Time","Taxo","Num")
Time <- spread(Time_long,key = 'Taxo',value = 'Num',fill=0)
Time$total <- rowSums(Time[, c(3:11)])
Time$long <- ifelse(Time$Collected_Time == '2021-Winter',111+2,
                    ifelse(Time$Collected_Time == '2022-Autumn',111.5+2+0.1,
                           ifelse(Time$Collected_Time == '2022-Spring',112+2+0.2,
                                  ifelse(Time$Collected_Time == '2022-Summer',112.5+2+0.3,
                                         ifelse(Time$Collected_Time == '2022-Winter',113+2+0.4,'0'
                    )))))
Time$lat <- ifelse(Time$Collected_Region == 'Anpu',22-0.8,
                    ifelse(Time$Collected_Region == 'Zhanjiang',21.5-1,0
                    ))
Time$long <- as.numeric (Time$long)
Time$lat <- as.numeric (Time$lat)
Time <- unique(subset(Time,Collected_Region == 'Anpu' | Collected_Region == 'Zhanjiang'))

Maps <- ggplot() +
  geom_sf(data=Guangdong,
          color = "#011910", linewidth=0.15)+
  annotation_scale(location = "tl",height = unit(0.25, "cm")) + #添加比例尺
  #geom_sf_text(data=Guangdong,aes(label = name_en, geometry = geometry), color = 'black', size=3)+
  geom_text(data=selected_wide,aes(label = Collected_Region, x = long-0.4,y=lat), color = 'black',
            size=6/lwd_pt)+

  geom_scatterpie_legend(log10(selected_wide$total)/10+0.1,x = 110,y=25,n = 3,
                         labeller=function(x) floor(10^((x-0.1)*10))) +
  theme_void()+
  theme(legend.position = c(0.9, 0.3),
        legend.direction = 'vertical',
        axis.title = element_blank(),
        plot.title = element_text(size = 7/lwd_pt),
        plot.background = element_rect(fill = "white", color = "white"),
  )+
  # geom_text(data=longlat3,aes(label = Collected_Region, x = long-0.2,y=lat), color = 'black', size=4)+
  geom_scatterpie(data=Time,
                  aes(x = long,y = lat,r=log10(total)/10+0.1),
                  cols=names,
                  alpha=.9)+
  geom_text(data=Time,aes(label = Collected_Time, x = long,y=20,angle=90), 
            color = 'black', size=6/lwd_pt)+
  annotate("segment",x=110.4107, y=21.25165,xend= 110.5+1, yend=21.5-1,size=1)+
  annotate("segment",x=110.5+1, y=21.5-1,xend= 110.7+2, yend=21.5-1,size=1)+
  annotate("segment",x=110.2779, y=21.57204,xend= 110.5+1,yend= 22-1,size=1)+
  annotate("segment",x=110.5+1, y=22-1,xend= 110.7+2,yend= 22-1,size=1)+
  annotate("segment",x=110.7+2, y=22-1+0.1,xend= 110.7+2,yend= 22-1-.1,size=1)+
  annotate("segment",x=110.7+2, y=21.5-1+.1,xend= 110.7+2,yend= 21.5-1-0.1,size=1)+
  geom_point(data=longlat3,aes(x=long,y=lat),size=3,color='#ECC97F')+
  geom_scatterpie(data=selected_wide,
                  aes(x = long,y = lat,r=log10(total)/10+0.1),
                  cols=names,
                  alpha=.9) +
  scale_fill_manual(values =color)
Maps
ggsave(filename = paste0(out_path,'/fig 1a.pdf',sep=''),  width = 121, height = 121,units = 'mm')
##绘制区域之间的柱状图
Region <- ggplot(selected_long, aes(x=Collected_Region,y = Num,fill=Taxo)) +
  geom_bar(stat = "identity",position = "fill",
           width = 0.8,colour = "black", size = 0.25,
           alpha = 1)+ coord_flip()+
  labs(x="",y="")+theme_bw()+scale_fill_manual(values =color)+
  ggtitle('Species composition (sites)')+ theme(legend.position = 'none')
Region
ggsave(filename = paste0(out_path,'/fig 1c.pdf',sep=''),  width = 90, height = 60,units = 'mm')
###绘制不同时间之间的柱状图
Time <- unique(Host_Taxo[,c("Region","Collected_Time","Collected_Region")])
selected <- merge(Host_long,Time)
selected <- subset(selected,Collected_Region=='Anpu'|Collected_Region=='Zhanjiang')
selected_long <- aggregate(selected$Num,
                           list(selected$Collected_Time,selected$Taxo),sum)
colnames(selected_long) <- c("Collected_Time","Taxo","Num")
Time <- ggplot(selected_long, aes(x=Collected_Time,y = Num,fill=Taxo)) +
  geom_bar(stat = "identity",position = "fill",
           width = 0.8,colour = "black", size = 0.25,
           alpha = 1)+ coord_flip()+
  labs(x="",y="")+theme_bw()+scale_fill_manual(values =color)+
  ggtitle('Species composition (season)')+ theme(legend.position = 'none')
Time
ggsave(filename = paste0(out_path,'/fig 1d.pdf',sep=''), Time, width = 90, height = 60,units = 'mm')
###绘制宿主进化树
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
## BiocManager::install("BiocUpgrade") ## you may need this
BiocManager::install("treeio")
install.packages("ape")
remotes::install_github("liamrevell/phytools")
library(ggtree)
library(ape)
library(phytools)
library(treeio)
name='COI_Ref_Query_Sequence_linsi_trimed.phy_phyml_tree.nwk'
###读取树
tree <- read.tree(paste(tree_path,name,sep='/'))
tree_rooted <- midpoint_root(tree)
ann <- read.csv(paste(data_path,'Ref_Host_Taxo.txt',sep='/'),sep='\t')
ann <- subset(ann,Ref_CDS %in% tree_rooted$tip.label)
Tree <- ggtree(tree_rooted) %<+% ann +
  geom_point2(aes(subset=label %in% unique(ann$Ref_CDS),color=Host_Taxo), size=3)+
  #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="red")+
  #geom_tiplab(align=TRUE, linetype='dashed', linesize=.3)+
  scale_x_continuous(breaks = c(seq(0,3.2,by=1),2),limits = c(-0.2,0.5))+
  #theme(legend.position=c(.2, .8))+
  geom_cladelab(node=133, label="Rattus", align=TRUE,offset.text=0)+
  geom_cladelab(node=173, label="Bandicota", align=TRUE,offset.text=0)+
  geom_cladelab(node=186, label="Crocidura", align=TRUE,offset.text=0)+
  # geom_cladelab(node=142, label="Rattus norvegicus", align=TRUE,offset.text=0.1)+
  geom_cladelab(node=180, label="Suncus", align=TRUE,offset.text=0)+
  geom_cladelab(node=130, label="Berylmys", align=TRUE,offset.text=0)+
  geom_cladelab(node=119, label="Mus", align=TRUE,offset.text=0)+
  geom_cladelab(node=107, label="Niviventer", align=TRUE,offset.text=0)+
  geom_cladelab(node=192, label="Sorex", align=TRUE,offset.text=0)+
  #geom_hilight(node=146,fill = "#789262",alpha = 0.4)+
  geom_treescale(x=0, y=0)+
  scale_color_manual(values=color)+theme(legend.position = 'none')
Tree
ggsave(filename = paste0(out_path,'/fig 1b.pdf',sep=''), plot = Tree, width = 64, height = 121,units = 'mm')

###对图形进行组合
install.packages("patchwork")
library(patchwork)

layout <- "
aaaabb
aaaabb
dddeee
"
Combine <- Maps+Tree+Region+Time+plot_layout(design = layout)+
  plot_annotation(tag_levels = "a")+ plot_layout(guides = 'collect') & theme(legend.position = 'none')

ggsave(filename = paste0(out_path,'/fig 1.pdf',sep=''), plot = Combine, width = 185, height = 161,units = 'mm')



