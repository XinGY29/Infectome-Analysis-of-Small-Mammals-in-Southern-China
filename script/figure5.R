###figure 5
#library
library(ggplot2)

# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")
##figure 5a
plot_data <- read.csv(paste0(data_path,'/fig5a.csv',sep=''),check.names = F)
plot_data <- plot_data[order(plot_data$dev_exp,decreasing = TRUE),]
legend_labels <- plot_data$var
names(legend_labels) <- plot_data$var
plot_data$var <- factor(plot_data$var,levels = c("Region of sample collected",
                                                 "Host species", "Time of sample collected", 
                                                 "Environment factors", "Unexplained"))

color = c("Region of sample collected" = "#DB432C",
          "Host species" = "#438870",
          "Time of sample collected" = "#838AAF",
          'Environment factors'='#C4B797',
          "Unexplained" = "gray")
dev_prop_plot = ggplot(plot_data) +
  geom_bar(aes(x="", y=-dev_exp, fill=var), color="black", stat="identity", width=1) +
  coord_polar("y", start= pi*0/180) +
  scale_fill_manual(values =color ,name="Variable", 
                    labels=legend_labels) + 
  theme(legend.position = "right", 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),  
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(filename = paste0(out_path,'/fig 5a.pdf',sep=''), plot = dev_prop_plot, width = 100, height = 45,units = 'mm')

##figure 5e
plot_data <- read.csv(paste0(data_path,'/fig5e.csv',sep=''),check.names = F)
plot_data <- plot_data[order(plot_data$dev_exp,decreasing = TRUE),]
legend_labels <- plot_data$var
names(legend_labels) <- plot_data$var
plot_data$var <- factor(plot_data$var,levels = c("Host species", "Region of sample collected",
                                                 "Environment factors","Time of sample collected",  "Unexplained"))

color = c("Region of sample collected" = "#DB432C",
          "Host species" = "#438870",
          "Time of sample collected" = "#838AAF",
          'Environment factors'='#C4B797',
          "Unexplained" = "gray")
dev_prop_plot = ggplot(plot_data) +
  geom_bar(aes(x="", y=-dev_exp, fill=var), color="black", stat="identity", width=1) +
  coord_polar("y", start= pi*0/180) +
  scale_fill_manual(values =color ,name="Variable", 
                    labels=legend_labels) + 
  theme(legend.position = "right", 
        legend.title = element_text(size=10),
        legend.text = element_text(size=8),  
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
ggsave(filename = paste0(out_path,'/fig 5e.pdf',sep=''), plot = dev_prop_plot, width = 100, height = 45,units = 'mm')

###figure 5b
pred_df_1 <- read.csv(paste0(data_path,'/fig5b.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$land <- factor(pred_df_1$land,levels = pred_df_1$land)
richness.time.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of pathogens \nper individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5b.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')

###figure 5c
pred_df_1 <- read.csv(paste0(data_path,'/fig5c.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$species <- factor(pred_df_1$species,levels = pred_df_1$species)
richness.time.plot = ggplot(aes(x=species, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of pathogens \nper individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5c.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')

###figure 5d
pred_df_1 <- read.csv(paste0(data_path,'/fig5d.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$land <- factor(pred_df_1$land,levels = pred_df_1$land)
richness.time.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of pathogens \nper individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5d.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')

###figure 5f
pred_df_1 <- read.csv(paste0(data_path,'/fig5f.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$land <- factor(pred_df_1$land,levels = pred_df_1$land)
richness.time.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of zoonotic pathogens \ninfect human per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5f.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')

###figure 5g
pred_df_1 <- read.csv(paste0(data_path,'/fig5g.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$land <- factor(pred_df_1$land,levels = pred_df_1$land)
richness.time.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of zoonotic pathogens \ninfect human per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5g.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')

###figure 5h
pred_df_1 <- read.csv(paste0(data_path,'/fig5h.csv',sep=''),check.names = F)
pred_df_1 <- pred_df_1[order(pred_df_1$mean,decreasing = TRUE),]
pred_df_1$time <- factor(pred_df_1$time,levels = pred_df_1$time)
richness.time.plot = ggplot(aes(x=time, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of zoonotic pathogens \ninfect human per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(filename = paste0(out_path,'/fig 5h.pdf',sep=''), plot = richness.time.plot, width = 45, height = 45,units = 'mm')


##figure 5i
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0,0.09,0.18), c("#339DB5","#F5C96B", "#C9352B"))
positivate.rate.table <- read.csv(paste0(data_path,'/fig5e1.positivate.rate.table.csv',sep=''),check.names = F,row.names = 1)
positivate.P.label <- read.csv(paste0(data_path,'/fig5e1.positivate.P.label.csv',sep=''),check.names = F,row.names = 1)
pathogen.names <- read.csv(paste0(data_path,'/fig5e1.pathogen.names.csv',sep=''),check.names = F,row.names = 1)
data <- read.csv(paste0(data_path,'/fig5e1.data.csv',sep=''),check.names = F,row.names = 1)
positivate.P.label[is.na(positivate.P.label)] <- '' 

pathogen.names <- pathogen.names$x
maxtri<-positivate.rate.table[,-1][pathogen.names,c("2021-Winter", "2022-Spring", "2022-Summer", "2022-Autumn", "2022-Winter", "residential aera", 
                                                    "agricultural aera", "Suncus murinus", "Rattus tanezumi", "Rattus norvegicus", "Bandicota indica", "Rattus losea", 
                                                    "Rattus andamanensis", "Berylmys sp.", "Mus caroli", "Niviventer lotipes", "Maoming", "Zhanjiang", 
                                                    "Anpu", "Heyuan", "Jieyang", "Shaoguang", "Foshan", "Shenzhen", 
                                                    "Yunfu")]*100
col_fun = colorRamp2(c(0,25,50), c("#339DB5","#F5C96B", "#C9352B"))
display_numbers = positivate.P.label[,-1][pathogen.names,c("2021-Winter", "2022-Spring", "2022-Summer", "2022-Autumn", "2022-Winter", "residential aera", 
                                                           "agricultural aera", "Suncus murinus", "Rattus tanezumi", "Rattus norvegicus", "Bandicota indica", "Rattus losea", 
                                                           "Rattus andamanensis", "Berylmys sp.", "Mus caroli", "Niviventer lotipes", "Maoming", "Zhanjiang", 
                                                           "Anpu", "Heyuan", "Jieyang", "Shaoguang", "Foshan", "Shenzhen", 
                                                           "Yunfu")]
rownames(data) <- data$Pathogen
row_Ann <- rowAnnotation(
  'Positive rate' = anno_barplot(data[pathogen.names,"rate"],gp = gpar(col = "#384259",fill="#384259"),
                                 width = unit(10,'mm')
  ),
  annotation_name_gp= gpar(fontsize = 5),simple_anno_size = unit(0.2, "cm")
)

pdf(paste0(out_path,'/fig 5i.pdf',sep=''), width = 8, height = 4)

hhh <- Heatmap(maxtri,name = "Positivate (%)",col = col_fun, cluster_rows=F,
               show_heatmap_legend = FALSE,
               width = ncol(maxtri)*unit(3, "mm"), height = nrow(maxtri)*unit(3, "mm"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(sprintf("%s", display_numbers[i, j]), x, y, gp = gpar(fontsize = 5))},
               rect_gp = gpar(col = 'gray', lwd = 0.1),row_names_side = c( "left"),
               row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5),
               column_split  = rep(c('Season','liveregion','host_Type','RTegion'),
                                   c(5,2,9,9)
               ),
               cluster_columns = F,right_annotation = row_Ann)

lgd_list = list(
  Legend(col_fun = col_fun, title = "Positivate rate(%)",labels_gp = gpar(fontsize = 5),
         title_gp = gpar(fontsize = 6),grid_width = unit(2, "mm"),title_position = "lefttop-rot")
)
draw(hhh,merge_legend = FALSE,annotation_legend_list = lgd_list)
dev.off()
