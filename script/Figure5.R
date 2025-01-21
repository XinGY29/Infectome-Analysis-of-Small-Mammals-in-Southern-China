####figure5
# R file path
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
##Path
data_path <- paste(script_dir, '..','raw_data', sep = "/")
out_path <- paste(script_dir, '..','output', sep = "/")
tree_path <- paste(script_dir, '..','tree', sep = "/")
##package
library(ggthemes)
library(tidyr)
library(pairwiseAdonis)
library(ggpubr)
library(patchwork)
library(cowplot)
library(ggplot2) 
library(grid)
library(gridExtra) 
library(dplyr)
library(RColorBrewer)
library(ComplexUpset)
library(ggsci)
library(dylyr)
library(vegan)
library(readr)
library(stats)
library(MuMIn)
library(purrr)
library(tidyverse)
###Model 1

setwd(data_path)
ecology_data <- read.csv("Env_information.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Individual_data <- read.csv("Individual_Information.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Individual_data$Hos_Taxo <- ifelse(Individual_data$Hos_Taxo=='Berylmys sp.','Berylmys bowersi',Individual_data$Hos_Taxo)
selected_data <- merge(ecology_data,Individual_data)
selected_data$Collected_Region <- ifelse(selected_data$Collected_Region == 'AP','Anpu',
                                         ifelse(selected_data$Collected_Region == 'ZJ','Zhanjiang',
                                                ifelse(selected_data$Collected_Region == 'DB','Maoming',
                                                       ifelse(selected_data$Collected_Region == 'HY','Heyuan',
                                                              ifelse(selected_data$Collected_Region == 'JY','Jieyang',
                                                                     ifelse(selected_data$Collected_Region == 'SG','Shaoguan',
                                                                            ifelse(selected_data$Collected_Region == 'SS','Foshan',
                                                                                   ifelse(selected_data$Collected_Region == 'SZ','Shenzhen',
                                                                                          ifelse(selected_data$Collected_Region == 'YF','Yunfu',''
                                                                                          )))))))))

# Function to selecte PCs with eigenvalues greater than 1
print_pcs_with_eigenvalue_above_1 <- function(pca_result) {
  eigenvalues <- pca_result$sdev^2
  var_explained <- summary(pca_result)$importance[2,] * 100  
  pcs_above_1 <- which(eigenvalues > 0.6)  
  if (length(pcs_above_1) > 0) {
    top_pcs_variances <- var_explained[pcs_above_1]
    print(top_pcs_variances)
    total_variance_explained <- sum(top_pcs_variances)
    cat("Total variance explained by PCs with eigenvalues > 1:", total_variance_explained, "%\n")
    return(pcs_above_1) } else {
      cat("No principal components have eigenvalues > 1.\n")
      return(NULL)}}

###clim_data_selection
data_clim <- selected_data[, c("Individual","aet","def","q","soil","PDSI","pet","ppt","srad","tmax","tmin","vap","vpd","ws")]
write.csv(data_clim,"data_clim.csv")
group_ids <- data_clim$Individual

pca_clim = data_clim %>%
  dplyr::select(-Individual) %>%
  scale() %>%
  prcomp()
clim_scores_with_ids <- cbind(Individual = group_ids, as.data.frame(pca_clim$x))
important_pcs <- print_pcs_with_eigenvalue_above_1(pca_clim)

if (!is.null(important_pcs)) {
  selected_scores <- pca_clim$x[, important_pcs, drop = FALSE]
  clim_selected_data <- cbind(Individual = group_ids, selected_scores)}

clim_selected_data <- as.data.frame(clim_selected_data)
clim_selected_data <- clim_selected_data %>%
  dplyr::rename(CPC1 = PC1, CPC2 = PC2, CPC3 = PC3)
clim_selected_data <- as.matrix(clim_selected_data)
###clim_PCA_plot
clim_pca_scores_with_location <- merge(clim_scores_with_ids, selected_data[, c("Individual", "Region")])
clim_PCA_plot <- ggplot(clim_pca_scores_with_location, aes(x = PC1, y = PC2, color = Region)) +
  geom_point(alpha = 0.8, size = 3) + 
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill=NA,linewidth = 1),  
        panel.grid.major = element_line(colour = "gray", linewidth = 0.6),  
        panel.grid.minor = element_line(colour = "gray", linewidth = 0.6)) + 
  labs(x = "Clim_PC1",
       y = "Clim_PC2",
       title = "PCA Plot for clim",
       color = "Region") +
  scale_color_ucscgb()
clim_PCA_plot
ggsave("clim_PCA_plot.pdf", plot = clim_PCA_plot, width = 5.5, height = 4.5)

###clim_loading
clim_loadings <- pca_clim$rotation
write.csv(clim_loadings, "clim_loadings.csv")

###land_data_selection
data_land <- selected_data[, c("Individual","land_used")]
write.csv(data_land,"data_land.csv")
group_ids <- data_land$Individual

###data_mergence
processed_data <- merge(selected_data, clim_selected_data)

library(MASS) 

###variable_categorization

###host_species
host_vars <- c("Hos_Taxo")

###climate_vars_and_vegetation_index_and_mammal_richness
environment_vars <- c("CPC1", "CPC2", "CPC3")

###land_use_and_population
anthrome_vars <- c("land_used", "NDVI", "Mammal_richness")

#supplementary_vars <- c("Longtitude", "Latitude","Collected_Region","Collected_Time")
#date_vars <- ("year", "season")
supplementary_vars <- c("Collected_Region","Collected_Time")

candidate_vars <- c(environment_vars, anthrome_vars, host_vars, supplementary_vars)
all_nona <- processed_data %>%
  dplyr::select(Collected_Time,Individual, any_of(candidate_vars)) %>%
  na.omit()
##Richness
library(vegan)
RMP_Table <- read.csv('F:\\Mouse_Result_24\\RPM_Table/Pathogens_Mammal-related_Wide_RPM_Table.csv',check.names = F)
RMP_Table <- subset(RMP_Table,select = -c(`Desulfovibrio piger`,`Oryza sativa`))
colnames(RMP_Table)[1] <- 'Individual'
RPM.table.long <- gather(RMP_Table,key = 'species',value = 'RPM',-'Individual')
Ann_Table <- read.csv('F:\\Mouse_Result_24\\RPM_Table/Pathongen_Taxon_Table_2.txt',sep='\t')
Ann_Table <- unique(Ann_Table)
Mammal.releted <- subset(Ann_Table,Host_Type == 'Mammal-related')
Mammal.releted <- Mammal.releted %>%
  filter(!(family %in% c('Picobirnaviridae','Anelloviridae')))
RPM.table.long <- merge(RPM.table.long,Mammal.releted)

time.series.wide <- spread(RPM.table.long[,c('Individual','species','RPM')],key = 'species',value = RPM,fill = 0)
rownames(time.series.wide) <- time.series.wide$Individual
time.series.wide <- time.series.wide[,-1]
Richness <- specnumber(t(time.series.wide), MARGIN = 2)#spe.rich =sobs
alpha <- as.data.frame(Richness)
alpha$Individual <- rownames(alpha)
#alpha <- read.csv('F:/Mouse_Result_24/RPM_Table/Individual_Mammal-related_Pathogen_Alpha.csv')

all_nona <- merge(all_nona,alpha)
all_nona$CPC1 <- as.numeric(all_nona$CPC1)
all_nona$CPC2 <- as.numeric(all_nona$CPC2)
all_nona$CPC3 <- as.numeric(all_nona$CPC3)
all_nona$land_used <- as.factor(all_nona$land_used)
all_nona$NDVI <- as.numeric(all_nona$NDVI)
all_nona$Hos_Taxo <- as.factor(all_nona$Hos_Taxo)
all_nona$Collected_Time <- as.factor(all_nona$Collected_Time)
all_nona$Collected_Region <- as.factor(all_nona$Collected_Region)
all_nona$Mammal_richness <- as.numeric(all_nona$Mammal_richness)
mean_virus_richness <- mean(all_nona$Richness)
variance_virus_richness <- var(all_nona$Richness)

write.csv(all_nona,"Richness_all_nona.csv")
enum_formula <- function(n) {
  comb <- combinat::combn(candidate_vars, n)
  if (any(is.null(dim(comb)))) {
    ret <- paste0("Richness~", paste0(comb, collapse = "+"))
  } else {
    ret <- apply(comb, 2, function(x) {
      paste0(
        "Richness~",
        paste0(x, collapse = "+")
      )
    })
  }
  return(ret)
}
###negative_binomial_model
all_fitted <- purrr::map(1:length(candidate_vars), enum_formula) %>%
  reduce(c) %>%
  purrr::map(~glm.nb(as.formula(.x), data = all_nona))

model_sel <- MuMIn::model.sel(all_fitted, rank=AIC)
model_index <- rownames(model_sel) %>% as.integer()

topN <- 10
var_explained <- all_fitted[model_index[1:topN]] %>%
  lapply(function(x) {
    drop_one <- drop1(x)
    drop_tbl <- drop_one %>%
      as_tibble() %>%
      mutate(term = rownames(drop_one), .before = 1) %>%
      mutate(dev_exp = (Deviance-x$deviance)/x$null.deviance) %>%
      dplyr::select(-Df, -Deviance, -AIC)
    drop_tbl[1, 2] <- (x$null.deviance - x$deviance) / x$null.deviance
    drop_tbl
  }) %>%
  reduce(full_join, by="term")

colnames(var_explained) = c("term", paste0("dev_exp_", 1:topN))
write.csv(var_explained, "var_explained_richness.csv")

nb_model <- MASS::glm.nb(Richness ~ CPC1 + CPC2 + CPC3 + Hos_Taxo + Collected_Region + NDVI + Mammal_richness +
                           Collected_Time,  data = all_nona)
#nb_model <- MASS::glm.nb(Richness ~ CPC3 + Collected_Region + Hos_Taxo + Collected_Time,  data = all_nona)
summary(nb_model)
x2.nb_model <- sum(residuals(nb_model, type = "pearson")^2)
x2.nb_model
qchisq(0.975, df.residual(nb_model))

poisson_model <- glm(Richness ~ CPC1 + CPC2 + CPC3 + Hos_Taxo + Collected_Region + NDVI + Mammal_richness +
                       Collected_Time,  data = all_nona, family = poisson)
summary(poisson_model)
x2.poisson_model <- sum(residuals(poisson_model, type = "pearson")^2)
x2.poisson_model
qchisq(0.975, df.residual(poisson_model))

###pro_calculate
sel_environment_vars <- c("CPC1","CPC2","CPC3","NDVI", "Mammal_richness")
f = paste0("Richness~", paste0(sel_environment_vars, collapse="+")) %>%
  as.formula() %>%
  glm.nb(data=all_nona)
dev_explained_environment_vars = (f$null.deviance-f$deviance) / f$null.deviance

f= paste0("Richness~", paste0(host_vars, collapse="+")) %>%
  as.formula() %>%
  glm.nb(data=all_nona)
dev_explained_host_vars = (f$null.deviance-f$deviance) / f$null.deviance
Time <- c("Collected_Time")
f= paste0("Richness~", paste0(Time, collapse="+")) %>%
  as.formula() %>%
  glm.nb(data=all_nona)
dev_explained_Time_vars = (f$null.deviance-f$deviance) / f$null.deviance
Region <- c("Collected_Region")
f= paste0("Richness~", paste0(Region, collapse="+")) %>%
  as.formula() %>%
  glm.nb(data=all_nona)
dev_explained_Region_vars = (f$null.deviance-f$deviance) / f$null.deviance

dev_prop = tibble(
  var = c("Environment factors","Host species", "Time of sample collected","Region of sample collected"),
  dev_exp = c(dev_explained_environment_vars, dev_explained_host_vars, 
              dev_explained_Time_vars,dev_explained_Region_vars)
)
dev_prop = rbind(dev_prop, c("var"="Unexplained",
                             "dev_exp"=1-sum(dev_prop$dev_exp)))
dev_prop$dev_exp = as.double(dev_prop$dev_exp)
write.csv(dev_prop, "dev_prop_richness.csv")

plot_data = dev_prop %>%
  arrange(dev_exp) %>%
  mutate(prop = dev_exp / sum(dev_exp) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop)
plot_data$var = factor(plot_data$var,levels = c('Region of sample collected','Host species','Time of sample collected',
                                                'Environment factors','Unexplained'))

category_percentages <- plot_data %>%
  group_by(var) %>%
  summarise(total_prop = sum(prop)) %>%
  arrange(desc(total_prop)) %>%
  mutate(legend_label = paste(var, sprintf("(%.1f%%)", total_prop)))
legend_labels <- category_percentages$legend_label
names(legend_labels) <- category_percentages$var
color = c("Region of sample collected" = "#DB432C",
          "Host species" = "#438870",
          "Time of sample collected" = "#838AAF",
          'Environment factors'='#C4B797',
          "Unexplained" = "gray")
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
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5 Richness of pie.pdf',sep=''),
  Richness.pie,
  width = 90,             # 宽
  height = 60,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
###CPC3
hhh <- all_nona[,c("Richness",'CPC3')]
colnames(hhh) <- c("Richness","Var")
CPC3.plot <- ggplot(data = hhh,aes(x=Var,y=Richness))+
  geom_point( size=2,alpha=0.5,color='gray')+
  geom_smooth(method = MASS::glm.nb,
              color='black',se = TRUE)+
  theme_bw()+xlab("CPC3")+
  ylab("Estimated richness of individuals")+
  theme(legend.position = "none", 
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),  
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size=5),axis.title = element_text(size=5))
CPC3.plot

###Collected_Time
Time_fitted = predict(nb_model, terms = "Collected_Time",type = "terms", se=T)

pred_df = tibble(Individual = all_nona$Individual, 
                 Collected_Time=all_nona$Collected_Time, 
                 sp_effect = Time_fitted$fit[,1],
                 se = Time_fitted$se.fit[,1]) %>%
  group_by(Collected_Time)
pred_df_1 <- aggregate(list(pred_df$sp_effect,pred_df$se), by = list(pred_df$Collected_Time), 
                       FUN = mean)
colnames(pred_df_1) <- c("land","sp_effect","se")
pred_df_1$land = factor(pred_df_1$land, levels=arrange(pred_df_1, desc(pred_df_1$sp_effect))$land)

richness.time.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of pathogen per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5 estimated richness of time.pdf',sep=''), 
       plot=richness.time.plot,width = 60,             # 宽
       height = 60,            # 高
       units = "mm",          # 单位
       dpi = 300,              # 分辨率DPI
       limitsize = FALSE)


###Collected_Region
Time_fitted = predict(nb_model, terms = "Collected_Region",type = "terms", se=T)

pred_df = tibble(Individual = all_nona$Individual, 
                 Collected_Time=all_nona$Collected_Region, 
                 sp_effect = Time_fitted$fit[,1],
                 se = Time_fitted$se.fit[,1]) %>%
  group_by(Collected_Time)
pred_df_1 <- aggregate(list(pred_df$sp_effect,pred_df$se), by = list(pred_df$Collected_Time), 
                       FUN = mean)
colnames(pred_df_1) <- c("land","sp_effect","se")
pred_df_1$land = factor(pred_df_1$land, levels=arrange(pred_df_1, desc(pred_df_1$sp_effect))$land)

richness.region.plot = ggplot(aes(x=land, y=sp_effect), data=pred_df_1) +
  geom_point(size=1.5) +
  geom_errorbar(aes(ymin=sp_effect-1.96*se, ymax=sp_effect+1.96*se, y=sp_effect), width=0, linewidth=1, alpha=0.5, color="blue") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme( axis.text.x = element_text(angle =30, hjust=1)) +
  ylab("Estimated effect size\n(number of pathogen per individual)")+
  xlab("")+
  theme( axis.text.x = element_text(angle =30, hjust=1,size=5),
         axis.text.y = element_text(size=5),
         axis.title.y=element_text(size=5))
ggsave(paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5 estimated richness of region.pdf',sep=''), 
       plot=richness.region.plot,width = 60,             # 宽
       height = 60,            # 高
       units = "mm",          # 单位
       dpi = 300,              # 分辨率DPI
       limitsize = FALSE)




##read RPM Table
RMP_Table <- read.csv('Pathogens_Mammal-related_Wide_RPM_Table.csv',check.names = F)
RMP_Table <- subset(RMP_Table,select = -c(`Desulfovibrio piger`,`Oryza sativa`))
colnames(RMP_Table)[1] <- 'Individual'
Ann_Table <- read.csv('Pathongen_Taxon_Table_2.txt',sep='\t')
Ann_Table <- unique(Ann_Table)
Ann_Table$Type <- ifelse(Ann_Table$Type == 'RNA_Virus','RNA Virus',
                         ifelse(Ann_Table$Type == 'Bac','Bacteria',
                                ifelse(Ann_Table$Type == 'Euk','Eukaryote',
                                       ifelse(Ann_Table$Type == 'DNA_Virus','DNA Virus',Ann_Table$Type))))
Mammal.releted <- subset(Ann_Table,Host_Type == 'Mammal-related')
Individual.infor <- read.csv('Full_Tissue_Individual.csv')
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
##colors
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color74 = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF",
            "#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
color20 <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
color37 = c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
Region <- unique(Individual.infor$Collected_Region)
Region_color <- color20[1:length(Region)]


###Site Sample Sum
RPM.table.long <- gather(RMP_Table,key = 'species',value = 'RPM',-'Individual')
RPM.table.long <- merge(RPM.table.long,Mammal.releted,all.x = TRUE)
RPM.table.long <- merge(RPM.table.long,Individual.infor,all.x = TRUE)
RPM.table.long$Site_Type <- paste(RPM.table.long$Collected_Region,RPM.table.long$Collected_Time,sep='_')

##Host constitute of echo Sample Site
###Mouse Sample distribution
Mouse.sample <- unique(RPM.table.long[,c('Individual','Hos_Taxo','Site_Type','Collected_Time','Collected_Region')])
Mouse.sample.sum <- aggregate(Mouse.sample$Individual,
                              list(Mouse.sample$Hos_Taxo,Mouse.sample$Site_Type,Mouse.sample$Collected_Region,
                                   Mouse.sample$Collected_Time),length)
colnames(Mouse.sample.sum) <- c('Hos_Taxo','Site_Type','Collected_Region','Collected_Time','Num')
###host distribution
host.total.plot <- ggplot(data = Mouse.sample.sum,aes(x=Site_Type,y=Num,fill=Hos_Taxo))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25)+
  scale_fill_manual(values = setNames(color20[1:length(unique(Mouse.sample.sum$Hos_Taxo))],unique(Mouse.sample.sum$Hos_Taxo)))+
  theme_few()+
  labs(x="")+
  theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
        axis.text.x=element_text(angle=90,hjust = 1))
time.host <- subset(Mouse.sample.sum,Collected_Region == 'Anpu' | Collected_Region=='Zhanjiang')
time.host.sum <- aggregate(time.host$Num,
                           list(time.host$Hos_Taxo,
                                time.host$Collected_Time),sum)
colnames(time.host.sum) <- c('Hos_Taxo','Collected_Time','Num')
host.time.plot <- ggplot(data = time.host.sum,aes(x=Collected_Time,y=Num,fill=Hos_Taxo))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25)+
  scale_fill_manual(values = setNames(color20[1:length(unique(Mouse.sample.sum$Hos_Taxo))],unique(Mouse.sample.sum$Hos_Taxo)))+
  theme_few()+
  labs(x="Time")+
  theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
        axis.text.x=element_text(angle=90,hjust = 1))
region.host <- subset(Mouse.sample.sum,Collected_Time == '2021-Winter' | 
                        Collected_Time=='2022-Winter' | 
                        Collected_Time == '2022-Autumn')

host.region.plot <- ggplot(data = region.host,aes(x=Site_Type,y=Num,fill=Hos_Taxo))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25)+
  scale_fill_manual(values = setNames(color20[1:length(unique(Mouse.sample.sum$Hos_Taxo))],unique(Mouse.sample.sum$Hos_Taxo)))+
  theme_few()+
  labs(x="Region")+
  theme(text=element_text(size= 6),legend.key.size = unit(3, 'mm'),
        axis.text.x=element_text(angle=90,hjust = 1))
host.region.plot
#host.plot.sum <- host.plot+host.region.plot+host.time.plot+plot_layout(guides = "collect")

##figure s Pathogens constitute of echo Sample Site
Region.data <- aggregate(RPM.table.long$RPM,
                         list(RPM.table.long$family,RPM.table.long$Collected_Region,
                              RPM.table.long$Type),mean)
colnames(Region.data) <- c('family','Region','Type','RPM')
dae.filll.color <- setNames(color37[1:length(unique(Region.data$family))],unique(Region.data$family))
Region.data$Type <- factor(Region.data$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
Region.plot <- ggplot(data = Region.data,aes(x=Region,y=log10(RPM+1),fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  facet_grid(Type ~ .)+ labs(x="")+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))
Region.plot

Selected <- subset(RPM.table.long,Collected_Region %in% c("Anpu","Zhanjiang"))
Time <- aggregate(Selected$RPM,list(Selected$family,Selected$Collected_Time,Selected$Type),sum)
colnames(Time) <- c("Family","Time","Type","RPM")
Time$Type <- factor(Time$Type,levels = c("RNA Virus","DNA Virus","Bacteria","Eukaryote"))
Time_plot <- ggplot(data = Time,aes(x=Time,y=log10(RPM+1),fill=Family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  facet_grid(Type ~ .)+ labs(x="")+
  theme(legend.position = 'none',
        text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))
Time_plot

##PCA/PCoA of Echo Sample Site
mean.rpm <- aggregate(RPM.table.long$RPM,list(RPM.table.long$species,RPM.table.long$Collected_Time,RPM.table.long$Type,
                                              RPM.table.long$Site_Type,RPM.table.long$Collected_Region),mean)
colnames(mean.rpm) <- c('species','Collected_Time','Pathogen_Type','Type','Collected_Region','RPM')
site.ann <- unique(mean.rpm[,c('Collected_Time','Type','Collected_Region')])
mean.rpm.wide <- spread(mean.rpm[,c('species','Type','RPM')],key = 'species',value = 'RPM',fill=0)
rownames(mean.rpm.wide) <- mean.rpm.wide$Type
mean.rpm.wide <- mean.rpm.wide[,-1]
pca <- prcomp(mean.rpm.wide)
pca$x
plot_data <- data.frame(pca$x)[1:2]
plot_data$Type <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCA1', 'PCA2')
plot_data <- merge(plot_data,site.ann,by = 'Type', all.x = TRUE)
Region<- unique(plot_data$Collected_Region)
Region_color <- colors[1:length(Region)]
plot_data$Collected_Time <- factor(plot_data$Collected_Time,levels =
                                     c("2021-Winter","2022-Spring","2022-Summer","2022-Autumn","2022-Winter"))
RT_PCA <- ggplot(data = plot_data, aes(x=PCA1, y=PCA2,color=Collected_Region,fill=Collected_Region,shape=Collected_Time)) +
  geom_point(alpha=.7, size=1)+
  scale_colour_manual(values = setNames(Region_color,Region))+
  scale_fill_manual(values = setNames(Region_color,Region))+
  scale_shape_manual(values=c("2021-Winter"=21,"2022-Autumn"=22,"2022-Spring"=23,
                              "2022-Summer"=24,"2022-Winter"=25))+
  theme_few()+coord_fixed(18)+
  theme(text=element_text(size= 6),legend.key.size = unit(1, 'mm'))
RT_PCA
## bar plot of echo Site
mean.rpm.dae <- aggregate(RPM.table.long$RPM,list(RPM.table.long$family,RPM.table.long$Collected_Time,RPM.table.long$Type,
                                                  RPM.table.long$Site_Type,RPM.table.long$Collected_Region),mean)
colnames(mean.rpm.dae) <- c('family','Collected_Time','Type','Site_Type','Collected_Region','RPM')
ggplot(data = mean.rpm.dae,aes(x=Site_Type,y=RPM,fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  facet_grid(Type ~ .)+ labs(x="")+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),
        axis.text.x=element_text(angle=60,hjust = 1))
dnavirus.mean <- subset(mean.rpm.dae,Type=='DNA Virus')
dnavirus.mean$family <- factor(dnavirus.mean$family,levels = unique(dnavirus.mean$family)) 
dnavirus.plot <- ggplot(data = dnavirus.mean,aes(x=Site_Type,y=RPM,fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  labs(x="")+coord_flip()+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),legend.location = 'top',
        axis.text.y=element_blank())
dnavirus.plot
rnavirus.mean <- subset(mean.rpm.dae,Type=='RNA Virus')
rnavirus.mean$family <- factor(rnavirus.mean$family,levels = unique(rnavirus.mean$family)) 
rnavirus.plot <- ggplot(data = rnavirus.mean,aes(x=Site_Type,y=RPM,fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  labs(x="")+coord_flip()+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),legend.location = 'top',
        axis.text.x=element_text(angle=90,hjust = 1))
rnavirus.plot
bac.mean <- subset(mean.rpm.dae,Type=='Bacteria')
bac.mean$family <- factor(bac.mean$family,levels = unique(bac.mean$family)) 
bac.plot <- ggplot(data = bac.mean,aes(x=Site_Type,y=RPM,fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  labs(x="")+coord_flip()+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),legend.location = 'top',
        axis.text.y=element_blank())
bac.plot
euk.mean <- subset(mean.rpm.dae,Type=='Eukaryote')
euk.mean$family <- factor(euk.mean$family,levels = unique(euk.mean$family)) 
euk.plot <- ggplot(data = euk.mean,aes(x=Site_Type,y=RPM,fill=family))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = dae.filll.color)+theme_few()+
  labs(x="")+coord_flip()+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),legend.location = 'top',
        axis.text.y=element_blank())+guides(color=guide_legend(ncol = 2,byrow = F))
euk.plot
relate.rpm.plot <- rnavirus.plot+dnavirus.plot+bac.plot+euk.plot + 
  plot_layout(ncol = 4, widths = c(1,1,1,1),guides = 'collect')
relate.rpm.plot
### Selected Time host
focus.host.list <- c('Bandicota indica','Rattus losea','Rattus norvegicus','Rattus tanezumi','Suncus murinus')
time.series <- subset(RPM.table.long,Collected_Region %in% c('Zhanjiang','Anpu') & Hos_Taxo %in% focus.host.list)
time.series$Type <- factor(time.series$Type,levels = c("RNA Virus",  "DNA Virus" ,"Eukaryote","Bacteria" )) 
time.series$Collected_Time <- factor(time.series$Collected_Time,levels = c("2021-Winter","2022-Spring","2022-Summer","2022-Autumn","2022-Winter")) 
##Selected Host
Host <- unique(time.series$Hos_Taxo)
Host.color <- colors[1:length(Host)]
host.time.series <- unique(time.series[,c('Individual','Collected_Time','Hos_Taxo')])
time.host.dis <- aggregate(time.series$Individual,
                           list(time.series$Collected_Time,time.series$Hos_Taxo),length)
colnames(time.host.dis) <- c('Time','Host','Num')
time.host.plot <- ggplot(data = time.host.dis,aes(x=Time,y=Num,fill=Host))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = setNames(Host.color,Host))+theme_few()+
  labs(x="")+
  theme(text=element_text(size= 6),legend.key.size = unit(2,'mm'),legend.location = 'top',
        axis.text.y=element_text(size=6,angle = 90))+guides(color=guide_legend(ncol = 2,byrow = F))
###绘制在不同时间共享的病毒情况 Venn
library(ComplexUpset)
library(ggsci)
pathogens.time <- aggregate(time.series$RPM,
                            list(time.series$species,time.series$Type,time.series$Collected_Time),mean)
colnames(pathogens.time) <- c('species','Type','Collected_Time','RPM')
pathogens.time <- subset(pathogens.time,RPM >0)
pathogens.time.wide <- spread(pathogens.time[,c('species','Type','Collected_Time','RPM')],key = Collected_Time,value = RPM,fill = 0)
time.nams <- c("2021-Winter","2022-Spring","2022-Summer","2022-Autumn","2022-Winter")
Type <- c("RNA Virus",  "DNA Virus" ,"Eukaryote","Bacteria" )
Type_color <- c('#DB432C','#438870','#838AAF','#C4B797')
lwd_pt <- .pt*72.27/96
time.upset.plot <- upset(pathogens.time.wide, time.nams,width_ratio=0.15,height_ratio=1,
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
                              legend.position = "none")),
                    set_sizes = (
                      upset_set_size()+geom_text(aes(label=..count..),size=4/lwd_pt, hjust=1.1, stat='count')+
                        ylab('Num of Pathogens')+
                        theme(axis.title.x = element_text(size=5),
                              axis.text = element_text(size=5),
                              text = element_text(size=5))
                    ),
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(text=element_text(size=5)),
                        'overall_sizes'=theme(axis.text.x=element_text(angle=90))
                      ))
)
time.upset.plot
##绘制不同类型的病毒随时间的变化规律
library(dylyr)
pathogens.time.wide <- pathogens.time.wide %>%  rowwise() %>% mutate(min = min(c(`2021-Winter`,`2022-Autumn`,
                                                          `2022-Spring`,`2022-Summer`,`2022-Winter`)))
pathogens.time.wide <- subset(pathogens.time.wide,min>0)
pathogens.time.wide <- subset(pathogens.time.wide,select = -c(min))
pathogens.time.long <- gather(pathogens.time.wide,key = "Time",value = RPM, -"Type", -'species')
pathogens.time.long$Time <- factor(pathogens.time.long$Time,levels = c("2021-Winter","2022-Spring","2022-Summer","2022-Autumn","2022-Winter"))
p10.sum.view <- ggplot(pathogens.time.long,aes(Time,log10(RPM+1),color=Time))+
  geom_line(aes(Time,log10(RPM+1),group=1),color="gray",size=1,linetype='solid')+labs(x="",y='log10(meanRPM+1)')+
  geom_point()+
  facet_wrap(vars(species),nrow =2) +
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),strip.text = element_text(size= 6),
        axis.text = element_text(size = 5,angle = 90))
p10.sum.view
###绘制不同科的病毒的时间变化规律
head(time.series)
times.pair <- list(c("2021-Winter","2022-Spring"),
                   c("2021-Winter","2022-Summer"),
                   c("2021-Winter","2022-Autumn"),
                   c("2021-Winter","2022-Winter"),
                   c("2022-Spring","2022-Summer"),
                   c("2022-Spring","2022-Autumn"),
                   c("2022-Spring","2022-Winter"),
                   c("2022-Summer","2022-Autumn"),
                   c("2022-Summer","2022-Winter"),
                   c("2022-Autumn","2022-Winter"))
time.names <- unique(time.series$Collected_Time)
time.names.color <- c('#55B7E6','#56BA77','#ED6E69','#9FAA3F','#BD79B6')
Pathogen.RPM.plot <- ggplot(subset(time.series,RPM >0),aes(Collected_Time,log10(RPM+1),color=Collected_Time))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
  stat_compare_means(method = 'kruskal.test',label='p.format')+
  # stat_compare_means(comparisons=times.pair,size=2,
  # label = "p.signif",method = 'wilcox.test')+
  theme_few()+facet_wrap(vars(Type),nrow =1)+
  scale_color_manual(values = setNames(time.names.color,time.names))+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),legend.position = 'bottom',
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))

Pathogen.dae.RPM.plot <- ggplot(subset(time.series,RPM >0),aes(Collected_Time,log10(RPM+1)))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
  stat_compare_means(comparisons=times.pair,size=2,
  label = "p.signif",method = 'wilcox.test')+
  theme_few()+facet_wrap(vars(family),nrow =6)+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
Pneumocystaceae.rpm.plot <- ggplot(subset(time.series,RPM >0 & family == 'Pneumocystaceae'),aes(Collected_Time,log10(RPM+1)))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='log10(RPM+1)')+ geom_jitter(size=0.5)+
  stat_compare_means(comparisons=list(c("2021-Winter","2022-Spring"),
                                      c("2022-Spring","2022-Summer"),
                                      c("2022-Spring","2022-Autumn"),
                                      c("2022-Spring","2022-Winter")
                                      ),size=2,
                     label = "p.signif",method = 'wilcox.test')+
  theme_few()+coord_fixed(1.7)+
  ggtitle('Pneumocystaceae')+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
Pneumocystaceae.rpm.plot
##挑选个别病毒描述时间变化规律
selected.pathogens <- unique(pathogens.time.long$species)
selected.pathogens <- c(
  'Norway rat pestivirus','Wenzhou virus','Orthohantavirus seoulense','Rat hepatitis E virus','Mischivirus E','Guangdong rosa-like virus','Betacoronavirus 1',
  'Guangdong rodent dependoparvovirus 1','Rat minute virus 2a','Porcine bocavirus',
  'Bartonella kosoyi','Chlamydia muridarum','Klebsiella variicola',
  'Angiostrongylus cantonensis','Brachylaima sp','Pneumocystis carinii')
time.series$species <- ifelse(time.series$species == 'Wnezhou virus','Wenzhou virus',time.series$species)
pathogen.focus.rpm <- subset(time.series,species %in% selected.pathogens)
pathogen.focus.rpm$species <- factor(pathogen.focus.rpm$species,levels = selected.pathogens)
pathogen.focus.rpm.mean <- aggregate(pathogen.focus.rpm$RPM,
                                     list(pathogen.focus.rpm$species,pathogen.focus.rpm$Hos_Taxo,pathogen.focus.rpm$Collected_Time),mean)
colnames(pathogen.focus.rpm.mean) <- c('species','Host','Collected_Time','RPM')
pathogen.focus.time.polt <- ggplot(subset(pathogen.focus.rpm.mean),aes(Collected_Time,log10(RPM+1)))+
  geom_line(aes(x = Collected_Time, y = log10(RPM+1), group=Host,color=Host),size=1,linetype='solid')+
  theme_few()+facet_wrap(vars(species),nrow =2,scales = "free",axis.labels = "all")+
  scale_color_manual(values=setNames(Host.color,Host))+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
pathogen.focus.time.polt

##Positivate Plot
##data = pathogen.focus.rpm
selected.focus1 <- aggregate(pathogen.focus.rpm$RPM, by = list(pathogen.focus.rpm$Individual,
                                                               pathogen.focus.rpm$species,
                                                               pathogen.focus.rpm$Collected_Time,
                                                               pathogen.focus.rpm$Hos_Taxo), 
                             FUN = sum)
colnames(selected.focus1) <- c('Individual','species','Time','Host','RPM')
selected.focus1$Value <- ifelse(selected.focus1$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Time + Value + species + Host, data=selected.focus1)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh_1 <- aggregate(hhhh$Freq, by = list(hhhh$Time,hhhh$species,hhhh$Host), FUN = sum)
colnames(hhhh_1) <- c('Time','species','Host','sum')
hhhhh <- merge(hhhh,hhhh_1)
time.posi.rate <- subset(hhhhh,Value == 'Positive')
time.posi.rate.selected <- subset(time.posi.rate,species %in% selected.pathogens)
time.posi.rate.selected$positive_rate <- time.posi.rate.selected$Freq/time.posi.rate.selected$sum*100

pathogen.focus.time.polt <- ggplot(subset(time.posi.rate.selected),aes(Time,positive_rate))+
  geom_line(aes(x = Time, y = positive_rate, group=Host,color=Host),size=1,linetype='solid')+
  theme_few()+facet_wrap(vars(species),nrow =2,scales = "free",axis.labels = "all")+
  scale_color_manual(values=setNames(Host.color,Host))+
  ylab('Positivate Rate (%)')+xlab('')+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
pathogen.focus.time.polt

##绘制所有病原变化规律，放置于附图 - 待定
##绘制不同个体里面的病原的Richness的信息
library(vegan)
pathogen.type <- unique(time.series[,c('species','Type')])
time.series.wide <- spread(time.series[,c('Individual','species','RPM')],key = 'species',value = RPM,fill = 0)
rownames(time.series.wide) <- time.series.wide$Individual
time.series.wide <- time.series.wide[,-1]
All <- specnumber(t(time.series.wide), MARGIN = 2)#spe.rich =sobs
DNA <- specnumber(t(time.series.wide[,subset(pathogen.type,Type == 'DNA Virus')$species]), MARGIN = 2)
RNA <- specnumber(t(time.series.wide[,subset(pathogen.type,Type == 'RNA Virus')$species]), MARGIN = 2)
BAC <- specnumber(t(time.series.wide[,subset(pathogen.type,Type == 'Bacteria')$species]), MARGIN = 2)
Euk <- specnumber(t(time.series.wide[,subset(pathogen.type,Type == 'Eukaryote')$species]), MARGIN = 2)

richness <- as.data.frame(cbind(All, DNA,RNA,BAC,Euk))
richness$Individual <- rownames(richness)
richness <- merge(richness,unique(time.series[,c('Individual','Collected_Time')]))
richness.long <- gather(richness,key = 'Type',value = 'Richness',-'Individual',-'Collected_Time')
richness.long$Type <- factor(richness.long$Type,levels = c('All','RNA','DNA','Euk','BAC'))
richness.plot <- ggplot(subset(richness.long,Type != 'All'),aes(Collected_Time,Richness,,color=Collected_Time))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness')+ geom_jitter(size=0.5)+
  stat_compare_means(method = 'kruskal.test',label='p.format')+
  # stat_compare_means(comparisons=times.pair,size=2,
  # label = "p.signif",method = 'wilcox.test')+
  theme_few()+facet_wrap(vars(Type),nrow =1)+
  scale_color_manual(values = setNames(time.names.color,time.names))+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),legend.position = 'bottom',
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
richness.plot
###病原在时间维度上的延续性
medium <- subset(time.series,RPM >0)
medium.patho <- unique(medium[,c('species','family','Type','Collected_Time')])
pathogens.time.dis <- aggregate(medium.patho$Collected_Time,
                                list(medium.patho$species,medium.patho$family,medium.patho$Type),length)
colnames(pathogens.time.dis) <- c('species','family','Type','Num')
pathogens.time.dis$Croee <- ifelse(pathogens.time.dis$Num > 1,'Cross Time pathogens','Single Time')
clade.sum <- aggregate(pathogens.time.dis$species,
                       list(pathogens.time.dis$Type,pathogens.time.dis$Croee),length)
colnames(clade.sum) <- c('Type','Cross','Num')
clade.plot <- ggplot(clade.sum, aes(x=Type,y = Num,fill=Cross)) +
  geom_bar(stat = "identity",
           width = 0.8, size = 0.25,
           alpha = 1)+ 
  labs(x="",y="Num of echo Type")+theme_few()+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))+coord_flip()
family.sum <- aggregate(pathogens.time.dis$species,
                       list(pathogens.time.dis$family,pathogens.time.dis$Croee),length)
colnames(family.sum) <- c('Type','Cross','Num')
family.plot <- ggplot(family.sum, aes(x=Type,y = Num,fill=Cross)) +
  geom_bar(stat = "identity",
           width = 0.8, size = 0.25,
           alpha = 1)+ 
  labs(x="",y="Num of echo Type")+theme_few()+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=60,hjust = 1))+coord_flip()
##绘制不Host同时空的PCA图象
Mouse.sample.sum.wide <- spread(Mouse.sample.sum[,c('Hos_Taxo','Site_Type','Num')],key = Hos_Taxo,value = Num,fill = 0)
rownames(Mouse.sample.sum.wide) <- Mouse.sample.sum.wide$Site_Type
Mouse.sample.sum.wide <- Mouse.sample.sum.wide[,-1]
host.pca <- prcomp(Mouse.sample.sum.wide)
host.pca$x
host.plot_data <- data.frame(host.pca$x)[1:2]
host.plot_data$Site_Type <- rownames(host.plot_data)
names(host.plot_data)[1:2] <- c('PCA1', 'PCA2')
host.plot_data <- merge(host.plot_data, unique(Mouse.sample.sum[,c('Site_Type','Collected_Region','Collected_Time')]), by = 'Site_Type', all.x = TRUE)

Host.PCA.plot <- ggplot(data = host.plot_data, aes(x=PCA1, y=PCA2,color=Collected_Region,fill=Collected_Region,shape=Collected_Time)) +
  geom_point(alpha=.7, size=1)+
  scale_colour_manual(values = setNames(Region_color,Region))+
  scale_fill_manual(values = setNames(Region_color,Region))+
  scale_shape_manual(values=c("2021-Winter"=21,"2022-Autumn"=22,"2022-Spring"=23,
                              "2022-Summer"=24,"2022-Winter"=25))+
  theme_few()+coord_fixed(1.4)+
  theme(text=element_text(size= 6),legend.key.size = unit(1, 'mm'))
##按照物种绘制不同时间区域的病毒的组成情况的PCoA
RPM.table.long.mean <- aggregate(RPM.table.long$RPM,
                                 list(RPM.table.long$species,RPM.table.long$Collected_Time,RPM.table.long$Hos_Taxo,
                                      RPM.table.long$Collected_Region),mean)
colnames(RPM.table.long.mean) <- c('species','Collected_Time','Hos_Taxo','Collected_Region','RPM') 
RPM.table.long.mean$Type <- paste(RPM.table.long.mean$Collected_Time,RPM.table.long.mean$Collected_Region,
                                  RPM.table.long.mean$Hos_Taxo,sep='-')
RPM.table.wide <- spread(RPM.table.long.mean[,c('Type','species','RPM')],key=species,value=RPM,fill=0)
rownames(RPM.table.wide)  <- RPM.table.wide$Type
RPM.table.wide <- RPM.table.wide[,-1]

df.dist = vegdist(RPM.table.wide,method='bray')
df.dist[is.na(df.dist)] <- 0
pcoa <- cmdscale(df.dist, k = (nrow(RPM.table.wide) - 1), eig = TRUE)
plot_data <- data.frame({pcoa$point})[1:2]
plot_data$Type <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')
eig = pcoa$eig
plot_data <- merge(plot_data, unique(RPM.table.long.mean[,c('Type','Collected_Time','Collected_Region','Hos_Taxo')]), 
                   by = 'Type', all.x = TRUE)
Region <- unique(plot_data$Collected_Region)
Region_color <- colors[1:length(Region)]
plot_data$Collected_Time <- factor(plot_data$Collected_Time,levels =c("2021-Winter","2022-Autumn","2022-Spring",
                                                  "2022-Summer","2022-Winter"))

host.time.region.PCoa <- ggplot(data = plot_data, aes(x=PCoA1, y=PCoA2,color=Collected_Region,fill=Collected_Region,shape=Hos_Taxo)) +
  geom_point(alpha=.7, size=1) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  scale_colour_manual(values = setNames(Region_color,Region))+
  scale_fill_manual(values = setNames(Region_color,Region))+
  scale_shape_manual(values=c("Bandicota indica" = 0,
                              "Berylmys sp." = 1,
                              "Mus caroli" = 2,
                              "Niviventer lotipes" = 5 ,
                              "Rattus andamanensis" = 6,
                              "Rattus losea" = 7,
                              "Rattus norvegicus" = 9,
                              "Rattus tanezumi" = 10,
                              "Suncus murinus" =  14))+
  #stat_ellipse(aes(fill = Collected_Region,color=Collected_Region), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  theme_few()+coord_fixed(1.2)+
  theme(text=element_text(size= 6),legend.key.size = unit(1, 'mm'))

##绘制不同区域的病毒多样性的
Region.Table1 <- subset(RPM.table.long,Collected_Time != '2022-Summer' )
Region.Table1 <- subset(RPM.table.long,Collected_Time != '2022-Spring' )

##绘制不同区域的Upset图象
RPM.region.long <- aggregate(RPM.table.long$RPM,
                               list(RPM.table.long$species,RPM.table.long$Type,RPM.table.long$Collected_Region
                                  ),mean)
colnames(RPM.region.long) <- c('species','Type','Collected_Region','RPM')
pathogens.region.wide <- spread(RPM.region.long[,c('species','Type','Collected_Region','RPM')],key = Collected_Region,value = RPM,fill = 0)
Region.names <- unique(RPM.region.long$Collected_Region)
region.upset.plot <- upset(pathogens.region.wide, Region.names,width_ratio=0.15,height_ratio=1,
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
                                   legend.position = "none")),
                         set_sizes = (
                           upset_set_size()+geom_text(aes(label=..count..),size=4/lwd_pt, hjust=1.1, stat='count')+
                             ylab('Num of Pathogens')+
                             theme(axis.title.x = element_text(size=5),
                                   axis.text = element_text(size=5),
                                   text = element_text(size=5))
                         ),
                         themes=upset_modify_themes(
                           list(
                             'intersections_matrix'=theme(text=element_text(size=5)),
                             'overall_sizes'=theme(axis.text.x=element_text(angle=90))
                           ))
)
###从不同的物种的角度考虑病原在不同的地区的变化规律 - 
region.compare <- aggregate(Region.Table1$RPM,
                            list(Region.Table1$Individual,Region.Table1$species,Region.Table1$Collected_Region,Region.Table1$Hos_Taxo
                            ),mean)
colnames(region.compare) <- c('Individual','species','Region','Host','RPM')
region.compare <- subset(region.compare,Host != 'Berylmys bowersi' & Host != "Mus caroli" &
                           Host != 'Niviventer lotipes')
region.compare.wide <- spread(region.compare[,c('Individual','species','RPM')],key = species,value=RPM)
rownames(region.compare.wide) <- region.compare.wide$Individual
region.compare.wide <- region.compare.wide[,-1]
Richness <- specnumber(t(region.compare.wide), MARGIN = 2)#spe.rich =sobs
Richness.region <- as.data.frame(Richness)
Richness.region$Individual <- rownames(Richness.region)
region.host <- unique(Richness.region$Host)
region.host.color <- color20[1:length(region.host)]
Richness.region <- merge(Richness.region,unique(region.compare[,c('Individual','Region','Host')]))
richness.region.plot1 <- ggplot(Richness.region,aes(Region,Richness,color=Host))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness')+ 
  scale_color_manual(values = setNames(region.host.color,region.host))+
  theme_few()+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
##分开绘制放在附图
pathogens.taxo <- unique(RPM.table.long[,c('species','Type')])
dna.names <- subset(pathogens.taxo,Type == 'DNA Virus')$species
rna.names <- subset(pathogens.taxo,Type == 'RNA Virus')$species
euk.names <- subset(pathogens.taxo,Type == 'Eukaryote')$species
bac.names <- subset(pathogens.taxo,Type == 'Bacteria')$species

RNA_Virus <- specnumber(t(region.compare.wide[,rna.names]), MARGIN = 2)#spe.rich =sobs
DNA_Virus <- specnumber(t(region.compare.wide[,dna.names]), MARGIN = 2)#spe.rich =sobs
Bacteria <- specnumber(t(region.compare.wide[,bac.names]), MARGIN = 2)#spe.rich =sobs
Eukaryote <- specnumber(t(region.compare.wide[,euk.names]), MARGIN = 2)#spe.rich =sobs

all.Shannon <- as.data.frame(cbind(RNA_Virus, DNA_Virus,Bacteria,Eukaryote))
all.Shannon$Individual <- rownames(all.Shannon)
all.Shannon <- merge(all.Shannon,unique(RPM.table.long[,c('Individual','Hos_Taxo','Collected_Region')]))
all.Shannon.long <- gather(all.Shannon,key = 'Pathogens_Type',value = 'Richness',
                           -'Individual',-'Hos_Taxo',-'Collected_Region')
region.type <- unique(all.Shannon.long$Collected_Region)
region.type.color <- colors[1:length(region.type)]
Region.Host.Richness <- ggplot(all.Shannon.long,aes(Collected_Region,Richness,color=Collected_Region))+
  geom_boxplot(outlier.shape = NA)+labs(x="",y='Richness')+ 
  scale_color_manual(values = setNames(region.type.color,region.type))+
  theme_few()+
  facet_grid(Pathogens_Type ~ Hos_Taxo)+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(angle = 90,hjust = 1))
##绘制不同类型的病毒在不同时间以及地区出现的数量规律-流向图 时间-病原-空间

##绘制病毒在不同地区不同物种之间的分布情况
region.compare <- aggregate(Region.Table1$RPM,
                              list(Region.Table1$species,Region.Table1$Collected_Region,Region.Table1$Hos_Taxo
                              ),mean)
colnames(region.compare) <- c('species','Region','Host','RPM')
region.compare <- subset(region.compare,Host != 'Berylmys bowersi' & Host != "Mus caroli" &
                           Host != 'Niviventer lotipes')
region.compare <- merge(region.compare,unique(Region.Table1[,c('species','Type')]))
region.compare$Type <- factor(region.compare$Type,levels=c("RNA Virus",  "DNA Virus" ,"Eukaryote","Bacteria" ))
region.compare <- region.compare[order(region.compare$Type,decreasing = TRUE),]

related.region.host.plot <- ggplot(data = region.compare,aes(x=species,y=log10(RPM+1),fill=Region))+
  geom_bar(stat = "identity",width = 0.8,size = 0.25,position = "fill"
  )+scale_fill_manual(values = setNames(Region_color,Region))+theme_few()+
  geom_point(aes(x=species,y=0,color=Type),shape=15)+
  scale_x_discrete(limits=factor(unique(region.compare$species)))+
  facet_grid( . ~ Host)+ labs(x="")+
  theme(text=element_text(size= 6),
        axis.text.x=element_text(angle=90,hjust = 1))+
  coord_flip()
###统计不同地区的病原的abundance+positivate
focus.pathogens <- data.frame(Genus=c('Angiostrongylus','Chlamydia','Pneumocystis','Pneumocystis','Pneumocystis','Pneumocystis','Bartonella','Bartonella','Rat hepatitis E virus','Betacoronavirus 1','Orthohantavirus seoulense'),
                              species=c('Angiostrongylus cantonensis','Chlamydia muridarum','Pneumocystis wakefieldiae','Pneumocystis carinii','Pneumocystis murina','Pneumocystis sus','Bartonella kosoyi','Bartonella sp','Rat hepatitis E virus','Betacoronavirus 1','Orthohantavirus seoulense'))
Region.Table1$species <- ifelse(Region.Table1$species == 'Wnezhou virus','Wenzhou virus',Region.Table1$species)
focus.pathogens <- c(
  'Norway rat pestivirus','Wenzhou virus','Orthohantavirus seoulense','Rat hepatitis E virus','Mischivirus E','Guangdong rosa-like virus','Betacoronavirus 1',
  'Guangdong rodent dependoparvovirus 1','Rat minute virus 2a','Porcine bocavirus',
  'Bartonella kosoyi','Chlamydia muridarum','Klebsiella variicola',
  'Angiostrongylus cantonensis','Brachylaima sp','Pneumocystis carinii')
selected.focus <- subset(Region.Table1,species %in% focus.pathogens)
abundance.focue <- subset(selected.focus,RPM > 0)
##Positivate
selected.focus1 <- aggregate(selected.focus$RPM, by = list(selected.focus$Individual,selected.focus$species,selected.focus$Collected_Region), 
                             FUN = sum)
colnames(selected.focus1) <- c('Individual','species','Region','RPM')
selected.focus1$Value <- ifelse(selected.focus1$RPM>0,'Positive','Negative')
mytable <- xtabs(~ Region + Value + species, data=selected.focus1)
hhhh <- ftable(mytable)
hhhh <- as.data.frame(hhhh)
hhhh_1 <- aggregate(hhhh$Freq, by = list(hhhh$Region,hhhh$species), FUN = sum)
colnames(hhhh_1) <- c("Region","species","Sum")
hhhh_1 <- subset(merge(hhhh,hhhh_1,all.x = TRUE),Value == "Positive")
hhhh_1 <- subset(hhhh_1,Freq >0 )
hhhh_1$Positiveta_Rate <- hhhh_1$Freq/hhhh_1$Sum*100
region.name <- unique(abundance.focue$Collected_Region)
region.name.color <- c('#43978F','#9EC4BE','#ABD0F1','#DCE9F4','#E56F5E','#F19685','#F6C957','#FFB77F','#FBE8D5')
abundance.focue$species <- factor(abundance.focue$species,levels = focus.pathogens)
hhhh_1$species <- factor(hhhh_1$species,levels = focus.pathogens)
Region.view <- ggplot()+
  geom_boxplot(data=abundance.focue,aes(Collected_Region,log10(RPM+1),color=Collected_Region),outlier.shape = NA)+
  labs(x="",y='Positivate rate(%) log10(RPM+1)')+ 
  geom_jitter(data=abundance.focue,aes(Collected_Region,log10(RPM+1),color=Collected_Region),size=0.5)+
  stat_compare_means(method = 'kruskal.test',label='p.format')+
  geom_bar(data = hhhh_1,aes(Region,-(Positiveta_Rate/10),fill=Region),stat = "identity")+
  scale_fill_manual(values = setNames(region.name.color,region.name))+
  scale_color_manual(values = setNames(region.name.color,region.name))+
  # stat_compare_means(comparisons=times.pair,size=2,
  # label = "p.signif",method = 'wilcox.test')+
  theme_few()+facet_wrap(vars(species),nrow =2)+
  theme(text=element_text(size= 5),legend.key.size = unit(3, 'mm'),
        axis.text = element_text(size = 5),axis.text.x = element_text(size = 5,angle=90,hjust=1,vjust=0))
Region.view
###Merge Figures
P1 <- plot_grid(Richness.pie,dev_prop_plot.sharedvirus, align = "h", axis = "bt",nrow = 1,rel_widths = c(1,1))+
  draw_plot_label(
    c("a", "b"),
    c(0, 1/2),
    c(1, 1),
    size = 8
  )
P1

P2 <- plot_grid(richness.region.plot,richness.time.plot,CPC3.plot,shared.virus.coi.dis, 
                align = "hv", axis = "bt",nrow = 1,rel_widths = c(1,1,1,1))+
  draw_plot_label(
    c("c", "d",'e','f'),
    c(0, 1/4,1/2,3/4),
    c(1, 1,1,1),
    size = 8
  )
P2
P3 <- plot_grid(Pathogen.RPM.plot,richness.plot+theme(legend.position = 'none'), align = "h",axis = "bt",ncol = 1,rel_heights = c(1,1))+
  draw_plot_label(
    c("g", "h"),
    c(0, 1/3),
    c(1, 1),
    size = 8
  )
P3
P4 <- plot_grid(pathogen.focus.time.polt+theme(legend.position = 'top'), align = "hv",
                axis = "bt",
                nrow = 1,rel_widths = c(4,5))+
  draw_plot_label(
    c("g"),
    c(0),
    c(1),
    size = 8
  )
P4
P5 <- plot_grid(Region.view+theme(legend.position = 'top'), align = "hv",
                axis = "bt",
                nrow = 1,rel_widths = c(4,5))+
  draw_plot_label(
    c("h"),
    c(0),
    c(1),
    size = 8
  )
P5

final <- plot_grid(NULL, NULL, P4,P5, align = "hv", axis = "bt",ncol = 1,rel_heights = c(1,1.2,1.3,1.2))
final
##save pictures
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5_V2.4.pdf',sep=''),
  final,
  width = 180,             # 宽
  height = 300,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5.1_hhh.pdf',sep=''),
  pathogen.focus.time.polt+theme(legend.position = 'top'),
  width = 120,             # 宽
  height = 90,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5.3_hhh.pdf',sep=''),
  P3,
  width = 90,             # 宽
  height = 120,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5.2_hhh.pdf',sep=''),
  P1+theme(legend.position = 'top'),
  width = 180,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5.4_hhh.pdf',sep=''),
  Region.view,
  width = 180,             # 宽
  height = 60,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5.6_hhh.pdf',sep=''),
  P2,scale = TRUE,
  width = 180,             # 宽
  height = 50,            # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.1 All.time.related.pathogens.pdf',sep=''),
  relate.rpm.plot,
  width = 180,             # 宽
  height = 120,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.2 Time.pathogens.relative.pdf',sep=''),
  Time_plot,
  width = 60,             # 宽
  height = 120,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.3 Pathogens Changes.pdf',sep=''),
  p10.sum.view,
  width = 180,             # 宽
  height = 80,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.4 Region.Pathogens.pdf',sep=''),
  Region.plot,
  width = 180,             # 宽
  height = 120,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.4 Time.rpm.compare.pdf',sep=''),
  Pathogen.dae.RPM.plot,
  width = 180,             # 宽
  height = 230,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.6 Host.region.richness.pdf',sep=''),
  Region.Host.Richness,
  width = 180,             # 宽
  height = 120,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)
ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure s5.5 Region.host.pathogens.relative.pdf',sep=''),
  related.region.host.plot,
  width = 180,             # 宽
  height = 150,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

ggsave(
  filename = paste('F:/Mouse_Result_24/Figure/TMP/','Figure 5 Test.pdf',sep=''),
  time.upset.plot,
  width = 180,             # 宽
  height = 150,           # 高
  units = "mm",          # 单位
  dpi = 300,              # 分辨率DPI
  limitsize = FALSE
)

