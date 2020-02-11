#Illumina Data analysis in R - EAT project

#_________________________________________________________________________________________________
#Oxygen consumption plot

oLight <- read.csv("~/Downloads/R-analysis/EAT_O2concentration-Light.csv")

oDark <- read.csv("~/Downloads/R-analysis/EAT_O2concentration-Dark.csv")

oLight[,1:2] <- NULL
oDark[,1:2] <- NULL

library(reshape2)
melt.oLight <-melt(oLight, id.vars = "Time",value.name = "Oxygen",variable.name = "Incubation")
melt.oDark <-melt(oDark, id.vars = "Time",value.name = "Oxygen",variable.name = "Incubation")

melt.oLight$LightSetting <- "Light"
melt.oDark$LightSetting <- "Dark"

melt.oDark$IncubationType <- "Seawater"
melt.oDark$IncubationType[334:666] <- "Fibre"
melt.oDark$IncubationType[667:999] <- "HDPE"
melt.oDark$IncubationType[1000:1332] <- "Tween"

melt.oLight$IncubationType <- "Seawater"
melt.oLight$IncubationType[334:666] <- "Fibre"
melt.oLight$IncubationType[667:999] <- "HDPE"
melt.oLight$IncubationType[1000:1332] <- "Tween"

#To make a mean line over the geom_point plot
library(dplyr) 
summaryLight <- melt.oLight %>% group_by(Time, LightSetting, IncubationType) %>% 
  summarise(Oxygen=mean(Oxygen))
summaryDark <- melt.oDark %>% group_by(Time, LightSetting, IncubationType) %>% 
  summarise(Oxygen=mean(Oxygen))

oTotal <- rbind(melt.oDark,melt.oLight)
summaryTotal <- rbind(summaryDark,summaryLight)

library(ggplot2)
ggplot(oTotal, aes(x = Time, y = Oxygen, color=LightSetting)) +
  geom_point(size=0.5) +
  facet_grid(IncubationType~.) +
  geom_line(data=summaryTotal, aes(x = Time, y = Oxygen), size = 0.8)


#_________________________________________________________________________________________________
#qPCR barplot

qPCR <- read.csv("~/Downloads/R-analysis/qPCR16S.csv")

qPCR[,5:9] <- NULL
qPCR <- qPCR[-c(5:8),]

library(ggplot2)
library(scales)
ggplot(qPCR, aes(x=IncubationType, y=X16S.copies.per.mL.incubation, col=LightSetting)) + 
  geom_boxplot() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  coord_flip()



#_________________________________________________________________________________________________
#Abundance OTUs from phylogenetic reconstruction

OTUtable<-read.csv("~/Downloads/R-analysis/EAT_otutable-taxonomy-HEATmapTOP10.csv",header=T)

NormOTUs <- data.frame(OTUtable)

NormOTUs$AquariumD1 <- rowMeans(NormOTUs[c('EAT01', 'EAT02', 'EAT03', 'EAT04')], na.rm=TRUE)
NormOTUs$AquariumD4 <- rowMeans(NormOTUs[c('EAT05', 'EAT07')], na.rm=TRUE)
NormOTUs$SWL <- rowMeans(NormOTUs[c('EAT09', 'EAT10','EAT11', 'EAT12','EAT13', 'EAT14')], na.rm=TRUE)
NormOTUs$FibreL <- rowMeans(NormOTUs[c('EAT15', 'EAT16','EAT17', 'EAT18','EAT19', 'EAT20')], na.rm=TRUE)
NormOTUs$HDPEL <- rowMeans(NormOTUs[c('EAT21', 'EAT22','EAT23', 'EAT24','EAT25', 'EAT26')], na.rm=TRUE)
NormOTUs$TweenL <- rowMeans(NormOTUs[c('EAT27', 'EAT28','EAT29', 'EAT30','EAT31', 'EAT32')], na.rm=TRUE)
NormOTUs$SWD <- rowMeans(NormOTUs[c('EAT33', 'EAT34','EAT35', 'EAT36','EAT38')], na.rm=TRUE)
NormOTUs$FibreD <- rowMeans(NormOTUs[c('EAT39', 'EAT40','EAT41', 'EAT42','EAT43', 'EAT44')], na.rm=TRUE)
NormOTUs$HDPED <- rowMeans(NormOTUs[c('EAT45', 'EAT46','EAT47', 'EAT48','EAT49', 'EAT50')], na.rm=TRUE)
NormOTUs$TweenD <- rowMeans(NormOTUs[c('EAT51', 'EAT52','EAT53','EAT55', 'EAT56')], na.rm=TRUE)

NormOTUs[,3:54] <- NULL
NormOTUs$AquariumD1 <- NULL
NormOTUs$AquariumD4 <- NULL
NormOTUs$Taxonomy <- NULL

NormOTUs[,2:9] <- round(NormOTUs[,2:9])

NormOTUs[NormOTUs < 10] <- 0

OTUs.interest <- c("Otu1","Otu12","Otu4","Otu34","Otu17","Otu35","Otu36","Otu28","Otu14","Otu63",
                   "Otu15","Otu6","Otu29","Otu2","Otu25","Otu9","Otu26","Otu70","Otu24","Otu77",
                   "Otu39","Otu10","Otu13","Otu5","Otu18","Otu22","Otu27","Otu21","Otu45","Otu8",
                   "Otu3","Otu20","Otu16","Otu49","Otu19","Otu23")

NormOTUs.OTU <- NormOTUs[NormOTUs$X.OTU.ID %in% OTUs.interest,]

library(reshape2)
NormOTUs.OTU.melt <-melt(NormOTUs.OTU, id.vars = "X.OTU.ID", 
                             value.name = "reads",variable.name = "treatments")

NormOTUs.OTU.melt$lightcondition <- rep("Light", times=nrow(NormOTUs.OTU.melt))
NormOTUs.OTU.melt[(nrow(NormOTUs.OTU)*4+1):(nrow(NormOTUs.OTU)*8),4] <- "Dark"
NormOTUs.OTU.melt$groups <- rep("SW", times=nrow(NormOTUs.OTU.melt))
NormOTUs.OTU.melt[c((nrow(NormOTUs.OTU)*1+1):(nrow(NormOTUs.OTU)*2),(nrow(NormOTUs.OTU)*5+1):(nrow(NormOTUs.OTU)*6)),5] <- "Fibre"
NormOTUs.OTU.melt[c((nrow(NormOTUs.OTU)*2+1):(nrow(NormOTUs.OTU)*3),(nrow(NormOTUs.OTU)*6+1):(nrow(NormOTUs.OTU)*7)),5] <- "HDPE"
NormOTUs.OTU.melt[c((nrow(NormOTUs.OTU)*3+1):(nrow(NormOTUs.OTU)*4),(nrow(NormOTUs.OTU)*7+1):(nrow(NormOTUs.OTU)*8)),5] <- "Tween"


library(ggplot2)
library(scales)

ggplot(NormOTUs.OTU.melt, aes(x = X.OTU.ID, y = reads, fill=groups)) +
  geom_bar(position="fill", stat = "identity") + 
  scale_y_continuous(name = "Reads") +
  coord_flip() +
  theme(axis.text = element_text(colour = "black", size=8), 
        axis.title.y = element_text(colour = "black", face = "bold"), 
        axis.title.x = element_text(colour = "black", face = "bold"), 
        plot.title = element_text(colour = "black", face = "bold"), 
        strip.text = element_text(colour = "white", face = "bold", size=12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
        strip.background = element_rect(fill = "black", colour = "black", size = 0.2),
        axis.line = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  #scale_fill_viridis(discrete = TRUE, option = "D") + 
  facet_grid(lightcondition~., scale="free", space="free_y")


#____________________________________________________________________________________________________________
#For NMDS plot
OTUtable<-read.csv("~/Downloads/R-analysis/EAT_otutable.csv",header=T, row.names = 1)

NormOTUs <- data.frame(OTUtable)

#to remove all OTUS which for all incubation are below 100
NormOTUs <- NormOTUs[apply(NormOTUs[c(1:52)], 1, function(x) !all(x<100)),]

NormOTUs$AquariumD1 <- rowMeans(NormOTUs[c('EAT01', 'EAT02', 'EAT03', 'EAT04')], na.rm=TRUE)
NormOTUs$AquariumD4 <- rowMeans(NormOTUs[c('EAT05', 'EAT07')], na.rm=TRUE)
NormOTUs$SW1L <- rowMeans(NormOTUs[c('EAT09', 'EAT10')], na.rm=TRUE)
NormOTUs$SW2L <- rowMeans(NormOTUs[c('EAT11', 'EAT12')], na.rm=TRUE)
NormOTUs$SW3L <- rowMeans(NormOTUs[c('EAT13', 'EAT14')], na.rm=TRUE)
NormOTUs$Fibre1L <- rowMeans(NormOTUs[c('EAT15', 'EAT16')], na.rm=TRUE)
NormOTUs$Fibre2L <- rowMeans(NormOTUs[c('EAT17', 'EAT18')], na.rm=TRUE)
NormOTUs$Fibre3L <- rowMeans(NormOTUs[c('EAT19', 'EAT20')], na.rm=TRUE)
NormOTUs$HDPE1L <- rowMeans(NormOTUs[c('EAT21', 'EAT22')], na.rm=TRUE)
NormOTUs$HDPE2L <- rowMeans(NormOTUs[c('EAT23', 'EAT24')], na.rm=TRUE)
NormOTUs$HDPE3L <- rowMeans(NormOTUs[c('EAT25', 'EAT26')], na.rm=TRUE)
NormOTUs$Tween1L <- rowMeans(NormOTUs[c('EAT27', 'EAT28')], na.rm=TRUE)
NormOTUs$Tween2L <- rowMeans(NormOTUs[c('EAT29', 'EAT30')], na.rm=TRUE)
NormOTUs$Tween3L <- rowMeans(NormOTUs[c('EAT31', 'EAT32')], na.rm=TRUE)
NormOTUs$SW1D <- rowMeans(NormOTUs[c('EAT33', 'EAT34')], na.rm=TRUE)
NormOTUs$SW2D <- rowMeans(NormOTUs[c('EAT35', 'EAT36')], na.rm=TRUE)
NormOTUs$SW3D <- NormOTUs$EAT38
NormOTUs$Fibre1D <- rowMeans(NormOTUs[c('EAT39', 'EAT40')], na.rm=TRUE)
NormOTUs$Fibre2D <- rowMeans(NormOTUs[c('EAT41', 'EAT42')], na.rm=TRUE)
NormOTUs$Fibre3D <- rowMeans(NormOTUs[c('EAT43', 'EAT44')], na.rm=TRUE)
NormOTUs$HDPE1D <- rowMeans(NormOTUs[c('EAT45', 'EAT46')], na.rm=TRUE)
NormOTUs$HDPE2D <- rowMeans(NormOTUs[c('EAT47', 'EAT48')], na.rm=TRUE)
NormOTUs$HDPE3D <- rowMeans(NormOTUs[c('EAT49', 'EAT50')], na.rm=TRUE)
NormOTUs$Tween1D <- rowMeans(NormOTUs[c('EAT51', 'EAT52')], na.rm=TRUE)
NormOTUs$Tween2D <- NormOTUs$EAT53
NormOTUs$Tween3D <- rowMeans(NormOTUs[c('EAT55', 'EAT56')], na.rm=TRUE)

NormOTUs[,1:52] <- NULL

#to transform --> inversion columns and lines
NormOTUs_t <- t(NormOTUs)

library(vegan)
NormOTUs_percent_t <- decostand(NormOTUs_t, "total") #VERY IMPORTANT TO TRANSFORM BEFORE DOING DECOSTAND, it acts on rows and not columns.

ord<-metaMDS(NormOTUs_percent_t)

#groups for plots
groups<-c("Aquarium","Aquarium","SWL","SWL","SWL","FibreL","FibreL","FibreL","HDPEL","HDPEL","HDPEL","TweenL","TweenL","TweenL",
          "SWD","SWD","SWD","FibreD","FibreD","FibreD","HDPED","HDPED","HDPED","TweenD","TweenD","TweenD")


#for plot: 
plot(ord, disp="sites", type = "n")
ordihull (ord, groups, col = "aquamarine")
ordispider(ord, groups, col="aquamarine", label = TRUE, size=8)
points(ord, disp="sites", pch= c(19), bg=as.numeric(groups), cex=1.3)

#Significance is the p value --> in my case: p=0.003
anosim(NormOTUs_percent_t, groups) #analysis of similarity 



#_____________________________________________________________________________________
#Heatmap
OTUtable1<-read.csv("~/Downloads/R-analysis/EAT_otutable-taxonomy-HEATmapTOP10.csv",header=T, row.names = 1)


#Just to check families
library(stringr)
NormOTUs_taxo <- str_split_fixed(OTUtable1$Taxonomy, ";", 6)

NormOTUs_tot <- cbind(NormOTUs_taxo, OTUtable1)

NormOTUs_tot$Taxonomy <- NULL

colnames(NormOTUs_tot)[1] <- "kingdoms"
colnames(NormOTUs_tot)[2] <- "phylums"
colnames(NormOTUs_tot)[3] <- "classes"
colnames(NormOTUs_tot)[4] <- "orders"
colnames(NormOTUs_tot)[5] <- "families"
colnames(NormOTUs_tot)[6] <- "genus-species"

NormOTUs <- data.frame(OTUtable1)

NormOTUs$AquariumD1 <- rowMeans(NormOTUs[c('EAT01', 'EAT02', 'EAT03', 'EAT04')], na.rm=TRUE)
NormOTUs$AquariumD4 <- rowMeans(NormOTUs[c('EAT05', 'EAT07')], na.rm=TRUE)
NormOTUs$SW1L <- rowMeans(NormOTUs[c('EAT09', 'EAT10')], na.rm=TRUE)
NormOTUs$SW2L <- rowMeans(NormOTUs[c('EAT11', 'EAT12')], na.rm=TRUE)
NormOTUs$SW3L <- rowMeans(NormOTUs[c('EAT13', 'EAT14')], na.rm=TRUE)
NormOTUs$Fibre1L <- rowMeans(NormOTUs[c('EAT15', 'EAT16')], na.rm=TRUE)
NormOTUs$Fibre2L <- rowMeans(NormOTUs[c('EAT17', 'EAT18')], na.rm=TRUE)
NormOTUs$Fibre3L <- rowMeans(NormOTUs[c('EAT19', 'EAT20')], na.rm=TRUE)
NormOTUs$HDPE1L <- rowMeans(NormOTUs[c('EAT21', 'EAT22')], na.rm=TRUE)
NormOTUs$HDPE2L <- rowMeans(NormOTUs[c('EAT23', 'EAT24')], na.rm=TRUE)
NormOTUs$HDPE3L <- rowMeans(NormOTUs[c('EAT25', 'EAT26')], na.rm=TRUE)
NormOTUs$Tween1L <- rowMeans(NormOTUs[c('EAT27', 'EAT28')], na.rm=TRUE)
NormOTUs$Tween2L <- rowMeans(NormOTUs[c('EAT29', 'EAT30')], na.rm=TRUE)
NormOTUs$Tween3L <- rowMeans(NormOTUs[c('EAT31', 'EAT32')], na.rm=TRUE)
NormOTUs$SW1D <- rowMeans(NormOTUs[c('EAT33', 'EAT34')], na.rm=TRUE)
NormOTUs$SW2D <- rowMeans(NormOTUs[c('EAT35', 'EAT36')], na.rm=TRUE)
NormOTUs$SW3D <- NormOTUs$EAT38
NormOTUs$Fibre1D <- rowMeans(NormOTUs[c('EAT39', 'EAT40')], na.rm=TRUE)
NormOTUs$Fibre2D <- rowMeans(NormOTUs[c('EAT41', 'EAT42')], na.rm=TRUE)
NormOTUs$Fibre3D <- rowMeans(NormOTUs[c('EAT43', 'EAT44')], na.rm=TRUE)
NormOTUs$HDPE1D <- rowMeans(NormOTUs[c('EAT45', 'EAT46')], na.rm=TRUE)
NormOTUs$HDPE2D <- rowMeans(NormOTUs[c('EAT47', 'EAT48')], na.rm=TRUE)
NormOTUs$HDPE3D <- rowMeans(NormOTUs[c('EAT49', 'EAT50')], na.rm=TRUE)
NormOTUs$Tween1D <- rowMeans(NormOTUs[c('EAT51', 'EAT52')], na.rm=TRUE)
NormOTUs$Tween2D <- NormOTUs$EAT53
NormOTUs$Tween3D <- rowMeans(NormOTUs[c('EAT55', 'EAT56')], na.rm=TRUE)

NormOTUs[,1:53] <- NULL

#create log transformation and check data
NormOTUs$AquariumD1 <- log10(NormOTUs$AquariumD1)
NormOTUs$AquariumD4 <- log10(NormOTUs$AquariumD4)
NormOTUs$SW1L <- log10(NormOTUs$SW1L)
NormOTUs$SW2L <- log10(NormOTUs$SW2L)
NormOTUs$SW3L <- log10(NormOTUs$SW3L)
NormOTUs$Fibre1L <- log10(NormOTUs$Fibre1L)
NormOTUs$Fibre2L <- log10(NormOTUs$Fibre2L)
NormOTUs$Fibre3L <- log10(NormOTUs$Fibre3L)
NormOTUs$HDPE1L <- log10(NormOTUs$HDPE1L)
NormOTUs$HDPE2L <- log10(NormOTUs$HDPE2L)
NormOTUs$HDPE3L <- log10(NormOTUs$HDPE3L)
NormOTUs$Tween1L <- log10(NormOTUs$Tween1L)
NormOTUs$Tween2L <- log10(NormOTUs$Tween2L)
NormOTUs$Tween3L <- log10(NormOTUs$Tween3L)
NormOTUs$SW1D <- log10(NormOTUs$SW1D)
NormOTUs$SW2D <- log10(NormOTUs$SW2D)
NormOTUs$SW3D <- log10(NormOTUs$SW3D)
NormOTUs$Fibre1D <- log10(NormOTUs$Fibre1D)
NormOTUs$Fibre2D <- log10(NormOTUs$Fibre2D)
NormOTUs$Fibre3D <- log10(NormOTUs$Fibre3D)
NormOTUs$HDPE1D <- log10(NormOTUs$HDPE1D)
NormOTUs$HDPE2D <- log10(NormOTUs$HDPE2D)
NormOTUs$HDPE3D <- log10(NormOTUs$HDPE3D)
NormOTUs$Tween1D <- log10(NormOTUs$Tween1D)
NormOTUs$Tween2D <- log10(NormOTUs$Tween2D)
NormOTUs$Tween3D <- log10(NormOTUs$Tween3D)

NormOTUs[NormOTUs == "-Inf"] <- 0
NormOTUs_t <- t(NormOTUs)

mydf <- data.frame(row.names = row.names(NormOTUs_t), Treatment = c("Aquarium","Aquarium","Seawater","Seawater","Seawater","Fibre",
                                                                            "Fibre","Fibre","HDPE","HDPE","HDPE","Tween","Tween","Tween",
                                                                            "Seawater","Seawater","Seawater","Fibre","Fibre","Fibre","HDPE",
                                                                            "HDPE","HDPE","Tween","Tween","Tween"),
                   LightSetting = c("Aquarium","Aquarium","Light","Light","Light","Light",
                                    "Light","Light","Light","Light","Light","Light","Light","Light",
                                    "Dark","Dark","Dark","Dark","Dark","Dark","Dark","Dark","Dark",
                                    "Dark","Dark","Dark"))

mydf2 <- data.frame(NormOTUs_tot[,1:2])
rows <- c(4,8,12,15,17,23,26,32,53,59,61,64
          ,76,89,93,98,100,111,112,119,133,136)
rows2 <- c(6,19,24,58,59,74,111,117,121,127,134)

mydf2$classes <- as.character(as.character(mydf2$classes))
mydf2[rows,2] <- NA
mydf2[is.na(mydf2$classes),2] <- "Unassigned"
mydf2[rows2,2] <- "Archaea"


mydf2[,1] <- NULL


library(RColorBrewer)
breaksList = seq(0, 20, by = 1)
library(pheatmap)
pheatmap(NormOTUs_t, color = colorRampPalette((brewer.pal(n = 9, name = "BuPu")))
         (length(breaksList)), border_color=NA,
         cluster_cols = TRUE, cluster_rows = TRUE,annotation_row = mydf, annotation_col = mydf2)


#_____________________________________________________________________________________
#Stats about OTU data

OTUtable1<-read.csv("~/Downloads/R-analysis/EAT_otutable-taxonomy.csv",header=T, row.names = 1)

#Sample Day4 not used - so removal
OTUtable1[,c(6,7)] <- NULL

#total number of reads
total_per_sample <- colSums(OTUtable1[,2:51])
sum(total_per_sample)

#second quality control
OTUtable1 <- OTUtable1[apply(OTUtable1[,2:51], 1, function(x) !all(x<10)),]

#number of unassigned OTUs
unassigned <- subset(OTUtable1, Taxonomy == "Unassigned")


#_____________________________________________________________________________________
#Arrangment of the date frames - Venn diagram

#Create a data frame from CSV file
OTUtable<-read.csv("~/Downloads/R-analysis/EAT_otutable-taxonomy-HEATmapTOP10.csv",header=T)

NormOTUs <- data.frame(OTUtable)

NormOTUs$AquariumD1 <- rowMeans(NormOTUs[c('EAT01', 'EAT02', 'EAT03', 'EAT04')], na.rm=TRUE)
NormOTUs$AquariumD4 <- rowMeans(NormOTUs[c('EAT05', 'EAT07')], na.rm=TRUE)
NormOTUs$SWL <- rowMeans(NormOTUs[c('EAT09', 'EAT10','EAT11', 'EAT12','EAT13', 'EAT14')], na.rm=TRUE)
NormOTUs$FibreL <- rowMeans(NormOTUs[c('EAT15', 'EAT16','EAT17', 'EAT18','EAT19', 'EAT20')], na.rm=TRUE)
NormOTUs$HDPEL <- rowMeans(NormOTUs[c('EAT21', 'EAT22','EAT23', 'EAT24','EAT25', 'EAT26')], na.rm=TRUE)
NormOTUs$TweenL <- rowMeans(NormOTUs[c('EAT27', 'EAT28','EAT29', 'EAT30','EAT31', 'EAT32')], na.rm=TRUE)
NormOTUs$SWD <- rowMeans(NormOTUs[c('EAT33', 'EAT34','EAT35', 'EAT36','EAT38')], na.rm=TRUE)
NormOTUs$FibreD <- rowMeans(NormOTUs[c('EAT39', 'EAT40','EAT41', 'EAT42','EAT43', 'EAT44')], na.rm=TRUE)
NormOTUs$HDPED <- rowMeans(NormOTUs[c('EAT45', 'EAT46','EAT47', 'EAT48','EAT49', 'EAT50')], na.rm=TRUE)
NormOTUs$TweenD <- rowMeans(NormOTUs[c('EAT51', 'EAT52','EAT53','EAT55', 'EAT56')], na.rm=TRUE)

NormOTUs[,c(3:54,56)] <- NULL

NormOTUs[,3:11] <- round(NormOTUs[,3:11])

#To split taxonomy into kingdom, phylum, class, family
library(stringr)
NormOTUs_taxo <- str_split_fixed(NormOTUs$Taxonomy, ";", 6)

#To merge the taxonomy split with the previous dataframe
NormOTUs_tot <- cbind(NormOTUs_taxo, NormOTUs)

#To remove Taxonomy column
NormOTUs_tot$Taxonomy <- NULL


#To rename column 1 to 5
colnames(NormOTUs_tot)[1] <- "kingdoms"
colnames(NormOTUs_tot)[2] <- "phylums"
colnames(NormOTUs_tot)[3] <- "classes"
colnames(NormOTUs_tot)[4] <- "orders"
colnames(NormOTUs_tot)[5] <- "families"
colnames(NormOTUs_tot)[6] <- "genus-species"


NormOTUs_tot_venn <- NormOTUs_tot

#replace all lower values than 10 with 0
NormOTUs_tot_venn[NormOTUs_tot_venn < 10] <- 0
NormOTUs_tot_venn[NormOTUs_tot_venn > 9] <- 1

NormOTUs_tot_venn <- NormOTUs_tot_venn[apply(NormOTUs_tot_venn[c(8:16)], 1, function(x) !all(x<1)),]

#sum columns in a new one
NormOTUs_tot_venn$Aquarium <- NormOTUs_tot_venn$AquariumD1
NormOTUs_tot_venn$SW <- NormOTUs_tot_venn$SWL + NormOTUs_tot_venn$SWD
NormOTUs_tot_venn$Fibre <- NormOTUs_tot_venn$FibreL + NormOTUs_tot_venn$FibreD
NormOTUs_tot_venn$HDPE <- NormOTUs_tot_venn$HDPEL + NormOTUs_tot_venn$HDPED
NormOTUs_tot_venn$Tween <- NormOTUs_tot_venn$TweenL + NormOTUs_tot_venn$TweenD

#Venn Diagram 
#source("http://www.bioconductor.org/biocLite.R")    
#biocLite("limma")
library(limma)

NormOTUs_venn <- NormOTUs_tot_venn[,c("Aquarium", "SW", "Fibre", "HDPE", "Tween")]

vennDiagram(vennCounts(NormOTUs_venn),
            counts.col=c("black"))


#_____________________________________________________________________________________
#Stack columns with ggplot2

OTUtable<-read.csv("~/Downloads/R-analysis/EAT_otutable-taxonomy.csv",header=T, row.names = 1)

NormOTUs <- data.frame(OTUtable)

NormOTUs$SW1L <- rowMeans(NormOTUs[c('EAT09', 'EAT10')], na.rm=TRUE)
NormOTUs$SW2L <- rowMeans(NormOTUs[c('EAT11', 'EAT12')], na.rm=TRUE)
NormOTUs$SW3L <- rowMeans(NormOTUs[c('EAT13', 'EAT14')], na.rm=TRUE)
NormOTUs$Fibre1L <- rowMeans(NormOTUs[c('EAT15', 'EAT16')], na.rm=TRUE)
NormOTUs$Fibre2L <- rowMeans(NormOTUs[c('EAT17', 'EAT18')], na.rm=TRUE)
NormOTUs$Fibre3L <- rowMeans(NormOTUs[c('EAT19', 'EAT20')], na.rm=TRUE)
NormOTUs$HDPE1L <- rowMeans(NormOTUs[c('EAT21', 'EAT22')], na.rm=TRUE)
NormOTUs$HDPE2L <- rowMeans(NormOTUs[c('EAT23', 'EAT24')], na.rm=TRUE)
NormOTUs$HDPE3L <- rowMeans(NormOTUs[c('EAT25', 'EAT26')], na.rm=TRUE)
NormOTUs$Tween1L <- rowMeans(NormOTUs[c('EAT27', 'EAT28')], na.rm=TRUE)
NormOTUs$Tween2L <- rowMeans(NormOTUs[c('EAT29', 'EAT30')], na.rm=TRUE)
NormOTUs$Tween3L <- rowMeans(NormOTUs[c('EAT31', 'EAT32')], na.rm=TRUE)
NormOTUs$SW1D <- rowMeans(NormOTUs[c('EAT33', 'EAT34')], na.rm=TRUE)
NormOTUs$SW2D <- rowMeans(NormOTUs[c('EAT35', 'EAT36')], na.rm=TRUE)
NormOTUs$SW3D <- NormOTUs$EAT38
NormOTUs$Fibre1D <- rowMeans(NormOTUs[c('EAT39', 'EAT40')], na.rm=TRUE)
NormOTUs$Fibre2D <- rowMeans(NormOTUs[c('EAT41', 'EAT42')], na.rm=TRUE)
NormOTUs$Fibre3D <- rowMeans(NormOTUs[c('EAT43', 'EAT44')], na.rm=TRUE)
NormOTUs$HDPE1D <- rowMeans(NormOTUs[c('EAT45', 'EAT46')], na.rm=TRUE)
NormOTUs$HDPE2D <- rowMeans(NormOTUs[c('EAT47', 'EAT48')], na.rm=TRUE)
NormOTUs$HDPE3D <- rowMeans(NormOTUs[c('EAT49', 'EAT50')], na.rm=TRUE)
NormOTUs$Tween1D <- rowMeans(NormOTUs[c('EAT51', 'EAT52')], na.rm=TRUE)
NormOTUs$Tween2D <- NormOTUs$EAT53
NormOTUs$Tween3D <- rowMeans(NormOTUs[c('EAT55', 'EAT56')], na.rm=TRUE)

NormOTUs[,8:53] <- NULL

#Round number to no decimal
NormOTUs[,2:31] <- round(NormOTUs[,2:31])

#To split taxonomy into kingdom, phylum, class, family
library(stringr)
NormOTUs_taxo <- str_split_fixed(NormOTUs$Taxonomy, ";", 6)

#To merge the taxonomy split with the previous dataframe
NormOTUs_tot <- cbind(NormOTUs_taxo, NormOTUs)

#To remove Taxonomy column
NormOTUs_tot$Taxonomy <- NULL

#To rename column 1 to 5
colnames(NormOTUs_tot)[1] <- "kingdoms"
colnames(NormOTUs_tot)[2] <- "phylums"
colnames(NormOTUs_tot)[3] <- "classes"
colnames(NormOTUs_tot)[4] <- "orders"
colnames(NormOTUs_tot)[5] <- "families"
colnames(NormOTUs_tot)[6] <- "genus-species"

NormOTUs_tot$families <- as.character(NormOTUs_tot$families)
NormOTUs_tot$families[NormOTUs_tot$kingdoms == "Unassigned"] <- "Unassigned"
NormOTUs_tot$families[NormOTUs_tot$families == ""] <- "Unclassified family"

NormOTUs_tot[,c(1:4,6)] <- NULL

# Aggregate
# use sum for column Y (the variable you want to aggregate according to X)
OTUtable1 <- cbind(aggregate(EAT01~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(EAT02 ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(EAT03 ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(EAT04 ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(EAT05 ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(EAT07 ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW1L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW2L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW3L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre1L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre2L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre3L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween1L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween2L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween3L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE1L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE2L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE3L ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW1D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW2D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(SW3D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre1D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre2D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Fibre3D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween1D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween2D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(Tween3D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE1D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE2D ~ families, data=NormOTUs_tot, FUN=sum),
                   aggregate(HDPE3D ~ families, data=NormOTUs_tot, FUN=sum))

OTUtable1[,3] <- NULL
OTUtable1[,4] <- NULL
OTUtable1[,5] <- NULL
OTUtable1[,6] <- NULL
OTUtable1[,7] <- NULL
OTUtable1[,8] <- NULL
OTUtable1[,9] <- NULL
OTUtable1[,10] <- NULL
OTUtable1[,11] <- NULL
OTUtable1[,12] <- NULL
OTUtable1[,13] <- NULL
OTUtable1[,14] <- NULL
OTUtable1[,15] <- NULL
OTUtable1[,16] <- NULL
OTUtable1[,17] <- NULL
OTUtable1[,18] <- NULL
OTUtable1[,19] <- NULL
OTUtable1[,20] <- NULL
OTUtable1[,21] <- NULL
OTUtable1[,22] <- NULL
OTUtable1[,23] <- NULL
OTUtable1[,24] <- NULL
OTUtable1[,25] <- NULL
OTUtable1[,26] <- NULL
OTUtable1[,27] <- NULL
OTUtable1[,28] <- NULL
OTUtable1[,29] <- NULL
OTUtable1[,30] <- NULL
OTUtable1[,31] <- NULL

OTUtable1[OTUtable1 < 10] <- 0
OTUtable1 <- OTUtable1[apply(OTUtable1[c(2:31)], 1, function(x) !all(x<1)),]

#to isolate the 20 most abundant families
OTUtable.compact <- OTUtable1[apply(OTUtable1[c(2:31)], 1, function(x) !all(x<1300)),]
#to isolate the other families
library(dplyr)
OTUtable.others <- OTUtable1 %>% anti_join(OTUtable.compact)

others <- colSums(OTUtable.others[,-1])
others <- c("others", others)

OTUtable.compact <- rbind(OTUtable.compact, others)

OTUtable.compact[2:31] <- lapply(OTUtable.compact[2:31], function(x) as.numeric(x))

summed <- colSums(OTUtable.compact[c(19,21:23,25:27),2:31])
OTUtable.compact <- OTUtable.compact[-c(19,21:23,25:27),]
summed <- c("others", summed)
OTUtable.compact <- rbind(OTUtable.compact, summed)

OTUtable.compact[2:31] <- lapply(OTUtable.compact[2:31], function(x) as.numeric(x))

#without others
OTUtable.compact <- OTUtable.compact[-21,]
OTUtable.compact <- OTUtable.compact[match(OTUtable.compact1$families, OTUtable.compact$families),]


#Turn your 'treatment' column into a character vector
OTUtable.compact$families <- as.character(OTUtable.compact$families)
#Then turn it back into a factor with the levels in the correct order
OTUtable.compact$families <- factor(OTUtable.compact$families, levels=unique(OTUtable.compact$families))

library(reshape2)
OTUtable.compact.melt <-melt(OTUtable.compact, id.vars = "families", 
                             value.name = "reads",variable.name = "treatments")

OTUtable.compact.melt$lightcondition <- rep("Aquarium", times=nrow(OTUtable.compact.melt))
OTUtable.compact.melt[(nrow(OTUtable.compact)*6+1):(nrow(OTUtable.compact)*18),4] <- "Light"
OTUtable.compact.melt[(nrow(OTUtable.compact)*18+1):(nrow(OTUtable.compact)*30),4] <- "Dark"
OTUtable.compact.melt$groups <- rep("AquariumDay1", times=nrow(OTUtable.compact.melt))
OTUtable.compact.melt[(nrow(OTUtable.compact)*4+1):(nrow(OTUtable.compact)*6),5] <- "AquariumDay4"
OTUtable.compact.melt[c((nrow(OTUtable.compact)*6+1):(nrow(OTUtable.compact)*9),(nrow(OTUtable.compact)*18+1):(nrow(OTUtable.compact)*21)),5] <- "SW"
OTUtable.compact.melt[c((nrow(OTUtable.compact)*9+1):(nrow(OTUtable.compact)*12),(nrow(OTUtable.compact)*21+1):(nrow(OTUtable.compact)*24)),5] <- "Fibre"
OTUtable.compact.melt[c((nrow(OTUtable.compact)*12+1):(nrow(OTUtable.compact)*15),(nrow(OTUtable.compact)*24+1):(nrow(OTUtable.compact)*27)),5] <- "Tween"
OTUtable.compact.melt[c((nrow(OTUtable.compact)*15+1):(nrow(OTUtable.compact)*18),(nrow(OTUtable.compact)*27+1):(nrow(OTUtable.compact)*30)),5] <- "HDPE"

OTUtable.compact.melt <- OTUtable.compact.melt[c(1:80, 121:600),]

library(ggplot2)
library(viridis)     

ggplot(OTUtable.compact.melt, aes(x = treatments, y = reads, fill = families)) +
  geom_bar(position = "fill", stat = "identity") + 
  scale_y_continuous(name = "Reads") +
  coord_flip() +
  theme(axis.text = element_text(colour = "black", size=8), 
        axis.title.y = element_text(colour = "black", face = "bold"), 
        axis.title.x = element_text(colour = "black", face = "bold"), 
        plot.title = element_text(colour = "black", face = "bold"), 
        strip.text = element_text(colour = "white", face = "bold", size=12),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.2),
        strip.background = element_rect(fill = "black", colour = "black", size = 0.2),
        axis.line = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right") +
  facet_grid(groups~., scale="free", space="free_y")
