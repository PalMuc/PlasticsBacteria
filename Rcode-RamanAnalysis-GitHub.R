#Raman spectroscopy of fibres - EAT project 

# plot 40 raman measurement fibres ----------------------------

pathfibre = "~/Downloads/R-analysis/Raman" #put in your actual path where the text files are saved
setwd(pathfibre) #this is to indicate which directory to use
list_fibre = list.files(path=pathfibre, pattern="*.txt", full.names=FALSE) #create list of text files
raman_fibre <- data.frame(raman_shift = seq(100, 1600, by=1500/59411)) #create dataframe with raman shift only

for (i in 1:40) {                         #create a look to integrate all raman text files in the one dataframe
  raman_datafibre <- scan(list_fibre[i])
  
  toDelete <- seq(1, length(raman_datafibre), 2)
  raman_datafibre <- raman_datafibre[-toDelete]
  
  raman_fibre <- cbind(raman_fibre, raman_datafibre)
}
colnames(raman_fibre) <- c("raman_shift", list_fibre) #to give the respective text file name to all columns

#To repeat 5x to reduce the number of measurements to 1858 measurements. 
toDelete2 <- seq(0, nrow(raman_fibre), 2)
last <- nrow(raman_fibre)/2 + 1
toDelete2 <- toDelete2[-last]
raman_fibre <- raman_fibre[-toDelete2,]
toDelete2 <- seq(0, nrow(raman_fibre), 2)
last <- nrow(raman_fibre)/2 + 1
toDelete2 <- toDelete2[-last]
raman_fibre <- raman_fibre[-toDelete2,]
toDelete2 <- seq(0, nrow(raman_fibre), 2)
last <- nrow(raman_fibre)/2 + 1
toDelete2 <- toDelete2[-last]
raman_fibre <- raman_fibre[-toDelete2,]
toDelete2 <- seq(0, nrow(raman_fibre), 2)
last <- nrow(raman_fibre)/2 + 1
toDelete2 <- toDelete2[-last]
raman_fibre <- raman_fibre[-toDelete2,]
toDelete2 <- seq(0, nrow(raman_fibre), 2)
last <- nrow(raman_fibre)/2 + 1
toDelete2 <- toDelete2[-last]
raman_fibre <- raman_fibre[-toDelete2,]



#Plot
library(ggplot2)
library(reshape2)
raman_fibre_40 <- raman_fibre

raman_fibre_40melt <- melt(raman_fibre_40,  id.vars = 'raman_shift', variable.name = 'fibre')

ggplot(raman_fibre_40melt, aes(x=raman_shift, y=value, color=fibre)) +
  geom_line(size=0.2) + scale_x_continuous(breaks=seq(100, 1600, 100)) +
  theme(legend.position="none") + ggtitle("fibre")


#Plot only a few
library(ggplot2)
library(reshape2)
raman_fibre_few <- raman_fibre[,c(1,3,25,40)]

raman_fibre_fewmelt <- melt(raman_fibre_few,  id.vars = 'raman_shift', variable.name = 'fibre')

ggplot(raman_fibre_fewmelt, aes(x=raman_shift, y=value, color=fibre)) +
  geom_line(size=0.2) + scale_x_continuous(breaks=seq(100, 1600, 100)) +
  theme(legend.position="left") + ggtitle("fibre")


#Fibers were manually separated, assigned and added into a tablw, than the Pie chart made
fibre_counts <- read.csv("~/Downloads/R-analysis/fibre_counts.csv", header = T)

slices <- fibre_counts$counts
lbls <- fibre_counts$types
pct <- round(fibre_counts$counts/sum(fibre_counts$counts)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Fibre-counts") 
