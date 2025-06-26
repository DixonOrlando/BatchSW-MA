###Generating schematic in R###

library(tidyverse)
library(ggplot2)
library(reshape2)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(Matrix)

SWdesmat <- function(Ts) {
  Xsw <- matrix(data=0, ncol = (Ts), nrow = (Ts-1))
  for(i in 1:(Ts-1)) {
    Xsw[i,(i+1):Ts] <- 1
  }
  return(Xsw)
}
#Batched SW schematic
batchSWscheme <- function(olap, S, Ts,  Kseq){
  #Creates a batched SW design schematic with
  #S batches of Ts-period stepped wedge designs, 
  #with overlap between batches of olap periods
  #and K clusters in each sequence
  
  g<- Ts- olap #gap between each batch and the next
  
  SWbatch <- matrix(data=NA, nrow=S*(Ts-1), ncol=Ts + (S-1)*(Ts-olap))  
  
  for(i in 1:S){
    SWbatch[((i-1)*(Ts-1) + 1):(i*(Ts-1)) , 
            ((i-1)*g +1):((i-1)*g + Ts)] <- SWdesmat(Ts)
    
  }
  
  return( SWbatch[sort(rep(1:nrow(SWbatch), Kseq)), ])
  
  
}

############Stepped wedge trial schematic############
Xdes00 <- SWdesmat(4)
melted_myx00 <- melt(Xdes00)
names(melted_myx00)[names(melted_myx00)=="Var1"] <- "Sequence"
names(melted_myx00)[names(melted_myx00)=="Var2"] <- "Period"

#Run this to get the plot.
ggplot(data =melted_myx00, aes(x=Period, y=Sequence, fill = factor(value))) + 
  geom_tile( colour = "grey50") +
  scale_y_reverse(breaks=c(1:9)) +
  scale_x_continuous(breaks=c(1:12)) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 5) +
  scale_fill_manual(values = color_palette1) + theme(legend.position="none")

############Batched stepped wedge trial schematic############
#Plot for design 1 in the main writing
Xdes00 <- batchSWscheme(olap=0, S=3, Ts=4,  Kseq=1)
melted_myx00 <- melt(Xdes00)
names(melted_myx00)[names(melted_myx00)=="Var1"] <- "Sequence"
names(melted_myx00)[names(melted_myx00)=="Var2"] <- "Period"
color_palette1 <-colorRampPalette(c( "white", "grey"))(2)
ggplot(data =melted_myx00, aes(x=Period, y=Sequence, fill = factor(value))) + 
  geom_tile( colour = "grey50") +
  scale_y_reverse(breaks=c(1:9)) +
  scale_x_continuous(breaks=c(1:12)) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 5) +
  scale_fill_manual(values = color_palette1, na.value = theme_grey()$panel.background$fill) + 
  theme(legend.position="none") +
  ggtitle("Design 1")


#Plot for design 2 in the main writing
Xdes00 <- batchSWscheme(olap=1, S=3, Ts=4,  Kseq=1)
melted_myx00 <- melt(Xdes00)
names(melted_myx00)[names(melted_myx00)=="Var1"] <- "Sequence"
names(melted_myx00)[names(melted_myx00)=="Var2"] <- "Period"
color_palette1 <-colorRampPalette(c( "white", "grey"))(2)
ggplot(data =melted_myx00, aes(x=Period, y=Sequence, fill = factor(value))) + 
  geom_tile( colour = "grey50") +
  scale_y_reverse(breaks=c(1:9)) +
  scale_x_continuous(breaks=c(1:12)) +
  theme(panel.grid.minor = element_blank()) +
  geom_text(aes(Period, Sequence, label = round(value,4)), color = "black", size = 5) +
  scale_fill_manual(values = color_palette1, na.value = theme_grey()$panel.background$fill) + 
  theme(legend.position="none") +
  ggtitle("Design 2")


