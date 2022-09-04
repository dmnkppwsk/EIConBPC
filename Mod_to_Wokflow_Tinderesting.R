#choose from which samples the BPC will be created (for example, choose experimental samples "E")
class_to_bpis <- "E"
to_bpis<-which(Experimental$class == class_to_bpis)


#Generate BPCs
for(i in to_bpis){
  bpis <- chromatogram(filterFile(xdata, i), aggregationFun = "max")
  nam <- paste("bpis", i, sep = "_")
  assign(nam, bpis)
}

#Generate additional metadata
Experimental2 <- Experimental
Experimental2$row_name <- rownames(Experimental2)

Experimental2$bpis_name <- ifelse(Experimental2$class == "E", Experimental2$bpis_name <- paste("bpis", Experimental2$row_name, sep = "_"), NA)
Experimental2$df_bpis_name <- ifelse(Experimental2$class == "E", Experimental2$df_bpis_name <- paste("df_bpis", Experimental2$row_name, sep = "_"), NA)  

#Generate and save data frames with BPC chromatogram data
for(i in to_bpis){
  nam_df <- paste("df_bpis", i, sep = "_")
  assign(nam_df, .ChromatogramsLongFormat(get(Experimental2$bpis_name[i])))
  save(nam_df, list=c(paste("df_bpis", i, sep = "_")), file=paste(nam_df, ".Rdata", sep=""))
}

#Generate and save more metadata
library(tidyverse)

short_Experimental<-Experimental2 %>% drop_na()
short_Experimental<-data_frame(short_Experimental$row_name, short_Experimental$class, short_Experimental$time, short_Experimental$Replicate, short_Experimental$bpis_name, short_Experimental$df_bpis_name)
short_Experimental2<-short_Experimental
short_Experimental3<-Experimental2 %>% drop_na()
short_Experimental3$df_filename<-paste(short_Experimental3$df_bpis_name, "Rdata", sep=".")


save(short_Experimental, file="short_Experimental.Rdata")
save(short_Experimental2, file="short_Experimental2.Rdata")
save(short_Experimental3, file="short_Experimental3.Rdata")


#Test without app
sn<-10
mzr<-755.24+c(-0.05, 0.05)
chr_sn <- chromatogram(filterFile(xdata, sn), mz = mzr, aggregationFun = "max")
df_EIC <- .ChromatogramsLongFormat(chr_sn)
dm2<-df_EIC
dm2[is.na(dm2<-df_EIC)] <- 0
df_EIC<-dm2
nam_df_bpis <- paste("df_bpis", sn, sep = "_")
target_df_bpis<-get(nam_df_bpis)
target_df_bpis$intensityEIC <- df_EIC$intensity

ggplot(data = target_df_bpis, aes(rtime/60)) + 
  geom_line(aes(y=intensity, colour="BPC"))+
  geom_line(aes(y=intensityEIC, colour="EIC"))+
  theme_classic() +
  labs(x = "Retention time", y = "Intensity")

####################################################################################################
###                                      EIConBPC                                                ###
####################################################################################################


#if needed
install.packages("DT")
install.packages("plotly")

runApp("EICconBPC.R")
