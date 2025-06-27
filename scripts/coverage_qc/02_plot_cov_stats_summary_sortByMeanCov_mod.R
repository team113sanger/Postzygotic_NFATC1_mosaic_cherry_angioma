#!/usr/bin/env Rscript 
## This script is used to plot the coverage stats summary for each sample and plot
# Get the working directory
here::i_am("scripts/coverage_qc/02_plot_cov_stats_summary_sortByMeanCov_mod.R")
# Load the required libraries
library(ggplot2)
library(here)
library(viridis)
library(dplyr)
library(reshape2)
library(cowplot)

progname<-"plot_cov_stats_summary_sortByMeanCov_mod.R"
# 1.Read the required information 
# Set he working directory
projdir<-here::here()
qcdir<-file.path(projdir, "results", "qc_plots", "depth")

setwd(projdir)
#Read the list of sample names within the project and folder
slist<-read.csv(file = file.path(projdir, "sample.list"), header = F, stringsAsFactors = F, sep="\t")
slist<-slist$V1

#Read the header for the coverage depth from the stats summary file 
header<-read.csv(file = file.path(projdir, "cov_stats_summary.tsv"), header = F, stringsAsFactors = F, sep="\t")
header<-header[1,]
data<-NULL
i<-NULL
for (i in 1:length(slist)){
    temp<-read.csv(file = file.path(projdir, paste(slist[i], ".covstats.tsv", sep="")), header = F, stringsAsFactors = F, sep="\t")
    data<-rbind(data, c(as.vector(t(temp[2,]))))  
}
data<-cbind(slist,data)
data<-as.data.frame(data)
colnames(data)<-as.vector(t(header[1,]))
data$`11+`<-as.numeric(data$`11+`)
data$`21+`<-as.numeric(data$`21+`)
data$`31+`<-as.numeric(data$`31+`)
data$`41+`<-as.numeric(data$`41+`)
data$`51+`<-as.numeric(data$`51+`)
data$`61+`<-as.numeric(data$`61+`)
data$`71+`<-as.numeric(data$`71+`)
data$`81+`<-as.numeric(data$`81+`)
data$`91+`<-as.numeric(data$`91+`)
data$`101+`<-as.numeric(data$`101+`)
data$`111+`<-as.numeric(data$`111+`)
data$`121+`<-as.numeric(data$`121+`)
data$`131+`<-as.numeric(data$`131+`)
data$`141+`<-as.numeric(data$`141+`)
data$`151+`<-as.numeric(data$`151+`)
data$Cov_Mean<-as.numeric(data$Cov_Mean)
#file="cov_stats_summary.tsv"
#data <- read.table(file, sep = "\t", header = T, check.names = F)
head(data)
mode(data)
class(data)
mode(data$Cov_Mean)
class(data)


data1 <- data %>% select(-Cov_Mean)
data1.melt <- melt(data1)

colnames(data1.melt) <- c("Sample", "Coverage", "Percent")

labels <- data1.melt %>% group_by(Sample) %>% filter(Percent > 80)  %>% top_n(n=-1, wt=Percent)

data2 <- data %>% select(Sample, Cov_Mean)
order <- data2 %>% arrange(Cov_Mean) %>% select(Sample)
order$Sample

data2.melt<-melt(data2)
colnames(data2.melt)<-c("Sample","Variable","Coverage")
data2.melt$Sample<-factor(data2.melt$Sample, levels=order$Sample)
data1.melt$Sample<-factor(data1.melt$Sample, levels=order$Sample)
str(data2.melt)
str(data1.melt)

plot<-ggplot(data1.melt, aes(Coverage, Sample, fill = Percent))+
  geom_tile()+
  scale_fill_viridis(option="turbo", direction=-1, begin=0.15, end=0.85)+
  theme(axis.text=element_text(size=11), axis.title=element_blank(), axis.ticks=element_blank())+
  scale_x_discrete(position = "top", expand=c(0,0))+
  geom_text(data=labels, aes(label=Percent))

quants <- as.numeric(quantile(data2$Cov_Mean ,c(.50, 0.25, .05, .01)))
quants1 <- which.max(sort(data2$Cov_Mean) > quants[1]) + 0.5
quants2 <- which.max(sort(data2$Cov_Mean) > quants[2]) + 0.5
quants3 <- which.max(sort(data2$Cov_Mean) > quants[3]) + 0.5
quants4 <- which.max(sort(data2$Cov_Mean) > quants[4]) + 0.5
quants.pos <- c(quants1, quants2, quants3, quants4)
quants
quants.pos

mean_cov_plot<-ggplot(data2.melt, aes(Variable, Sample, fill=Coverage)) + geom_tile() + 
  scale_fill_viridis(option="turbo", direction=-1, begin=0.15, end=0.85) + 
  scale_x_discrete(position = "top", expand=c(0,0)) + 
  geom_text(aes(label=Coverage)) + 
  theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.text=element_text(size=12))

png(file.path(projdir, "summary_cov_stats_ordered.png"), width=800, height=1000)
plot_grid(plot, mean_cov_plot, align = "h", rel_widths = c(5,1))
dev.off()


#Get the samples which had a coverage >20X across 80% of the baits positions 
#GEt all thje information for 21+
samp_20cov<-data1.melt[data1.melt$Coverage=="21+", ]

#Add the information 
samp_20cov_final<-samp_20cov[ samp_20cov$Percent>79.99, ]
write.table(samp_20cov_final, file = file.path(projdir, "Final_samples_passed_20cov.txt"), quote = F, col.names = T, row.names = F, sep="\t")

