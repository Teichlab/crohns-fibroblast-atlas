library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(summarytools)
library(rstatix)


### importing data from Flowjo analysis ###
flowjo <- read.csv("20269122_Fib_SN_Neutrophil_activation.csv", header=T, stringsAsFactors = T)
view(dfSummary(flowjo))
flowjo$treatment_cytokines<-factor(flowjo$treatment_cytokines, levels=c("unstimulated","stimulated"))
flowjo$treatment_epigenetic<-factor(flowjo$treatment_epigenetic, levels=c("DMSO","BEL"))
flowjo$supernatant_dilution <-factor(flowjo$supernatant_dilution , levels=c("1in2","1in4","1in8","1in16"))
flowjo$treatment_combined<-factor(flowjo$treatment_combined, levels=c("NA_NA_NA","control_medium_NA_NA","fibroblast_SN_unstimulated_DMSO","fibroblast_SN_stimulated_DMSO","fibroblast_SN_stimulated_BEL"))
flowjo.tidy<-gather(flowjo, stat, value, -sample,-timepoint,-supernatant_dilution, -replicate,-treatment,-treatment_cytokines,-treatment_epigenetic,-treatment_combined)
flowjo.tidy$stat<-as.factor(flowjo.tidy$stat)

stat.test<-flowjo.tidy %>% 
  filter(timepoint != "BL") %>% 
  filter (timepoint %in% c("4h",NA)) %>%
  filter (supernatant_dilution %in% c("1in2",NA)) %>% 
  group_by(stat) %>%
  pairwise_t_test(value~treatment_combined,ref.group = "fibroblast_SN_stimulated_DMSO", p.adjust.method = "fdr")
View(stat.test)

#significant differences after 4h in the folliwing stats:
stat_selection_mfi <- c(
                    "granulocytes.Single.Cells.Live.CD14.CD16....Geometric.Mean..CD16.",
                    "granulocytes.Single.Cells.Live.CD14.CD16....Geometric.Mean..CXCR1.",
                    "granulocytes.Single.Cells.Live.CD14.CD16....Geometric.Mean..CXCR2.")
stat_selection_proportion <- c("granulocytes.Single.Cells.Live.CD16.CXCR1hi...Freq..of.Parent....",
                               "granulocytes.Single.Cells.Live.CD16.CXCR1lo...Freq..of.Parent....",
                               "granulocytes.Single.Cells.Live.CD16.CXCR2hi...Freq..of.Parent....",
                               "granulocytes.Single.Cells.Live.CD16.CXCR2lo...Freq..of.Parent....",
                               "granulocytes.Single.Cells.Live.CD16hi...Freq..of.Parent....",
                               "granulocytes.Single.Cells.Live.CD16lo...Freq..of.Parent....")
#plotting to confirm
flowjo.tidy %>% filter(timepoint != "BL") %>% 
  filter (timepoint %in% c("30min",NA)) %>%
  filter (supernatant_dilution %in% c("1in2",NA)) %>%
  filter (stat %in% stat_selection_mfi) %>%
  ggplot(aes(x=treatment_combined, y=value, colour=treatment_combined))+ 
  geom_violin()+
  geom_point(position=position_jitterdodge(),size=1)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, color = "black", position = position_dodge(width = .25))+
  xlab("treatment")+ylab("gMFI")+
  facet_wrap(.~stat,scales = "free")+
  #stat_compare_means(method= "wilcox.test")+
  #scale_fill_manual(values=c("#FF1500","#F7A41D","#363B8E","#FFFFFF"))+
  #geom_text_repel(aes(label=sample))+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4,colour="black"),
        axis.text.y = element_text(colour="black"))

flowjo.tidy %>% filter(timepoint != "BL") %>% 
  filter (timepoint %in% c("30min",NA)) %>%
  filter (supernatant_dilution %in% c("1in2",NA)) %>%
  filter (stat %in% stat_selection_proportion) %>%
  ggplot(aes(x=treatment_combined, y=value, colour=treatment_combined))+ 
  geom_violin()+
  geom_point(position=position_jitterdodge(),size=1)+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 1, color = "black", position = position_dodge(width = .25))+
  xlab("treatment")+ylab("gMFI")+
  facet_wrap(.~stat,scales = "free",ncol=2)+
  #stat_compare_means(method= "wilcox.test")+
  #scale_fill_manual(values=c("#FF1500","#F7A41D","#363B8E","#FFFFFF"))+
  #geom_text_repel(aes(label=sample))+
  theme(text=element_text(size=12,colour="black",face="bold"),
        axis.line = element_line(colour="black"),
        axis.ticks.length =unit(.25,"cm"),
        axis.ticks=element_line(colour="black"),
        panel.grid.major=element_line(colour="white"),
        panel.background=element_rect(fill="white"),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.4,colour="black"),
        axis.text.y = element_text(colour="black"))

