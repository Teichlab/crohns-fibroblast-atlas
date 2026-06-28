library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)


### 260108-LuminexResults-LXSAHM-24 #L163957 ###

#loading Luminex calculated test results#
luminex<-read.csv("260108_LuminexResults_cleaned.csv", header=T, stringsAsFactors = T)
luminex$replicate<-factor(luminex$replicate)
luminex$treatment<-factor(luminex$treatment, levels = c("untreated","cytokines","cytokines_belinostat","cytokines_gskj4"))
luminex.tidy<-gather(luminex, target,concentration, -c(location,sample_name,dilution,replicate,treatment))

luminex.tidy %>% filter (dilution != "NA") %>%
  ggplot(aes(x=treatment, y=concentration, colour=treatment))+ 
  geom_violin()+
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "crossbar", 
               width = 0.25,
               position = position_dodge(width = .25))+
  geom_point(position=position_jitterdodge(),size=1)+
  xlab("")+ylab("CXCL5 (pg/ml)")+labs(shape="epimodulator_treatment")+
  facet_wrap(dilution~target,scales = "free")+
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

stat.test<-luminex.tidy %>% filter (dilution == "1in2") %>% 
  group_by(target) %>%
  anova_test(concentration~treatment)
View(stat.test)

stat.test<-luminex.tidy %>% filter (dilution == "1in2") %>% 
  group_by(target) %>%
  pairwise_t_test(concentration~treatment, p.adjust.method = "fdr")
View(stat.test)
