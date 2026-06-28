library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(rstatix)
library(summarytools)
library(ggcorrplot)

replace_factor_na <- function(x){
  x <- as.character(x)
  x <- if_else(is.na(x), "n", x)
  x <- as.factor(x)
}

### Analysis of Pathology Scoring on Crohn's specimen ###

#loading results#
pathreports<-read.csv("260122_Term_Scoring_ResearchCollections_CaseControl.csv", header=T, stringsAsFactors = T, na.strings = c("","NA","na"))
#checking for duplicate cases and removing#
cases<-as.data.frame(colnames(pathreports))
cases<-cases %>% separate (`colnames(pathreports)`,c("gi","path","case"),"_")
cases$case[70:128]<-cases$path[70:128]
which(duplicated(cases$case))
dup_cases<-c(101,102, 103, 104, 105, 106, 107, 108, 112, 115, 116, 120, 123, 124, 126, 128)
pathreports_dedup<-pathreports[,-dup_cases]
pathreports<-pathreports_dedup %>% mutate_if (is.factor, replace_factor_na)
pathreports<-t(pathreports)
colnames(pathreports)<-pathreports[1,]
pathreports<-pathreports[-1,]

view(dfSummary(pathreports))

#model matrix to convert categorical data: negative correlation points towards mutually exclusive; postive towards congruent)
model.matrix(~0+., data=as.data.frame(pathreports)) %>% 
  cor(use="pairwise.complete.obs") %>% 
  ggcorrplot(method = "circle",
             show.diag=FALSE, 
             type="lower", 
             lab=TRUE, 
             lab_size=2, 
             ggtheme = theme_bw())

#visually validating strong correlations#
pathreports %>% ggplot () + 
  geom_jitter (aes(x=submucosal_fibrosis, y=submucosal_adipose_change), width=0.2, height=0.2)
pathreports %>% ggplot () + 
  geom_jitter (aes(x=crypt_abscess, y=cryptitis), width=0.2, height=0.2)
pathreports %>% ggplot () + 
  geom_jitter (aes(x=muscular_thickening, y=neuronal_changes), width=0.2, height=0.2)
pathreports %>% ggplot () + 
  geom_jitter (aes(x=granuloma, y=fissure), width=0.2, height=0.2)
pathreports %>% ggplot () + 
  geom_jitter (aes(x=lymphoid_aggregate, y=lymphoid_follicle), width=0.2, height=0.2)

#feature-level stats: gathering terms into categories and counting occurences per case#
features<-as.data.frame(pathreports)
fibrosis<-c("transmural_fibrosis","submucosal_fibrosis","muscularis_propria_fibrosis")
muscular_change<-c("muscular_thickening","muscularis_mucosa_thickening","muscularis_propria_thickening")
acute_inflammation<-c("cryptitis","crypt_abscess","erosion","ulcer") 
chronic_inflammation<-c("lymphoid_aggregate","lymphoid_follicle","granuloma")
all_features<-c("transmural_fibrosis","submucosal_fibrosis","muscularis_propria_fibrosis","muscular_thickening","muscularis_mucosa_thickening","muscularis_propria_thickening","cryptitis","crypt_abscess","erosion","ulcer","lymphoid_aggregate","lymphoid_follicle","granuloma")

features$nfeatures_total<-rowSums((features[,all_features]=="y"))
features$nfeatures_fibrosis<-rowSums(features[,fibrosis]=="y")
features$nfeatures_muscular_change<-rowSums(features[,muscular_change]=="y")
features$nfeatures_acute_inflammation<-rowSums(features[,acute_inflammation]=="y")
features$nfeatures_chronic_inflammation<-rowSums(features[,chronic_inflammation]=="y")

features_table<-features[,1:23]
ggplot(features_table, aes(x="", fill=factor(ulcer))) +
  geom_bar(width = 1) + 
  coord_polar("y", start=0) + 
  scale_fill_manual(values = c("grey","black"))+
  theme_void()
  
features_subset<-features[,24:28]
features_subset$case<-row.names(features_subset)
features_subset.tidy <- features_subset %>% gather("features","count", -case)

features_subset.tidy %>% 
  filter(features != "nfeatures_total") %>%
  ggplot(aes(x=reorder(case,count), y=count, fill=features)) + 
  geom_bar(position="stack", stat="identity")+
  #facet_wrap(.~features)+
  theme(text=element_text(colour="black",face="bold"),
      axis.line = element_line(colour="black"),
      axis.ticks.length =unit(.25,"cm"),
      axis.ticks=element_line(colour="black"),
      panel.grid.major=element_line(colour="white"),
      panel.background=element_rect(fill="white"),
      axis.text.x = element_text(size=8, angle = 90,hjust=1,vjust=0.4,colour="black"),
      axis.text.y = element_text(size=14, colour="black"))

