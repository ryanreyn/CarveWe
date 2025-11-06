#This file contains the relevant functions for supporting the main CarveWe scripts.
#Load all necessary packages
library(gplots)
library(tidyverse)
library(patchwork)
library(vegan)
library(multcompView)

#We need to do some bootstrapping analysis per a request from Cameron. I think what we need to do is bootstrap the genomes in each cluster and then compute sums and relative abundances -> get a CI on relative abundance per cluster per region -> using functions from Ben Temperton
bootstrap_sum_pi <- function(data, func, size=1000){
  bs_sample = sample(data, size=size, replace=TRUE)
  do.call(func, list(x=bs_sample))
}

generate_bootstrap_replicates <- function(data, func, n=1000, size=1000){
  replicate(n, do.call(bootstrap_sum_pi, list(data=data, func=func, size=size)))
}

#We need to perform many iterations of clustering and optimize according to an objective function, so we will use this pre-built helper function to look at the distance between cluster means
clusterMeanDist <- function(clusters){
  cluster.means = c()
  
  for(c in unique(clusters)){
    temp.members <- which(clusters == c)
    
    if(length(temp.members) > 1){
      temp.dist <- som.dist[temp.members,]
      temp.dist <- temp.dist[,temp.members]
      cluster.means <- append(cluster.means, mean(temp.dist))
    }else(cluster.means <- append(cluster.means,0))
  }
  
  return(mean(cluster.means))
  
}

#For code cleanliness and reproducibility we are functionalizing our calls for comparative values between the experimental data and growth tests using models built on the experimental organismal genoems
compute.comparisons <- function(data){
  only_agreements <- data%>%
    filter(agreement_status %in% c("Perfect Agreement","Single Rescue","Double Rescue"))
  
  accuracy <- data%>%
    summarise(accuracy=sum(agreement_status%in%c("Perfect Agreement","Single Rescue","Double Rescue"))/n())
  
  true_positives <- only_agreements%>%
    filter(`Growth/No Growth`==1)%>%
    summarise(count=n())%>%
    as.matrix()
  
  false_positives <- data%>%
    filter(agreement_status=="Model Growth, Exp. No Growth")%>%
    summarise(count=n())%>%
    as.matrix()
  
  false_negatives <- data%>%
    filter(agreement_status%in%c("No Growth, Pathway Present","No Growth, No Pathway"))%>%
    summarise(count=n())%>%
    as.matrix()
  
  precision <- true_positives/(true_positives+false_positives)
  
  recall <- true_positives/(true_positives+false_negatives)
  
  return(list(accuracy, precision, recall))
}

#We perform a variety of statistical significance tests using Tukey's HSD metric and ANOVA models so we will functionalize it here
ANOVA.test <- function(data, vars){
  col1 <- data%>%
    select(vars[1])%>%
    as.matrix()
  col2 <- data%>%
    select(vars[2])%>%
    as.matrix()
  model=lm(col1 ~ col2)
  ANOVA=aov(model)
  
  tukey <- TukeyHSD(x=ANOVA,"col2",conf.level=0.95)
  Tukey.levels <- tukey[["col2"]][,4]
  tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  tukey.labels$treatment <- rownames(tukey.labels)
  tukey.labels <- tukey.labels[order(tukey.labels$treatment),]
  numeric.labels <- match(tukey.labels[,1],unique(tukey.labels[,1]))
  return(list(data, tukey.labels, numeric.labels, tukey))
}