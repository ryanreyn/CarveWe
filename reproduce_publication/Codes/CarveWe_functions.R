#This file contains the relevant functions for supporting the main CarveWe scripts.
#Load all necessary packages
library(gplots)
library(tidyverse)
library(patchwork)
library(vegan)
library(multcompView)
library(rlang)
library(reshape2)

#We need to do some bootstrapping analysis per a request from Cameron. I think what we need to do is bootstrap the genomes in each cluster and then compute sums and relative abundances -> get a CI on relative abundance per cluster per region -> using functions from Ben Temperton
bootstrap_sum_pi <- function(data, func, size=1000){
  bs_sample = sample(data, size=size, replace=TRUE)
  do.call(func, list(x=bs_sample))
}

generate_bootstrap_replicates <- function(data, func, n=1000, size=1000){
  replicate(n, do.call(bootstrap_sum_pi, list(data=data, func=func, size=size)))
}

#We need to perform many iterations of clustering and optimize according to an objective function, so we will use this pre-built helper function to look at the distance between cluster means
clusterMeanDist <- function(clusters, som.dist){
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

kruskal_test <- function(data, vars){
  vars_1 = vars[1]
  vars_2 = vars[2]
  col1 <- data%>%
    select(vars[1])%>%
    as.matrix()
  col2 <- data%>%
    select(vars[2])%>%
    as.matrix()%>%
    as.numeric()
  kruskal.list <- list()
  for (i in 1:max(col2)){
    kruskal.list[[i]] <- data%>%
      mutate(across(all_of(vars_2),as.numeric))%>%
      filter(!!sym(vars_2)==i)%>%
      select({{ vars_1 }})%>%
      as.matrix()
  }
  kruskal.out <- kruskal.test(kruskal.list)
  
  #Run a post-hoc pairwise Wilcox rank sum test
  wilcox.out<-pairwise.wilcox.test(col1, col2, p.adjust.method = "bonferroni")
  p.mat<-wilcox.out$p.value
  p.df<-melt(p.mat, na.rm = TRUE)
  p.convert<-p.df$value%>%as.matrix()%>%as.vector()
  names(p.convert)<-paste0(p.df$Var1, "-", p.df$Var2)
  
  wilcox.labels <- data.frame(multcompLetters(p.convert)['Letters'])
  wilcox.labels$treatment <- rownames(wilcox.labels)
  wilcox.labels <- wilcox.labels[order(wilcox.labels$treatment),]
  numeric.labels <- match(wilcox.labels[,1],unique(wilcox.labels[,1]))
    
  return(list(kruskal.out, wilcox.out, wilcox.labels, numeric.labels))
}

generate.soms <- function(data, grid.info){
  som.model <- som(data,grid=grid.info,rlen=1500,alpha=c(0.025,0.01))
  
  #Comment this line out if you plan to generate a map de novo
  # load("Publication_Data/final-kohonen-som.RData")
  
  #This has been an issue but we need to set the specific seed present *after* running the original SOM model to be able to reproduce the exact results
  # .Random.seed <- good.seed
  
  #Plot the training progress to see and also to ensure right parameters were used
  # plot(som.model,type="changes")
  
  carve.som <- som.model$codes[[1]]
  som.dist <- as.matrix(dist(carve.som))
  
  return(list(som.model, carve.som, som.dist))
}

generate.clusters <- function(data, num_clusters=8){
  carve.som <- data$codes[[1]]
  som.dist <- as.matrix(dist(carve.som))
  try.k <- 2:20
  
  #There's a nebulous issue with defining seeds while generating a markdown file, so at the moment we are loading a pre-saved data object containing the relevant cluster definitions using the correct seeds. You can uncomment this if you wish to run the chunks directly.
  
  cluster.dist.eval <- as.data.frame(matrix(ncol = 3, nrow = (length(try.k))))
  colnames(cluster.dist.eval) <- c('k', 'kmeans', 'hclust')
  
  for(i in 1:length(try.k)) {
    cluster.dist.eval[i, 'k'] <- try.k[i]
    cluster.dist.eval[i, 'kmeans'] <- clusterMeanDist(kmeans(carve.som, centers = try.k[i], iter.max = 20)$cluster, som.dist)
    cluster.dist.eval[i, 'hclust'] <- clusterMeanDist(cutree(hclust(vegdist(carve.som,method="euclidean")), k = try.k[i]), som.dist)
  }
  
  #Plot the kmeans distances as a function of the number of clusters for both clustering methodologies
  dist.plot <- plot(cluster.dist.eval[, 'kmeans'] ~ try.k,
       type = 'l', xlab = "Number of clusters (k)", ylab = "Intracluster distance", cex.axis = 1.5, cex.lab = 1.5)
  
  lines(cluster.dist.eval[, 'hclust'] ~ try.k,
        col = 'red')
  
  legend('topright',
         inset = c(-0.35,-0.05),
         legend = c('k-means', 'hierarchical'),
         col = c('black', 'red'),
         cex = 1.5,
         lty = c(1, 1))
  
  
  cluster.tries <- list()
  
  for(k in num_clusters){
    ## k-means clustering
    
    som.cluster.k <- kmeans(carve.som, centers = k, iter.max = 10000, nstart = 10)$cluster     # k-means
    
    ## hierarchical clustering
    
    som.dist <- dist(carve.som) # hierarchical, step 1
    som.cluster.h <- cutree(hclust(som.dist), k = k) # hierarchical, step 2
    
    ## capture outputs
    
    cluster.tries[[paste0('som.cluster.k.', k)]] <- som.cluster.k
    cluster.tries[[paste0('som.cluster.h.', k)]] <- som.cluster.h
  }
  
  obs_density <- data.frame(node=som.model$unit.classif)%>%group_by(node)%>%summarise(count=n())
  # ggplot(obs_density,aes(x=count))+geom_histogram()+custom_theme
  
  #Want to find the average variance and number of models represented in each cluster
  clusters <- som.cluster.k
  ids <- som.model$unit.classif
  sample_clusters <- som.cluster.k[som.model$unit.classif]
  
  na.nodes <- which(c(1:gridsize^2)%in%unique(obs_density$node)==FALSE)
  som.cluster.k[na.nodes] <- NA
  
  all_vars <- c()
  for (i in 1:k){
    curr_clust=i
    
    curr_models <- which(clusters==curr_clust)
    
    if (length(curr_models)>1){
      curr_vars <- apply(metabolite_matrix[curr_models,],2,sd)%>%mean()
      all_vars <- c(all_vars,curr_vars)
    } else {
      all_vars <- c(all_vars,0)
    }
  }
  
  return(list(som.cluster.k, sample_clusters, dist.plot))
}

process.GRUMP.data <- function(raw_data){
  data <- raw_data%>%
    separate(.,taxonomy,sep=";",into=c("domain","phylum","class","order","family","genus","species"))
  
  total_abund <- data%>%
    select(where(is.double))%>%
    as.matrix()%>%
    rowSums()
  
  data <- data%>%
    mutate(total_abundance = total_abund)
  return(data)
}

#This function takes the raw ASV abundance data and companion BLAST alignment data to genomes/clusters and uses them to create relative abundance dataframes per sampling station
compute.ASV.abundances <- function(blast_data, cluster_data, asv_data){
  merge_data <- left_join(blast_data, cluster_data, by=join_by(genome==value))
  
  piv.asv <- asv_data%>%
    select(OTU_ID, -total_abundance, where(is.double))%>%
    pivot_longer(-OTU_ID, names_to = "Stations", values_to = "ASV_abundance")
  
  merge.asv <- left_join(piv.asv, select(merge_data, c(qseqid, clusters)), 
                         by = join_by(OTU_ID == qseqid),
                         relationship = "many-to-many")
  
  final.asv <- merge.asv%>%
    group_by(Stations, clusters)%>%
    summarise(asv_sum = sum(ASV_abundance),
              .groups = 'drop')%>%
    group_by(Stations)%>%
    mutate(perc = asv_sum/sum(asv_sum))
  
  perc.data <- final.asv%>%
    filter(!is.na(clusters))%>%
    group_by(Stations)%>%
    summarise(total_perc = sum(perc))%>%
    ungroup()%>%
    summarise(min = min(total_perc),
              mean = mean(total_perc),
              median = median(total_perc),
              max = max(total_perc))
  
  return(list(final.asv, perc.data))
}

#This function takes some relative abundance data and metadata in and creates dataframes/plots analyzing the patterns in relative abundance
analyze.rel.abundances <- function(data, metadata, use.metadata = TRUE){
  #Merge the relative abundance data with the metadata for Longhurst Province info
  if (use.metadata == TRUE){
    NA.perc.data <- left_join(data, select(metadata, c(Cruise,Longhurst_Long)), 
                              by=join_by(Stations==Cruise))%>%
      group_by(Stations, clusters)%>%
      mutate(Longhurst_Long = str_split(Longhurst_Long,"\t")[[1]][2])%>%
      distinct(.keep_all=TRUE)
  } else {
    NA.perc.data <- data
  }
  
  #Filter out the "total abundance" group
  NA.perc.data <- NA.perc.data%>%
    filter(Stations != "total_abundance")
  
  #Filter out the NA clusters (unmapped abundance) and recompute relative percentages while also shortening the Longhurst names
  perc.data <- NA.perc.data%>%
    filter(!is.na(clusters))%>%
    mutate(clusters = as.factor(clusters))%>%
    group_by(Stations)%>%
    mutate(perc = asv_sum/sum(asv_sum))
  
  #Keep the NA fraction and shorten the Longhurst names
  NA.perc.data <- NA.perc.data%>%
    mutate(clusters = as.factor(clusters))
  
  #Pull out the NA fractions and merge into dataframe to be able to rearrange the stations
  NA_frac <- NA.perc.data%>%
    filter(is.na(clusters))%>%
    select(Stations, perc)%>%
    rename(NA_frac=perc) 
  
  NA.perc.data <- NA.perc.data%>%
    ungroup()%>%
    left_join(., NA_frac, by="Stations")
  
  if (use.metadata == TRUE){
    ordered.NA.perc.data <- NA.perc.data%>%
      arrange(Longhurst_Long, NA_frac)
  }
  
  #Plot the relative abundances when excluding the unmapped fraction (possibly add a colorbar of the unmapped fraction at some point)
  if (use.metadata == TRUE){
  perc.plot <- ggplot(perc.data,aes(x=Stations,y=perc,fill=clusters))+
    geom_point(aes(x=Stations,y=1.005,color=Longhurst_Long),size=2,alpha=0.5)+
    geom_bar(position="stack",stat="identity",width=1)+
    coord_cartesian(ylim=c(0,1.01),expand=FALSE)+
    scale_fill_manual(values=c("#cc5a43","#7764cb","#75ab3d","#c062bb","#54a77a","#c75980","#b68f40","#678ecd"))+
    theme_bw()+
    labs(x="",fill="Cluster")+
    theme(axis.text.x=element_text(angle=45,size=4.5,hjust=1,vjust=1.1),axis.ticks.x=element_blank(),panel.grid=element_blank())+
    guides(color=guide_legend(), fill=guide_legend())+
    custom_theme
  
  #Raw abundance plot including the NA fraction to see how much of the community is unmapped
  NA.perc.plot <- ggplot(NA.perc.data,aes(x=Stations,y=perc,fill=clusters))+
    geom_point(aes(x=Stations,y=1.005,color=Longhurst_Long),size=2,alpha=0.5)+
    geom_bar(position="stack",stat="identity",width=1)+
    coord_cartesian(ylim=c(0,1.01),expand=FALSE)+
    scale_fill_manual(values=c("#cc5a43","#7764cb","#75ab3d","#c062bb","#54a77a","#c75980","#b68f40","#678ecd"), na.value = "gray80")+
    theme_bw()+
    labs(x="",fill="Cluster")+
    theme(axis.text.x=element_text(angle=45,size=4.5,hjust=1,vjust=1.1),axis.ticks.x=element_blank(),panel.grid=element_blank())+
    guides(color=guide_legend(), fill=guide_legend())+
    custom_theme
  
  ordered.NA.perc.plot <- ggplot(ordered.NA.perc.data,aes(x=fct_inorder(Stations),y=perc,fill=clusters))+
    geom_point(aes(x=fct_inorder(Stations),y=1.005,color=Longhurst_Long),size=2,alpha=0.5)+
    geom_bar(position="stack",stat="identity",width=1)+
    coord_cartesian(ylim=c(0,1.01),expand=FALSE)+
    scale_fill_manual(values=c("#cc5a43","#7764cb","#75ab3d","#c062bb","#54a77a","#c75980","#b68f40","#678ecd"), na.value = "gray80")+
    theme_bw()+
    labs(x="",fill="Cluster")+
    theme(axis.text.x=element_text(angle=45,size=4.5,hjust=1,vjust=1.1),axis.ticks.x=element_blank(),panel.grid=element_blank())+
    guides(color=guide_legend(), fill=guide_legend())+
    custom_theme
  } else {
    perc.plot <- ggplot(perc.data,aes(x=Stations,y=perc,fill=clusters))+
      geom_bar(position="stack",stat="identity",width=1)+
      coord_cartesian(ylim=c(0,1.01),expand=FALSE)+
      scale_fill_manual(values=c("#cc5a43","#7764cb","#75ab3d","#c062bb","#54a77a","#c75980","#b68f40","#678ecd"), na.value = "gray80")+
      theme_bw()+
      labs(x="",fill="Cluster")+
      theme(axis.text.x=element_text(angle=45,size=4.5,hjust=1,vjust=1.1),axis.ticks.x=element_blank(),panel.grid=element_blank())+
      guides(color=guide_legend(), fill=guide_legend())+
      custom_theme
    
    NA.perc.plot <- ggplot(NA.perc.data,aes(x=Stations,y=perc,fill=clusters))+
      geom_bar(position="stack",stat="identity",width=1)+
      coord_cartesian(ylim=c(0,1.01),expand=FALSE)+
      scale_fill_manual(values=c("#cc5a43","#7764cb","#75ab3d","#c062bb","#54a77a","#c75980","#b68f40","#678ecd"), na.value = "gray80")+
      theme_bw()+
      labs(x="",fill="Cluster")+
      theme(axis.text.x=element_text(angle=45,size=4.5,hjust=1,vjust=1.1),axis.ticks.x=element_blank(),panel.grid=element_blank())+
      guides(color=guide_legend(), fill=guide_legend())+
      custom_theme
  }
  
  if (use.metadata == TRUE){
    return(list(perc.data, NA.perc.data, perc.plot, NA.perc.plot, ordered.NA.perc.plot))
  } else {
    return(list(perc.data, NA.perc.data, perc.plot, NA.perc.plot))
  }
}

