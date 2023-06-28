#Load in packages
library(SOMbrero)
library(gplots)
library(tidyverse)
library(patchwork)
library(purrr)
library(ragg)

#Set a random seed
set.seed(123)

#Load in the necessary data files and reformat as necessary
metabolite_dictionary<-read.csv("Data/all_max_flux_met_names_diamond.csv")
metabolite_info<-read.csv("Data/all_max_flux_diamond_transposed.csv")
file_headers<-read.csv("Data/all_max_flux_headers_diamond.csv")

#Prepare matrix for input into SOMs
metabolite_matrix<-as.matrix(metabolite_info[,-1])
dimnames(metabolite_matrix)<-list(as.matrix(file_headers[,2]),as.matrix(metabolite_dictionary[,2]))

#Read in information of the reaction frequency for the genomes to subset a "high quality" set of genomes
rxn_freq_info<-read.csv("Data/rxn_info.csv")
thresh<-0.5

high.quality.genomes<-rxn_freq_info%>%filter(mean_freq>thresh)%>%select(X)%>%as.matrix()
subset_rows<-which(rownames(metabolite_matrix)%in%high.quality.genomes)
metabolite_matrix<-metabolite_matrix[subset_rows,]

#Need to remove any columns with no flux in any models
zerofluxes<-which(colSums(metabolite_matrix)==0)
if (is_empty(zerofluxes)==FALSE){
  metabolite_matrix<-metabolite_matrix[,-zerofluxes]
}

#Need to scale the data before input into SOM algorithm
scaled_matrix<-scale(metabolite_matrix)

#We would like to remove certain metabolites we don't think are super interesting and have been shown to drive SOM patterns
remove_metab<-which(colnames(scaled_matrix) %in% c("Manganese","Zinc","Co2+","Calcium","Chloride","Potassium","Copper","Magnesium"))
scaled_matrix<-scaled_matrix[,-remove_metab]

#Run SOMs, for the full scaled dataset we have selected a 9x8 grid but other grids are being tested
gridrows=10; gridcols=10; gridsize=gridrows*gridcols
initSOM(dimension=c(gridrows,gridcols),maxit=2000,nb.save=10,scaling="none")
carve.som<-trainSOM(x.data=scaled_matrix)

#Analyze SOMs
#Use superClass with a chosen k value (manually selected) and then index it to get the clusters for each genome
num_clusters=6
carve.sc<-superClass(carve.som,k=num_clusters)

clusters<-carve.sc$cluster
ids<-carve.sc$som$clustering
sample_clusters<-clusters[ids]

#We also want to examine the recipes to find genomes where all or nearly all recipes are identical
genome_IDs<-unique(rownames(scaled_matrix))
dupe_count<-0
dupe_IDs<-c()
for (i in 1:length(genome_IDs)){
  curr_id<-genome_IDs[i]
  curr_indices<-which(rownames(scaled_matrix) %in% curr_id)
  sub_mat<-scaled_matrix[curr_indices,]
  dupe_check<-unique(sub_mat)
  if (dim(dupe_check)[1]==1){
    print(paste("Dupe found in genome",curr_id))
    dupe_IDs<-cbind(dupe_IDs,curr_id)
    dupe_count<-dupe_count+1
  }
}

#Try running a PCA on the scaled data to then overlay with the SOM clustering (not super informative)
data_pca<-prcomp(scaled_matrix)
points<-data_pca$x%>%as_tibble()%>%mutate(.,clusters=sample_clusters)
points$clusters<-as_factor(points$clusters)

ggplot(points,aes(x=PC1,y=PC2,colour=clusters))+geom_point()+theme_bw()

#Now try reformulating the dataframe by metabolite and look at summary stats
scaled_df<-as_tibble(scaled_matrix)%>%mutate(.,clusters=sample_clusters,genomes=rownames(scaled_matrix))%>%
  pivot_longer(-c(clusters,genomes),names_to="metabolite",values_to="scaled_value")
scaled_df$clusters<-as_factor(scaled_df$clusters)

summarize_df<-scaled_df%>%group_by(clusters,metabolite)%>%summarize(mean=mean(scaled_value,na.rm=TRUE),std=sd(scaled_value,na.rm=TRUE))

ggplot(summarize_df,aes(x=metabolite,y=mean,colour=clusters))+geom_point()+theme_bw()

#Repeating the above process but thresholding out metabolites
threshold_df<-summarize_df%>%filter(abs(mean)>0.1)
ggplot(threshold_df,aes(x=metabolite,y=mean,colour=clusters))+geom_point()+theme_bw()+theme(axis.text.x=element_blank())+ylab("Scaled Mean Flux")
ggplot(threshold_df,aes(x=mean,y=std,colour=clusters))+geom_point()+theme_bw()+xlab("Scaled Mean Flux")+ylab("Scaled Standard Deviation")

#Try the same procedure but for the unscaled dataset
raw_df<-as_tibble(metabolite_matrix)%>%mutate(.,clusters=sample_clusters,genomes=rownames(scaled_matrix))%>%
  pivot_longer(-c(clusters,genomes),names_to="metabolite",values_to="scaled_value")
raw_df$clusters<-as_factor(raw_df$clusters)

summarize_df<-raw_df%>%group_by(clusters,metabolite)%>%summarize(mean=mean(scaled_value,na.rm=TRUE),std=sd(scaled_value,na.rm=TRUE))

threshold_df<-summarize_df%>%filter(abs(mean)>10)

ggplot(threshold_df,aes(x=metabolite,y=mean,colour=clusters))+geom_point()+theme_bw()
ggplot(threshold_df,aes(x=mean,y=std,colour=clusters))+geom_point()+theme_bw()

#---------------------------------------------------------------------------------------------

#We want to compute several metrics on the data itself to determine need for thresholding and excluding certain metabolites

#First we want to look at the average coefficients of variance per metabolite per genome (one value per genome)
#Use the genome IDs to subset all models for a genome
genome_IDs<-unique(rownames(metabolite_matrix))
average_coeffsvar<-vector(mode="integer",length=length(genome_IDs))
for (i in 1:length(genome_IDs)){
  curr_rows<-which(rownames(metabolite_matrix)==genome_IDs[i])
  curr_matrix<-metabolite_matrix[curr_rows,]
  #Need to remove metabolites that aren't present in any recipes
  zerocols<-which(colSums(curr_matrix)==0)
  curr_matrix<-curr_matrix[,-zerocols]
  
  #Now compute the coefficient of variance (mean/sd), need to also exclude zeros here in each column
  curr_means<-apply(curr_matrix,2,function(v) v[which(v>0)]%>%mean())
  curr_sds<-apply(curr_matrix,2,function(v) v[which(v>0)]%>%sd())
  #Any NA values reflect metabolites with only one non-zero flux so we will modify those to 0
  curr_sds[which(is.na(curr_sds)==TRUE)]<-0
  
  coeffs_var<-curr_sds/curr_means
  average_coeffsvar[i]<-mean(coeffs_var)
}

#Convert to a dataframe and plot as a histogram
coeffs_df<-as_tibble(average_coeffsvar)
ggplot(coeffs_df,aes(x=value))+geom_histogram(color="blue",fill="white")+labs(x="Average Coefficient of Variance",y="Number of Genomes")+theme_bw()

#We want to look at the evenness of model distribution across the clusters to look for genomes to remove
ids<-carve.sc$som$clustering
sample_clusters<-clusters[ids]
counts_table<-vector(mode="integer",length=60)
for (i in 1:length(genome_IDs)){
  curr_rows<-which(rownames(scaled_matrix)==genome_IDs[i])
  curr_distr<-sample_clusters[curr_rows]%>%table()%>%unique()
  counts_table[curr_distr]<-counts_table[curr_distr]+1
}
evenness_df<-as_tibble(counts_table)
ggplot(evenness_df,aes(x=c(1:60),y=value))+
  geom_point()+
  labs(x="Number of replicates in cluster",y="Number of Genomes")+
  theme_bw()

#We also want to identify for different levels of K how much genome replicates are being split into multiple nodes
sizes<-c(2:20)
# genome_IDs<-unique(file_headers$model_name)
consensus<-matrix(data=0,nrow=length(sizes),ncol=1)
breadth<-matrix(data=0,nrow=length(sizes),ncol=1)
evenness<-matrix(data=0,nrow=length(sizes),ncol=1)
for (i in 1:length(sizes)){
  curr_size=sizes[i]
  curr.sc<-superClass(carve.som,k=curr_size)
  clusters<-curr.sc$cluster
  ids<-curr.sc$som$clustering
  sample_clusters<-clusters[ids]
  
  tmp.cons<-matrix(data=0,nrow=length(genome_IDs),ncol=1)
  tmp.breadth<-matrix(data=0,nrow=length(genome_IDs),ncol=1)
  tmp.evenness<-matrix(data=0,nrow=length(genome_IDs),ncol=1)
  for (j in 1:length(genome_IDs)){
    curr_genome<-genome_IDs[j]
    #curr.cons<-(which(file_headers$model_name %in% curr_genome)%>%sample_clusters[.]%>%table()%>%max())/60
    #curr.breadth<-(which(file_headers$model_name %in% curr_genome)%>%sample_clusters[.]%>%unique()%>%length())/60
    curr_cluster<-which(rownames(scaled_matrix) %in% curr_genome)%>%sample_clusters[.]
    curr.cons<-(curr_cluster%>%table()%>%max())/60
    curr.breadth<-(curr_cluster%>%unique()%>%length())/60
    if (length(table(curr_cluster))==1){
      curr.evenness<-1
    } else {
      curr.evenness<-(max(table(curr_cluster))-min(table(curr_cluster)))/60
    }
    
    tmp.cons[j]<-curr.cons
    tmp.breadth[j]<-curr.breadth
    tmp.evenness[j]<-curr.evenness
  }
  consensus[i]<-mean(tmp.cons)
  breadth[i]<-mean(tmp.breadth)
  evenness[i]<-mean(tmp.evenness)
}
consensus_df<-cbind(sizes,consensus,breadth,evenness)%>%as.data.frame()
colnames(consensus_df)<-c("Cluster Size","Consensus","Breadth","Evenness")
consensus_df<-consensus_df%>%pivot_longer(-c(`Cluster Size`),names_to="Metric",values_to="Value")

ggplot(consensus_df,aes(x=`Cluster Size`,y=Value,color=Metric))+geom_line()+ylim(0,1)+theme_bw()

#New code c/o Chase that can generate some new types of plots similar to those in the kohonen package
df <- carve.som$prototypes %>% as.data.frame()

df$nodes <- rownames(df) 
df$cluster <- carve.sc$cluster

df$row_names <- rep(1:gridcols,times = gridrows)
df$column_names <- rep(1:gridrows,each = gridcols)

ggplot(df[,c(tmp,554:557)],
       aes(x = row_names, y = column_names, fill = `L-Phenylalanine`, color = as.factor(cluster))) +
  geom_point(size = 8, pch = 21, stroke = 2) +
  scale_fill_gradient(low = "white", high = "black") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom", legend.key = element_blank(),
        legend.box = "vertical") +
  labs(color = "Cluster") +
  coord_equal()


test_df <- df[,c(1:10,554:557)]

test_df <- test_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                    names_to = "metabolite", values_to = "scaled_flux") %>%
  group_by(nodes) %>% mutate(min = min(scaled_flux),
                             max = max(scaled_flux),
                             polar_size = ((scaled_flux-min)/(max-min)))

test_df$metabolite <- factor(test_df$metabolite, levels = unique(test_df$metabolite))

plot_df = test_df %>%
  nest(-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = metabolite, y = polar_size, fill = metabolite) +
                       geom_bar(stat = "identity", show.legend = FALSE, width = 1) +
                       coord_polar() +
                       scale_fill_manual(values = c("#c35ea5",
                                                    "#55b14b",
                                                    "#9a5ed1",
                                                    "#8b9c41",
                                                    "#7474c4",
                                                    "#c18c41",
                                                    "#619fd8",
                                                    "#cb5537",
                                                    "#4aab8b",
                                                    "#ca586f"))))

plot_df$row_names <- rep(1:gridcols,times = gridrows)
plot_df$column_names <- rep(1:gridrows,each = gridcols)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.5, xmax = ..1+0.5,
                                               ymin = ..2-0.5, ymax = ..2+0.5))) %>%
  pull(annotation)


out_plot <- ggplot() + 
  geom_point(data = df[,554:557], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 25,stroke=2) +
  geom_col(data = test_df,
           aes(0,0, fill = metabolite)) +
  coord_cartesian(expand = 1.4) +
  plot_annotations + coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = ggplot2::margin(5,5,5,5)) +
  scale_fill_manual(values = c("#c35ea5",
                               "#55b14b",
                               "#9a5ed1",
                               "#8b9c41",
                               "#7474c4",
                               "#c18c41",
                               "#619fd8",
                               "#cb5537",
                               "#4aab8b",
                               "#ca586f")) + labs(color="Cluster",fill = "Metabolite") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

# agg_png(filename = "./som_wedges.png", width = 8, height = 12, units = "in", res = 400)
# plot(out_plot)
# dev.off()

#Code to examine the phylogeny of SOM clusters
#Load phylogeny and match it to genome IDs
phylo<-read.csv("full-data_order-classifications.tsv",sep="\t",header=FALSE)
phylo_loc<-match(genome_IDs,phylo[,1])%>%.[-which(is.na(.))]
phylo_IDs<-phylo[phylo_loc,1]
phylo_classifications<-phylo[phylo_loc,2]

#Regenerate a superClass object of the correct value K
carve.sc<-superClass(carve.som,k=num_clusters)
clusters<-carve.sc$cluster
ids<-carve.sc$som$clustering
sample_clusters<-clusters[ids]

#Get the bulk estimates of what gridpoints and clusters the phylogeny place into
indices<-which(file_headers[,2] %in% phylo_IDs)
phylo_clusters<-sample_clusters[indices]
phylo_grid_ids<-ids[indices]

#Pull the dataframe information pertinent to ordering the nodes and clusters
info_df<-df[,521:524]

#Now we want to construct a data frame where the phylogeny is broken out into columns for "wedgie" plots that we can combine with info df
#We will do this by creating a matrix of grid points by cluster and fill it out with relative amounts of each taxonomic id
phylo_mat<-matrix(data=0,nrow=gridsize,ncol=15)
for (i in 1:length(unique(ids))){
  names_dict<-unique(phylo[,2])
  curr_IDs<-which(phylo_grid_ids==i)%>%names()%>%unique()
  phylo_values<-which(phylo[,1] %in% curr_IDs)%>%phylo[.,2]%>%table()%>%as.matrix()%>%t()
  name_order<-match(colnames(phylo_values),names_dict)
  phylo_mat[i,name_order]<-phylo_values
}
colnames(phylo_mat)<-names_dict
phylo_df<-as.data.frame(phylo_mat)%>%mutate(.,nodes=c(1:gridsize))%>%merge(.,info_df,by="nodes")

#Now we want to plot the relative abundances of phylo for each point
phylo_df <- phylo_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                    names_to = "taxonomic id", values_to = "abundance") %>%
  group_by(nodes) %>% mutate(min = min(abundance),
                             max = max(abundance),
                             polar_size = ((abundance-min)/(max-min)))

phylo_df$`taxonomic id` <- factor(phylo_df$`taxonomic id`, levels = unique(phylo_df$`taxonomic id`))

plot_df = phylo_df %>%
  nest(-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = `taxonomic id`, y = polar_size, fill = `taxonomic id`) +
                       geom_bar(stat = "identity", show.legend = FALSE, width = 1) +
                       coord_polar() +
                       scale_fill_manual(values = c("#c56863",
                                                    "#64b850",
                                                    "#c353bb",
                                                    "#b7b441",
                                                    "#7067d7",
                                                    "#de923e",
                                                    "#598dcc",
                                                    "#cc4e32",
                                                    "#53c1a9",
                                                    "#d14074",
                                                    "#3e8851",
                                                    "#9871bc",
                                                    "#758539",
                                                    "#c16e9c",
                                                    "#a4793f"))))

plot_df$row_names <- rep(1:gridcols,times = gridrows)
plot_df$column_names <- rep(1:gridrows,each = gridcols)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.5, xmax = ..1+0.5,
                                               ymin = ..2-0.5, ymax = ..2+0.5))) %>%
  pull(annotation)

out_plot <- ggplot() + 
  geom_point(data = df[,521:524], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 25,stroke=2) +
  geom_col(data = phylo_df,
           aes(0,0, fill = `taxonomic id`)) +
  coord_cartesian(expand = 1.4) +
  plot_annotations + coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.margin = ggplot2::margin(5,5,5,5)) +
  scale_fill_manual(values = c("#c56863",
                               "#64b850",
                               "#c353bb",
                               "#b7b441",
                               "#7067d7",
                               "#de923e",
                               "#598dcc",
                               "#cc4e32",
                               "#53c1a9",
                               "#d14074",
                               "#3e8851",
                               "#9871bc",
                               "#758539",
                               "#c16e9c",
                               "#a4793f")) + labs(color="Cluster",fill = "Taxonomic ID") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))
