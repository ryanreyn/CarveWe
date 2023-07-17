#Load in packages
library(SOMbrero)
library(gplots)
library(tidyverse)
library(patchwork)
library(purrr)
library(ragg)

##### Read in Data ####

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
metabolite_matrix<-metabolite_matrix[,-remove_metab]; scaled_matrix<-scaled_matrix[,-remove_metab]

#Run SOMs, for the full scaled dataset we have selected a 9x8 grid but other grids are being tested
gridrows=10; gridcols=10; gridsize=gridrows*gridcols
initSOM(dimension=c(gridrows,gridcols),maxit=2000,nb.save=10,scaling="none")
set.seed(123)
carve.som<-trainSOM(x.data=scaled_matrix)

#Analyze SOMs
#Use superClass with a chosen k value (manually selected) and then index it to get the clusters for each genome
num_clusters=7
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

#Determine the majority cluster to which each genome belongs
majority_clust<-c()
set.seed(123)
for (i in 1:length(genome_IDs)){
  curr_id<-genome_IDs[i]
  curr_indices<-which(rownames(scaled_matrix) %in% curr_id)
  curr_clust<-sample_clusters[curr_indices]%>%table()
  
  major<-which(curr_clust %in% max(curr_clust)) %>% names(curr_clust)[.] %>% as.numeric()
  #There are a handful of models (14 in the current construction, 10x10 grid for the >50% genomes) that have two SOM clusters that tie for majority so we are randomly selecting in those cases
  if (length(major)==1){
    majority_clust<-c(majority_clust,major)
  } else {
    #Functionality to visualize the genomes with ties
    print(paste("This genome has",length(major),"SOM clusters with",max(curr_clust),"models each."))
    #print(curr_clust)
    
    majority_clust<-c(majority_clust,sample(major,1))
  }
}
cluster_df<-as_tibble(genome_IDs)%>%mutate(.,clusters=majority_clust)

#We want to repeat this process for the majority SOM grid point which is more specific than SOM cluster
majority_gridpoint<-c()
set.seed(123)
for (i in 1:length(genome_IDs)){
  curr_id<-genome_IDs[i]
  curr_indices<-which(rownames(scaled_matrix) %in% curr_id)
  curr_point<-ids[curr_indices]%>%table()
  
  major<-which(curr_point %in% max(curr_point)) %>% names(curr_point)[.] %>% as.numeric()
  #There are a handful of models (14 in the current construction, 10x10 grid for the >50% genomes) that have two SOM clusters that tie for majority so we are randomly selecting in those cases
  if (length(major)==1){
    majority_gridpoint<-c(majority_gridpoint,major)
  } else {
    #Functionality to visualize the genomes with ties
    print(paste("This genome has",length(major),"SOM grid points with",max(curr_point),"models each."))
    #print(curr_clust)
    
    majority_gridpoint<-c(majority_gridpoint,sample(major,1))
  }
}
gridpoint_df<-as_tibble(genome_IDs)%>%mutate(.,gridpoints=majority_gridpoint)
gridpoint_table<-table(gridpoint_df$gridpoints)

#### Growth Data ####

#Plot a grid map with genomes as density, growth rates as a heatmap and colored borders by cluster
growth_dat<-read.csv("Data/growth_violin_data.csv")

df <- carve.som$prototypes %>% as.data.frame()
classifiers_df<-data.frame(nodes=rownames(df))%>%mutate(.,growth=0,`Number Genomes`=0)
classifiers_df$`Number Genomes`[names(gridpoint_table)%>%as.numeric()]<-gridpoint_table
for (i in 1:length(unique(gridpoint_df$gridpoints))){
  curr_point<-unique(gridpoint_df$gridpoints)[i]
  point_indices<-which(gridpoint_df$gridpoints %in% curr_point)
  point_matches<-match(gridpoint_df$value[point_indices],growth_dat$X,nomatch=0)
  point_matches<-which(point_matches>0)%>%point_matches[.]
  if (length(point_matches)>1){
    classifiers_df$growth[curr_point]<-mean(growth_dat$dCUB[point_matches])
  } else {
    classifiers_df$growth[curr_point]<-growth_dat$dCUB[point_matches]
  }
}
#Change any zero entries to NA
zero_loc<-which(classifiers_df$growth==0)
classifiers_df$growth[zero_loc]<-NA

classifiers_df$cluster <- carve.sc$cluster

classifiers_df$row_names <- rep(1:gridcols,times = gridrows)
classifiers_df$column_names <- rep(1:gridrows,each = gridcols)

classifiers_df<-classifiers_df%>%select(.,c(nodes,cluster,row_names,column_names),everything())

#Create a dataframe object that will draw the boxes around the cluster nodes
cluster_squares_df<-data.frame(cluster=c(1:num_clusters))
cluster_squares_df<-mutate(cluster_squares_df,xmin=0,xmax=0,ymin=0,ymax=0)
for (i in 1:num_clusters){
  curr_clust<-i
  curr_points<-filter(classifiers_df,cluster==i)%>%select(row_names,column_names)
  xmin=min(curr_points$column_names -0.4); xmax=max(curr_points$column_names +0.4); ymin=min(curr_points$row_names -0.4); ymax=max(curr_points$row_names +0.4)
  cluster_squares_df[i,c("xmin","xmax","ymin","ymax")]<-c(xmin,xmax,ymin,ymax)
}

ggplot(classifiers_df,
       aes(x = column_names, y = row_names)) +
  geom_point(aes(size = `Number Genomes`,fill=growth),pch=21,stroke=1) +
  scale_fill_gradient2(low = "red", high = "blue",mid="white",midpoint=-0.08,na.value="grey50",breaks=c(-0.3,-0.2,-0.1),limits=c(-0.3,-0.1)) +
  scale_size_continuous(range=c(4,16),breaks = c(0,150,300))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom", legend.key = element_blank(),
        legend.box = "vertical")+
  guides(color=guide_legend(override.aes=list(size=5)))+
  labs(size="Number of Genomes",fill="dCUB value") +
  coord_equal()

asinh_trans = function(){
  scales::trans_new("asinh",
                    transform = asinh,
                    inverse = sinh)
  
}

ggplot() +
  geom_rect(data=cluster_squares_df,mapping=aes(xmin=xmin, xmax = xmax, ymin = ymin, ymax= ymax,color=as.factor(cluster)),fill= "white", lwd = 2) +
  geom_point(data= classifiers_df,aes(x = column_names, y = row_names, size = `Number Genomes`,fill=growth),pch=21,stroke=1) +
  scale_fill_gradientn(colors = c("red","white","blue"), values = scales::rescale(c(-0.25,-0.09833212,-0.08)),
                       guide = "colorbar", limits=c(-0.25,-0.08), 
                       breaks = c(-0.25,-0.15,-0.08), labels = c("<-0.25","-0.15",">-0.08"),
                       oob = scales::squish, trans = "asinh") +
  scale_size_continuous(range=c(4,14),breaks = c(0,150,300))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom", legend.key = element_blank(),
        legend.box = "vertical")+
  guides(color=guide_legend(override.aes=list(size=5)))+
  labs(size="Number of Genomes",fill="dCUB value",color="Cluster")+
  coord_equal()

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

#We will now use manual classifications of key metabolites as identified within the thresholded group to coarse-grain our resolution of this issue for better data visualization
#Read in classifications
class_file<-read.csv("Data/classified_metabolites.csv")
class_df<-select(class_file,c(Metabolite,Higher.Level.Classification))

matching_classifiers<-which(class_df$Metabolite %in% unique(threshold_df$metabolite))
class_df<-class_df[matching_classifiers,]
classifications<-class_df$Higher.Level.Classification

#Now we want to match these to the metabolites in our threshold df and re-plot
metab_matches<-match(threshold_df$metabolite,class_df$Metabolite)
threshold_df$metab_class<-class_df$Higher.Level.Classification[metab_matches]

ggplot(threshold_df,aes(x=metab_class,y=mean,colour=clusters))+geom_violin(aes(x=metab_class,fill=clusters))+theme_bw()

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

#----------------------------------------------------------------------------------------------#
#### Metabolite Classes ####


#Version of plots designed to do metabolite classes
#Create the prototypes dataframe, then modify it down to metabolite classes by averaging
df <- carve.som$prototypes %>% as.data.frame()
classifiers_df<-data.frame(nodes=rownames(df))
for (i in 1:length(unique(class_df$Higher.Level.Classification))){
  curr_class<-unique(class_df$Higher.Level.Classification)[i]
  metab_indices<-which(class_df$Higher.Level.Classification %in% curr_class)
  metab_matches<-match(class_df$Metabolite[metab_indices],colnames(df))
  if (length(metab_matches)>1){
    classifiers_df<-mutate(classifiers_df,!!curr_class:=rowSums(df[,metab_matches]))
  } else {
    classifiers_df<-mutate(classifiers_df,!!curr_class:=sum(df[,metab_matches]))
  }
}

classifiers_df$cluster <- carve.sc$cluster

classifiers_df$row_names <- rep(1:gridcols,times = gridrows)
classifiers_df$column_names <- rep(1:gridrows,each = gridcols)

min_vals<-apply(classifiers_df[5:18],2,min)

classifiers_df<-classifiers_df%>%select(.,c(nodes,cluster,row_names,column_names),everything())

ggplot(classifiers_df,
       aes(x = column_names, y = row_names, fill = `Peptides`, color = as.factor(cluster))) +
  geom_tile(size = 0.25) +
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

#Create the wedge plots for every point in the SOM cluster with wedges for each metabolite class
pivot_df<-classifiers_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                    names_to = "metab_class", values_to = "scaled_flux") %>%
  group_by(nodes) %>% mutate(min = min(scaled_flux),
                             max = max(scaled_flux),
                             polar_size = ((scaled_flux)/(max-min)))

pivot_df$metabolite <- factor(pivot_df$metab_class, levels = unique(pivot_df$metab_class))

plot_df = pivot_df %>%
  nest(-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = metab_class, y = polar_size, fill = metab_class) +
                       geom_bar(stat = "identity", show.legend = FALSE, width = 1) +
                       coord_polar() +
                       scale_fill_manual(values = c("#7663cf",
                                                    "#8db43b",
                                                    "#c851b1",
                                                    "#55bc63",
                                                    "#d54767",
                                                    "#4fbeae",
                                                    "#d15236",
                                                    "#588dcc",
                                                    "#d2a23a",
                                                    "#a879bf",
                                                    "#458149",
                                                    "#bd6484",
                                                    "#89883f",
                                                    "#ba7648"))))

plot_df$row_names <- rep(1:gridcols,times = gridrows)
plot_df$column_names <- rep(1:gridrows,each = gridcols)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.5, xmax = ..1+0.5,
                                               ymin = ..2-0.5, ymax = ..2+0.5))) %>%
  pull(annotation)


out_plot <- ggplot() + 
  geom_point(data = classifiers_df[,1:4], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 20,stroke=2) +
  geom_col(data = pivot_df,
           aes(0,0, fill = metab_class)) +
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
  scale_fill_manual(values = c("#7663cf",
                               "#8db43b",
                               "#c851b1",
                               "#55bc63",
                               "#d54767",
                               "#4fbeae",
                               "#d15236",
                               "#588dcc",
                               "#d2a23a",
                               "#a879bf",
                               "#458149",
                               "#bd6484",
                               "#89883f",
                               "#ba7648")) + labs(color="Cluster",fill = "Metabolite Class") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

#We also want to try to make a version where we only display the average relative flux values for each 
#Create a mean df from the classifiers df
all_means<-c()
for (i in 1:num_clusters){
  curr_means<-filter(classifiers_df,cluster==i)%>%select(!c(nodes,cluster,row_names,column_names))%>%colMeans(.)
  all_means<-rbind(all_means,curr_means)
}
mean_df<-as.data.frame(all_means,row.names = FALSE)

mean_df$nodes<-c(1:num_clusters)
mean_df$row_names<-c(1,2,1,2,1,2,3)
mean_df$column_names<-c(1,1,2,2,3,3,3)
mean_df$cluster<-c(1:num_clusters)

#Optional block of code for if you want to threshold

mean_df<-mean_df%>%select(.,c(nodes,cluster,row_names,column_names),everything())

pivot_df<-mean_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                          names_to = "metab_class", values_to = "scaled_flux") %>%
  group_by(nodes) %>% mutate(min = min(scaled_flux),
                             max = max(scaled_flux),
                             polar_size = ((scaled_flux)/(max-min)))

pivot_df$metabolite <- factor(pivot_df$metab_class, levels = unique(pivot_df$metab_class))

plot_df = pivot_df %>%
  nest(data=-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = metab_class, y = polar_size, fill = metab_class) +
                       geom_bar(stat = "identity", show.legend = FALSE, width = 1) +
                       coord_polar() +
                       scale_fill_manual(values = c("#7663cf",
                                                    "#8db43b",
                                                    "#c851b1",
                                                    "#55bc63",
                                                    "#d54767",
                                                    "#4fbeae",
                                                    "#d15236",
                                                    "#588dcc",
                                                    "#d2a23a",
                                                    "#a879bf",
                                                    "#458149",
                                                    "#bd6484",
                                                    "#89883f",
                                                    "#ba7648"))))

plot_df$row_names <- c(1,2,1,2,1,2,3)
plot_df$column_names <- c(1,1,2,2,3,3,3)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.5, xmax = ..1+0.5,
                                               ymin = ..2-0.5, ymax = ..2+0.5))) %>%
  pull(annotation)


out_plot <- ggplot() + 
  geom_point(data = mean_df[,1:4], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 40,stroke=2) +
  geom_col(data = pivot_df,
           aes(0,0, fill = metab_class)) +
  coord_cartesian(expand = 10) +
  expand_limits(x=5,y=5) +
  plot_annotations + coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = ggplot2::margin(5,5,5,5)) +
  scale_fill_manual(values = c("#7663cf",
                               "#8db43b",
                               "#c851b1",
                               "#55bc63",
                               "#d54767",
                               "#4fbeae",
                               "#d15236",
                               "#588dcc",
                               "#d2a23a",
                               "#a879bf",
                               "#458149",
                               "#bd6484",
                               "#89883f",
                               "#ba7648")) + labs(color="Cluster",fill = "Metabolite Class") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5),color=guide_legend(title.position="top",override.aes=list(size=20)))

#-------------------------------------------------------------------------------------------------------#
#### Quartile Plots ####

#We now want to attempt a different way of plotting the data by metabolite class where we bin everything into quartiles and apply a set of discrete values to the wedges (either high, low-high, high-low, or low) based on quartiles
metabolite_df<-as.data.frame(metabolite_matrix)
classifiers_df<-data.frame(Model=rownames(metabolite_df),cluster=sample_clusters)
for (i in 1:length(unique(class_df$Higher.Level.Classification))){
  curr_class<-unique(class_df$Higher.Level.Classification)[i]
  metab_indices<-which(class_df$Higher.Level.Classification %in% curr_class)
  metab_matches<-match(class_df$Metabolite[metab_indices],colnames(metabolite_matrix))
  if (length(metab_matches)>1){
    classifiers_df<-mutate(classifiers_df,!!curr_class:=rowMeans(metabolite_df[,metab_matches]))
  } else {
    classifiers_df<-mutate(classifiers_df,!!curr_class:=metabolite_df[,metab_matches])
  }
}
#Add a numeric tag so all models from same ensemble have same tag
num_genomes<-dim(classifiers_df)[1]/60
genomes<-rep(c(1:num_genomes),each=60)
classifiers_df<-mutate(classifiers_df,Genome=genomes,nodes=ids)%>%select(Model,Genome,nodes,cluster,everything())
data_classifiers_df<-classifiers_df

#Now average each ensemble together by genome
# mean_df<-data.frame(Genome=unique(rownames(metabolite_matrix)),Model=unique(classifiers_df$Genome))
# mean_df[,colnames(classifiers_df)[c(-1:-3)]]<-NA
# for (i in 1:num_genomes){
#   curr_means<-filter(classifiers_df,Genome==i)%>%select(-c(Model,Genome,cluster))%>%colMeans(.)
#   mean_df[i,c(-1:-2)]<-curr_means
# }

#We now want to add either the majority cluster or grid point value to the dataframe
#Starting with the majority cluster
# mean_df<-mutate(mean_df,cluster=majority_clust)%>%select(Genome,Model,cluster,everything())

#Averaging values by SOM grid point for all models (no reduction to genomes)
mean_df<-data.frame(nodes=rownames(carve.som$prototypes))
mean_df[,colnames(classifiers_df)[c(-1:-4)]]<-NA
for (i in 1:gridsize){
  curr_means<-filter(classifiers_df,nodes==i)%>%select(-c(Genome,Model,cluster,nodes))%>%colMeans(.)
  mean_df[i,-1]<-curr_means
}

#Next we need to unscale the SOM prototypes to assign them to bins. NOTE: there is a slight different between the original data and the unscaled data because of rounding errors (probably)
col_means<-colMeans(metabolite_matrix)
col_sds<-apply(metabolite_matrix,2,sd)

df<-data.frame(carve.som$prototypes)
for (i in 1:dim(df)[1]){
  df[i,]<-(df[i,]*col_sds)+col_means
}

#Now we will compute quartiles for each metabolite class on the full data
# data_df<-select(data_classifiers_df,-c(Genome,Model,cluster,nodes))
# quartile_mat<-matrix(data=0,nrow=5,ncol=dim(data_df)[2])
# for (i in 1:dim(quartile_mat)[2]){
#   quartile_mat[,i]<-quantile(data_df[,i])
# }

#The unscaled prototypes then need to be averaged by metabolite class
classifiers_df<-data.frame(nodes=rownames(carve.som$prototypes))
for (i in 1:length(unique(class_df$Higher.Level.Classification))){
  curr_class<-unique(class_df$Higher.Level.Classification)[i]
  metab_indices<-which(class_df$Higher.Level.Classification %in% curr_class)
  metab_matches<-match(class_df$Metabolite[metab_indices],colnames(metabolite_matrix))
  if (length(metab_matches)>1){
    classifiers_df<-mutate(classifiers_df,!!curr_class:=rowMeans(df[,metab_matches]))
  } else {
    classifiers_df<-mutate(classifiers_df,!!curr_class:=(df[,metab_matches]))
  }
}

#Now we can compute the quartiles on the SOM prototypes and assign values to each SOM prototype for each metabolite class
quant_df<-select(classifiers_df,-nodes)
quartile_mat<-matrix(data=0,nrow=5,ncol=dim(quant_df)[2])
for (i in 1:dim(quartile_mat)[2]){
  quartile_mat[,i]<-quantile(quant_df[,i])
}

quartile_values<-c(0,0.33,0.66,1)
quartile_df<-data.frame(nodes=rownames(carve.som$prototypes))
quartile_df[,colnames(classifiers_df)[-1]]<-NA
for (i in 2:dim(classifiers_df)[2]){
  for (j in 1:dim(classifiers_df)[1]){
    quartile_check<-which(quartile_mat[,i-1]>=classifiers_df[j,i])%>%min()
    if (quartile_check<5){
      quartile_df[j,i]<-quartile_values[quartile_check]
    } else {
      quartile_df[j,i]<-quartile_values[quartile_check-1]
    }
  }
}

#Plot the quartile information
quartile_df$cluster <- carve.sc$cluster

quartile_df$row_names <- rep(1:gridcols,times = gridrows)
quartile_df$column_names <- rep(1:gridrows,each = gridcols)

quartile_df<-quartile_df%>%select(.,c(nodes,cluster,row_names,column_names),everything())

pivot_df<-quartile_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                          names_to = "metab_class", values_to = "flux quartiles") %>%
  group_by(nodes)

pivot_df$metabolite <- factor(pivot_df$metab_class, levels = unique(pivot_df$metab_class))

plot_df = pivot_df %>%
  nest(data=-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = metab_class, y = `flux quartiles`, fill = metab_class) +
                       geom_bar(stat = "identity", show.legend = FALSE, width = 1) +
                       coord_polar() +
                       scale_fill_manual(values = c("#7663cf",
                                                    "#8db43b",
                                                    "#c851b1",
                                                    "#55bc63",
                                                    "#d54767",
                                                    "#4fbeae",
                                                    "#d15236",
                                                    "#588dcc",
                                                    "#d2a23a",
                                                    "#a879bf",
                                                    "#458149",
                                                    "#bd6484",
                                                    "#89883f"))))

plot_df$row_names <- rep(1:gridcols,times = gridrows)
plot_df$column_names <- rep(1:gridrows,each = gridcols)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.45, xmax = ..1+0.45,
                                               ymin = ..2-0.45, ymax = ..2+0.45))) %>%
  pull(annotation)


out_plot <- ggplot() + 
  geom_point(data = quartile_df[,1:4], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 20,stroke=2) +
  geom_col(data = pivot_df,
           aes(0,0, fill = metab_class)) +
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
  scale_fill_manual(values = c("#7663cf",
                               "#8db43b",
                               "#c851b1",
                               "#55bc63",
                               "#d54767",
                               "#4fbeae",
                               "#d15236",
                               "#588dcc",
                               "#d2a23a",
                               "#a879bf",
                               "#458149",
                               "#bd6484",
                               "#89883f")) + labs(color="Cluster",fill = "Metabolite Class") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

#---------------------------

#Function for saving a plot nicely
# agg_png(filename = "Figures/prelim_som-wedges_by-class_averaged.png", width =250, height = 250, units = "mm", res = 400)
# plot(out_plot)
# dev.off()

#-------------------------------------------------------------------------------------------------#
#### Legacy Plots ####

#New code c/o Chase that can generate some new types of plots similar to those in the kohonen package -> this version is for individual metabolites
df <- carve.som$prototypes %>% as.data.frame()

df$nodes <- rownames(df) 
df$cluster <- carve.sc$cluster

df$row_names <- rep(1:gridcols,times = gridrows)
df$column_names <- rep(1:gridrows,each = gridcols)

ggplot(df[,c(tmp,541:544)],
       aes(x = row_names, y = column_names, fill = `N-Acetylneuraminate`, color = as.factor(cluster))) +
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


test_df <- df[,c(1:10,541:544)]

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
                       scale_fill_manual(values = c("#7663cf",
                                                    "#8db43b",
                                                    "#c851b1",
                                                    "#55bc63",
                                                    "#d54767",
                                                    "#4fbeae",
                                                    "#d15236",
                                                    "#588dcc",
                                                    "#d2a23a",
                                                    "#a879bf",
                                                    "#458149",
                                                    "#bd6484",
                                                    "#89883f",
                                                    "#ba7648"))))

plot_df$row_names <- rep(1:gridcols,times = gridrows)
plot_df$column_names <- rep(1:gridrows,each = gridcols)

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.5, xmax = ..1+0.5,
                                               ymin = ..2-0.5, ymax = ..2+0.5))) %>%
  pull(annotation)


out_plot <- ggplot() + 
  geom_point(data = df[,541:544], aes(x = column_names, y = row_names,color=as.factor(cluster)),
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
  scale_fill_manual(values = c("#7663cf",
                               "#8db43b",
                               "#c851b1",
                               "#55bc63",
                               "#d54767",
                               "#4fbeae",
                               "#d15236",
                               "#588dcc",
                               "#d2a23a",
                               "#a879bf",
                               "#458149",
                               "#bd6484",
                               "#89883f",
                               "#ba7648")) + labs(color="Cluster",fill = "Metabolite") + 
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5))

# agg_png(filename = "./som_wedges.png", width = 8, height = 12, units = "in", res = 400)
# plot(out_plot)
# dev.off()

#-------------------------------------------------------------------------------------------------#
#### Phylogeny Plots ####

#Code to examine the phylogeny of SOM clusters
#Load phylogeny and match it to genome IDs
phylo<-read.csv("Data/derep_all-GTDB-phylo.tsv",sep="\t",header=TRUE)

split_phylo<-separate(phylo,GTDB.Taxonomy,sep=";",into=c("domain","phylum","class","order","family","genus","species"))
groups<-split_phylo%>%group_by(order)%>%summarise(count=n())%>%arrange(desc(count))

groups$type<-"Other"
groups$type[1:15]<-groups$order[1:15]
groups$type<-gsub("o__","",groups$type)
split_phylo$Genome.prokka<-gsub("-prokka","",split_phylo$Genome.prokka)
split_phylo$group<-groups$type[match(split_phylo$order,groups$order)]

#Isolate the correct subset of genomes and reduce the dataframe to match
split_phylo_hq<-split_phylo[split_phylo$Genome.prokka %in% genome_IDs,]%>%select(Genome.prokka,group)

#Regenerate a superClass object of the correct value K
carve.sc<-superClass(carve.som,k=num_clusters)
clusters<-carve.sc$cluster
ids<-carve.sc$som$clustering
sample_clusters<-clusters[ids]

#Use the majority clustering previously calculated to assign a cluster to each genome
split_phylo_hq<-mutate(split_phylo_hq,cluster=majority_clust)
phylo_mat<-matrix(data=0,nrow=num_clusters,ncol=16)
colnames(phylo_mat)<-unique(split_phylo_hq$group)
for (i in 1:num_clusters){
  curr_cluster<-split_phylo_hq%>%select(cluster,group)%>%filter(cluster==i)
  curr_phylo<-table(curr_cluster$group)
  phylo_cols<-match(names(curr_phylo),colnames(phylo_mat))
  phylo_mat[i,phylo_cols]<-curr_phylo
}
phylo_df<-as.data.frame(phylo_mat,row.names = FALSE)%>%mutate(.,cluster=c(1:7))%>%select(cluster,everything())

# ggplot(filter(phylo_df,cluster==1),aes(x="",y=value,fill=group))+
#   geom_col(color="black")+
#   geom_text(aes(x=1.35,label=value),position=position_stack(vjust=0.5))+
#   coord_polar(theta="y")+
#   theme_void()

#Now we want to plot all the pie charts together
phylo_df$nodes<-c(1:num_clusters)
phylo_df$row_names<-c(1,2,1,2,1,2,3)
phylo_df$column_names<-c(1,1,2,2,3,3,3)
phylo_df$cluster<-c(1:num_clusters)

#Optional block of code for if you want to threshold

phylo_df<-phylo_df%>%select(.,c(nodes,cluster,row_names,column_names),everything())

pivot_df<-phylo_df %>% pivot_longer(-c(nodes, cluster, row_names, column_names),
                                   names_to = "Group", values_to = "Number of Genomes") %>%
  group_by(nodes)

pivot_df$Order <- factor(pivot_df$Group, levels = unique(pivot_df$Group))

plot_df = pivot_df %>%
  nest(data=-nodes) %>%
  mutate(plot = map2(data, nodes,
                     ~ ggplot(.x) + theme_void() +
                       aes(x = "", y = `Number of Genomes`, fill = Order)+
                       geom_col(color="black",show.legend=FALSE)+
                         coord_polar(theta="y")+
                       scale_fill_manual(values=c("#b8617c",
                                                  "#63b750",
                                                  "#895bc9",
                                                  "#b2b53b",
                                                  "#c84ca3",
                                                  "#3d854f",
                                                  "#d74164",
                                                  "#54bf9f",
                                                  "#cf5230",
                                                  "#54acd8",
                                                  "#d48e36",
                                                  "#5f7ac7",
                                                  "#6e772c",
                                                  "#bf83c9",
                                                  "#bca262",
                                                  "#b7694b"))))

pivot_df <- pivot_df %>% group_by(cluster) %>% mutate(perc = `Number of Genomes`/sum(`Number of Genomes`))

order_df <- pivot_df %>% group_by(Order) %>% summarise(mean = mean(perc)) %>%
  arrange(mean)

pivot_df$Order <- factor(pivot_df$Order, levels = order_df$Order)

ggplot(pivot_df, aes(x = as.factor(cluster), y = perc, fill = Order)) +
  geom_bar(stat = "identity") + coord_cartesian(expand = FALSE) + labs(x = "Cluster", y = "Proportion")+
  scale_fill_manual(values = c("#b8617c",
                               "#63b750",
                               "#895bc9",
                               "#b2b53b",
                               "#c84ca3",
                               "#3d854f",
                               "#d74164",
                               "#54bf9f",
                               "#cf5230",
                               "#54acd8",
                               "#d48e36",
                               "#5f7ac7",
                               "#6e772c",
                               "#bf83c9",
                               "#bca262",
                               "#b7694b")) 


plot_df$row_names <- (c(1,2,1,2,1,2,3)-1)*2+0.55
plot_df$column_names <- (c(1,1,2,2,3,3,3)-1)*2+0.25
plot_df$cluster<-phylo_df$cluster
plot_df<-select(plot_df,nodes,cluster,row_names,column_names,everything())

plot_annotations = plot_df %>% 
  mutate(annotation = pmap(list(column_names, row_names, plot),
                           ~ annotation_custom(ggplotGrob(..3),
                                               xmin = ..1-0.9, xmax = ..1+0.9,
                                               ymin = ..2-0.9, ymax = ..2+0.9))) %>%
  pull(annotation)


out_plot <- ggplot() +
  geom_col(data = pivot_df,
           aes(0,0, fill = Order))+
  geom_point(data = plot_df[,1:4], aes(x = column_names, y = row_names,color=as.factor(cluster)),
             pch = 21, size = 82.5,stroke=2)+
  coord_cartesian(expand = 10) +
  expand_limits(x=5,y=5) +
  plot_annotations + coord_equal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.margin = ggplot2::margin(5,5,5,5))+ 
        labs(color="Cluster")+
  guides(fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5),color=guide_legend(title.position="top",override.aes=list(size=20)))+
  scale_fill_manual(values = c("#b8617c",
                               "#63b750",
                               "#895bc9",
                               "#b2b53b",
                               "#c84ca3",
                               "#3d854f",
                               "#d74164",
                               "#54bf9f",
                               "#cf5230",
                               "#54acd8",
                               "#d48e36",
                               "#5f7ac7",
                               "#6e772c",
                               "#bf83c9",
                               "#bca262",
                               "#b7694b")) 


#Get the bulk estimates of what gridpoints and clusters the phylogeny place into
# indices<-which(file_headers[,2] %in% split_phylo_hq$Genome.prokka)
# phylo_clusters<-sample_clusters[indices]
# phylo_grid_ids<-ids[indices]

#Pull the dataframe information pertinent to ordering the nodes and clusters


# info_df<-df[,521:524]

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
