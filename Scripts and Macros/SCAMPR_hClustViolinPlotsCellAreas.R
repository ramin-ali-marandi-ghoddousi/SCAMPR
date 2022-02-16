####---------- SCAMPR_hClustViolinPlotsCellAreas ---------#####

##### Use this script to cluster cells into groups based on HiPlex RNAScope gene expression data and generate gene expression violin plots and cell area boxplots.
##### alimaran@usc.edu, 01-17-2021 

#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 60   
#*2) Define the gene columns that you want to remove on line 73 
#*3) Use correlational instead of euclidean distance matrix for clustering (see note in line 97)

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

# Clear environment
rm(list = ls()) 

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")

## Choose hierarchical clustering method: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
hc_method <- "ward.D"

## Choose and store the  number of clusters desired for dataset.  
# If unsure, use NbClust() package for determining optimal number of clusters for your data.
predicted_clusters <- 6

## Choose the 2 clusters for which violin plots will be generated
clusterA <- 1
clusterB <- 2


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages(c( "dplyr", "stringr", "dendextend", "colorspace", "reshape2", "ggplot2"))

library(dplyr) # For piping commands
library(stringr) # For manipulating column (gene) names
library(dendextend) # For assigning correct labels to each cluster
library(colorspace)  # For color schemes
library(reshape2) # For reshaping dataframe to long format for easy plotting in ggplot
library(ggplot2) # For plotting


#************* DATA WRANGLING AND NORMALIZATION ************#

## Check structure and dimensions of matrix
head(Raw.Data)
dim(Raw.Data)

## Convert Area Fraction gene expression to gene expression in pixels or um. To do this, multiply the Area Fraction columns by the ROI_Area column 
#and divide by 100. 
head(Raw.Data) # Determine which columns contain gene expression information.
geneColumns <- which(!(colnames(Raw.Data) %in% c("X", "Animal_ID", "Image_ID", "ROI_Area")))
geneIndices <- geneColumns # If desired, indices of the columns that represent your genes can be saved manually (example: geneIndices <- 5:16)
Raw.Data[,geneIndices] <- Raw.Data[,geneIndices]*(Raw.Data$ROI_Area)/100 # Calculate gene expression in pixels or um (will depend on your image scale setting in ImageJ).
head(Raw.Data)

## If desired, convert gene names to first letter capital (Mouse gene syntax)
newnames <- str_to_title(colnames(Raw.Data[geneIndices]))
colnames(Raw.Data)[geneIndices] <- newnames
head(Raw.Data)

## If desired, remove cells that have zero gene expression across all assayed genes.
nrow(Raw.Data)
zeroExpressors <- which(rowSums(Raw.Data[geneIndices])==0) # Get row indices for cells that don't express any genes.
Raw.Data <- Raw.Data[-zeroExpressors,]# Remove zero-expressing cells
nrow(Raw.Data)

# ## If desired, remove columns that represent genes that you do not want to analyze ("Cckar", "Piezo1", "Npy2r")
# genesToRemove <- c("Cckar", "Piezo1", "Npy2r")
# removeIndices <- which(colnames(Raw.Data) %in% genesToRemove)
# Raw.Data <- subset(Raw.Data, select = -c(removeIndices))
# head(Raw.Data)
# geneIndices <- min(geneIndices):(max(geneIndices)-length(genesToRemove)) # Define new gene indices

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,geneIndices] <- log1p(Normal.Data[,geneIndices]) 
range(Normal.Data[,geneIndices])

## Save ROI_Area column
ROI_Area <- Normal.Data$ROI_Area

## Store total number of cells as a variable
total_cells <- nrow(Normal.Data)  


#************* HEIRARCHICAL CLUSTERING ***************#

## Perform hierarchical clustering on gene expression data.  
# You can use get_dist(method = "pearson) or get_dist(method = "spearman) instead of dist() to get correlation distance measures.
hc <- Normal.Data[, geneIndices]    %>%   # Select columns with gene expression data
  dist()                        %>% # Compute euclidean distance/dissimilarity matrix or correlation matrix
  hclust(method = hc_method)   # Compute hierarchical clusters

## Identify which cells fall into which cluster and add this to data table. k = the number of clusters.
Cluster_ID <- cutree(hc, k = predicted_clusters)

## Append Cluster ID column to Normal.Data matrix for downstream plotting, and see how many cells fall into each cluster.
Normal.Data <- cbind(Normal.Data,Cluster_ID) 
table(Normal.Data$Cluster_ID) # Determine number of cells that fall into each cluster

## Generate dendrogram to determine which cluster maps to which label.  This is needed for color and cluster label consistency.
dend <- as.dendrogram(hc) # Generate dendrogram
dend.cutree <- cutree(dend, k = predicted_clusters, order_clusters_as_data = FALSE) # Cut dendrogam into k clusters.  Order data by the order of the labels in the dendrogram
idx <- order(as.numeric(names(dend.cutree))) # Grab the index of each cell, then sort indices based on cell number
dend.cutree <- dend.cutree[idx] # Reorder the dendrogram based on sorted indices (not labels)
tbl <- table(Cluster_ID, dend.cutree) # Match the number of cells per cluster to the number of cells per label
lbls <- apply(tbl,2,which.max) # Match cluster number to label number


#************* GENE EXPRESSION VIOLIN PLOTS FOR ALL CLUSTERS OR CLUSTERS OF CHOICE ***************#

## Subset full data matrix to create new data matrix for cells that fall into 2 clsters of choice (chosen in line 34-35)
SelectedClusters <- rbind(Normal.Data[Normal.Data$Cluster_ID == clusterA,] ,
                          Normal.Data[Normal.Data$Cluster_ID == clusterB,])

## Create color map for all clusters
cluster_color_map <- data.frame(Cluster_ID = as.numeric(lbls), 
                                Cluster_Color = c(rainbow_hcl(predicted_clusters)))

cluster_color_map <- cluster_color_map[order(cluster_color_map$Cluster_ID),2]

## Reshape Raw and Normalized datasets to long format, making it compatible with ggplot() for downstream plotting.
Genes <- colnames(Raw.Data[geneIndices])  # Get list column names representing each gene
Genes

SelectedClusters <- reshape2::melt(SelectedClusters,
                         id.vars = c("X", "Image_ID","ROI_Area", "Cluster_ID"),
                         measure.vars = Genes,
                         variable.name = "Gene",
                         value.name = "Expression")

Allclusters <- reshape2::melt(Normal.Data,
                    id.vars = c("X", "Image_ID","ROI_Area", "Cluster_ID"),
                    measure.vars = Genes,
                    variable.name = "Gene",
                    value.name = "Expression")

## Plot gene expression for clusters of choice
ggplot(SelectedClusters, aes(x = Gene, y=Expression, fill = as.factor(Cluster_ID)))+
  geom_violin(alpha = 0.8, scale="width")+
  geom_point(position = position_jitterdodge(dodge.width = .9), size = .2, alpha = .1)+
  theme_bw()+
  ggtitle(paste("Cluster", clusterA,"and",clusterB,"Gene Expression"))+
  xlab("")+
  ylab("LogNormal Gene Expression")+
  guides(fill=guide_legend(title="Cluster\nNumber"))+
  scale_fill_manual(values = c(cluster_color_map[as.numeric(levels(factor(SelectedClusters$Cluster_ID))[1])],
                               cluster_color_map[as.numeric(levels(factor(SelectedClusters$Cluster_ID))[2])]))+
  theme(plot.title = element_text(hjust = 0.5, size=28, face="bold"))+
  theme(axis.text.x=element_text(angle=90, vjust = .5, size = 20, face = "italic"), 
        axis.text.y=element_text(angle=90, vjust = .5, size = 20),
        axis.title.y=element_text(angle=90, vjust = 2.5, size = 22, face = "bold"),
        legend.title = element_text(size=21, face = "bold"),
        legend.text = element_text(size = 20))+
  theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"))


## Plot gene expression by cluster.  Group by cluster (1st graph) or by gene (2nd graph)
ggplot(Allclusters, aes(x=Gene,y=Expression, fill = Gene))+
  geom_violin(alpha = 0.4,width = .8, scale="width")+
  geom_point(position = "jitter",  size = .1, alpha = .1)+
  facet_wrap(~as.factor(Cluster_ID), scale = "free")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, vjust = .5), legend.position = "none")+
  ggtitle("Gene Expression by Cluster")+
  ylab("LogNormal Gene Expression")

ggplot(Allclusters, aes(x= as.factor(Cluster_ID),y=Expression, fill = as.factor(Cluster_ID)))+
  geom_violin(alpha = 0.4,width = .8, scale="width")+
  geom_point(position = "jitter",  size = .1, alpha = .1)+
  facet_wrap(~Gene, scale = "free")+
  theme(axis.text.x=element_text(vjust = .5), legend.position = "none")+
  ylab("LogNormal Gene Expression")+
  ggtitle("Gene Expression by Cluster")+
  ylab("LogNormal Gene Expression")+
  xlab("Cluster")


#************* CELL AREA BOXPLOTS FOR ALL CLUSTERS ***************#

# Generate boxplots

par(cex=1.5, font = 2)
boxplot(Normal.Data$ROI_Area~Normal.Data$Cluster_ID, ylab="ROI Area px^2", xlab = "Cluster ID", col=cluster_color_map)
title("Cell Size of Clusters", line =.75)


# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
