####---------- SCAMPR_hClustHeatmapTopomaps ---------#####

##### Use this script to cluster cells into groups based on HiPlex RNAScope gene expression data and generate single-cell, gene expression heat and topographical maps.
##### alimaran@usc.edu, 01-17-2021 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 64   
#*2) Define the gene columns that you want to remove on line 80 
#*3) Define the border color for cells in the topomap on line 237 and 289
#*4) Use correlational instead of euclidean distance matrix for clustering (see note in line 98)

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

# Clear environment
rm(list = ls()) 

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")

## Select the section that you want to map by specifying it's Image ID (can also specify on line 189)
unique(Raw.Data$Image_ID) # List image IDs for all images.
Section <- "18-56~s1-Sec4"

## Choose hierarchical clustering method: ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
hc_method <- "ward.D"

## Choose and store the  number of clusters desired for dataset.  
# If unsure, use NbClust() package for determining optimal number of clusters for your data.
predicted_clusters <- 6


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages(c("factoextra", "dplyr", "stringr", "dendextend", "colorspace", "gplots", "RImageJROI" ))

library(stringr) # For manipulating column (gene) names
library(dplyr) # For piping commands
library(dendextend) # For assigning correct labels to each cluster
library(colorspace)  # For color schemes
library(factoextra) # For visualizing clusters
library(gplots) # For generating heatmaps
library(RImageJROI) # Used to load ImageJ ROIs into R
library(plotrix) # For adding a color legend


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

## Store total number of cells as a variable
total_cells <- nrow(Normal.Data)  


#************* HEIRARCHICAL CLUSTERING ***************#

## Perform hierarchical clustering using Euclidean distances on gene expression data.  
# NOTE: You can use get_dist(method = "pearson) or get_dist(method = "spearman) instead of dist() to get correlation distance measures.
hc <- Normal.Data[, geneIndices]    %>%   # Select columns with gene expression data
  dist()                        %>% # Compute euclidean distance/dissimilarity matrix or correlation matrix
  hclust(method = hc_method)   # Compute hierarchical clusters

## Identify which cells fall into which cluster and add this to data table. k = the number of clusters.
Cluster_ID <- cutree(hc, k = predicted_clusters)

## Append Cluster ID column to Normal.Data and Raw.Data matrices for downstream plotting, and see how many cells fall into each cluster.
Raw.Data <- cbind(Raw.Data,Cluster_ID) 
Normal.Data <- cbind(Normal.Data,Cluster_ID) 
table(Cluster_ID) # Display number of cells that fall into each cluster

## Generate dendrogram to determine which cluster maps to which label.  This is needed for color and cluster label consistency.
dend <- as.dendrogram(hc) # Generate dendrogram
dend.cutree <- cutree(dend, k = predicted_clusters, order_clusters_as_data = FALSE) # Cut dendrogam into k clusters.  Order data by the order of the labels in the dendrogram
idx <- order(as.numeric(names(dend.cutree))) # Grab the index of each cell, then sort indices based on cell number
dend.cutree <- dend.cutree[idx] # Reorder the dendrogram based on sorted indices (not labels)
tbl <- table(Cluster_ID, dend.cutree) # Match the number of cells per cluster to the number of cells per label
lbls <- apply(tbl,2,which.max) # Match cluster number to label number

## Assign colors to each cluster
dend <- color_branches(dend, k = predicted_clusters,   
                       col= rainbow_hcl(predicted_clusters), 
                       groupLabels = lbls) 
dend <- rotate(dend, 1:total_cells)   # Rotate dendogram and order based on total cells
dend <- hang.dendrogram(dend,hang_height=100000)   # Hang dendogram generated from hclust

## Plot dendrogram to see if it looks as expected.
par(cex=1.2, font = 1.5)
plot(dend, 
     main = "Clustered Dendogram", 
     horiz =  TRUE)


#************* GENERATE HEATMAP (WITH DENDROGRAM) OF CLUSTERS ************#

## Create vector of values to create separators for the clusters in the heatmap
separators <-  rev(as.data.frame(table(Normal.Data$Cluster_ID))[c(as.numeric(lbls)),2])
separators <- cumsum(separators)[-length(separators)] # Get the size of each cluster and sum consecutively to get separator locations on heatmap
colnames(Normal.Data)

## If desired, choose the order of genes on x-axis
dput(colnames(Normal.Data[geneIndices])) # Get list of genes seperated by commas

gene_list <- as.character(c("Pvalb", "Gpr65", "Ntsr1", "Glp1r", "Piezo2", "Htr3b",      
                            "Scn10a", "Trpv1", "Npy2r", "Cckar","Piezo1", "Scn1a")) # Reorder genes

## If desired, convert gene names to italics in heatmap graphics
make_italics <- function(x) {
  as.expression(lapply(x, function(y) bquote(italic(.(y)))))
}

make_italics(gene_list)

## Plot dendrogram heatmap.  Setting Rowv to the variable containing the dendrogram created above
# will order the cells based on the cluster that they fall into and will plot the dendrogram on the y-axis.
# See https://www.rdocumentation.org/packages/gplots/versions/3.1.1/topics/heatmap.2 for more information on heatmap.2() arguments.
gplots::heatmap.2(as.matrix(Normal.Data[,gene_list]), 
                  main = "Nodose Clusters",
                  dendrogram = "row",
                  scale = "column",
                  Rowv = dend,
                  Colv = FALSE, # this is to make sure the columns are not ordered
                  labCol = make_italics(gene_list),
                  trace="none",          
                  margins =c(7,.8),
                  keysize = .6,
                  lhei=c(1,6), 
                  lwid=c(1,2),
                  key.title = NA,
                  key.ylab = NULL,
                  rowsep = separators,
                  sepcolor ="azure",
                  sepwidth=c(0,16),
                  key.xlab = "Gene Expression Z-Score",
                  key.par=list(mgp=c(1.15, .5, 0),
                               mar=c(3, 3, 1, 0)),
                  density.info = "none",
                  symm=F,symkey=F,symbreaks=T,
                  col = redgreen(5) # Set color scheme and number of color breaks here
)



#************* SPATIAL VISUALIZATION - CLUSTER IDENTITY OF EACH CELL MAPPED ONTO CELLULAR ROIS ***************#

## If dataframe contains expression data for multiple sections/images, specify the section that you want to map using
# the Image_ID column, then subset into new data frame
head(Normal.Data)
dim(Normal.Data)
# Section <- "18-56~s1-Sec4" 
Data.Subset <- Normal.Data[grep(Section, Normal.Data$Image_ID), ] # Note: Partial matches to "Section" will also be included
dim(Data.Subset)

## Import cellular ROIs from the ImageJ ROI.zip file that corresponds the section chosen in the above step.
roi_list <- read.ijzip(file.choose())

## If desired, import hand-drawn annotations from the ImageJ ROI.zip file that corresponds to the section chosen in the above step.
annotations <- read.ijzip(file.choose())
length(annotations)

## ImageJ and R have different y-axis orientation.  Flip y-coordinate for cellular and annotation ROIs.
for (i in 1:length(roi_list)){
  roi_list[[i]]$coords[,2] <- roi_list[[i]]$coords[,2]*-1
  roi_list[[i]]$yrange <- roi_list[[i]]$yrange*-1
}


for (i in 1:length(annotations)){
  annotations[[i]]$coords[,2] <- annotations[[i]]$coords[,2]*-1
  annotations[[i]]$yrange <- annotations[[i]]$yrange*-1
}

## Generate a colormap by assigning a color to each cluster.  You can change color type by changing palette = .
# See https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html for more information on color palettes.
cluster_color_map <- data.frame(Cluster_ID = as.numeric(lbls), 
                                Cluster_Color = c(rainbow_hcl(predicted_clusters)))

cluster_color_map

## Plot all cell outlines to generate plot area for spatial topomap
dev.off()
par(mar=c(4,4,3,4.5), xaxt="n", yaxt = "n")
plot(roi_list, add = FALSE, asp = 1, col="white") # col = "white" so only titles will show up.  If you want to see the outlines
# of the zero-expressing cells that were removed in line 71 in the final plot, have col = "color" (ex. col = "blue").
par(cex=1, font = 2)
title(main = paste("All Clusters Normal Euclidean"), line = 1)
par(cex=1.3, font = 2)
title(ylab = "Dorsal-Ventral Axis", line=0.5)
title(xlab = "Rostral-Caudal Axis", line=0.5)

## Fill all cells with color corresponding to cluster ID.  You can have the borders be a specified color of interest
# or the same color as the rest of the cells using border = . (this may take a few minutes)
for (row in 1:nrow(Data.Subset)) {
  current_cell_id <- Data.Subset[row, "X"]
  current_cluster_id <- Data.Subset[row, "Cluster_ID"]
  current_cluster_color <- cluster_color_map[which(cluster_color_map$Cluster_ID == current_cluster_id), "Cluster_Color"]
  polygon(roi_list[[current_cell_id]][["coords"]][,"x"], roi_list[[current_cell_id]][["coords"]][,"y"],
          border = current_cluster_color, col = current_cluster_color)
}

## Add annotations (specify color, line-type, and line-width using lty, lwd, and col)
for (i in 1:length(annotations)) {
  current_line_id <- Data.Subset[i, "X"]
  lines(annotations[[current_line_id]][["coords"]][,"x"], annotations[[current_line_id]][["coords"]][,"y"],
        lty= "dashed", lwd = 3, col = "black")
}

## Add color legend
breaks = as.numeric(lbls)
colors = rainbow_hcl(predicted_clusters)
color.legend(xl = par("usr")[2], 
             yb = par("usr")[3]*0.7, 
             xr = par("usr")[2]*1.05, 
             yt = par("usr")[3]*0.3, # xl, yb, xr, and yt values determine the size of the legend
             rev(breaks),rev(colors),
             gradient="y",
             align = "rb") # "align = " determines the location of the legend

## Create separate cluster topology maps for each individual cluster.
# NOTE: Un-comment lines 300-305 if you want to add manually drawn annotations

cluster = seq(1,predicted_clusters,1) #Create a vector of cluster numbers.  You can also set this to a cluster number

for(cluster in 1:length(cluster)) {
  
  # create color map for single clusters
  ncluster = length(lbls)
  cluster_index <- which(as.numeric(lbls)==cluster)
  other_indices <- which(as.numeric(lbls)!=cluster)
  colormap <- replace(as.numeric(lbls),other_indices, "white")# replace other clusters with "white"
  colormap <- replace(colormap,cluster_index, rainbow_hcl(ncluster)[cluster_index]) # replace cluster with rainbow_hcl color corresponding to its index
  
  cluster_color_map <- data.frame(Cluster_ID = as.numeric(lbls),
                                  Cluster_Color = colormap)

  # Plot all cells within ijzip list
  par(cex=1.1, font = 2)
  par(mar=c(5,6,5,4),xaxt="n", yaxt = "n")
  plot(roi_list, add = FALSE,asp = 1, col="white");
  title(ylab = "Rostral-Caudal Axis", line=4 )
  title(main = paste("Cluster",cluster, "(Normalized, Euclidean)"), line = 1)
  title(xlab = "Dorsal-Ventral Axis", line = 3)
  
  # Plot all cells, colored based on cluster
  for (row in 1:nrow(Data.Subset)) {
    current_cell_id <- Data.Subset[row, "X"]
    current_cluster_id <- Data.Subset[row, "Cluster_ID"]
    current_cluster_color <- cluster_color_map[which(cluster_color_map$Cluster_ID == current_cluster_id), "Cluster_Color"]
    polygon(roi_list[[current_cell_id]][["coords"]][,"x"], roi_list[[current_cell_id]][["coords"]][,"y"],
            border = current_cluster_color, col = current_cluster_color)
  }
  
  ## Add color legend
  breaks = as.numeric(lbls)
  colors = rainbow_hcl(predicted_clusters)
  color.legend(xl = par("usr")[2], 
               yb = par("usr")[3]*0.7, 
               xr = par("usr")[2]*1.05, 
               yt = par("usr")[3]*0.3, # xl, yb, xr, and yt values determine the size of the legend
               rev(breaks),rev(colors),
               gradient="y",
               align = "rb") # "align = " determines the location of the legend
  
  # # Add annotation lines
  # for (i in 1:length(annotations)) {
  #   current_line_id <- Data.Subset[i, "X"]
  #   lines(annotations[[current_line_id]][["coords"]][,"x"], annotations[[current_line_id]][["coords"]][,"y"],
  #         lty= "dashed", lwd = 3, col = "black")
  # }

}


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
