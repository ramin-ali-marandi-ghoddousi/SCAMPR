####---------- SCAMPR_ExpressionTopomaps ---------#####

##### Use this script to generate single-cell, gene expression topographical maps for HiPlex RNAScope gene expression data
##### alimaran@usc.edu, 01-07-2021 

                
#************* OPTIONAL USER INPUTS ************#
                
#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 58   
#*2) Define the gene columns that you want to remove on line 74

#**NOTE 2: Lines that start with "If desired" can be skipped                

                                
#*************REQUIRED USER INPUTS ************#
                
# Clear environment
rm(list = ls()) 
                
## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")
                
## Select the section that you want to map by specifying it's Image ID (can also specify on line 95)
unique(Raw.Data$Image_ID) # List image IDs for all images.
Section <- "18-56~s1-Sec4" 
                
## Select the gene whose expression you want to map (can also specify on line 121)
Gene = "Pvalb"


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages(c("reshape2", "ggplot2", "dplyr", "stringr", "pacman", "colorspace", "plotrix", "RImageJROI"))

library(reshape2)# For reshaping dataframe to long format for easy plotting in ggplot
library(ggplot2)# For plotting
library(dplyr) # For getting means by group
library(stringr) # For manipulating column (gene) names
library(pacman) # For hierarchical clustering and loading/unloading packages
library(colorspace) # For color schemes
library(plotrix) # For adding a color legend
library(RImageJROI) # Used to load ImageJ ROIs into R


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

## If desired, remove columns that represent genes that you do not want to analyze ("Cckar", "Piezo1", "Npy2r") prior to normalization
genesToRemove <- c("Cckar", "Piezo1", "Npy2r")
removeIndices <- which(colnames(Raw.Data) %in% genesToRemove)
Raw.Data <- subset(Raw.Data, select = -c(removeIndices))
head(Raw.Data)
geneIndices <- min(geneIndices):(max(geneIndices)-length(genesToRemove)) # Define new gene indices

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,geneIndices] <- log1p(Normal.Data[,geneIndices]) 
range(Normal.Data[,geneIndices])

## Save ROI_Area column
ROI_Area <- Normal.Data$ROI_Area


#************* SPATIAL VISUALIZATION - LOG NORMALIZED GENE EXPRESSION LEVELS MAPPED ONTO CELLULAR ROIS ***************#

## If dataframe contains expression data for multiple sections/images, specify the section that you want to map using
# the Image_ID column, then subset into new data frame
head(Normal.Data)
dim(Normal.Data)
#Section <- "18-56~s1-Sec4" 
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


## Choose gene that you want to plot and assign each cell a bin label (0 to range) based on expression level of that gene.
# Add bin information as column to dataframe.
#Gene = "Pvalb"

start <- floor(range(Data.Subset[Gene]))[1] # round down the lowest expression value and save
end <- ceiling(range(Data.Subset[Gene]))[2] # round up the highest expression value and save
data.range = end-start # determine modified range of values using start and end

# Assign each cell into a bin
Data.Subset$categories = cut(as.numeric(unlist(Data.Subset[Gene])),
                             breaks = seq(start, end, length.out = data.range + 1),
                             labels = seq(0,data.range-1,1),
                             include.lowest = T)
head(Data.Subset)
table(Data.Subset$categories)


## Generate a colormap by assigning a color to each bin label.  You can change color type by changing palette = .
# See https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html for more information on color palettes.
regression_color_map <- data.frame(categories = seq(0,data.range-1,1), 
                                   categories_color = rev(sequential_hcl(data.range, palette = "Heat")))
regression_color_map

## Plot all cell outlines to generate plot area for spatial topomap
dev.off()
par(mar=c(4,4,3,4.5), xaxt="n", yaxt = "n")
plot(roi_list, add = FALSE, asp = 1, col="white") # col = "white" so only titles will show up.  If you want to see the outlines

# of the zero-expressing cells that were removed in line 66 in the final plot, have col = "color" (ex. col = "blue").
par(cex=1, font = 2)
title(main = bquote(bold(paste("V1 Cortex ", bolditalic(.(Gene)), " Expression"))), line=0.75)
par(cex=1.3, font = 2)
title(ylab = "Dorsal-Ventral Axis", line=0.5)
title(xlab = "Rostral-Caudal Axis", line=0.5)

## Fill all cells with color corresponding to bin value.  You can have the borders be a specified color of interest
# or the same color as the rest of the cells using border = .
for (row in 1:nrow(Data.Subset)) {
  current_cell_id <- Data.Subset[row, "X"]
  current_cluster_id <- Data.Subset[row, "categories"]
  current_cluster_color <- regression_color_map[which(regression_color_map$categories == current_cluster_id), "categories_color"]
  polygon(roi_list[[current_cell_id]][["coords"]][,"x"], roi_list[[current_cell_id]][["coords"]][,"y"],
          border = current_cluster_color, col = current_cluster_color)
}

## Add annotations (specify color, line-type, and line-width using lty, lwd, and col)
for (i in 1:length(annotations)) {
  current_line_id <- Data.Subset[i, "X"]
  lines(annotations[[current_line_id]][["coords"]][,"x"], annotations[[current_line_id]][["coords"]][,"y"],
        lty= "dashed", lwd = 4, col = "green")
}

## Add color legend
breaks = round(seq(start, end, length.out = data.range+1))
colors = rev(sequential_hcl(data.range+1, palette = "Heat"))
color.legend(xl = par("usr")[2], 
             yb = par("usr")[3]*0.7, 
             xr = par("usr")[2]*1.05, 
             yt = par("usr")[3]*0.3, # xl, yb, xr, and yt values determine the size of the legend
             rev(breaks),rev(colors),
             gradient="y",
             align = "rb") # "align = " determines the location of the legend

#************* CLEAN UP ****************#

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)

