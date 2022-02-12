####---------- SCAMPR_CorrelationCirclePlots ---------#####

##### Use this script to generate single-cell, pairwise correlation matrices for HiPlex RNAScope gene expression data
##### alimaran@usc.edu, 01-21-2021 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 48   
#*2) Define the gene columns that you want to remove on line 64 

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

# Clear environment
rm(list = ls()) 

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages(c("stringr","corrplot","colorspace", "tcltk", "pacman"))

library(stringr) # For manipulating column (gene) names
library(corrplot) # For plotting correlation matrices
library(colorspace) # For color schemes
library(tcltk) # For choosing a directory using Windows file explorer
library(pacman) # For hierarchical clustering


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

## If desired, remove columns that represent genes that you do not want to analyze ("Cckar", "Piezo1", "Npy2r")
genesToRemove <- c("Cckar", "Piezo1", "Npy2r")
removeIndices <- which(colnames(Raw.Data) %in% genesToRemove)
Raw.Data <- subset(Raw.Data, select = -c(removeIndices))
head(Raw.Data)
geneIndices <- min(geneIndices):(max(geneIndices)-length(genesToRemove)) # Define new gene indices

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log.
# Include ROI_AREA column
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)] <- log1p(Normal.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)]) 
range(Normal.Data[,geneIndices])

## Scale raw (non-normalized) dataset if you wish to used scaled values for correlations.
Scaled.Raw.Data <- Raw.Data
head(Scaled.Raw.Data)
Scaled.Raw.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)] <- scale(Scaled.Raw.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)], center = F, scale = T)
range(Scaled.Raw.Data[,geneIndices])

#************* CALCULATED AND PLOT PAIRWISE CORRELATIONS ************#

## Calculate pairwise pearson correlation coefficients for normalized and scaled data
Normal.Data.cor <- cor(Normal.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)], method = "pearson") 
Scaled.Raw.Data.cor <- cor(Scaled.Raw.Data[,c(which(colnames(Normal.Data) == "ROI_Area"),geneIndices)], method = "pearson")

## Plot correlation matrices for normalized and scaled data.  Save as .tif file using tiff() and dev.off() for high resolution images.
setwd(tk_choose.dir()) # Set working directory.  Plot will be saved here.

## Normalized data plot
tiff(filename = "NormalDataCorPlot.tif", width = 4000, height = 4000, res = 600) # Set plot parameters and plot name.

corrplot(Normal.Data.cor, method = "circle", type = "upper", 
         col = diverging_hcl(10, palette = "Blue-Red 3"),
         tl.col = "black",
         tl.srt = 45,
         #tl.pos = "n",
         order = "hclust",
         hclust.method = "ward.D2",
         title = "LogNormal Correlations Nodose",
         mar = c(10,0,3,0),
         outline = T)

dev.off() # Run dev.off() to save plot in your working directory.

## Scaled data plot
tiff(filename = "ScaledDataCorPlot.tif", width = 4000, height = 4000, res = 600) # Set plot parameters and plot name.

corrplot(Scaled.Raw.Data.cor, method = "circle", type = "upper", 
         col = diverging_hcl(10, palette = "Blue-Red 3"),
         tl.col = "black",
         tl.srt = 45,
         #tl.pos = "n",
         order = "hclust",
         title = "Scaled Correlations Nodose",
         mar = c(10,0,3,0),
         outline = T)

dev.off()# Run dev.off() to save plot in your working directory.

# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)



