####---------- SCAMPR_LoessCoexpressionTopomaps ---------#####

##### Use this script to generate single-cell, gene co-expression topographical maps for HiPlex RNAScope gene expression data
##### alimaran@usc.edu, 01-07-2021 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 59   
#*2) Define the gene columns that you want to remove on line 75 
#*3) Define the border color for cells in the topomap on line 193

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

# Clear environment
rm(list = ls()) 

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")

## Select the section that you want to map by specifying it's Image ID (can also specify on line 95)
unique(Raw.Data$Image_ID) # List image IDs for all images.
Section <- "18-56~s1-Sec4" 

## Select the two genes whose co-expression you want to map 
colnames(Raw.Data) # List genes by listing column headers.
gene1 <- "Scn10a" # Case sensitive
gene2 <- "Trpv1" # Case sensitive


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages(c( "ggplot2", "stringr", "pacman", "colorspace", "plotrix", "RImageJROI"))

library(ggplot2)# For plotting 
library(stringr) # For manipulating column (gene) names
library(colorspace) # For color schemes
library(plotrix) # For adding a color legend
library(RImageJROI) # Used to load ImageJ ROIs into R
library(pacman) # For unloading packages


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

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,geneIndices] <- log1p(Normal.Data[,geneIndices]) 
range(Normal.Data[,geneIndices])


#************* COMPUTE AND PLOT LOESS REGRESSION SCATTERPLOT FOR TWO GENES ************#

## Compute local polynomial regression on log-normalized data
fit <- loess(data = Normal.Data, formula = Normal.Data[,gene1]~Normal.Data[,gene2]) # fit the model
Normal.Data$predicted <- predict(fit)   # Add the predicted values as column
Normal.Data$Residuals <- residuals(fit) # Add the residual values as column

## Bin each cell cell based on its residual size (distance from the loess line)
range(Normal.Data$Residuals) # find the range of your residuals
start <- floor(range(Normal.Data$Residuals))[1] # save the lowest residual
end <- ceiling(range(Normal.Data$Residuals))[2] # save the highest residual
data.range <- abs(start) + end

## Plot scatterplot with LOESS curve. Use color and size parameters in geom_point() and to assign color and size to 
# each point (cell) based on it's distance from the LOESS curve and the color gradient map created in scale_color_gradientn(). 

ggplot(Normal.Data, aes(x = Normal.Data[,gene2] , y = Normal.Data[,gene1])) +
  geom_smooth(method = "loess", se = FALSE, color = "green") +        # regression line  
  geom_segment(aes(xend = get(gene2), yend = predicted), alpha = .1) +  # draw line from points to regression line
  geom_point(aes( color=Residuals, alpha=abs(Residuals)), size = 4) +  # assign size and color of the points
  scale_color_gradientn(colors = diverging_hcl(data.range, palette = "Blue-Red"),
                        breaks = round(seq(start, end, length.out = data.range + 2)))+   # color of the points mapped to residual size - grey smaller, blue larger negative, red larger positive
  guides(size = "none", alpha = "none") +   # Size legend removed
  labs(title = bquote(bold("Nodose"~bolditalic(.(gene1))~"and"~bolditalic(.(gene2))~"Log-Normalized Coexpression")))+
  labs(x = bquote(bold("Log-Normalized"~bolditalic(.(gene2))~"Expression")))+
  labs(y = bquote(bold("Log-Normalized"~bolditalic(.(gene1))~"Expression")))+
  geom_point(aes(y = predicted), shape = 1) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.text.x = element_text( size = 14),
        axis.text.y = element_text( size = 14),
        axis.title.x = element_text( size = 18, face = "bold", vjust = -1.5),
        axis.title.y = element_text(size = 18, face = "bold", vjust = 3),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 16, face = "bold"))+
  theme(plot.margin = unit(c(1,1,1,1), "cm"))


#************* SPATIAL VISUALIZATION - COEXPRESSION LEVELS (RESIDUALS) MAPPED ONTO CELLULAR ROIS ***************#

## If dataframe contains expression data for multiple sections/images, specify the section that you want to map using
# the Image_ID column, then subset into new data frame
head(Normal.Data)
dim(Normal.Data)
#Section <- "18-56~s1-Sec4 
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


## Assign each cell into a bin based on it's residual size and add as column to the dataframe.
range(Data.Subset$Residuals) # find the range of your residuals
start <- floor(range(Data.Subset$Residuals))[1] # round down and save the lowest residual
end <- ceiling(range(Data.Subset$Residuals))[2] # round up and save the highest residual
data.range <- end-start # determine modified range of residuals using start and end

# Assign cells into bins
Data.Subset$categories = cut(Data.Subset$Residuals,
                             breaks = seq(start, end, length.out = data.range),
                             labels = seq(1,data.range-1,1))
head(Data.Subset)
table(Data.Subset$categories)

## Generate a colormap by assigning a color to each bin label.  You can change color type by changing palette = .
# See https://colorspace.r-forge.r-project.org/articles/hcl_palettes.html for more information on color palettes.
regression_color_map <- data.frame(categories = seq(1,data.range,1), 
                                   categories_color = diverging_hcl(data.range, palette = "Blue-Red"))
regression_color_map


## Plot all cell outlines to generate plot area for spatial topomap
dev.off()
par(mar=c(4,4,3,4.5), xaxt="n", yaxt = "n")
plot(roi_list, add = FALSE, asp = 1, col="white") # col = "white" so only titles will show up.  If you want to see the outlines
# of the zero-expressing cells that were removed in line 65 in the final plot, have col = "color" (ex. col = "blue").
par(cex=1, font = 2)
title(main = bquote(bold("Nodose"~bolditalic(.(gene1))~"and"~bolditalic(.(gene2))~"Log-Normalized Coexpression")),line = 0.75)
par(cex=1.3, font = 2)
title(ylab = "Dorsal-Ventral Axis", line=0.5)
title(xlab = "Rostral-Caudal Axis", line=0.5)

## Fill all cells with color corresponding to bin (residual) value.  You can have the borders be a specified color of interest
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
colors = diverging_hcl(data.range+1, palette = "Blue-Red")
color.legend(xl = par("usr")[2], 
             yb = par("usr")[3]*0.7, 
             xr = par("usr")[2]*1.05, 
             yt = par("usr")[3]*0.3, # xl, yb, xr, and yt values determine the size of the legend
             rev(breaks),rev(colors),
             gradient="y",
             align = "rb") # "align = " determines the location of the legend


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
