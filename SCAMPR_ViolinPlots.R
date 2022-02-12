####---------- SCAMPR_ViolinPlots ---------#####

#### Use this script to generate violin plots for HiPlex RNAScope gene expression data
#### alimaran@usc.edu, 01-03-2021
                 
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
                 
## Install packages, and load libraries

#install.packages(c("reshape2", "ggplot2", "dplyr", "stringr", "pacman", "colorspace", "plotrix"))
library(reshape2)# For reshaping dataframe to long format for easy plotting in ggplot
library(ggplot2)# For plotting
library(dplyr) # For getting means by group
library(stringr) # For manipulating column (gene) names
library(pacman) # For hierarchical clustering
library(colorspace) # For color schemes
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

## Reshape Normalized datasets to long format, making it compatible with ggplot() for downstream plotting.
Genes <- colnames(Normal.Data[geneIndices]) # Get list column names representing each gene
Genes
Normal.long <- melt(Normal.Data,
                id.vars = c("X", "Image_ID","ROI_Area"),
                measure.vars = Genes,
                variable.name = "Gene",
                value.name = "Expression")

## If desired, remove cells that have an expression of 0 for each particular gene.  Because the data is in long format,
# this will only remove the cell from the dataset in a gene-specific manner.  So if Cell X has 0 expression for Pvalb
# but an expression of 2 for Ntsr1, it will only be removed from the PVALB expressing subset.
Normal.long.no.zeros<-subset(Normal.long, Normal.long[,"Expression"]>0)


#************* PLOT! ****************#

## Generate violin plots for data with and without zero-expressing cells.
ggplot(Normal.long, aes(x = reorder(Gene,-Expression, median), y = Expression, fill = Gene)) +
  theme_bw() +
  geom_violin(alpha = 0.4, width = .8, trim = T, scale = "width") +
  geom_point(aes(fill = Gene), size = 0.1, shape = 21, alpha = 0.07,
             position=position_jitter(width=.35)) +
  labs(title = "Gene Expression Visual Cortex",
       x = NULL,
       y = "Normalized Expression") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1.1, size = 14, face = "italic"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18, vjust=2, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, face = "bold"),
        legend.position = "")+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))

ggplot(Normal.long.no.zeros, aes(x = reorder(Gene,-Expression, median), y = Expression, fill = Gene)) +
  theme_bw() +
  geom_violin(alpha = 0.4, width = .8, trim = T, scale = "width") +
  geom_point(aes(fill = Gene), size = 0.1, shape = 21, alpha = 0.07,
             position=position_jitter(width=.35)) +
  labs(title = "Gene Expression Visual Cortex",
       x = NULL,
       y = "Normal Expression No Zeros") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1.1, size = 14, face = "italic"),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18, vjust=2, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, face = "bold"),
        legend.position = "")+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))


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
