####---------- SCAMPR_CorrelationRegressionComparisonsMatrix ---------#####

##### Use this script to compare pairwise gene co-expression patterns from Hiplex RNAscope between two experimental groups
##### alimaran@usc.edu, 01-26-2021 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 77  
#*2) Define the gene columns that you want to remove on line 93

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

# Clear environment
rm(list = ls()) 

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")

## Select names for two experimental groups and identify which animal IDs fall into each group.
Group1name = "CAU" # Choose title for animals that fall into group 1
Group2name = "ELS" # Choose title for animals that fall into group 2

unique(Raw.Data$Animal_ID) # List all Animal IDs in dataset
Group1 = c("18-55", "18-56") # Animal IDs that fall into group 1
Group2 = c("18-43", "18-44") # Animal IDs that fall into group 2

## Plot Settings: Specify methods for smoothing and correlation calculations and the plotting colors for each experimental group.
smooth.method = "lm" # Alternatives are "glm", "gam", and "loess".
cor.method = "pearson" # Some alternatives are  "spearman" and "kendall".
group1col = "blue"
group2col = "red"


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages("mgsub", "dplyr", "stringr", "ggplot2", "GGally", "ggpubr", "pacman")

library(mgsub) # For adding group information
library(dplyr) # For reordering columns in dataframe
library(stringr) # For manipulating column (gene) names
library(ggplot2) # For plotting
library(GGally) # For generating correlation matrix
library(ggpubr) # For generating regression slopes
library(pacman) # For unloading libraries


#************* DATA WRANGLING AND NORMALIZATION ************#

## Check structure and dimensions of matrix
head(Raw.Data)
dim(Raw.Data)

# Add group information
Raw.Data$Group <- Raw.Data$Animal_ID
head(Raw.Data)

Raw.Data$Group <- mgsub(string = Raw.Data$Group, pattern = Group1, replacement = rep(Group1name,length(Group2))) # Sub group1 name for group 1 animal id's
Raw.Data$Group <- mgsub(string = Raw.Data$Group, pattern = Group2, replacement = rep(Group2name,length(Group2))) # Sub group2 name for group 2 animal id's

head(Raw.Data)
Raw.Data <- Raw.Data %>% relocate(Group,.after=Image_ID) # If desired, move group column to proceed Image_ID column
head(Raw.Data)

table(Raw.Data$Animal_ID, Raw.Data$Group)

## Convert Area Fraction gene expression to gene expression in pixels or um. To do this, multiply the Area Fraction columns by the ROI_Area column 
#and divide by 100. 
head(Raw.Data) # Determine which columns contain gene expression information.
geneColumns <- which(!(colnames(Raw.Data) %in% c("X", "Animal_ID", "Image_ID", "Group", "ROI_Area")))
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

# ## If desired, remove columns that represent genes that you do not want to analyze ("Cckar", "Piezo1", "Npy2r") prior to normalization
# genesToRemove <- c("Cckar", "Piezo1", "Npy2r")
# removeIndices <- which(colnames(Raw.Data) %in% genesToRemove)
# Raw.Data <- subset(Raw.Data, select = -c(removeIndices))
# head(Raw.Data)
# geneIndices <- min(geneIndices):(max(geneIndices)-length(genesToRemove)) # Define new gene indices

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log.
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,geneIndices] <- log1p(Normal.Data[,geneIndices]) 
head(Normal.Data)
range(Normal.Data[,geneIndices])


#************* GENERATING REGRESSION/CORRELATION PLOTS ************#

## Create function to generate individual regression/smoothing plots that will end up in the matrix
upperFn <- function(data, mapping, emap=NULL, method = smooth.method, ...) {
  mapping <- ggplot2:::new_aes( c(mapping, emap))
  p <- ggplot(data = data, mapping = mapping) +
    geom_smooth(method = method, ...) +
    stat_regline_equation(label.y = c(max(range(Normal.Data[geneIndices])*.9),max(range(Normal.Data[geneIndices]))*.75),label.x=c(0,0), aes(label = ..eq.label..), parse = T, size=4) +
    scale_color_manual(values = c(group1col,group2col))+
    xlim(range(Normal.Data[geneIndices]))+
    ylim(range(Normal.Data[geneIndices]))+
    theme_classic()# to get the white background and prominent axis
    
  p
}

## Generate and store correlation/regression plot matrix using above functions and and correlations calculated by ggally_cor().
plot1 <- ggpairs(Normal.Data, columns=geneIndices, mapping = aes(colour = Group),axisLabels = "none",
                 
                 upper = list(continuous = wrap(upperFn, method = smooth.method,
                                                fullrange=F, se=FALSE,
                                                emap=aes(color=Group), size=1)),
                 
                 diag = list(continuous = function(data, mapping, ...) {
                   ggally_diagAxis(data = data, mapping = mapping, gridLabelSize=1, displayGrid = FALSE, size = 4) + 
                     scale_colour_manual(values = c("black", "black"))+
                     theme_classic()+
                     theme(text = element_text(face = "italic"),
                           axis.ticks.x=element_blank(),
                           axis.ticks.y=element_blank(),
                           panel.border = element_rect(colour = "black", fill=NA, size=3.5))}),
                 
                 lower = list(continuous = function(data, mapping, ...) {
                   ggally_cor(data = data, mapping = mapping, size = 4, stars = F, method = cor.method) + 
                     scale_colour_manual(values = c(group1col, group2col))})
                 
                 
)

## Print plot and save at 1500 x 1500 dimensions for best text-to-data proportions.
print(plot1 + theme(strip.placement = "outside",
                    strip.text = element_text(size = 22, face="italic")))


#************* CLEAN UP ************#

## Clear environment
rm(list = ls()) 

## Clear packages
p_unload(all)  # Remove all add-ons

## Clear console
cat("\014")  # ctrl+L

## Clear mind :)
