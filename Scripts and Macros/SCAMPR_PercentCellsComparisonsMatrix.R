####---------- SCAMPR_PercentCellsComparisonsMatrix ---------#####

##### Use this script to compare the per-animal percentages of cells expressing each gene between experimental groups
##### alimaran@usc.edu, 02-02-2022 


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

## Plot Settings: Set statistical test for comparing means and the plotting colors for each experimental group.
group1col = "blue"
group2col = "red"
statTest = "t.test" # Alternatives are "wilcox.test" for two groups. Can use "aov", "anova",or "kruskal.test" for multiple groups.


#************* REQUIRED PACKAGES ************#

## Install packages and load libraries

#install.packages("mgsub", dplyr", "stringr", "reshape2", data.able", "ggplot2", "tidyverse", "ggpubr",  "pacman")

library(mgsub) # For adding group information
library(dplyr) # For calculating group means
library(stringr) # For manipulating column (gene) names
library(reshape2) # For converting data into long format
library(data.table) # To convert data.frame to data.table
library(ggplot2) # For generating plots
library(ggpubr) # For generating group statistics on plots
library(pacman) # For unloading libraries


#************* DATA WRANGLING AND NORMALIZATION ************#

## Check structure and dimensions of matrix
head(Raw.Data)
dim(Raw.Data)

# Add group information
Raw.Data$Group <- Raw.Data$Animal_ID
head(Raw.Data)

Raw.Data$Group <- mgsub(string = Raw.Data$Group, pattern = Group1, replacement = rep(Group1name,length(Group1))) # Sub group1 name for group 1 animal id's
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


#************* CALCULATE PERCENTAGE OF CELLS EXPRESSION EACH GENE ************#

## Generate vector of gene names in alphabetical order.
geneNames <- sort(colnames(Normal.Data[geneIndices]))

## Reshape data to long format using vector of gene names
Normal.Data.long <- reshape2::melt(Normal.Data,
                         id.vars = c("X", "Animal_ID", "Image_ID", "Group", "ROI_Area"),
                         measure.vars = geneNames,
                         value.name = "Expression",
                         variable.name = "Gene")

head(Normal.Data.long)

## Calculate the number of cells in each animal that express each gene and save as a table.
expressionCounts <- Normal.Data.long %>%
  group_by(Group, Animal_ID, Gene) %>%
  filter(Expression > 0) %>%
  summarize(counts = n())
expressionCounts

## Calculate the total number of cells in each animal and save as a table.
totalCells <- Normal.Data.long %>%
  group_by(Group, Animal_ID, Gene) %>%
  summarize(total_cells = n())
totalCells

## Merge the expression cell counts and total cell counts tables and use to generate new column with percentages of cells in each animal that express each gene.
expressionPercentages <- cbind(expressionCounts, totalCells["total_cells"])
expressionPercentages$Percentages <- expressionPercentages$counts/expressionPercentages$total_cells
expressionPercentages

## Convert new dataframe to long format for plotting
expressionPercentages <- reshape2::melt(expressionPercentages,
                               id.vars = c("Group","Animal_ID", "counts", "total_cells", "Percentages"),
                               measure.vars = "Gene",
                               value.name = "Gene")

#************* SET INDIVIDUAL Y-AXES AND PLOT PERCENT OF CELLS EXPRESSING EACH GENE ************#

## Generate data table of means y-axis plot minimums and maximums for each gene
expressionPercentages <- data.table(expressionPercentages)
expressionPercentages[,y_min := min(Percentages*0.75), by = Gene]
expressionPercentages[,y_max := max(Percentages*1.25), by = Gene]
head(expressionPercentages)

## Plot percent of cells expressing each gene as boxplots!
ggplot(expressionPercentages, aes(x = Group, y = Percentages, fill = Group))+
  geom_boxplot(position = position_dodge(1),alpha = 0.6, outlier.color = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 1), size = 2, alpha = .8)+
  scale_fill_manual(values = c(group1col, group2col))+
  stat_compare_means(method = statTest, label.y.npc = 0.9, size = 4.5)+
  theme_bw()+
  ggtitle("% of All Cells Expressing Gene")+
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.title.x = element_text(size=16, face = "bold"),
        axis.title.y = element_text(size=16, face = "bold"),
        axis.text.x = element_text(size=13, face = "bold"),
        axis.text.y = element_text(size=13, face = "bold"),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15, face = "bold"),
        strip.text = element_text(size = 16, face = "italic"))+
  facet_wrap(~Gene, scales = "free_y")+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))

## Calculate mean % of cells expressing each gene per animal.
expressionPercentagesMeans <- expressionPercentages %>%
  group_by(Group, Gene) %>%
  summarize(Means = mean(Percentages))

print(expressionPercentagesMeans, n=nrow(expressionPercentagesMeans))


#************* CLEAN UP ****************#

## Clear environment
rm(list = ls()) 

## Clear packages
p_unload(all)  # Remove all add-ons

## Clear plots
dev.off()  # But only if there IS a plot

## Clear console
cat("\014")  # ctrl+L

## Clear mind :)
