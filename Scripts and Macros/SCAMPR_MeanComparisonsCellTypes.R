####---------- SCAMPR_MeanComparisonsCellTypes ---------#####

##### Use this script to compare the per-animal mean expression of select genes from Hiplex RNAscope between two experimental groups
##### alimaran@usc.edu, 02-02-2022 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 85   
#*2) Define the expression cutoff of a gene that will be used to define a cell-type (lines 111-112)

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

## Select cell-type marker genes whose expression levels will be used to subset data.
cellType1_gene <- "Scn1a" # Case sensitive
cellType2_gene <- "Scn10a" # Case sensitive

## Choose genes whose means will be plotted in each cell type
Gene1 <- "Pvalb" # For cell type 1. Case sensitive
Gene2 <- "Ntsr1" # For cell type 2. Case sensitive

## Plot Settings: Set statistical test for comparing means and the plotting colors for each experimental group.
group1col = "blue"
group2col = "red"
statTest = "t.test" # Alternatives are "wilcox.test" for two groups. Can use "aov", "anova",or "kruskal.test" for multiple groups.


#************* REQUIRED PACKAGES ************#

## Install packages, and load libraries

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

## Normalization of dataset. Add pseudocount of 1 to each gene expression value and normalize using natural log.
Normal.Data <- Raw.Data # generate new data table to store values
Normal.Data[,geneIndices] <- log1p(Normal.Data[,geneIndices]) 
head(Normal.Data)
range(Normal.Data[,geneIndices])


#************* GENERATE CELL-TYPE SUBSETS BASED ON EXPRESSION LEVELS OF SPECIFIC GENES ************#

## Select the  minimum expression level of a gene required for a cell to be included in the new subset.
gene1min <- 0 # If 0, all cells that have > 0 expression of cell-type 1 marker gene will be included in subset
gene2min <- 0

## Generate new dataframes with subsets of cells that express cell-type 1 marker gene or cell-type 2 marker gene at a desired level
cellType1_subset <- Normal.Data %>% 
  filter(get(cellType1_gene) > gene1min)

nrow(Normal.Data)
nrow(cellType1_subset)

cellType2_subset <- Normal.Data %>% 
  filter(get(cellType2_gene) > gene2min)

nrow(Normal.Data)
nrow(cellType2_subset)

## Generate new dataframes with subsets of cells that DON'T express cell-type 1 and cell-type 2 marker genes
non_cellType1_subset <- Normal.Data %>% 
  filter(get(cellType1_gene) == 0)

nrow(Normal.Data)
nrow(non_cellType1_subset)

non_cellType2_subset <- Normal.Data %>% 
  filter(get(cellType2_gene) == 0)

nrow(Normal.Data)
nrow(non_cellType2_subset)


#************* CALCULATE GROUP MEANS ************#

## Generate vector of gene names in alphabetical order.
geneNames <- sort(colnames(Normal.Data[geneIndices]))

## Reshape dataframes into long format.
cellType1_subset_long <- reshape2::melt(cellType1_subset,
                          id.vars = c("X", "Image_ID","ROI_Area","Group","Animal_ID"),
                          measure.vars = geneNames,
                          value.name = "Expression",
                          variable.name = "Gene")

cellType2_subset_long <- reshape2::melt(cellType2_subset,
                                    id.vars = c("X", "Image_ID","ROI_Area","Group","Animal_ID"),
                                    measure.vars = geneNames,
                                    value.name = "Expression",
                                    variable.name = "Gene")

non_cellType1_subset_long <- reshape2::melt(non_cellType1_subset,
                                    id.vars = c("X", "Image_ID","ROI_Area","Group","Animal_ID"),
                                    measure.vars = geneNames,
                                    value.name = "Expression",
                                    variable.name = "Gene")

non_cellType2_subset_long <- reshape2::melt(non_cellType2_subset,
                                    id.vars = c("X", "Image_ID","ROI_Area","Group","Animal_ID"),
                                    measure.vars = geneNames,
                                    value.name = "Expression",
                                    variable.name = "Gene")


## Calculate mean expression of Gene 1 and Gene 2 in the different cell types
cellType1_subset_means <- cellType1_subset_long  %>%
  group_by(Group, Animal_ID, Gene) %>%
  filter(Gene %in% Gene1 )%>%
  summarize(Means = mean(Expression))
head(cellType1_subset_means)

cellType2_subset_means <- cellType2_subset_long  %>%
  group_by(Group, Animal_ID, Gene) %>%
  filter(Gene %in% Gene2)%>%
  summarize(Means = mean(Expression))
head(cellType2_subset_means)

non_cellType1_subset_means <- non_cellType1_subset_long  %>%
  group_by(Group, Animal_ID, Gene) %>%
  filter(Gene %in% Gene1)%>%
  summarize(Means = mean(Expression))
head(non_cellType1_subset_means)

non_cellType2_subset_means <- non_cellType2_subset_long  %>%
  group_by(Group, Animal_ID, Gene) %>%
  filter(Gene %in% Gene2)%>%
  summarize(Means = mean(Expression))
head(non_cellType2_subset_means)

#************* PLOT GROUP MEANS ************#

## Convert means dataframe into a data table and generate y-axis plot minimums and maximums for each gene
cellType1_subset_means <- data.table(cellType1_subset_means)
cellType1_subset_means[,y_min := min(Means*0.75), by = Gene]
cellType1_subset_means[,y_max := max(Means*1.25), by = Gene]

cellType2_subset_means <- data.table(cellType2_subset_means)
cellType2_subset_means[,y_min := min(Means*0.75), by = Gene]
cellType2_subset_means[,y_max := max(Means*1.25), by = Gene]

non_cellType1_subset_means <- data.table(non_cellType1_subset_means)
non_cellType1_subset_means[,y_min := min(Means*0.75), by = Gene]
non_cellType1_subset_means[,y_max := max(Means*1.25), by = Gene]

non_cellType2_subset_means <- data.table(non_cellType2_subset_means)
non_cellType2_subset_means[,y_min := min(Means*0.75), by = Gene]
non_cellType2_subset_means[,y_max := max(Means*1.25), by = Gene]

## Plot means of selected genes in cell type 1 and cell type 2 
ggplot(non_cellType2_subset_means, aes(x = Group, y=Means, fill = Group))+
  geom_boxplot(position = position_dodge(1),alpha = 0.6, outlier.color = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 1), size = 2, alpha = .8)+
  scale_fill_manual(values = c(group1col, group2col))+
  stat_compare_means(method = statTest, label.y.npc = .9,label.x = 1.25, size = 4.5)+
  theme_bw()+
  labs(title = bquote(atop(bold("Mean"~bolditalic(.(Gene1))~"Expression"), 
                           bold("in"~bolditalic(.(cellType1_gene)*"+")~"Cells"))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15),
        strip.text = element_text(size = 15, face = "italic"))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))

ggplot(cellType2_subset_means, aes(x = Group, y=Means, fill = Group))+
  geom_boxplot(position = position_dodge(1),alpha = 0.6, outlier.color = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 1), size = 2, alpha = .8)+
  scale_fill_manual(values = c(group1col, group2col))+
  stat_compare_means(method = statTest, label.y.npc = .9,label.x = 1.25, size = 4.5)+
  theme_bw()+
  labs(title = bquote(atop(bold("Mean"~bolditalic(.(Gene2))~"Expression"), 
                           bold("in"~bolditalic(.(cellType2_gene)*"+")~"Cells"))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15),
        strip.text = element_text(size = 15, face = "italic"))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))


## Plot means of selected genes in cells that ARE NOT cell type 1 and cell type 2 
ggplot(non_cellType1_subset_means, aes(x = Group, y=Means, fill = Group))+
  geom_boxplot(position = position_dodge(1),alpha = 0.6, outlier.color = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 1), size = 2, alpha = .8)+
  scale_fill_manual(values = c(group1col, group2col))+
  stat_compare_means(method = statTest, label.y.npc = .9,label.x = 1.25, size = 4.5)+
  theme_bw()+
  labs(title = bquote(atop(bold("Mean"~bolditalic(.(Gene1))~"Expression"), 
                           bold("in"~bolditalic(.(cellType1_gene)*"-")~"Cells"))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15),
        strip.text = element_text(size = 15, face = "italic"))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))

ggplot(non_cellType2_subset_means, aes(x = Group, y=Means, fill = Group))+
  geom_boxplot(position = position_dodge(1),alpha = 0.6, outlier.color = NA)+
  geom_point(position = position_jitterdodge(dodge.width = 1), size = 2, alpha = .8)+
  scale_fill_manual(values = c(group1col, group2col))+
  stat_compare_means(method = statTest, label.y.npc = .9,label.x = 1.25, size = 4.5)+
  theme_bw()+
  labs(title = bquote(atop(bold("Mean"~bolditalic(.(Gene2))~"Expression"), 
                           bold("in"~bolditalic(.(cellType2_gene)*"-")~"Cells"))))+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15),
        strip.text = element_text(size = 15, face = "italic"))+
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))


#************* CLEAN UP ****************#

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L
