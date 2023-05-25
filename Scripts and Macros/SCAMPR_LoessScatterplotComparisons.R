####---------- SCAMPR_LoessScatterplotComparisons ---------#####

##### Use this script to compare the per-animal, pairwise gene co-expression patterns between experimental groups
##### alimaran@usc.edu, 02-02-2022 


#************* OPTIONAL USER INPUTS ************#

#**NOTE 1: Optional user inputs in this script are as follows:
#*1) Manually define the column indices of genes in the dataset on line 81   
#*2) Define the expression cutoff of a gene that will be used to define a cell-type (lines 160-161)

#**NOTE 2: Lines that start with "If desired" can be skipped                


#*************REQUIRED USER INPUTS ************#

## Import gene-expression-matrix .csv file from SCAMPR:AreaFraction FIJI ImageJ Macro
Raw.Data <-  read.csv(file.choose(),header = TRUE, sep = ",", dec = ".")

## Select names for two experimental groups and identify which animal IDs fall into each group.
Group1name = "CAU" # Choose title for animals that fall into group 1
Group2name = "ELS" # Choose title for animals that fall into group 2

unique(Raw.Data$Animal_ID) # List all Animal IDs in dataset
Group1 = c("18-55", "18-56") # Animal IDs that fall into group 1
Group2 = c("18-43", "18-44") # Animal IDs that fall into group 2

## Select gene pair of interest for pairwise comparison. 
cell_type1_gene <- "Scn1a" # Save cell-type 1 marker gene
sub1_gene2 <- "Pvalb" # Save gene that you want to plot in cell-type 1

## Select a second gene pair for for pairwise comparisons. (optional)
cell_type2_gene <- "Scn10a" # Save cell-type 1 marker gene
sub2_gene2 <- "Ntsr1" # Save gene that you want to plot in cell-type 2

## Plot Settings: Set the smoothing method, correlation method, and the plotting colors for each experimental group.
group1col = "blue" # Choose dot/line color for group 1
group2col = "red" # Choose dot/line color for group 2
smooth.method = "loess" # Alternatives are "lm', 'glm", and "gam".
cor.method = "pearson" # Some alternatives are  "spearman" and "kendall".


#************* REQUIRED PACKAGES ************#

## Install packages, and load libraries

# install.packages("mgsub", "dplyr", "stringr", "ggplot2", "ggpubr", "pacman")

library(mgsub) # For adding group information
library(dplyr) # For calculating group means
library(stringr) # For manipulating column (gene) names
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


#************* GENERATE SUBSETS OF DATA BASED ON GENE PAIRS OF INTEREST ************#

## Save the indices of the columns that DON'T represent your genes.
nonGeneIndices <- 1:(min(geneIndices)-1) 

## Generate new data frames with just your gene pairs of interest
subset1 <- cbind(Normal.Data[nonGeneIndices],
                 cell_type1_gene = Normal.Data[cell_type1_gene],
                 sub1_gene2 = Normal.Data[sub1_gene2])
head(subset1)

subset2 <- cbind(Normal.Data[nonGeneIndices],
                 cell_type2_gene = Normal.Data[cell_type2_gene],
                 sub2_gene2 = Normal.Data[sub2_gene2])
head(subset2)


#************* PLOT LOESS OR LM LINES AND CORRELATION COEFFICIENTS ************#

## Plot smoothed scatterplot!
ggplot(subset1, aes_string(x = cell_type1_gene, y = sub1_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub1_gene2))~"and"~bolditalic(.(cell_type1_gene))), 
                           bold("Coexpression"))))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))

ggplot(subset2, aes_string(x = cell_type2_gene, y = sub2_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub2_gene2))~"and"~bolditalic(.(cell_type2_gene))), 
                           bold("Coexpression"))))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))


#************* GENERATE SUBSETS BASED ON EXPRESSION LEVELS OF GENE 1 OR GENE 2 AND RE-PLOT ************#

## For each subset, discard cells that have < X expression for BOTH genes.  Determine how many cells are left.

X1 <- 0 # Choose normalized gene expression cutoff for subset 1.
X2 <- 0 # Choose normalized gene expression cutoff for subset 2.

subset1_bothGenes_nozeros <- subset1[(subset1[cell_type1_gene] > X1)& (subset1[sub1_gene2] > X1),]
subset2_bothGenes_nozeros <- subset2[(subset2[cell_type2_gene] > X2)& (subset2[sub2_gene2] > X2),]

nrow(subset1)
nrow(subset1_bothGenes_nozeros)
nrow(subset2)
nrow(subset2_bothGenes_nozeros)

## For each subset, discard cells that have < X expression for your cell-type marker gene and any expression level for the other gene.  Determine how many cells are left.
subset1_gene1_nozeros <- subset1[(subset1[cell_type1_gene] > X1),]
subset2_gene1_nozeros <- subset2[(subset2[cell_type2_gene] > X2),]

nrow(subset1)
nrow(subset2_gene1_nozeros)
nrow(subset2)
nrow(subset1_gene1_nozeros)

#************* PLOT LOESS OR LM LINES AND CORRELATION COEFFICIENTS FOR NEW SUBSETS ************#

ggplot(subset1_bothGenes_nozeros, aes_string(x = cell_type1_gene, y = sub1_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub1_gene2))~"and"~bolditalic(.(cell_type1_gene))), bold("Coexpression (Zeros Excluded)"))))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))

ggplot(subset2_bothGenes_nozeros, aes_string(x = cell_type2_gene, y = sub2_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub2_gene2))~"and"~bolditalic(.(cell_type2_gene))), bold("Coexpression (Zeros Excluded)"))))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))

ggplot(subset1_gene1_nozeros, aes_string(x = cell_type1_gene, y = sub1_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub1_gene2))~"and"~bolditalic(.(cell_type1_gene))), 
                           bold("Coexpression in")~bolditalic(.(cell_type1_gene)*"+")~bold("Cells"))))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))

ggplot(subset2_gene1_nozeros, aes_string(x = cell_type2_gene, y = sub2_gene2, color = "Group"))+
  geom_point(alpha = .7, shape = 21, size = 3)+
  geom_smooth(method = smooth.method, size = 1.5 )+
  scale_color_manual(values=c(group1col, group2col))+
  stat_cor(method = cor.method, digits = 3, size = 6)+
  labs(title = bquote(atop(bold(bold(.(Group1name))~"vs"~bold(.(Group2name))~bolditalic(.(sub2_gene2))~"and"~bolditalic(.(cell_type2_gene))), 
                           bold("Coexpression in")~bolditalic(.(cell_type2_gene)*"+")~bold("Cells"))))+  
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        axis.title.x = element_text(size=18, face = "italic"),
        axis.title.y = element_text(size=18, face = "italic"),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=15))

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
