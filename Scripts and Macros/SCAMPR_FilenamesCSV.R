####---------- SCAMPR_FilenamesCSV ---------#####

##### Use this script to generate a blank key CSV file composed of TIFF image names.  This key will be used with the SCAMPR:AreaFraction ImageJ Macro.  
##### alimaran@usc.edu, 01-26-2021 

#************* REQUIRED PACKAGES ************#

## Clear Environment, install packages, and load libraries
rm(list = ls()) 
#install.packages("tcltk", "pacman")

library(tcltk) # For choosing a directory using Windows file explorer
library(pacman) # For loading/unloading packages

#************* GENERATE AND SAVE IMAGE NAMES ************#

## Save the path where your images are located
image.dir <- tk_choose.dir()

## Get list of file names from chosen directory
file.names <- list.files(path = image.dir,pattern = ".tif", full.names = FALSE, recursive = TRUE)
file.names <- strsplit(file.names, split = "/")
file.names <- sapply(file.names, function(x) x[length(file.names[[1]])])

## If desired, remove DAPI and HuCD file names from key.
remove <- c("HuCD", "DAPI")
file.names <- file.names[-grep(paste0(remove, collapse = "|"), file.names)]

## Choose a file name and save location for your key 
save.dir = tk_choose.dir()
key.name = "keyFinal"

## Save file names into a .csv file
write.table(x = cbind(file.names, 0,0,0), 
            file = paste(save.dir,"/",key.name,".csv", sep = ""),
            row.names = FALSE,
            col.names = c("Image Title","RBR","Threshold Minimum", "Average Particle Size"),
            sep = ",")
          

# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear packages
p_unload(all)  # Remove all add-ons


# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
