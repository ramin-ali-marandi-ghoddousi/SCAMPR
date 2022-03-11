# SCAMPR - **S**ingle-**C**ell **A**utomated **M**ultiplex **P**ipeline for **R**NA

SCAMPR utilizes ImageJ macros and R scripts for quantification, spatial representation, and comparative analysis of single cell _in situ_ hybridization data.  We recommend first 
running the sample images and ROIs (provided) through the pipeline prior to running your own. 

![image](https://user-images.githubusercontent.com/64667688/154179103-586c8f7d-3bc2-44eb-81bc-4ba6cda1a6f8.png)<br><br><br>

## Required Software

1. [Fiji ImageJ](https://imagej.net/software/fiji/) and plugins:
   - SCF plugin -- To convert cellpose masks to ImageJ ROIs. To dowload plugin in ImageJ:  
      1. Select _Help --> Update... --> Manage Update Sites_.  
      2. Check the box next to SCF MPI CBG.
      3. Apply changes and restart ImageJ.
3. [R Studio Version 1.4.1717 or greater](https://www.rstudio.com/products/rstudio/)
4. [R Version 1.4.0 or greater](https://www.r-project.org/)
5. [MouseLand Cellpose](https://github.com/MouseLand/cellpose) (Download not required.  Google Colab notebook provided in Step 3)
6. [XQuartz](https://www.xquartz.org/) (Only required for macOS)<br><br>

## 1. Generate flattened, maximun intensity z-projections

Open each microscopy file in ImageJ and generate flattened image:  _Image > Stacks > Z Project..._.  Select  _Maximum Intensity_ for projection type.  Save flattened image as TIFF file.<br><br><br>
##### &#x2757;&#x2757;_IMPORTANT CONSIDERATIONS_&#x2757;&#x2757;

- These flattened TIFF files are fed into the SCAMPR_AreaFraction.ijm macro in step 4 to quantify mRNA expression.  For this macro to run without throwing errors, the flattened TIFF image files for each tissue section must be stored in seperate folders.  Example: The flattened TIFF images from the 3-4 rounds of imaging for the first tissue section on a slide should be stored in one folder and the flattened TIFF images from the 3-4 rounds of imaging for the second tissue section should be stored in a seperate folder (see sample data and image below). Everything to the left of the " _ " in the image file name must &#x1F536; EXACTLY &#x1F536; match the name of its containing folder.<br><br><br>
![image](https://user-images.githubusercontent.com/64667688/157580492-5a1b7517-a712-406c-9cef-cebbe5e5d8a8.png)<br><br><br>
- The following portions of the image file names (in red) &#x1F536;**_cannot_**&#x1F536;  be changed:<br><br>

![image](https://user-images.githubusercontent.com/64667688/153970966-b7e56ef6-a5f7-4a39-b9b0-2746af59ecba.png)<br><br><br>

## 2. Register TIFF microscopy files (optional)

If you only performed single-round _in situ_ hybridization (e.g Multiplex RNAscope), you can skip this step.  Conversely, if you performed multi-round _in situ_ hybridization (e.g HiPlex RNAscope), utilize the [ACD BIO Registration Software](https://acdbio.com/store/rnascope-hiplex-image-registration-software.html) or ImageJ to register the TIFF image files from the multiple rounds of imaging.<br><br><br>

To register images in ImageJ, use the _Register Virtual Stack Slices_ and _Transform Virtual Stack Slices_ (_Plugins --> Transform_) plugins. The following protocol outlines how to do this for one tissue section (unregistered sample images are provided in the "Registration Sample" folder):
1. Place all flattened TIFF images into a folder named Input_Images and organize them into seperate subfolders folders in the following manner:
   -  Place all DAPI images from each microscopy round into one folder.
   -  Place all other mRNA and HuCD images from all other microscopy rounds into folders corresponding to their microscopy round (i.e genes 1/2/3/4 from round 1 into a folder named Round1, genes 5/6/7/8 from microscopy round 2 into a folder named Round2, etc.)
2. Create an Ouput_Images folder and generate subfolders where the registered DAPI and round-specific mRNA/HuCD images will be saved (i.e DAPI_Images, Round1, Round2, Round3, etc).
3. Create a Tranformation_files folder and generate round-specific subforlders. 
4. Run the _Register Virtual Stack Slices_ plugin (_Plugins --> Registration_)



## 3. Generate Cellular ROIs

Using the HuC/D or DAPI images that correspond to each sample, run [cellpose](https://github.com/MouseLand/cellpose) locally on your device or by utilizing a modified google colab notebook to segment each image into individual cells.  The colab notbeook is written by Pradeep Rajasekhar from the Monash Institute of Pharmaceutical Sciences and was modified for simplicity and to save all masks as a zip file. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/11m15HUl_StmZUiZjQ2QiaEKaToWLYKn2#scrollTo=Drx_6GbEHIOj)  

The output will be cellular mask files that can then be converted into ImageJ ROIs using the SCF plugin in ImageJ.

To convert masks to ROIs, open a mask file in ImageJ and select _SCF > Segmentation > LabelMap to ROI Manager (2D)_ to get ROIs.  Select all ROIs in the ROI Manager and save. &#x2757;**_Ensure that each ROI folder matches the corresponding image title (see sample data and image below)._**

![image](https://user-images.githubusercontent.com/64667688/154147059-49c96f2a-ba9c-44aa-b73d-7608e8d1bba9.png)<br><br><br>


## 4. Quantification of Gene Expression in ImageJ

In this step, a key file will be generated to help calculate the expression of each gene in each cellular ROI.  Gene expression is calculated as an area-fraction (the percent of a cellular ROI that contains fluorescent signal).  These expression values can be used as-is during the data visualization and analyses in Step 4, or can be converted into a value representing the number of pixels or ums that are positive for a fluorescent signal.  This conversion will take place in the R scripts used for visualization and analysis.

#### Required Data

1. Flattened TIFF image files (from Step 1 or 2).
2. ROIs for each TIFF image file (from Step 3).
3. CSV Key file.
    - To generate this key, run the _SCAMPR_FilenamesCSV.R_ script to generate a Key template in CSV format.  Fill in the template by opening each TIFF file in ImageJ and noting:
    1. The ideal Rolling-ball-radius: _Process > Subtract_ Background. Hit Okay.
    2. The upper threshold value: _Ctrl + Shift + T_. Use lower adjustement bar to change value.
    3. Average particle size (optional): _Analyze > Set Measurements_ -- Check Area and Limit to Threshold. Draw a circle around an isolated mRNA particle using the Oval                selection tool and type 'M' on your keyboard.  The area will give you the number of pixels with signal.

##### &#x2757;&#x2757;_IMPORTANT CONSIDERATIONS_&#x2757;&#x2757;
The names of each TIFF image file should &#x1F536; EXACTLY &#x1F536; match the names in the Image Title column in the CSV Key file.  Using the _SCAMPR_FilenamesCSV.R_ to pre-populate the CSV Key file should help ensure that this is the case.


#### Required Code:

1. SCAMPR_FilenamesCSV.R
2. SCAMPR_AreaFraction.ijm


#### Outputs:

|Key.csv template generated from SCAMPR_FilenamesCSV.R (replace 0’s with real values)| Gene X Cell count matrix file generated from SCAMPR_AreaFraction.ijm|
|  --- | ---  |
| ![image](https://user-images.githubusercontent.com/64667688/153974233-da44ba8b-849d-4ad5-9349-8174d70dea51.png) | ![image](https://user-images.githubusercontent.com/64667688/153974260-4e6878dd-d3fa-47f6-8a54-767485b760c8.png) |

<br><br>
## 5. Visualization and Analysis

The R scripts in this section can be used to visualize the level of gene expression in specific samples, the level of gene co-expression in these samples, and can be used to perform gene expression comparisions between two experimental groups on a broad and cell-type specific manner.   Example outputs are provided for each R script below using the sample data.  

### Violin Plots and Expression Spatial Topomaps

Generate normalized gene expression violin plots and map these gene expression values back onto the tissue.

#### Required Data

1. Gene X Cell Count Matrix from Step 4
2. ROIs from Step 3

#### Required Code:

1. SCAMPR_ViolinPlots.R
2. SCAMPR_ExpressionTopomaps.R

#### Outputs:

|SCAMPR_ViolinPlots.R| SCAMPR_ExpressionTopomaps.R| 
|  --- | ---  |
| ![image](https://user-images.githubusercontent.com/64667688/153976379-5d37f1a0-0fe6-433a-80b3-33c367b859a7.png) | ![image](https://user-images.githubusercontent.com/64667688/153976355-e65d5950-e214-436b-9071-d90723b09d97.png) |

<br><br>
### Heirarchical Clustering and Cluster Spatial Topomaps

Perform heirarchical clustering on data, compare gene expression and cell size across different clustering, and map clusters back onto the tissue.

#### Required Data

1. Gene X Cell Count Matrix from Step 4
2. ROIs from Step 3

#### Required Code:

1. SCAMPR_hClustHeatmapTopomaps.R
2. SCAMPR_hClustViolinPlotsCellAreas.R

#### Outputs:

|SCAMPR_hClustHeatmapTopomaps.R | SCAMPR_hClustViolinPlotsCellAreas.R  | 
|  --- | --- |
| ![image](https://user-images.githubusercontent.com/64667688/154167325-3287b5c6-7a53-4f80-aa93-78d26c804bfa.png) | ![image](https://user-images.githubusercontent.com/64667688/154167175-bcea1e2a-5e75-4943-a1fd-ee39490c91dc.png) |

<br><br>
### Co-expression Correlation Plots and Co-expression Spatial Topomaps

Plot pairwise co-expressoin circle plots, co-expression scatterplots with smoothing lines, and map co-expression patterns back onto tissue.

#### Required Data

1. Gene X Cell Count Matrix from Step 4
2. ROIs from Step 3

#### Required Code:

1. SCAMPR_CorrelationCirclePlots.R
2. SCAMPR_LoessCoexpressionTopomaps.R

#### Outputs:

|SCAMPR_CorrelationCirclePlots.R | SCAMPR_LoessCoexpressionTopomaps.R  | 
|  --- | --- |
| ![image](https://user-images.githubusercontent.com/64667688/153977754-f0c75825-a8b8-4557-9cdb-15c137c5d998.png) | ![image](https://user-images.githubusercontent.com/64667688/153977445-316b70ae-a94a-4039-9d16-bc3bdb6cf70f.png) |

<br><br>
### Comparative (between-group) Analysis

Compare pairwise gene co-expression patterns between two experimental groups.  You can compare pairwise correlation, mean gene expression, and the percent of cells expressing each gene.

#### Required Data

1. Gene X Cell Count Matrix from Step 4

#### Required Code:

1. SCAMPR_CorrelationRegressionComparisonsMatrix.R
2. SCAMPR_MeanComparisonsMatrix.R
3. SCAMPR_PercentCellsComparisonsMatrix.R

#### Outputs:

| SCAMPR_CorrelationRegressionComparisonsMatrix.R | SCAMPR_MeanComparisonsMatrix.R  |  SCAMPR_MeanComparisonsMatrix.R |
|  --- |  ---  |  ---  |
| ![image](https://user-images.githubusercontent.com/64667688/154169686-020cd78e-1ad2-44d6-9d7c-468dc172d219.png) | ![image](https://user-images.githubusercontent.com/64667688/154169702-a34404c7-61de-4000-8676-34c565f8b4ba.png) | ![image](https://user-images.githubusercontent.com/64667688/154169720-0fb05f72-41ea-4c52-a8f4-652db1750476.png) |

<br><br>
### Cell-Type Specific Comparative (between-group) Analysis

Generate cell-type specific subsets of data, then compare pairwise correlation, mean gene expression, and the percent of cells expressing each gene in each cell type between experimental groups.

#### Required Data

1. Gene X Cell Count Matrix from Step 4

#### Required Code:

1. SCAMPR_LoessScatterplotComparisons.R
2. SCAMPR_MeanComparisonsCellTypes.R
3. SCAMPR_PercentCellsComparisonsCellTypes.R

#### Outputs:

| SCAMPR_LoessScatterplotComparisons.R | SCAMPR_MeanComparisonsCellTypes.R  |  SCAMPR_PercentCellsComparisonsCellTypes.R |
|  --- |  ---  |  ---  |
| ![image](https://user-images.githubusercontent.com/64667688/154170184-cc3ee472-3b42-49f2-b519-bd056605078b.png) | ![image](https://user-images.githubusercontent.com/64667688/154169047-9d29d1f6-3ba8-4104-beea-28bb3571e5c3.png) | ![image](https://user-images.githubusercontent.com/64667688/154172852-a16a13ee-8e06-45a1-8325-5bbbf72c151b.png) |



