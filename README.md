# SCAMPR - **S**ingle-**C**ell **A**utomated **M**ultiplex **P**ipeline for **R**NA  [![DOI](https://zenodo.org/badge/443220284.svg)](https://zenodo.org/badge/latestdoi/443220284)

Please cite the following publication if utilizing SCAMPR: 
[Ali Marandi Ghoddousi et al., 2022, Cell Reports Methods 2, 100316 October 24, 2022. https://doi.org/10.1016/j.crmeth.2022.100316](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(22)00206-5#%20) 

LAST UPDATED: 09/26/2022

SCAMPR utilizes ImageJ macros and R scripts for quantification, spatial representation, and comparative analysis of single cell _in situ_ hybridization data.  We recommend first 
running the sample images and ROIs (provided) through the pipeline prior to running your own. 

![image](https://user-images.githubusercontent.com/64667688/154179103-586c8f7d-3bc2-44eb-81bc-4ba6cda1a6f8.png)<br><br><br>

## Required Software and GitHub Instructions

1. [Fiji ImageJ](https://imagej.net/software/fiji/)
3. [R Studio Version 1.4.1717 or greater](https://www.rstudio.com/products/rstudio/)
4. [R Version 1.4.0 or greater](https://www.r-project.org/)
5. [MouseLand Cellpose](https://github.com/MouseLand/cellpose) (Download not required.  Google Colab notebook provided in Step 3)
6. [XQuartz](https://www.xquartz.org/) (Only required for macOS)

The best way to run SCAMPR is to download the entire repository as a ZIP file (_Code --> Dowload ZIP_) so that you have access to all of the scripts and macros as they were written.  Then, follow the instructions in steps 1-4 below to quantify mRNA signal and scroll through the different subheadings of step 5 to run your analysis and spatial mapping of choice.<br><br>

![image](https://user-images.githubusercontent.com/64667688/158484529-261429a2-f207-4eb9-a3e2-8892a6989e4e.png)<br><br>


## 1. Generate flattened, maximum intensity z-projections

Open each microscopy file in ImageJ and generate flattened images:  _Image --> Stacks --> Z Project..._.  Select  _Maximum Intensity_ for projection type.  Split the different channels in the stack into separate images (_Image --> Color --> Split Channels_) and save flattened images as 8-bit (_Image --> Type --> 8-bit_) TIFF files.<br><br><br>
##### &#x2757;&#x2757;_IMPORTANT CONSIDERATIONS_&#x2757;&#x2757;

- These flattened TIFF files are fed into the SCAMPR_AreaFraction.ijm macro in step 4 to quantify mRNA expression.  For this macro to run without throwing errors, the flattened TIFF image files for each tissue section must be stored in separate folders.  Example: The flattened TIFF images from the 3-4 rounds of imaging for the first tissue section on a slide should be stored in one folder and the flattened TIFF images from the 3-4 rounds of imaging for the second tissue section should be stored in a separate folder (see sample data and image below). Everything to the left of the " _ " in the image file name must &#x1F536; EXACTLY &#x1F536; match the name of its containing folder.<br><br><br>
![image](https://user-images.githubusercontent.com/64667688/157580492-5a1b7517-a712-406c-9cef-cebbe5e5d8a8.png)<br><br><br>
- The following portions of the image file names (in red) &#x1F536;**_cannot_**&#x1F536;  be changed:<br><br>

![image](https://user-images.githubusercontent.com/64667688/153970966-b7e56ef6-a5f7-4a39-b9b0-2746af59ecba.png)<br><br><br>

## 2. Register TIFF microscopy files (optional)

If you only performed single-round _in situ_ hybridization (e.g Multiplex RNAscope), you can skip this step.  Conversely, if you performed multi-round _in situ_ hybridization (e.g HiPlex RNAscope), utilize the [ACD BIO Registration Software](https://acdbio.com/store/rnascope-hiplex-image-registration-software.html) or ImageJ to register the TIFF image files from the multiple rounds of imaging.

_**Note:** We have noticed on MacOS machines that the images produced by ACD's image registration software are not directly compatible with Cellpose. To fix this issue, open the registered images in ImageJ and re-save them to a new folder._<br><br>

To register images in ImageJ, use the _Register Virtual Stack Slices_ and _Transform Virtual Stack Slices_ plugins. The _Register Virtual Stack Slices_ plugin will use one of your DAPI images as an anchor and will register the DAPI images from the other rounds of imaging onto the anchor image.  It will also generate a transformation file for each DAPI registration that can then be applied to the other non-DAPI images for that round.

The following protocol outlines how to do this for one tissue section.  A completed example is provided in the _Sample_Registration_Completed_ folder to help organize your folders and subfolders.  In addition, sample images are provided in the _Sample_Registration_Practice_ folder and can be run through the below pipeline prior to running your own data.  The contents of the _Sample_Registration_Practice_ folder after it has gone through the pipeline should match the contents of the _Sample_Registration_Practice_ folder:
1. Create an Input_Images folder. Create a subfolder corresponding to the name of the tissue section (see "Important Considerations" in step 1 for naming conventions). Place all flattened TIFF images into this subfolder.  Create a seperate subfolder with just the DAPI images from all rounds.
2. Create an Ouput_Images folder.  Create a subfolder corresponding to the name of the tissue section (see "Important Considerations" in step 1 for naming conventions). This is where your final registered images will be saved.
3. Create a Tranformation_files folder. Create a subfolder corresponding to the name of the tissue section (see "Important Considerations" in step 1 for naming conventions). This is where your transformation files from the registration will be saved, copied, and renamed.
4. Run the _Register Virtual Stack Slices_ plugin (_Plugins --> Registration_)
   - Select folder with all of your DAPI input images
   - Select folder where you want your registered DAPI images to be saved
   - Choose "Affine" for Feature extraction model and "Elastic" for Registration Model.
   - Select "Save Transforms" and click OK.
   - Select folder where transformation matrices will be saved.
   - Select one of your DAPI input images.  All other images will be registered to this DAPI image.
   - An image stack will pop up with all of your registered DAPI images.  Scroll through to ensure that registration worked properly.
5. Place each transformation file into its appropriate subfolder.  **_Make copies_** of each transformation file in each subfolder until you have the same number of transformation files as you do images for that round.  For example, if you have 3 mRNA images that correspond to the T0-dapi image, make sure that you make 3 copies of the T0-dapi transformation file (see image below).
6. Rename each copy of the transformation files to match the image name that it corresponds to.  For example, if the T0-dapi file corresponds to the T1-met, T2-met, and T3-rorb mRNA channels, rename the three copies of the T0-dapi transformation file to T1-met, T2-met, and T3-rorb (see image below).<br><br>  
<br><br>
![image](https://github.com/ramin-ali-marandi-ghoddousi/SCAMPR/assets/64667688/54e57ddb-7de2-4828-ada7-37dbcd134200)<br><br>

7. Run the _Transform Virtual Stack Slices_ plugin (_Plugins --> Transform_)
   -  Choose folder containing input images for imaging round.
   -  Choose output folder for registered images.
   -  Choose folder containing copied and renamed transformation files.
   -  Check _interpolate_ and run.
8. Crop the image stack that pops up by using the box tool and the dapi images to generate a crop box and typing the keyboard shortcut (Ctrl+Shift+X).  Convert the stack into individual images and save.
9. If necessary, change each individual image titles to match the required naming conventions.  See "Important Considerations" in Step 1 for file and folder naming conventions.<br><br><br>



## 3. Generate Cellular ROIs

Using the HuC/D or DAPI images that correspond to each sample, run [cellpose](https://github.com/MouseLand/cellpose) locally on your device or by utilizing a modified google colab notebook to segment each image into individual cells.  The colab notbeook is written by Pradeep Rajasekhar from the Monash Institute of Pharmaceutical Sciences and was modified for simplicity and to save all masks as a zip file. [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1Su8AKlmJoCvX4fa7BYFlSMAqPTex1Oku#scrollTo=Nb90LCrotIx4)

The output will be cellular mask files that can then be converted into ImageJ ROIs using cellpose's [imagej_roi_converter.py script](https://github.com/MouseLand/cellpose/blob/master/imagej_roi_converter.py).  We have included a copy of this script in the SCAMPR "Scripts and Macros" folder (_Cellpose_Outline_to_ROI_Converter.py_).

To convert Cellpose outputs to ROIs, open the HuCD TIFF file as well as the _Cellpose_Outline_to_ROI_Converter.py_ script in ImageJ.  Click "Run" on the macro window, navigate to the .zip folder containing your Cellpose outputs, and select the outline.txt output file corresponding to the appropriate HuCD TIFF image file.  Select all ROIs in the ROI Manager and save. &#x2757;**_Ensure that each ROI folder matches the corresponding image title (see sample data and image below)._**

![image](https://user-images.githubusercontent.com/64667688/154147059-49c96f2a-ba9c-44aa-b73d-7608e8d1bba9.png)<br><br><br>


## 4. Quantification of Gene Expression in ImageJ

In this step, a key file will be generated to help calculate the expression of each gene in each cellular ROI.  Gene expression is calculated as an area-fraction (the percent of a cellular ROI that contains fluorescent signal).  These expression values can be used as-is during the data visualization and analyses in Step 4, or can be converted into a value representing the number of pixels or ums that are positive for a fluorescent signal.  This conversion will take place in the R scripts used for visualization and analysis.

_Note:  We highly recommend that the user run the macro on the Sample Data first to gain familiarity with the process before proceeding with their own datasest._

#### Required Data

1. Flattened TIFF image files (from Step 1 or 2).
2. ROIs for each TIFF image file (from Step 3).
3. CSV Key file.
    
    This Key file can be generated in two ways: 1) Manually, 2) or in a semiautomated manner in ImageJ.  Manual generation of the key can lead to user-related bias, which can be mitigated with proper blinding.  The semiautomated method is less prone to user error but may be less accurate for low-intensity images in some datasets.  Descriptions on how to generate the key for both methods are outlined below:
    <br>
    
   ----------------------- Semiautomated Key Generation and Gene Expression Quantification -----------------------
    
   To generate this key, run SCAMPR_AreaFraction_SemiAuto.ijm. &#x2757;_**Make sure to remove all non-mRNA images (DAPI,HuCD,etc) from the Image folder prior to running the macro**_&#x2757;
   
   There are some options that can be changed and sections that require user input:
   
   1. The default rbr for processing is set to 1 and can be changed in line 10 of the code ("ExerRBR = ").
   2. The image with the highest average intensity will be identified (Representative Image).  The script will ask the user to choose the location of that image.
   3. Two copies of Representative Image will be opened by the script.  The one is white is the ground truth, and the one in red is the image that the user will used to adjust and identify the ideal threshold.  Place both images side-by-side and use the zoom tool to zoom in ALL THE WAY onto one cell in both images.  Then use the threshold adjustment bar to find the ideal threshold.  Remember this number. Click OK.
   ![image](https://user-images.githubusercontent.com/64667688/179127400-1be7261b-f311-4a66-80ce-d35df9bd05bb.png)

   4. Input the threshold identified in the prior step into the box.
   5. SCAMPR quantifies gene expression as the percentage of pixels in a cell that are positive for signal.  This is later converted to the number of pixels in each cell that are positive for signal.  SCAMPR also allows for the option of estimating the numper of particles in each cell by dividing the number of pixels by an average particle size. The default average particle size is set to 1 for ALL images by default.  If you want to properly estimate the number of particles in each cell, find the Key.csv file, open each individual image in a **NEW** imageJ window, remove background (_Process > Subtract_ Background, rbr = 1) and threshold (_Ctrl + Shift + T_) using the corrected threshold value in the key, zoom in on 2-3 random "spots", average the number of pixels comprising the spots, and input this number into the Key.csv file. You can average calculate the number of pixels in a spot visually, or by circling the spot using the freehand selection tool and going to _Analyze > Measure_ or hitting "m" on your keyboard. Click OK in the macro pop-up window to resume with the anlysis.
   6. After the script runs, to ensure that proper threshold values were generated, test the values from the Key.csv file on some of the low signal intensity images towards the top of the file.  Make sure to remove background prior to checking thresholding values.
   
   No further step is required as this script will generate the Key, then use that key to quantify signal for each ROI in each image.  The output will be two csv files with %Area or particle count gene expression values for each gene in each cell.
   <br>
    
   ----------------------- Manual Key Generation and Gene Expression Quantification -----------------------
    
    To generate this key, run the _SCAMPR_FilenamesCSV.R_ script in R Studio to generate a template for the key in CSV format.  Fill in the template by opening each TIFF file in ImageJ and noting:
    1. The ideal Rolling-ball-radius (rbr): _Process > Subtract_ Background. Choose a pixel value and hit check _Preview_.  Click _OK_ once you have your value.  _**Note:** We have found that a low rbr of 1 works well across all images.  Once you have selected an rbr for a few of the images the same value can be filled for the others_
    2. The minimum threshold value: _Ctrl + Shift + T_. Check _Dark background_ and use the lower adjustement bar to change value. Write down ideal minimum threshold value and close image without applying or saving.
    3. Average particle size (optional): _Analyze > Set Measurements_ -- Check Area and Limit to Threshold. Draw a circle around an isolated mRNA particle using the Oval                selection tool and type 'M' on your keyboard.  The area will give you the number of pixels with signal.

   
##### &#x2757;&#x2757;_IMPORTANT CONSIDERATIONS FOR MANUAL KEY GENERATION_&#x2757;&#x2757;
The names of each TIFF image file should &#x1F536; EXACTLY &#x1F536; match the names in the Image Title column in the CSV Key file.  Using the _SCAMPR_FilenamesCSV.R_ to pre-populate the CSV Key file should help ensure that this is the case.
Once you have the flattened TIFF image files, the ROIs for each image file, and the CSV key file, drag the SCAMPR_AreaFraction.ijm (for manually generated Key)  into FIJI ImageJ and run it to calculate gene expression.

#### Required Code:

1. SCAMPR_FilenamesCSV.R
2. SCAMPR_AreaFraction.ijm or SCAMPR_AreaFraction_SemiAuto.ijm


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
### Hierarchical Clustering and Cluster Spatial Topomaps

Perform hierarchical clustering on data, compare gene expression and cell size across different clusters, and map clusters back onto the tissue.

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
| ![image](https://user-images.githubusercontent.com/64667688/154170184-cc3ee472-3b42-49f2-b519-bd056605078b.png) | ![image](https://user-images.githubusercontent.com/64667688/154169047-9d29d1f6-3ba8-4104-beea-28bb3571e5c3.png) | ![image](https://user-images.githubusercontent.com/64667688/159094409-088efbeb-f709-4627-a351-9e1e3a62c7fd.png) |



