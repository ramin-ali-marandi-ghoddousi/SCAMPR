/******* MAX LOG-MEANS SEMIAUTOMATED THRESHOLDING METHOD *******/
setBatchMode(true); //	Runs macro without opening any image windows.

// Version 8: 
//// 1. Uses rolling-ball-radius of 5 for background removal to generate the image intensity histograms. Can change this value in line 9.
//// 2. Uses rolling-ball-radius of 1 for background removal when quantifying signal. Can change this value in line 10.
//// 3. Uses the image with the highest LOGMEANS intensity to find "Representative Image".

HistRBR = 5;
ExperRBR = 1;

//	************************************************************* UPLOAD DATA AND SET DIRECTORIES ***********************************************************

showMessage("Choose Image Folder.","Choose folder containing all of the IMAGE subfolders.");
path_imagefolder=getDirectory("Choose a data folder"); // Variable that stores path to the folder containing image subfolders.
subDir = getFileList(path_imagefolder); // Variable that stores the names of each individual image subfolder.

showMessage("Choose folder containing ROIs.");
path_roifolder = getDirectory("");	//	Variable that stores path to folder containing ROIs.

showMessage("Choose folder to save RESULTS files.");
path_processedfiles = getDirectory("");	//	Variable that stores path to results folder where files will be saved.

showMessage("Choose folder to save KEY file with image-specific thresholds.");
path_keyfile = getDirectory(""); //	Variable that stores path to the Key.csv file containing the thresholds


// *********************************************** GENERATE KEY WITH OPTIMAL THRESHOLD VALUES FOR EACH IMAGE *************************************************

// Peform background subtraction and calculate mean intensity for each image. Rolling ball radius for this step is set in line 9.
run("Set Measurements...", "mean display redirect=None decimal=3");

for(i = 0; i < subDir.length; i++){
	files = getFileList(path_imagefolder + subDir[i]);
	for(j = 0; j < files.length; j++){
		open(path_imagefolder + subDir[i] + files[j]);
		run("Subtract Background...", "sliding disable rolling="+HistRBR);
		run("Measure");
		close();
	}
}

// Add Image Names to the results table

ImageFiles = newArray();
for (i=0; i<nResults; i++) {
   files= getResultLabel(i);
   ImageFiles = Array.concat(ImageFiles,files);
}

IJ.renameResults("Results");
for (i = 0; i < nResults; i++) {
    imageNames = ImageFiles[i];
    setResult("Image Title", i, imageNames);
}
updateResults();

// Remove "Label" column

Image_Title = Table.getColumn("Image Title")
Mean = Table.getColumn("Mean");

IJ.renameResults("Results");
close("*");

Table.create("Results");
Table.setColumn("Image Title", Image_Title);
Table.setColumn("Mean", Mean);
//Table.update;

//Reorder the columns in the table
IJ.renameResults("Results");
for (i = 0; i < nResults; i++) {
    names = getResult("Mean", i);
    setResult("Means", i, names);
}
Table.deleteColumn("Mean");
updateResults();

IJ.renameResults("Results");
for (i = 0; i < nResults; i++) {
    names = getResult("Means", i);
    setResult("Mean", i, names);
}
Table.deleteColumn("Means");
updateResults();

// Calculate the log of the mean intensity for each image and add it as a logMeans column in results table
IJ.renameResults("Results"); // otherwise below does not work...
for (row=0; row<nResults; row++) {
	Mean = getResult("Mean", row);
    setResult("logMeans", row, log(Mean+1));
}
updateResults();

// Sort images in results table by their LogMean intensity and calculate the max of the logMean intensities.
sortedTable = Table.sort("logMeans");
logMeans = Table.getColumn("logMeans");

logMeans = Array.sort(logMeans);

maxValue = logMeans[logMeans.length-1];

for (i=0; i<logMeans.length; i++){
	
	if(logMeans[i] == maxValue){
	index = i+1;
	}
}


// Divide the logMean intensity of each image in the results table by the max logMean intensity to generate a correction factor column in the results table.
IJ.renameResults("Results"); // otherwise below does not work...
for (row=0; row<nResults; row++) {
	logMean = getResult("logMeans", row);
    setResult("CorrectionFactor", row, logMean/maxValue);
}
updateResults();

//Array.print(logMeans);
//print("Max_Index =", index);
print("Max =",maxValue);

// Identify the image that has the highest logMean intensity (Representative Image).
RepImage = getResultString("Image Title", index-1);
print("Representative Image=", RepImage);

// Open Representative Image, remove background, and threshold.
showMessage("Choose folder.","Choose folder containing the Representative Image:" + RepImage);
RepImagePath=getDirectory("Choose a data folder"); //Get Rep Image directory

setBatchMode(false);
open(RepImagePath + RepImage);
selectWindow(RepImage);
run("Enhance Contrast", "saturated=0.35");
run("Brightness/Contrast...");
run("Subtract Background...", "rolling=1 sliding disable");
open(RepImagePath + RepImage);
run("Subtract Background...", "rolling=1 sliding disable");
run("Threshold...");
setAutoThreshold("RenyiEntropy dark");
call("ij.plugin.frame.ThresholdAdjuster.setMode", "Red");

// Manually determine threshold for representative image
waitForUser("1. Place the ground truth image (white) and the image that you will use to determine the optimal threshold (red) side-by-side. \n                                           2. Zoom in on the same a single cell in both images.\n                           3. Use upper slide bar to choose optimal threshold value.  Input value in next screen.");


close("*");
Dialog.create("Choose Threshold");
Dialog.addNumber("Optimal threshold value", 0);
Dialog.show();
Rep_Threshold = Dialog.getNumber()
print(Rep_Threshold);

setBatchMode(true);

// Add a corrected threshold value for each image by multiplying the representative threshold by the image-specific correction factor
for (row=0; row<nResults; row++) {
	CorrectedMean = getResult("CorrectionFactor", row);
    setResult("Corrected Threshold", row, round(CorrectedMean * Rep_Threshold));
}
updateResults();

// Add a default average particle size to the results table.
partSize = 1;
for (row=0; row<nResults; row++) {
    setResult("Average Particle Size", row, partSize);
}
updateResults();

// Save the results table as a "Key.csv" file.
saveAs("Results", path_keyfile + "Key.csv");
selectWindow("Results"); 
run("Close");

// Add image-specific average particle size if needed
waitForUser("                                                The average particle size for each image is set to 1.\n                                                                        ------------------------                                                                        \n                                   Do you want to manually specify an average particle size for each image? \n                                                                        ------------------------                                                \nIf yes, locate Key.csv file, note the image-specific particle size, save the Key.csv file. Click 'OK' to continue");



// **************************************************  QUANTIFY SIGNAL IN EACH IMAGE **********************************************************

// Process each image using the thresholds that were determined using the semi-automated method that were saved in the Key.csv file.  
processNestedFiles(path_imagefolder);	//	Pass the main image folder containing image subfolders through the processNestedFiles function (see Lines #44 - #67).

showMessage("Almost done!  Creating summary .csv files.  Click OK to proceed.");
createSummary(path_processedfiles);	//	Pass the main results folder containing results subfolders through the createSummary function (see Lines #143 - #222).

combineCSVs(path_processedfiles);	//	Create a final .csv file for combined Area Fraction Results and combined Area Division Results (see lines #225 - #253).

showMessage("                                                     Macro finished. Results are saved in the designated folder. :) \nPlease open a random subset of images and manually verify the semi-automated threshold values in the Key.csv file.")	//	Open dialog box that indicates to user that the macro has finished running.


//	********************************************************************* FUNCTIONS ****************************************************************************

function checkKeyValue(header) {	//	Creates a function that checks the .csv key for a specific average particle size, RBR, and threshold minimum.

	open(path_keyfile+ "Key.csv");
	for (f = 0; f < Table.size; f++) {
		key_row = Table.getString("Image Title", f);
		if(key_row == tif_title) {
			key_value = Table.get(header, f);
		}
	}
	close("*.csv");
	return key_value;
}


function processNestedFiles(dir) {	//	Creates a recursive processNestedFiles function that the main image folder path can be passed through. If an image subfolder is returned, then all image files within that subfolder are passed through recursively and encounter the "else" commands.
	
	files_currentdir = getFileList(dir);	//	Variable that stores an array of the folders/files in the current directory that is passed through the function.
	
	for (a = 0; a < files_currentdir.length; a++) {	
			
		if (endsWith(files_currentdir[a], "/")) {
			path_resultsdir = path_processedfiles + files_currentdir[a];	//	Variable that stores the path of a new results folder associated with an image subfolder.
			File.makeDirectory(path_resultsdir);	//	Create a new results subfolder within the main processed files folder.
			processNestedFiles("" + dir + files_currentdir[a]);	//	For each image subfolder, pass the files within it through the processNestedFiles function.
			
			path_logfile = path_resultsdir + "Analysis Variables Log.txt";
			selectWindow("Log");
			saveAs("text", path_logfile);
			print("\\Clear");
		}
		else {
			if (!matches(files_currentdir[a], ".*T0.*")) {
				path_tifimage = dir + files_currentdir[a];
				processImageFile(path_tifimage);	//	For each .tif file in this subfolder path, call the function processImageFile (see Lines # - #).
			}
		}
	}
}


function processImageFile(path_tifimage) {	//	Creates a processImageFile function that processes all .tif image paths that are passed through it.
	
	if (endsWith(path_tifimage, ".tif")) {

		//	Open the .tif image.
		open(path_tifimage);
		tif_title = getTitle();
		tif_titlearray = split(tif_title, "_");
		tif_titleshort = substring(tif_title,0,indexOf(tif_title,".tif"));
		tif_target = substring(tif_titlearray[1],0,indexOf(tif_titlearray[1],".tif"));
		current_subfolder = tif_titlearray[0];

		//	Set desired output measurements - ROI area, ROI location (centroid), and ROI area fraction.
		run("Set Measurements...", "area centroid area_fraction redirect=None decimal=3");

		//	Subtract background from .tif file using rolling-ball radius set in .csv key.
		radius = ExperRBR;
		run("Subtract Background...", "rolling=radius");

		//	Apply thresholding on .tif image using the lower (min) threshold value set in .csv key		
		min_threshold = checkKeyValue("Corrected Threshold");
		max_threshold = 255;
		setThreshold(min_threshold, 255);	//	Set lower threshold to value set in .csv key
		setOption("BlackBackground", true);	//	Assume image has a black background.
		run("Convert to Mask");	//	Apply threshold parameters and create a binary image.

		//	Open ROI Manager and show all ROIs on .tif image.
		path_roifile = path_roifolder + "ROI_" + current_subfolder + ".zip";
		run("ROI Manager...");
		roiManager("open", path_roifile);
		roiManager("Show All without labels");
		roiCount = roiManager("count");	//	Variable that stores total # of ROIs open in ROI Manager.
		
		//	Measure ROI area, ROI x- & y-location (centroid), and ROI area fraction.
		roiManager("select", Array.getSequence(roiCount));
		roiManager("multi-measure measure_all append");

		//	Add (Thresholded Area / Particle Area) measurements into results table.
		TAdivPA_results = newArray(nResults);
		partarea = checkKeyValue("Average Particle Size");
		for (e = 0; e < nResults; e++) {
			area_result = getResult("Area", e);	//	Area of one cell.
			areafraction_percent = getResult("%Area", e);	//	Area fraction (%) in one cell.
			areafraction_proportion = areafraction_percent / 100;	//	Area fraction (proportion) in one cell.
			TAdivPA_results[e] = (area_result * areafraction_proportion) / partarea;	//	Thresholded area divided by particle area in .csv key.
			setResult("(ThresArea / PartArea) Count", e, TAdivPA_results[e]);	//	Creates a new column in results table and populates it with (Thresholded Area / Particle Area) results.
		}

		//	Save measurements as .csv file.
		results_filename = "Results_" + tif_titleshort + ".csv";	// Variable that stores complete results file name.
		path_results = path_resultsdir + results_filename;	//	Variable that stores path to results file associated with this .tif image.
		saveAs("results", path_results);

		//	Print all variables associated with current .tif image into log.
		print("");
		print("Image: " + tif_titlearray[0]);
		print("Target: " + tif_target);
		print("Average Particle Area: " + partarea);
		print("Rolling-Ball Radius: " + radius);
		print("Lower Threshold (Minimum): " + min_threshold);
		print("Upper Threshold (Maximum): " + max_threshold);
		print("");

		//	Clear results, clears ROI manager, and closes .tif image.
		run("Clear Results");
		roiManager("reset");
		selectWindow(tif_title);
		close();
		
	}
}


function createSummary(dir) {	//	

	files_currentdir = getFileList(dir);	//	Variable that stores an array of the folders/files in the current directory that is passed through the function.
	variables_array = newArray("%Area", "(ThresArea / PartArea) Count");
	variables_filenames = newArray("Area-Fraction-Summary_", "Area-Division-Count-Summary_");
	variables_foldernames = newArray("Area Fraction Summaries", "Area Division Summaries");
	
	for (l = 0; l < variables_foldernames.length; l++) {
		path_summaryresults = path_processedfiles + variables_foldernames[l];
		File.makeDirectory(path_summaryresults);
	}
				
	for (b = 0; b < files_currentdir.length; b++) {

		if (endsWith(files_currentdir[b], "/") && !endsWith(files_currentdir[b], "Summaries") && !startsWith(files_currentdir[b], "Cyto")) {
			current_imageID = File.getName(files_currentdir[b]);
			for (k = 0; k < variables_array.length; k++) {
				createSummary("" + dir + files_currentdir[b]);	//	For each results subfolder, pass the files within it through the createSummary function.
				results_filename = variables_filenames[k] + current_imageID + ".csv";
				path_resultsdir = path_processedfiles + variables_foldernames[k] + "/";
				path_results = path_resultsdir + results_filename;
				saveAs("Results", path_results);
				run("Clear Results");
			}
		}
		else {
			if (startsWith(files_currentdir[b], "Results_")) {
				path_resultsfile = dir + files_currentdir[b];

				// Open the .csv result file associated with one gene target.
				open(path_resultsfile);
				results_title = getInfo("window.title");
				results_titlearray = split(results_title, "_");
				results_titlearray = substring(results_titlearray[2],0);
				results_titlearray = split(results_titlearray, "-");
				results_target = substring(results_titlearray[1],0,indexOf(results_titlearray[1],".csv"));


				// Save Image ID , Animal ID, and ROI Area columns into the results.
					imageID_array = newArray(Table.size);
					animalID_array = newArray(Table.size);
					ROIarea_array = newArray(Table.size);
					
					for (q = 0; q < Table.size; q++) {
						current_animalID = current_imageID;
						current_animalID = split(current_animalID, "~");
						current_animalID = substring(current_animalID[0],0);
						animalID_array[q] = current_animalID;
					}
					for (p = 0; p < Table.size; p++) {
						setResult("Animal_ID", p, animalID_array[p]);
					}
					for (g = 0; g < Table.size; g++) {					
						imageID_array[g] = current_imageID;
					}
					for (h = 0; h < Table.size; h++) {
						setResult("Image_ID", h, imageID_array[h]);
					}
					for (i = 0; i < Table.size; i++) {
						ROIarea_array[i] = Table.get("Area", i);
					}
					for (j = 0; j < Table.size; j++) {
						setResult("ROI_Area", j, ROIarea_array[j]);
					}
					
				
				// Save the variable data (either Area Fraction data or Area Division data) into the results.
				variable_results = newArray(Table.size);
				for (c = 0; c < Table.size; c++) {
					variable_results[c] = Table.get(variables_array[k], c);
				}

				for (d = 0; d < Table.size; d++) {
					setResult(results_target, d, variable_results[d]);
				}
				close("*.csv");
			}
		}
	}
}


function combineCSVs(dir) {
	
	variables_foldernames = newArray("Area Fraction Summaries", "Area Division Summaries");
	combined_filenames = newArray("All-Area-Fraction-Summary-Data.csv", "All-Area-Division-Summary-Data.csv");
	flines = newArray();
	
	for (m = 0; m < variables_foldernames.length; m++) {
		path_summaryresults = path_processedfiles + variables_foldernames[m];
		all_summaryresults = getFileList(path_summaryresults);
		
		for (n = 0; n < all_summaryresults.length; n++) {
			file_name = File.getName(all_summaryresults[n]);
			str = File.openAsString(path_summaryresults + "/" + file_name);
			
			if (n==0){
				flines=str;
			}
			else{
				substr = substring(flines, indexOf(flines, '\n')+1, lengthOf(flines));
				flines = str + substr;
				
			}
			
			}
			File.saveString(flines,path_processedfiles + combined_filenames[m]);
			print(flines);
		}

}


