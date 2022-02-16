 //	********************************************************* MACRO FOR QUATIFYING RNASCOPE SIGNAL IN IMAGEJ ********************************************************************************
 //	******************************************************* alimaran@usc.edu and vmm6714@gmail.com, 02/02/2022 ******************************************************************************

//	************************************************************* UPLOAD DATA AND SET DIRECTORIES ***********************************************************

setBatchMode(true);	//	Runs macro without opening any image windows.

showMessage("Choose folder containing IMAGE subfolders.");
path_imagefolder = getDirectory("");	//	Variable that stores path to the folder containing image subfolders.

showMessage("Choose folder containing ROIs.");
path_roifolder = getDirectory("");	//	Variable that stores path to folder containing ROIs.

showMessage("Choose folder to save RESULTS files.");
path_processedfiles = getDirectory("");	//	Variable that stores path to results folder where files will be saved.

showMessage("Choose .csv containing KEY for average particle size, RBR, and minimum threshold."); 
path_keyfile = File.openDialog("");

processNestedFiles(path_imagefolder);	//	Pass the main image folder containing image subfolders through the processNestedFiles function (see Lines #44 - #67).

showMessage("Almost done!  Creating summary .csv files.  Click OK to proceed.");
createSummary(path_processedfiles);	//	Pass the main results folder containing results subfolders through the createSummary function (see Lines #143 - #222).

combineCSVs(path_processedfiles);	//	Create a final .csv file for combined Area Fraction Results and combined Area Division Results (see lines #225 - #253).

showMessage("Macro finished. Results are saved in the designated folder. :)")	//	Open dialog box that indicates to user that the macro has finished running.


//	********************************************************************* FUNCTIONS ****************************************************************************

function checkKeyValue(header) {	//	Creates a function that checks the .csv key for a specific average particle size, RBR, and threshold minimum.

	open(path_keyfile);
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
		radius = checkKeyValue("RBR");
		run("Subtract Background...", "rolling=radius");

		//	Apply thresholding on .tif image using the lower (min) threshold value set in .csv key		
		min_threshold = checkKeyValue("Threshold Minimum");
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


