/*===============================================================================================
 * Macro        : QIBC (Quantitative Image-Based Cytometry) Analysis Pipeline
 * Description  : Automated pipeline for analyzing both single-plane images and Z-stacks. 
 *                For Z-stacks, it dynamically generates Z-projections (Sum for robust 
 *                quantification, Max for visualization), while single planes are 
 *                processed directly. The macro segments nuclei using Cellpose (DAPI 
 *                channel), filters them by size and borders using MorphoLibJ, applies 
 *                background subtraction, and classifies cells into four populations 
 *                (Double+, Double-, EdU+, H2Ax+) based on user-defined intensity thresholds.
 *
 * Reference    : Sara Martín-Vírgala et al. ,Slow RNAPII elongation enhances naive pluripotency 
 * 				  rewiring while maintaining high replication fork speed.Sci. Adv.12,eadz6211(2026).
 * 				  DOI:10.1126/sciadv.adz6211
 *
 * Facility     : Advanced Light Microscopy Facility (SMOA), CBM
 * Authors      : Carlos Gallego-Garcia
 * Date Created : 2025-07-23
 * Version      : 1.1
 * Last Updated : 2026-10-09
 * 
 * Plugins required:
 *  - PTBIOP (Cellpose wrapper for ImageJ)
 *  - IJPB plugins (MorpholibJ for label processing)
 * 
 *=============================================================================================*/
 
 // Workspace initialization: Clears previous images, ROIs, and memory
	workspaceInitialization();
	
// Parameter Menu: GUI for user-defined variables
	Dialog.create("Parameter Menu");
	
// Image input directory
	Dialog.addMessage("Input Data", 14, "#8A2BE2");
	Dialog.addDirectory("Input images folder:", ""); 
	
// Cellpose environment settings (Requires local Conda environment)
	Dialog.addMessage("Cellpose Settings", 14, "#1E90FF");
	Dialog.addString("Cellpose environment path:", ""); Dialog.addToSameRow();
	Dialog.addNumber("Cellpose Diameter", 50);
	
// Nuclei selection criteria (Size filtering)
	Dialog.addMessage("Nuclei Area Selection", 14, "#32CD32");
	Dialog.addNumber("Minimum nuclei size (px²):", 250); 
	Dialog.addToSameRow();
	Dialog.addNumber("Maximum nuclei size (px²):", 3000);
	
// Intensity thresholds for population classification
	Dialog.addMessage("Intensity levels", 14, "#DAA520");
	Dialog.addNumber("EdU Intensity Threshold", 2300); Dialog.addToSameRow();
	Dialog.addNumber("H2Ax Intensity Threshold", 2100);
	
// Background Rolling Ball radius for shading correction
	Dialog.addMessage("Background Subtraction", 14, "#FF0000");
	Dialog.addNumber("EdU Background radius", 25); Dialog.addToSameRow();
	Dialog.addNumber("H2Ax Background radius", 40);
	Dialog.show();
	
	// Retrieve variables from the dialog
	inputPath 		= Dialog.getString();
	outputPath 		= inputPath+File.separator+"Resultados"+File.separator; File.makeDirectory(outputPath);
	cellposePath 	= Dialog.getString();
	cellposeDiam 	= Dialog.getNumber();
	minCell 		= Dialog.getNumber();
	maxCell 		= Dialog.getNumber();
	eduThres 		= Dialog.getNumber();
	h2axThres 		= Dialog.getNumber();
	eduBkg 			= Dialog.getNumber();
	h2axBkg 		= Dialog.getNumber();
	
// Images Processing Loop
	row=0; 
	Table.create("Final_Data"); // Table for population counts
	Table.create("Mean_Intensities"); // Table for raw single-cell data
	
	imgList = getFileList(inputPath);
	
	for (i = 0; i < imgList.length; i++) {
	// Filter for supported biological image formats
		if (endsWith(imgList[i], ".nd2") || endsWith(imgList[i], ".czi") || endsWith(imgList[i], ".lsm") || endsWith(imgList[i], ".tiff")){
			
		// Open image using Bio-Formats to handle metadata correctly
			run("Bio-Formats Importer", "open=["+inputPath+imgList[i]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			rawID = getImageID();
			imgName = File.nameWithoutExtension; rename(imgName);
			
			getDimensions(width, height, channels, slices, frames);
			
			quantImgName = imgName; // By default, a single plane image will be analyzed
			
		// Z-Projection handling: Sum for robust quantification, Max for qualitative output
			if (slices > 1) {
			// 1. MAX projection for better visualization
				run("Z Project...", "projection=[Max Intensity]");
				rawID = getImageID(); // reassign the rawID from the raw image to the MAX projection
				saveAs("tiff", outputPath+imgName+"_MAX_projection.tiff");
				rename("MAX_"+imgName);
				
			// 2. Crear SUM para cuantificación
				selectImage(imgName);
				run("Z Project...", "projection=[Sum Slices]");
				quantImgName = "SUM_"+imgName; // Reassign the quantification image variable from the raw to the SUM projection
				
			// 3. Closes the raw image to save RAM
				selectImage(imgName); 
				close(); 
			}
			
		// Background Subtraction on quantitative image (SUM o Single Plane)
			selectImage(quantImgName);
			Stack.setChannel(2); // EdU Channel
			run("Subtract Background...", "rolling="+eduBkg+" slice");
			Stack.setChannel(3); // H2Ax Channel
			run("Subtract Background...", "rolling="+h2axBkg+" slice");
			
		// Nuclei Segmentation using Cellpose on DAPI (Channel 1)
			selectImage(rawID); run("Duplicate...", "title=Dapi duplicate channels=1"); run("Tile");
			run("Cellpose ...", "env_path="+cellposePath+" env_type=conda model=cpsam model_path=path\\to\\own_cellpose_model diameter="+cellposeDiam+" ch1=0 ch2=0 additional_flags=[--use_gpu, --flow_threshold, 0.9, --cellprob_threshold, 0.4]");
			
		// Filter ROIs (remove borders and size outliers)
			nRois = labelFiltering(minCell, maxCell, outputPath, imgName);
			
		// Measure Intensity for each channel on the quantitative image
			selectImage(quantImgName); Stack.setChannel(1); // DAPI
			roiManager("measure");
			dapiMean = Table.getColumn("RawIntDen");
			close("Results");
			
			selectImage(quantImgName); Stack.setChannel(2); // EdU
			roiManager("measure");
			eduMean = Table.getColumn("Mean");
			close("Results");
			
			selectImage(quantImgName); Stack.setChannel(3); // H2Ax
			roiManager("measure");
			h2aMean = Table.getColumn("Mean");
			close("Results");
			
		// Append measurements to the raw data table
			selectWindow("Mean_Intensities");
			Table.setColumn(imgName+"_Dapi_Total_Int", dapiMean);
			Table.setColumn(imgName+"_EdU_Int", eduMean);
			Table.setColumn(imgName+"_H2Ax_Int", h2aMean);
			Table.update;
			
		// Population Classification Arrays
			dobleNeg = newArray(); // EdU- / H2Ax-
			doblePos = newArray(); // EdU+ / H2Ax+
			eduOnly  = newArray(); // EdU+ / H2Ax-
			h2aOnly  = newArray(); // EdU- / H2Ax+
			
		// Classify each nucleus based on user-defined thresholds
			for (a = 0; a < eduMean.length; a++) {
				if (eduMean[a] < eduThres && h2aMean[a] < h2axThres) { dobleNeg = Array.concat(dobleNeg,a); }
				else if (eduMean[a] >= eduThres && h2aMean[a] >= h2axThres) { doblePos = Array.concat(doblePos,a); }
				else if (eduMean[a] >= eduThres && h2aMean[a] < h2axThres) { eduOnly = Array.concat(eduOnly,a); }
				else if (eduMean[a] < eduThres && h2aMean[a] >= h2axThres) { h2aOnly = Array.concat(h2aOnly,a); }
			}
			
		// Data Saving: Summary statistics per image
			selectWindow("Final_Data");
			Table.set("Label", row, imgName);
			Table.set("Total Cells", row, nRois);
			
		// Absolute counts
			Table.set("EdU + / H2Ax +", row, doblePos.length);
			Table.set("EdU + / H2Ax -", row, eduOnly.length);
			Table.set("EdU - / H2Ax +", row, h2aOnly.length);
			Table.set("EdU - / H2Ax -", row, dobleNeg.length);
			
			// Percentages relative to total cells
			Table.set("EdU + / H2Ax + (%)", row, doblePos.length/nRois*100);
			Table.set("EdU + / H2Ax - (%)", row, eduOnly.length/nRois*100);
			Table.set("EdU - / H2Ax + (%)", row, h2aOnly.length/nRois*100);
			Table.set("EdU - / H2Ax - (%)", row, dobleNeg.length/nRois*100);
			
			Table.update; row++;
			
		// Visualization: Generate and save an overlay map of the classified populations
			selectImage(rawID);
			Property.set("CompositeProjection", "Sum");
			Stack.setDisplayMode("composite");
			run("RGB Color"); rgbID = getImageID();
			close(rawID); selectImage(rgbID);
			roiManager("show all without labels"); roiManager("Set Line Width", 2);
			
		// Color code ROIs by population
			roiManager("select", doblePos); roiManager("Set Color", "green");
			roiManager("select", dobleNeg); roiManager("Set Color", "red");
			roiManager("select", eduOnly); roiManager("Set Color", "blue");
			roiManager("select", h2aOnly); roiManager("Set Color", "orange"); 
			
			roiManager("show none"); roiManager("show all"); run("Flatten");
			saveAs("Tiff", outputPath+imgName+"_Classification_Map.tiff");
					
		// Cleanup for next iteration
			run("Collect Garbage"); run("Close All"); roiManager("reset");
		}
	}
	
// Export final data tables to Excel-compatible format
	selectWindow("Final_Data");
	saveAs("Results", outputPath+"Resultados_Contajes.xls");
	run("Close");
	
	selectWindow("Mean_Intensities");
	saveAs("Results", outputPath+"Resultados_Intensidad.xls");
	run("Close");
	
	
/*===========================================================================
 *	FUNCTIONS
 *=========================================================================*/

	function workspaceInitialization(){
	// Prepares the ImageJ workspace by closing open images, windows, and resetting the ROI manager. Sets measurement parameters.
		if (nImages > 0) { run("Close All"); } 
		roiManager("reset"); 
		selectWindow("ROI Manager"); run("Close"); 
		close("Log"); close("Results"); 
		run("Collect Garbage");	
		run("Bio-Formats Macro Extensions"); 
		
		if (isOpen("Final_Data")) { 
			selectWindow("Final_Data"); Table.deleteRows(0, 0, "Final_Data"); 
		} else { 
			Table.create("Final_Data"); 
		}
	// Configure global measurements (Mean for intensities, RawIntDen for DAPI)
		run("Set Measurements...", "area mean min integrated redirect=None decimal=3");
	}
	
	function labelFiltering(minCell, maxCell, output_folder, imgName) {
	// Processes Cellpose label output using MorphoLibJ. Removes labels touching the image borders and filters out ROIs
	// outside the user-defined size range (minCell, maxCell). Saves the final valid ROIs to a zip file.
		
	// Handle border killing depending on how Cellpose named the output window
		if (isOpen("Label Image")) { selectImage("Label Image"); run("Kill Borders"); } 
		else if (isOpen("Dapi-cellpose")) { selectImage("Dapi-cellpose"); rename("Label Image"); run("Remove Border Labels", "left right top bottom"); }
		close("Label Image");
		
	// Apply MorphoLibJ size filters
		selectImage("Label-killBorders"); 
		run("Label Size Filtering", "operation=Greater_Than size="+minCell); 
		close("Label-killBorders");
		
		selectImage("Label-killBorders-sizeFilt"); 
		run("Label Size Filtering", "operation=Lower_Than size="+maxCell); 
		close("Label-killBorders-sizeFilt");
		
	// Convert final filtered label map to ImageJ ROIs
		run("Label image to ROIs"); 
		close("Label-killBorders-sizeFilt-sizeFilt");
		
		nRois = roiManager("count");
		roiManager("save", output_folder+imgName+"_ROIset_nuclei.zip");
		
		return nRois;
	}