/*
Started on June, 2022

This script aims to analyze confocal images of C2C12 cells treated with BCAA and EPS in terms of LD, PLIN5, PLIN2, PGC1a and PCK1 contents and their colocalization. 
Aditionally, lipid droplet size is also measured. 
In each image, the measurements are preformed for both myoblasts and myotubes (marked with MF-20)

code by Vasco Fachada

Principal Investigator Heikki Kainulainen
*/

while (nImages>0) { 
		selectImage(nImages); 
        close();    				  
	} 
	
list = getList("window.titles");
	for (i=0; i<list.length; i++){
    	winame = list[i];
      	selectWindow(winame);
      	print(winame);
     	run("Close");
     }

start = getTime();

//Finding the files
dir = getDirectory("Select the 'READY' directory containing your processed .tif files");
files = getFileList(dir);

setBatchMode(true);
dataTypes = newArray("tubeData.tif", "blastData.tif");
datafile = dir+"FINAL_COLOC_ANALYSIS.csv";
if (File.exists(datafile)){
	File.delete(datafile);	
}

run("Bio-Formats Macro Extensions");
for(s=0; s<files.length; s++){
	if (substring(files[s], lengthOf(files[s])-13, lengthOf(files[s])) == "blastData.tif"){
		if (substring(files[s], 0, 1) == "1"){
			marker4 = "PLIN2";
		}else if (substring(files[s], 0, 1) == "2"){
			marker4 = "PGC1a";
		}else{
			marker4 = "PCK1";
		}

		if (substring(files[s], 3, 6) == "EPS"){
			EPS = "yes";
		}else{
			EPS = "no";
		}

		if (substring(files[s], 7, 9) == "08"){
			BCAA = "yes";
		}else{
			BCAA = "no";
		}
		
		for(dt=0; dt<dataTypes.length; dt++){
			curFile = substring(files[s], 0, 18)+dataTypes[dt];
			print(curFile);
			id = dir+curFile;
			
			if (substring(dataTypes[dt], 0, 1) == "t"){
				cellType = "tubes";
			}else{
				cellType = "blasts";
			}
			
			print("4th marker is: "+marker4);
			print("Cell Type is: "+cellType);
			print("EPS: "+EPS);
			print("BCAA: "+BCAA);
			openFile();
			prepareImgs();
			densiValues = density();
			colocValues = coloc();
			writeResults(densiValues, colocValues);
		}
		//waitForUser;
		while (nImages>0) { 
				selectImage(nImages); 
				close();   
		} 
		call("java.lang.System.gc");
	}
}
end = getTime();
//saveAs("Text", savingDir+File.separator+"analysis_Log.txt");
Dialog.create("DONE!");
Dialog.addMessage("The analysis was completed within "+(end-start)/1000+" seconds ("+(end-start)/1000/60+" minutes).\n \nI can tell you, these numbers were hard to crunch, hope you get a Nobel for this ;)");
Dialog.show();


function openFile(){
	//OPENING FILE
	run("Bio-Formats Importer", "open=["+id+"] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename("main");
	getVoxelSize(pWidth, pHeight, depth, unit); //detecting the voxel size of the current raw file, so the final image also has the same scale.
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Properties...", "channels="+channels+" slices="+slices+" frames="+frames+" unit="+unit+" pixel_width="+pWidth+" pixel_height="+pHeight+" voxel_depth="+depth+" global");
	run("Set Measurements...", "area mean modal min integrated redirect=None decimal=3");
}

function prepareImgs(){
	//Preparing and segmenting MF-20 for myotubes
	selectWindow("main");
	setSlice(5);
	run("Duplicate...", "title="+cellType+"Nuclei");
	run("HiLo");
	setAutoThreshold("Otsu dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	
	run("Fill Holes");
	run("Despeckle");
	run("Remove Outliers...", "radius=30 threshold=50 which=Dark");
	run("Dilate");
	run("Dilate");
	run("Fill Holes");
	
	if (cellType == "tubes"){
		selectWindow("main");
		setSlice(2);
		run("Duplicate...", "title=tubeBin");
		run("HiLo");
		setAutoThreshold("Percentile dark");
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Remove Outliers...", "radius=10 threshold=50 which=Bright");
		run("Fill Holes");
		imageCalculator("Subtract create", "tubeBin", "tubesNuclei");
		selectWindow("Result of tubeBin");
		rename("tubesCytosol");
	}else{
		imageCalculator("Add create", "tubeBin","blastsNuclei");
		selectWindow("Result of tubeBin");
		run("Invert");
		rename("blastsCytosol");
	}

	selectWindow("main");
	setSlice(1);
	run("Duplicate...", "title="+cellType+"grayPLIN5");
	selectWindow("main");
	setSlice(3);
	run("Duplicate...", "title="+cellType+"grayLDs");
	selectWindow("main");
	setSlice(4);
	run("Duplicate...", "title="+cellType+"gray"+marker4);
	close("main");

	//prepare gray image from full PLIN5-LDs coloc in order to use later with PGC1-alpha
	run("Colocalization Threshold", "channel_1="+cellType+"grayPLIN5 channel_2="+cellType+"grayLDs use=[Channel 1] channel=[Red : Green] show include set mander's mander's_0 number %");	
	selectWindow("Colocalized Pixel Map RGB Image");
	run("Split Channels");
	selectWindow("Colocalized Pixel Map RGB Image (blue)");
	rename(cellType+"grayPLIN5_LDs");
	if (isOpen("Results") == true){
		selectWindow("Results");
		run("Close");
	}
			
}


function density(){	
	selectWindow(cellType+"Cytosol");
	run("Create Selection");
	LDsigC = measure_density("LDs");	
	PLIN5sigC = measure_density("PLIN5");	
	mk4sigC = measure_density(marker4);
	
    selectWindow(cellType+"Nuclei");
	run("Create Selection");
	LDsigN = measure_density("LDs");	
	PLIN5sigN = measure_density("PLIN5");	
	mk4sigN = measure_density(marker4);	

	densiValues = newArray(LDsigC, PLIN5sigC, mk4sigC, LDsigN, PLIN5sigN, mk4sigN); 
	return densiValues;

	function measure_density(marker){
		selectWindow(cellType+"gray"+marker);
		run("Restore Selection");
		List.setMeasurements;
		return List.getValue("Mean");
	}
}
			


function coloc(){
	
	PLIN5LDc = measure_coloc("Cytosol", "PLIN5", "LDs");
	PLIN5LDn = measure_coloc("Nuclei", "PLIN5", "LDs");
	LDPLIN5c = measure_coloc("Cytosol", "LDs", "PLIN5");
	LDPLIN5n = measure_coloc("Nuclei", "LDs", "PLIN5");
	
	mk4LDc = measure_coloc("Cytosol", marker4, "LDs");
	mk4LDn = measure_coloc("Nuclei", marker4, "LDs");
	LDmk4c = measure_coloc("Cytosol", "LDs", marker4);
	LDmk4n = measure_coloc("Nuclei", "LDs", marker4);
	
	PLIN5mk4c = measure_coloc("Cytosol", "PLIN5", marker4);
	PLIN5mk4n = measure_coloc("Nuclei", "PLIN5", marker4);
	mk4PLIN5c = measure_coloc("Cytosol", marker4, "PLIN5");
	mk4PLIN5n = measure_coloc("Nuclei", marker4, "PLIN5");

	if (marker4 == "PGC1a"){
		pgcPLIN5LDc = measure_coloc("Cytosol", marker4, "PLIN5_LDs");
		pgcPLIN5LDn = measure_coloc("Nuclei", marker4, "PLIN5_LDs");
		PLIN5LDpgcc = measure_coloc("Cytosol", "PLIN5_LDs", marker4);
		PLIN5LDpgcn = measure_coloc("Nuclei", "PLIN5_LDs", marker4);		
	}else{
		pgcPLIN5LDc = "NaN";
		pgcPLIN5LDn = "NaN";
		PLIN5LDpgcc = "NaN";
		PLIN5LDpgcn = "NaN";
	}
	
	colocValues = newArray(PLIN5LDc, PLIN5LDn, LDPLIN5c, LDPLIN5n, mk4LDc, mk4LDn, LDmk4c, LDmk4n, PLIN5mk4c, PLIN5mk4n, mk4PLIN5c, mk4PLIN5n, pgcPLIN5LDc, pgcPLIN5LDn, PLIN5LDpgcc, PLIN5LDpgcn);
	return colocValues;
	
	function measure_coloc(cc, mrkA, mrkB){
		selectWindow(cellType+cc);
		run("Create Selection");
		print("Performing intensity correlation analysis between "+mrkA+" and "+mrkB+" in "+cellType+cc);
		//waitForUser;
		selectWindow(cellType+"gray"+mrkA);
		run("Restore Selection");
		
		run("Colocalization Threshold", "channel_1="+cellType+"gray"+mrkA+" channel_2="+cellType+"gray"+mrkB+" use=[Channel 1] channel=[Red : Green] include set mander's mander's_0 number %");
		setOption("ScaleConversions", true);
		
		if (isOpen("Results") == false){
			tM1 = "NaN";
		}else{
			selectWindow("Results");
			wholeStr=replace(getInfo("window.contents"), ',', '.');
			lines=split(replace(wholeStr, '%', ''), "\n");
			values=split(lines[2], "\t");
			tM1 = parseFloat(values[5]);
			//tM2 = parseFloat(values[6]);
			selectWindow("Results"); 
			run("Close");
		} 
		return tM1
	}
}



function writeResults(densiValues, colocValues){
	PLIN5LDc = colocValues[0];
	PLIN5LDn = colocValues[1];
	LDPLIN5c = colocValues[2];
	LDPLIN5n = colocValues[3];
	mk4LDc = colocValues[4];
	mk4LDn = colocValues[5];
	LDmk4c = colocValues[6];
	LDmk4n = colocValues[7];
	PLIN5mk4c = colocValues[8];
	PLIN5mk4n = colocValues[9];
	mk4PLIN5c = colocValues[10];
	mk4PLIN5n = colocValues[11];
	pgcPLIN5LDc = colocValues[12];
	pgcPLIN5LDn = colocValues[13];
	PLIN5LDpgcc = colocValues[14];
	PLIN5LDpgcn = colocValues[15];

	LDsigC = densiValues[0];
	PLIN5sigC = densiValues[1]; 
	mk4sigC = densiValues[2]; 
	LDsigN = densiValues[3]; 
	PLIN5sigN = densiValues[4]; 
	mk4sigN = densiValues[5];
	
	//Create datafile in case there is none
	if (!File.exists(datafile)){
    	// Create a header
    	File.append("fileName,BCAA,EPS,cell type,marker4," 
    	+ "LD_signal_fractionc,LD_signal_fractionn,PLIN5_signal_fractionc,PLIN5_signal_fractionn,mrk4_signal_fractionc,mrk4_signal_fractionn,"
    	+ "PLIN5.vs.LDc,PLIN5.vs.LDn,LD.vs.PLIN5c,LD.vs.PLIN5n,"
    	+ "mrk4.vs.LDc,mrk4.vs.LDn,LD.vs.mrk4c,LD.vs.mrk4n,"
    	+ "PLIN5.vs.mrk4c,PLIN5.vs.mrk4n,mrk4.vs.PLIN5c,mrk4.vs.PLIN5n,PGC1a.vs.PLIN5-LDc,PGC1a.vs.PLIN5-LDn,PLIN5-LD.vs.PGC1ac,PLIN5-LD.vs.PGC1an", datafile);
	}


    // insert data into table
	File.append(curFile+","+BCAA+","+EPS+","+cellType+","+marker4+","
	+LDsigC+","+LDsigN+","+PLIN5sigC+","+PLIN5sigN+","+mk4sigC+","+mk4sigN+","
	+PLIN5LDc+","+PLIN5LDn+","+LDPLIN5c+","+LDPLIN5n+","
	+mk4LDc+","+mk4LDn+","+LDmk4c+","+LDmk4n+","
	+PLIN5mk4c+","+PLIN5mk4n+","+mk4PLIN5c+","+mk4PLIN5n+","+pgcPLIN5LDc+","+pgcPLIN5LDn+","+PLIN5LDpgcc+","+PLIN5LDpgcn, datafile);
	
	//saveData();
}