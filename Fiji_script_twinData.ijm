//Study concerning...
//This script aims at batch processing 3 sometimes 4 channel data, where the first channel->sarcolemma, the second channel->neutral lipids, the 3rd channel->PLIN marker (2 or 5) and sometimes a 4th channel-> fiber type (I or II).
//gray value density, binary density, size, colocalization and fiber type are to be measured
//Publication main author Vasco Fachada
//Code by Vasco Fachada
denVar = "Mean";
//denVar = "IntDen";

while (nImages>0) { 
		selectImage(nImages); 
        close();    				  
	} 
	
list = getList("window.titles");
	for (i=0; i<list.length; i++){
    	winame = list[i];
      	selectWindow(winame);
      	//print(winame);
     	run("Close");
     }

start = getTime();
//Finding the files
dir = getDirectory("Select the directory containing all your twin data");
sFolders = getFileList(dir);
//Creating folder where ready images and raw data file will be be saved
//savingDir = dir+File.separator+"results";
//File.makeDirectory(savingDir);

setBatchMode(true);
newImage("new", "8-bit black", 1024, 1024, 1);
stacks = newArray("LD_stack_PLIN2", "LDbin_stack_PLIN2", "LD_stack_PLIN5", "LDbin_stack_PLIN5");
for (st=0; st<stacks.length; st++) {
	saveAs("tiff", dir+File.separator+stacks[st]+".tif");
}
close();


run("Bio-Formats Macro Extensions");
for(s=0; s<sFolders.length; s++) {
	
	if (substring(sFolders[s], 0,1) !="x"){
		subjFold = dir+sFolders[s];	
		curSubj = sFolders[s];
		vFolders = getFileList(subjFold);	
		for(v=0; v<vFolders.length; v++) {
		
			if (substring(vFolders[v], 5) =="_LD_20x1/"){
				varFold = subjFold+File.separator+vFolders[v];
				plinVar = substring(vFolders[v], 0, 5);
				datafile = subjFold+File.separator+plinVar+"_datafile.csv";
				if (File.exists(datafile)){
					File.delete(datafile);	
				}
				files = getFileList(varFold);			
				for(f=0; f<files.length; f++) {
					if (endsWith(files[f], ".tif")){
						id = varFold+File.separator+files[f];
						Ext.setId(id);
						filename=substring(files[f], 0, lengthOf(files[f])-4);
						print('opening file: '+filename);
						print('Analyzing '+plinVar+' in subject '+curSubj);

						openFile();
						//lipids();
						prepareImg();
						runCells();
						
						while (nImages>0) { 
							selectImage(nImages); 
        					close();   
						} 
						call("java.lang.System.gc");
					}
				}
			}
		}
	}
}

//makeMontage();

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


function prepareImg(){
	run("8-bit");
	run("Duplicate...", "title=main_membrane");
	setAutoThreshold("Percentile");
	run("Convert to Mask");
	run("Analyze Particles...", "size=500.00-Infinity circularity=0.50-1.00 add");

	selectWindow("main");
	run("Duplicate...", "title=gray_main duplicate channels=2-3");
	run("Subtract Background...", "rolling=10 stack");
	setSlice(2);
	run("Duplicate...", "title=grayPLIN");
	selectWindow("gray_main");
	setSlice(1);
	run("Duplicate...", "title=grayLD");	
	run("Duplicate...", "title=binLD duplicate channels=2");
	run("Convert to Mask", "method=Moments background=Dark");
	run("Watershed");
	
	run("Select All");
	run("Copy");
	open(dir+File.separator+"LDbin_stack_"+plinVar+".tif");
	run("Add Slice");
	run("Paste");
	run("Set Label...", "label="+curSubj);
	saveAs("tiff", dir+File.separator+"LDbin_stack_"+plinVar+".tif");
	close();
	selectWindow("grayLD");
	run("Select All");
	run("Copy");
	open(dir+File.separator+"LD_stack_"+plinVar+".tif");
	run("Add Slice");
	run("Paste");
	run("Set Label...", "label="+curSubj);
	saveAs("tiff", dir+File.separator+"LD_stack_"+plinVar+".tif");
	close();
		
	selectWindow("main");
	setSlice(4);
	run("Duplicate...", "title=FT");
	setAutoThreshold("Percentile dark");
	run("Convert to Mask");
}

function runCells(){
	for (i=0; i<roiManager("count"); i++){
		cellVars = detectFiberType();
		LDVars = measureLDs();
		PLINdensity = measure_PLIN();
		colVars = colocalization();
		
		writeResults(cellVars[1], LDVars[0], LDVars[1], cellVars[0], LDVars[2], PLINdensity, colVars[0], colVars[1], colVars[2], colVars[3]);
	}
	roiManager("Delete");
    roiManager("Delete");
}


function detectFiberType(){
	selectWindow("FT");
	roiManager("select", i);
	List.setMeasurements;
	cellArea = List.getValue("Area");
	if (List.getValue("Mean") > 150){
		ft = "Type I";
	}else{
		ft = "Type II";
	}
	run("Clear Results");

	cellVars = newArray(cellArea, ft);
	return cellVars;
}


function measureLDs(){
	selectWindow("grayLD");
	roiManager("select", i);
    List.setMeasurements;
	LDdensity = List.getValue(denVar);
    
	selectWindow("binLD");
	roiManager("select", i);
    run("Analyze Particles...", "  circularity=0.0-1.00 show=Nothing summarize");
    selectWindow("Summary");
	IJ.renameResults("Summary","Results");
	LDcount = getResult("Count", 0);
	LDsize = getResult("Average Size", 0);

	selectWindow("Results"); 
	run("Close"); 

	LDVars = newArray(LDdensity, LDcount, LDsize);
	return LDVars;
}


function measure_PLIN(){
	selectWindow("grayPLIN");
	roiManager("select", i);
    List.setMeasurements;
	PLINdensity = List.getValue(denVar);

	run("Clear Results"); 
	return PLINdensity;
}
			


function colocalization(){
	selectWindow("grayLD");
	roiManager("select", i);
	run("Colocalization Threshold", "channel_1=grayLD channel_2=grayPLIN use=[Channel 1] channel=[Red : Green] set mander's mander's_0 number %");
	setOption("ScaleConversions", true);
	setOption("ScaleConversions", true);

	wholeStr=replace(getInfo("window.contents"), ',', '.');
	lines=split(replace(wholeStr, '%', ''), "\n");
	values=split(lines[2], "\t");
	M1 = parseFloat(values[3]);
	M2 = parseFloat(values[4]);
	tM1 = parseFloat(values[5]);
	tM2 = parseFloat(values[6]);
	//VolPerc = parseFloat(values[7]);

	selectWindow("Results"); 
	run("Close"); 

	colVars = newArray(M1, M2, tM1, tM2);
	return colVars;
}




function writeResults(fiberType, LDdensity, LDcount, cellArea, LDsize, PLINdensity, M1, M2, tM1, tM2){
	
	//Create datafile in case there is none
	if (!File.exists(datafile)){
    	// Create a header
    	File.append("Cell_ID,PLIN_stained,FiberType,"
    	+"LD_signal_density,LD_#_density,LD_AvgSize,LD#,PLIN_signal_density,"
    	+"LD_PLIN_Manders,PLIN_LD_Manders,LD_PLIN_tManders,PLIN_LD_tManders", datafile);
	}


    // insert data into table
	File.append(Roi.getName+","+plinVar+","+fiberType+","
	+LDdensity+","+LDcount/cellArea+","+LDsize+","+LDcount+","+PLINdensity+","
	+M1+","+M2+","+tM1+","+tM2+",", datafile);
	
	//saveData();
}


function makeMontage(){
	for (st=0; st<stacks.length; st++) {
		open(dir+File.separator+stacks[st]+".tif");
		run("Delete Slice");
		run("Fire");
		run("Make Montage...", "columns="+nSlices/2+" rows="+nSlices/4+" scale=1 font=24 label");
		saveAs("tiff", dir+File.separator+"montage_"+stacks[st]+".tif");
		close();
	}
}
