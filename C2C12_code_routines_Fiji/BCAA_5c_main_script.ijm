/*
Started on August, 2020

This script aims to analyze confocal images of C2C12 cells treated with BCAA and EPS in terms of LD, PLIN5, PLIN2, PGC1a and PCK1 contents and their colocalization. 
Aditionally, lipid droplet size is also measured. 
In each image, the measurements are preformed for both myoblasts and myotubes (marked with MF-20)

code by Vasco Fachada

Principal Investigator Heikki Kainulainen
*/


//Finding the .czi files
dir = getDirectory("Select a directory containing your CZI files you want to analyze");
files = getFileList(dir);
start = getTime();

setBatchMode(true);
run("Misc...", "divide=Infinity hide run");

k=0;
n=0;
pixelTotal = 0;
sizeLoaded = 0;

run("Bio-Formats Macro Extensions");
for(f=0; f<files.length; f++) {

	if(endsWith(files[f], ".czi")) { 
		
		//Clearing previous images 
		while (nImages>0) { 
			selectImage(nImages); 
    		close();
		}
		
		k++;
		id = dir+files[f];
		
		//defining tyoe of file
		magnif = substring(id, indexOf(id, "BCAA_")+5, lengthOf(id)-4);
		if (magnif == "63x"){
			thickness = 0.8;
		}
		else{
			thickness = 1.8;
		}
		if (substring(id, indexOf(id, "BCAA_")-6, indexOf(id, "BCAA_")-3) == 'CTR'){
			EPS = 'no';
		}
		else{
			EPS = 'yes';
		}
		if (substring(id, indexOf(id, "BCAA_")-2, indexOf(id, "BCAA_")) == 'no'){
			BCAA = 'no';
		}
		else{
			BCAA = '08';
		}
		if (substring(id, indexOf(id, "BCAA_")-9, indexOf(id, "BCAA_")-8) == '1'){
			guest_marker = 'PLIN2';
		}
		else if (substring(id, indexOf(id, "BCAA_")-9, indexOf(id, "BCAA_")-8) == '2'){
			guest_marker = 'PGC1a';
		}
		else{
			guest_marker = 'PCK1';
		}
		print(magnif);
		print(EPS);
		print(BCAA);
		print(guest_marker);
		Ext.setId(id);
		Ext.getSeriesCount(seriesCount);
		n+=seriesCount;

		//Creating folder where ready images and raw data file will be be saved
		savingDir = dir+File.separator+"deconvolved"+File.separator+"READY";

		//READING METADATA
		run("Bio-Formats Importer", "open=["+id+"] color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack series_"+(1));
		fullMeta=getMetadata("Info");
	

		//Looking and determining number of channels in data
		channels_nmr = 	parseInt(substring(fullMeta, indexOf(fullMeta, "SizeC =")+8, indexOf(fullMeta, "SizeC =")+9)); 
		channels = newArray(channels_nmr);
		colors = newArray(channels_nmr);
		markers = newArray(channels_nmr);		
				
		for (ii=0; ii<channels_nmr; ii++){
			channels[ii] = substring(fullMeta, indexOf(fullMeta, "EmissionWavelength #"+ii+1+" =")+24, indexOf(fullMeta, "EmissionWavelength #"+ii+1+" =")+27);
			//print(channels[ii]);
			if (channels[ii] == "435") markers[ii] = "DAPI";
			if (channels[ii] == "580") markers[ii] = "PLIN5";
			if (channels[ii] == "585") markers[ii] = "LD540";
			if (channels[ii] == "572") markers[ii] = guest_marker;
			if (channels[ii] == "593") markers[ii] = "MF-20";
			//print(markers[ii]);
		}
		//waitForUser;
		getVoxelSize(width, height, depth, unit); //detecting the voxel size of the current raw file, so the final image also has the same scale.
		close();

		if (!File.exists(savingDir)){			
			File.makeDirectory(savingDir);
		}
		
		//OPENING TILES from .czi file orderly. For analysis (parallel macro) and stiching

		run("Bio-Formats Importer", "open=["+id+"] color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack series_");

		fullName	= getTitle();
		dirName 	= substring(fullName, 0,lastIndexOf(fullName, ".czi"));
		fileName 	= substring(fullName, 0, lengthOf(fullName)-4);

		id_= dir+File.separator+"deconvolved"+File.separator+substring(files[f], 0, lengthOf(files[f])-4)+".tif";
		print("listen!!!!: "+id_);
		run("Bio-Formats Importer", "open=["+id_+"] color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack series_"); //these two lines are opening and running only deconvoluted data
		
		rename("main");
		run("Properties...", "channels=5 slices=1 frames=1 unit="+unit+" pixel_width="+width+" pixel_height="+height+" voxel_depth=1 global");
		run("Duplicate...", "title=visual duplicate");
		Stack.setDisplayMode("color");
		//Stack.setDisplayMode("composite");
		//run("Stack to RGB");
		saveAs("tiff", savingDir+File.separator+fileName+".tif");
		close(fileName+".tif");
		close("visual");
		
		imageProcessing();
		print("cgheckpoint");
		getData();		
		saveImgs();								
	}	
}	  				
					

run("Misc...", "divide=Infinity run");
run("Set Scale...", "distance=1 known=1 pixel=1 unit=microns");//removes global calibration, necessary to analyze different magnification images in one run
//Clearing previous images 
while (nImages>0) { 
	selectImage(nImages); 
    close();
}

end = getTime();
Dialog.create("Here you go");
Dialog.addMessage("Done! \n...within "+(end-start)/1000+" seconds ("+(end-start)/1000/60+" minutes).\n \nYou know, that's Super Computer kind of stuff...clap clap!");
Dialog.show();

function imageProcessing(){
	//Preparing DAPI image for weka segmentation
	open("D:/VASCO/bcaa_final_testing/regular_BCAAEPS/repeat/"+magnif+"/temp_dapi/"+fileName+".tif_class.tif"); //check this path if the dataset changes!!! "repeat" or not!
	run("Size...", "width=2048 height=2048 depth=1 constrain average interpolation=Bilinear");
	setAutoThreshold("Default");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	outlier_rad = 50;
	if (magnif == "20x"){
		//run("Invert");
		outlier_rad = 25;
	}
	run("Watershed");
	run("Remove Outliers...", "radius="+outlier_rad+" threshold=50 which=Dark");
	rename("weka_DAPI_ready");

	//Preparing and segmenting MF-20 for myotubes
	selectWindow("main");
	setSlice(2);
	run("Duplicate...", "title=tubes_ready");
	run("HiLo");
	setAutoThreshold("Triangle dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	if (magnif == "63x"){
		run("Fill Holes");
		run("Despeckle");
		run("Remove Outliers...", "radius=30 threshold=50 which=Dark");
		run("Dilate");
		run("Dilate");
		run("Fill Holes");
	}else{
		run("Despeckle");
		run("Remove Outliers...", "radius=20 threshold=50 which=Dark");
		run("Dilate");
		run("Remove Outliers...", "radius=25 threshold=50 which=Bright");		
	}
	
	//Create mask with nuclei only from myoblasts
	imageCalculator("Subtract create", "weka_DAPI_ready","tubes_ready");
	run("Remove Outliers...", "radius=5 threshold=50 which=Dark");
	rename("blast_nuclei");

	//Both blast_nuclei and tubes_ready NEED TO BE A STACK OF 5!!
	run("Duplicate...", "title=blast_nuclei_");//first the blasts
	run("Duplicate...", "title=blast_nuclei__");
	run("Duplicate...", "title=blast_nuclei___");
	run("Duplicate...", "title=blast_nuclei___");
	run("Images to Stack", "method=[Copy (center)] name=Stack title=blast_nuclei");
	rename("blast_nuclei_stack");

	selectWindow("tubes_ready");
	run("Duplicate...", "title=tubes_ready_");//then the tubes
	run("Duplicate...", "title=tubes_ready__");
	run("Duplicate...", "title=tubes_ready___");
	run("Duplicate...", "title=tubes_ready____");
	run("Images to Stack", "method=[Copy (center)] name=Stack title=tubes_ready");
	rename("tubes_ready_stack");
	

	
	//Subtract tubes and blast_nuclei masks from gray scale data in "main"
	selectWindow("main");
	run("Duplicate...", "title=blast_data duplicate");
	run("Invert", "stack");
	run("Duplicate...", "title=tube_data duplicate");
	
	imageCalculator("Subtract create stack", "tubes_ready_stack","tube_data");
	selectWindow("Result of tubes_ready_stack");
	//run("Invert", "stack");
	rename("tubes_data_ready");
	
	selectWindow("tubes_ready_stack");
	run("Invert", "stack");	
	imageCalculator("Subtract create stack", "tubes_ready_stack","blast_data");//first the blasts
	selectWindow("Result of tubes_ready_stack");
	//run("Invert", "stack");
	rename("blast_data_ready");
}





function getData(){
	print("Starting to measure and collect data...");
	//First the basics. Get the areas of the respective signals for normalizations
	selectWindow("tubes_ready_stack");
	run("Create Selection");
	run("Make Inverse");
	List.setMeasurements;
	tubes_vol = List.getValue ("Area")*thickness;

	selectWindow("blast_nuclei_stack");
	run("Create Selection");
	List.setMeasurements;
	blasts_vol = List.getValue ("Area")*thickness;
	//Count nuclei number
	run("Analyze Particles...", "clear");
	nucNum=nResults();
	close();
	
	//MEASURING AND CREATING VARIABLES
	//measure signal
	selectWindow("tubes_data_ready");
	setSlice(1);
	List.setMeasurements;
	tubPLIN5sigfr = List.getValue ("Mean") / tubes_vol;
	setSlice(3);
	List.setMeasurements;
	tubLDsigfr = List.getValue ("Mean") / tubes_vol;
	setSlice(4);
	List.setMeasurements;
	tubGuestsigfr = List.getValue ("Mean") / tubes_vol;
	selectWindow("blast_data_ready");
	setSlice(1);
	List.setMeasurements;
	blaPLIN5sigfr = List.getValue ("Mean") / blasts_vol;
	setSlice(3);
	List.setMeasurements;
	blaLDsigfr = List.getValue ("Mean") / blasts_vol;
	setSlice(4);
	List.setMeasurements;
	blaGuestsigfr = List.getValue ("Mean") / blasts_vol;

	//COLOCALIZATION ANALYSIS
	selectWindow("tubes_data_ready");
	setSlice(1);
	run("Duplicate...", "title=tubes_PLIN5");
	selectWindow("tubes_data_ready");
	setSlice(3);
	run("Duplicate...", "title=tubesLD");
	selectWindow("tubes_data_ready");
	setSlice(4);
	run("Duplicate...", "title=tubes_"+guest_marker+"");
	run("Colocalization Threshold", "channel_1=tubesLD channel_2=tubes_PLIN5 use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"tubeLD_PLIN5_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close"); 
	run("Colocalization Threshold", "channel_1=tubesLD channel_2=tubes_"+guest_marker+" use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"tubeLD_"+guest_marker+"_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close");
	run("Colocalization Threshold", "channel_1=tubes_PLIN5 channel_2=tubes_"+guest_marker+" use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"tubePLIN5_"+guest_marker+"_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close");

	selectWindow("blast_data_ready");
	setSlice(1);
	run("Duplicate...", "title=blasts_PLIN5");
	selectWindow("blast_data_ready");
	setSlice(3);
	run("Duplicate...", "title=blastsLD");
	selectWindow("blast_data_ready");
	setSlice(4);
	run("Duplicate...", "title=blasts_"+guest_marker+"");
	run("Colocalization Threshold", "channel_1=blastsLD channel_2=blasts_PLIN5 use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"blastLD_PLIN5_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close"); 
	run("Colocalization Threshold", "channel_1=blastsLD channel_2=blasts_"+guest_marker+" use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"blastLD_"+guest_marker+"_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close");
	run("Colocalization Threshold", "channel_1=blasts_PLIN5 channel_2=blasts_"+guest_marker+" use=None channel=[Red : Blue] show show include");
	selectWindow("Results");
	saveAs("Text", savingDir+File.separator+fileName+"blastPLIN5_"+guest_marker+"_coloc"+magnif+".tsv");
	selectWindow("Results"); 
	run("Close"); 

	//measure LDs size and number
	selectWindow("tubesLD");//first in myotubes
	setAutoThreshold("MaxEntropy dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Invert");
	run("Analyze Particles...", "  show=Nothing display clear");
	tubLDnum = nResults();
	size_count = 0;
	for(ld=0; ld<nResults()-1; ld++) {
		size_count = size_count + getResult("Area",ld);
	}
	tubLDarea = size_count/nResults();
	tubLDdiam = 2*sqrt(tubLDarea/3.14159);
	tubLDnum_per_mm = (tubLDnum / tubes_vol) * 1000;
	selectWindow("Results"); 
	run("Close");
	
	
	selectWindow("blastsLD");//then in myoblasts
	setAutoThreshold("MaxEntropy dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Invert");
	run("Analyze Particles...", "  show=Nothing display clear");
	blastLDnum = nResults();
	size_count = 0;
	for(ld=0; ld<nResults()-1; ld++) {
		size_count = size_count + getResult("Area",ld);
	}
	blastLDarea = size_count/nResults();
	blaLDdiam = 2*sqrt(blastLDarea/3.14159);
	blaLDnum_per_mm = (blastLDnum / blasts_vol) * 1000;//per cubic mm
	selectWindow("Results"); 
	run("Close");


/*
print(fileName); 
print(BCAA); 
print(EPS); 
print(tubLDsigfr); 
print(tubPLIN5sigfr); 
print(tubGuestsigfr); 
print(tubLDnum_per_mm); 
print(tubLDdiam); 
print(blaLDsigfr); 
print(blaPLIN5sigfr); 
print(blaGuestsigfr); 
print(blaLDnum_per_mm);
print(blaLDdiam);
*/		
	data_file = savingDir+File.separator+"data_table_"+magnif+"_"+guest_marker+".txt";
	if (!File.exists(data_file)){
    	// Create a header
    	File.append("fileName\tBCAA\tEPS\tTube_LD_sigfrac\tTube_PLIN5_sigfrac\tTube_"+guest_marker+"_sigfrac\tTube_LD_#\tTube_LD_#per_mm3\tTube_LD_diam\tBlast_LD_sigfrac\tBlast_PLIN5_sigfrac\tBlast_"+guest_marker+"_sigfrac\tBlast_LD_#\tBlast_LD_#per_mm3\tBlast_LD_diam", data_file);

	}
	File.append(fileName+ "\t" + BCAA + "\t" + EPS + "\t" + tubLDsigfr + "\t" + tubPLIN5sigfr + "\t" + tubGuestsigfr + "\t" + tubLDnum+ "\t" + tubLDnum_per_mm + "\t" + tubLDdiam + "\t" + blaLDsigfr + "\t" + blaPLIN5sigfr + "\t" + blaGuestsigfr + "\t" + blastLDnum + "\t" + blaLDnum_per_mm + "\t" + blaLDdiam , data_file);							
}

function saveImgs(){
	//save image only with myotube data
	selectWindow("weka_DAPI_ready");
	run("Invert");
	selectWindow("tubes_data_ready");
	setSlice(5);
	run("Invert");
	imageCalculator("Add create", "weka_DAPI_ready","tubes_data_ready");
	selectWindow("Result of weka_DAPI_ready");
	run("Select All");
	run("Copy");
	selectWindow("tubes_data_ready");
	run("Paste");
	run("Invert");
	saveAs("tiff", savingDir+File.separator+fileName+"_tubeData.tif");
	selectWindow("Result of weka_DAPI_ready");
	close();
	
	//save image only with myblast data
	selectWindow("blast_data_ready");
	setSlice(5);
	run("Invert");
	imageCalculator("Add create", "weka_DAPI_ready","blast_data_ready");
	selectWindow("Result of weka_DAPI_ready");
	run("Select All");
	run("Copy");
	selectWindow("blast_data_ready");
	run("Paste");
	run("Invert");	
	saveAs("tiff", savingDir+File.separator+fileName+"_blastData.tif");

	//save image with segmented LDs over MF-20 enhanced signal
	selectWindow("blastsLD");
	run("Red");
	selectWindow("tubesLD");
	run("Red");
	selectWindow("main");
	setSlice(2);
	run("Duplicate...", "title=MF-20");
	//run("Enhance Contrast", "saturated=0.35");
	//run("Apply LUT");
	run("Merge Channels...", "c1=tubesLD c2=MF-20 c7=blastsLD keep");
	saveAs("tiff", savingDir+File.separator+fileName+"_binLD.tif");
}
