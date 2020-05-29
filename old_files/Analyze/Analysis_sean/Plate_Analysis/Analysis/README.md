# plate_analysis

##Overview

This code assists you in extracting color data from the wells of a microplate and
with plotting/analyzing the data. To read the plate, you run the python file
**FindWellPositionsImproved.py** and to analyze the resulting data you run the
 python files **PlotColorTimeSeries** and **PlotLollipopGraph.py**

##Read Plate

###What You Need to Run the Code

To extract data from the microplate, run FindWellPositions.py with a single
argument specifying the name of a folder where you have an input file called
**FindWellPositions.inp**. This folder be located at the location
**../inputs/** relative to the directory from which you run this function.

**FindWellPositions.inp** must define four strings: **img_info**, 
**save_location**, **folder_location**, and **reoriented_plates_location**. An
example line for each and what they mean is shown below:

**img_info = ../inputs/ElectronUptake/plates.csv**

The **plates** csv file.

**save_location = ../inputs/ElectronUptake**

The folder location where you want to save the output files.

**folder_location = ../../../../../Photographic Data/Farshid**

The folder where you find folders containing images of a plate taken at
different times.

**reoriented_plates_location =None**

Ignore this. This was an attempt to use code to make sure all the plates
were oriented in the same direction. There's probably a better way of doing this
anyway. I'll delete this completely once I get the chance to check if removing
it would break anything.

###CSV File Inputs

**plates**

This is a csv file that has three columns. The headings are: "PlateFolder",
"WellInfo", and "ImageType". The "PlateFolder" tells you the folder within
folder_location where you can find the plate images. "WellInfo" points you to
the **well_info** csv that tells you what genes are within each well. 
"ImageType" just says the format of the image; for example: '.jpg' 
(the '.' is important)

**well_info**

Need columns labelled, "Test Well" and "Gene Name" (in that order). The "Test
Well" column identifies the well with the row being specified by a capital letter
combined with a two digit number. For example, the top left row would be "A01".

"Gene Name" is just a a string that identifies the mutant in the well. The same
mutant should be labelled the same even if it shows up in multiple wells.

###Using the GUI

Whe the GUI pops up

(a) Enter correct number of rows and columns of the plate.  
(b) Click the leftmost point on the top-left well, there is a colored dot left to indicate the selection.  
    &nbsp; &nbsp; &nbsp; To reselect the point, click the [Undo Well Define] button.  
(c) Click the rightmost point on the bottom-right well, same as step (b)  
(d) Click the [Test Well Location] button, black circles will show on the graph to indicate the found well positions.  
    &nbsp; &nbsp; &nbsp; To reset the found well positions, click the [Undo Well Define] button, and repeat the step (a) - (d)  
    &nbsp; &nbsp; &nbsp; To observe the RGB value of specific well, click the corresponding black circle.   
(e) After verifying the found well positions, click the [Approve Well Locations] button to process all the images in the directory and output the image information as xml and csv.

###Outputs

The code will ouptut an xml file and a csv file and put them in the
specified save folder location.

For the csv file: column A (starting
at row 2) lists the gene names. The rows alternate between the mean and the
standard deviation "yellow" color measurement for the genes. The columns (starting
at column B) shows the value at that time with the first row showing the time
stamp at which this value was measured.

The xml file has an entry for each well that lists the mean and median color
measurements (for blue, green, red, and "yellow") as well as the well it was
found and all the information listed in the **well_info** csv for the well. There
is other stuff, but the only important things are the color measurements, the
well location, and the gene name specified in **well_info**.


##Analyze Plate

###PlotColorTimeSeries.py

This function brings up a GUI that will allow you to plot the color
change of multiple mutants of interest on the same plot and for a desired amount
of time.

####Running the Code

Run PlotColorTimeSeries.py with a single
argument specifying the name of a folder where you have an input file called
**PlotColorTimeSeries.inp**. This folder be located at the location
**../inputs/** relative to the directory from which you run this function.

**PlotColorTimeSeries.inp** should have a single line specifying the location of
an xml file generated from **FindWellPositionsImproved.py**. For example:

color_info_file = ../inputs/ElectronUptake/ElectronUptake_colors.xml

####Using the GUI

(a) Select target gene(s) in the upper box.  
(b) Select information wanted for the plot in the bottom box.  
(c) Enter [Start Time] and [End Time] in the corresponding field. If no values is entered, the plot contain the whole time series.  
(d) Click the [Plot Stuff] button, the plot will pop up. To save or modify the plot, use the button shown in the bottom of the figure window.  

###PlotLollipopGraph.py

This file allows you to plot "lollipop" time series plots side by side
for the gene of interest. You input the times of interest in the input file
and running the file brings up a GUI that allows you to select the genes you
want to plot. You must select at least two genes.

####Running the Code

Run PlotLollipopGraph.py with a single
argument specifying the name of a folder where you have an input file called
**PlotLollipopGraph.inp**. This folder be located at the location
**../inputs/** relative to the directory from which you run this function.

**PlotLollipopGraph.inp** must define three strings: "color_info_file", 
"times", and "save_location". For example, your file could look like

**color_info_file = ../inputs/ElectronUptake/ElectronUptake_colors.xml**

Same input as for PlotColorTimeSeries.py

**times = 0,100,200,300**

Comma separated list of time in minutes

**save_location = ../inputs/ElectronUptake/ElectronUptakeLollipop.png**

Output file

####Using the GUI
Just select the genes in the drop down and press "Make Lollipop Graphs"



