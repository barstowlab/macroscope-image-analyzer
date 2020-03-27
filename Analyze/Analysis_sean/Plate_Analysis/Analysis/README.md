# plate_analysis

1. Run FindWellPositionsImporved.py with target directory. In the GUI,  
    (a) Enter correct number of rows and columns of the plate.  
    (b) Click the leftmost point on the top-left well, there is a colored dot left to indicate the selection.  
        &nbsp; &nbsp; &nbsp; To reselect the point, click the [Undo Well Define] button.  
    (c) Click the rightmost point on the bottom-right well, same as step (b)  
    (d) Click the [Test Well Location] button, black circles will show on the graph to indicate the found well positions.  
        &nbsp; &nbsp; &nbsp; To reset the found well positions, click the [Undo Well Define] button, and repeat the step (a) - (d)  
        &nbsp; &nbsp; &nbsp; To observe the RGB value of specific well, click the corresponding black circle.   
    (e) After verifying the found well positions, click the [Approve Well Locations] button to process all the images in the directory and output the image information as xml and csv.
        
2. To plot the image information, run PlotColorTimeSeries.py with target directory. In the GUI,  
    (a) Select target gene(s) in the upper box.  
    (b) Select information wanted for the plot in the bottom box.  
    (c) Enter [Start Time] and [End Time] in the corresponding field. If no values is entered, the plot contain the whole process.  
    (d) Click the [Plot Stuff] button, the plot will pop up. To save or modify the plot, use the button shown in the bottom of the figure window.  


