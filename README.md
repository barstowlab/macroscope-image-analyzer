# Macroscope Image Analyzer

This program was designed in the Barstow Lab at Cornell University to increase the throughput of well image analysis. The program performs four separate functions: generate barcodes, organize by barcodes, collect data, and analyze data. Generate Barcodes is intended to aid the creation of the barcodes themselves to be placed on 96-well plates to make data collection more robust. Organizing by Barcodes is intended to take a set of images from many plates, and to organize them into independent folders per plate renaming the files to give more information about the images including the time and plate number. Collecting Data analyzes images from a plate and creates an XML and CSV file which contains the data collected at each of the time intervals from which the photos were taken of the plate for reds, greens, blues, and yellows. Finally Analyze Plates takes the XML file made from Collect Data and generates a series of graphs to display the data in a straightforward and easy-to-understand manner.

## Getting Started

The Well Analyzer was written in python3. While python2 should also be backwards compatible, it is not guaranteed to function. To install python via Conda go [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Prerequisites

This program requires zbar, a package that is standard with Windows installations but needs to be installed manually on a Mac. Below is the command needed to install zbar.

```
brew install zbar
```

 If your mac does not have homebrew then run this second command first:

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```

### Installing

The ```requriements``` file in the base directory conatains all of the necessary steps for installing dependencies. To install them all, simply run

```
pip3 install -r requirements.txt
```

### Dependencies

#####Note: These are all automatically installed during the previous Installing step.

OpenCV-python: This package is used for reading in images of the plates into the program.

pyzbar: This package is used for reading information from barcodes.

Pillow: This package is used for image analysis.

numpy: This package is used in mathematical operations.

matplotlib: This package is used in plotting graphs.

## Usage

To run the code, open the base directory containing ```WellAnalyzer.py``` and the ```helper_functions.py``` files. Then simply run:

```python3 WellAnalyzer.py```

The program will then open and prompt the user for which function to run. 

### Generate Barcodes

This step is used for generating barcode numbers to be placed on the plates

#### Inputs:

All of the inputs are prompted by the GUI. First select if the plate is an assay or storage plate which changes how the barcodes are generated. Next choose a prefix and optional postfix to the barcodes. Next include a file name, note the name must end in ```.csv```. Next choose the max rows and columns; this is determined by what sheets are being used for printing the barcodes. Next choose whether to pad columns or rows. Last choose a first barcode number and a final barcode number. If the range exceeds that which could fit in the max rows and columns, the rest will not be used. Lastly select a destination where you want the CSV file to be saved, and click ```Generate Barcodes```

#### Outputs:

In the destination selected by ```Select Destination``` a new CSV file will be made with the barcode text. To make the barcodes themselves, use a barcode font.

### Organize by Barcodes

This step reads through the images in a selected folder, and copies them with new names based on the time the photo was taken and the barcode into a new directory. This step is used for organizing the photographic data based on the barcodes within the images themselves.

#### Inputs:

First select if the plate is an assay plate or a storage plate which changes how the renaming occurs. Next select a destination directory where the new images will be saved. Note that each plate will be categorized into their own folders within the destination directory based on their barcodes. Next select a source directory where all the images are stored that are to be organized. Lastly click ```Organize by Barcodes``` to complete the operation.

#### Outputs:

The program reads through all of the image in the source directory and finds the barcodes in the images. If a barcode is found, for each barcode a new folder will be made in the destination directory based on all but the last two digits of the barcode (those are reserved for checking the correctness of the barcode reading). In the respective folders, the images will be copied over with new names of form XXXX-TTTTT where the XXXX corresponds to the barcode number and the TTTTT to the timestamp of the image. 

### Collect Data

This step is imperative to drawing data from the images. It takes a folder with only a set of images from a single plate that you want to collect data from. Make sure that the well "A01" is in the top left corner and that all of the plates have the exact same placement within the image. Note that sometimes jpg viewing programs will show images oriented differently then the underlying file. Make sure the file itself is oriented properly. The program then reads through the images and produces an XML and a CSV file with the data. The XML file is used for the Analysis step, and the CSV file is a human-readable format. Finally a new jpg is made which shows how the program placed wells within the image.

#### Inputs:

The program prompts the user for three inputs. First a folder location where the a set of images over a course of time from a single plate are stored. Second is an assay plate info CSV which is the mapping file from wells to gene names. Finally it asks for a save location which is the directory to where the collected data will be written to. The assay plate info CSV must follow the structure of the example one given in the ```test_data``` folder mentioned in testing. Once the program is running, it will pull up an image of the first plate in the series and request well definitions. Click first on the far left side of well A01 and then second on the far right side of well H12. It is important to click directly on those positions. Next click test well locations and make sure that the wells are the desired wells. 


#### Outputs:

The program produces three outputs all to the save location. First is an XML file which contains all of the data collected from the images organized by well. Currently the program collects: mean-red, median-red, mean-blue, median-blue, mean-green, median-green, mean-yellow, and median-yellow.  To include more colors, the code itself must be changed. The XML file is used only for the Analyze step and should not be altered. The second output is  CSV file. This file contains the data from the images. The top row is the times. The subsequent row is a single color and a single gene and the intensity values of the colors at those times. This file is human-readable. Lastly a jpg file is made which shows the placement of the wells within the image.


### Analyze Data

This is the last, and most elucidating step. It takes the data collected and generated in previous steps, and displays them in a format where the user can see all 96 wells at once over a time period. It also allows the user to select genes which are of interest, and graph them all together or save them out to a CSV file for other data analysis. The program can handle analyzing different plates and selecting hits across multiple plates.

#### Inputs / Outputs / Usage:

The program is broken into three parts. The left is data input. The middle is data visualization, and the right is data output. First select a source folder in the top left. This should be the folder that contains the XML files made from the Collect Data stage. Next select one of the XML files from the box above the button making it highlighted blue. Next select a color from the list below for graphing. Below that is two boxes where the user can choose to manually select a min and max time for the graphs. Radio buttons below dictate whether user inputted times will be used, or if the max time frame will be used via default bounds. Finally the update graph button calls the functions to generate all of the graphs which will appear in the middle section. The graphs correspond to the wells in the plate with the top left being A01. The bottom of the middle section shows what the current graph is depicting. On the right side, if "View Graph" is selected, then when you click on one of the subgraphs in the middle, a larger version of it will open up. If "Select Hits" is selected, then clicking on a subgraph will highlight it yellow and add it to a list on the right hand side. The list shows the gene name and the plate from which it came. Selecting one of the hits from the list, then clicking remove hit removes it from the list. You can also simply reclick on the hit gene subgraph. Save Hits prompts the user to choose a location to save the hits in a CSV file too. Finally Plot Hits will plot the data of all the hits. It uses the color determined by the current set of graphs in the middle, thus to change color, you must select a new color on the left, click and click update graph updating the central graphs.

## Testing

A folder named ```test_data``` has been included in the home directory. This folder includes 11 images of a single plate over a period of time. The images come from a real assay performed by Buz Barstow. The folder also includes an accompanying CSV file which has a mapping from the wells to the gene names within those wells. This folder is intended to be an example of the data needed for the program to work. Use it your first time running the code so that you can understand the end-to-end process.


## Authors

* **Lucas Jurgensen**
* **Sean Medin**
* **Buz Barstow**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details




