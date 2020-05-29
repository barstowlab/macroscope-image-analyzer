# ------------------------------------------------------------------------------------------------ #
# FindWellPositionsImproved.py
# Author: Sean Medin
# The purpose of this code is to identify wells in a series of plate images and write their positions on the plate and
# their color change over time in an xml file
#
# Note: borrows liberally from archive/FindWellPositionsAndMeasureColorsOnAssayPlate.py
# Created: November 1st 2019 by Sean Medin
# Last Updated: 11-17-19 by Sean Medin
# ------------------------------------------------------------------------------------------------ #

from utils.input import get_input
import sys

import csv
from utils.find_wells_helper import *

# Some configuration information
new_pixel_width = 800

scaledMarkFilePostfix = r"_scaled_marked.png"
wellColorsPostfix = r"_colors.xml"

# Parses input file
argv = sys.argv
inputParameters = ["img_info", "folder_location", "save_location", "reoriented_plates_location"]
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../inputs/' + proj_name + '/'
file_name = file_intro + 'FindWellPositionsImproved.inp'
inputParameterValues = get_input(file_name, inputParameters)

# This argument gives the path to a csv file with information about the plates to analyze
# The csv file has three columns: "PlateFolder", "WellInfo", and "ImageType"
# "PlateFolder": each entry in this column is a folder where we have a time series of a single plate; the folder should
# be written as relative to the "folder_location" variable seen below
# "WellInfo": is a relative path to a csv file containing information about what's in each well. Wells should be in
# Capital letter = row, number = column format (ex: A7)
# "ImageType" gives the type of images found in PlateFolder (ex: ".jpg", the "." is necessary)
img_info = inputParameterValues["img_info"]

# the location to save the generated xml file and sample positioning of the plates
save_location = inputParameterValues["save_location"]
print(save_location)

# location where to get the "PlateFolder"s
folder_location = inputParameterValues["folder_location"]

# If this is "None", then no need to reorient plates
reoriented_plates_location = inputParameterValues["reoriented_plates_location"]

# get all plate paths and well information
well_info = {}
plate_names = []
img_types = {}
with open(img_info, newline='') as csvfile:
    plate_reader = csv.reader(csvfile, delimiter=',')
    for plate in plate_reader:
        if plate[0] == "PlateFolder":
            continue
        plate_names.append(plate[0])
        if plate[1] != "None":
            well_info[plate[0]] = parse_well_info(plate[1])
        img_types[plate[0]] = plate[2]

# makes sure all plates are facing the same direction (and resaves them all in another folder)
if reoriented_plates_location != "None":
    orient_images(plate_names, img_types, folder_location, reoriented_plates_location)
    folder_location = reoriented_plates_location

# Finds well positions on all plates and saves scaled images, well coordinates, and colors inside the wells
print("plate_names", plate_names, "\n")
print("img_types", img_types, "\n")
print("folder_location", folder_location, "\n")
print("well_info", well_info, "\n")
print("new_pixel_width", new_pixel_width, "\n")
print("save_location", save_location, "\n")
print("scaledMarkFilePostfix", scaledMarkFilePostfix, "\n")
print("wellColorsPostfix", wellColorsPostfix, "\n")



find_wells(plate_names, img_types, folder_location, well_info, new_pixel_width, save_location,
           scaledMarkFilePostfix, wellColorsPostfix)
