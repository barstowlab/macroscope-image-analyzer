# Requires python2

# Read Barcodes and Organize Files

from macroscopeUtils4B import RenameAndCopyImageFile
from specutils9 import GenerateFileList
from tkinter import filedialog
import tkinter as tk

#TODO - add capability to select files and not hardcode
root = tk.Tk()
root.withdraw()
print("select source directory")
sourceBaseDir = filedialog.askdirectory()
print("select destination directory")
destBaseDir = filedialog.askdirectory()



# systemBaseDir = "/media/psf/Dropbox for Business/"
# # systemBaseDir = "/Users/buz/Dropbox (BarstowLab)/"


# sourceBaseDir = systemBaseDir \
# + "/BarstowLab Shared Folder/Data/Macroscope Photos/2015-11-13 - Library Replication for AQDS Reduction Assay/Photographs Prior to AQDS Intro (47-78)/"

# destBaseDir = systemBaseDir \
# + "/BarstowLab Shared Folder/" \
# + "Analysis/2015-11-13 - Library Replication for AQDS Reduction Assay/" \
# + "Photographs Prior to AQDS Intro (47-78)/"


fileList = GenerateFileList(directory=sourceBaseDir, regex=".*\.jpg")

for file in fileList:
	RenameAndCopyImageFile(destBaseDir, sourceBaseDir + '/' + file, \
	diagnostics=True, preflight=False, organizeIntoDirectories=False)
