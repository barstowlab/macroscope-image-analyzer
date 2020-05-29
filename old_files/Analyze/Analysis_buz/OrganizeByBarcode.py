# Read Barcodes and Organize Files

from macroscopeUtils4B import RenameAndCopyImageFile
from specutils9 import GenerateFileList


systemBaseDir = "/media/psf/Dropbox for Business/"
# systemBaseDir = "/Users/buz/Dropbox (BarstowLab)/"


sourceBaseDir = systemBaseDir \
+ "/BarstowLab Shared Folder/Data/Macroscope Photos/" \
+ "2015-11-13 - Library Replication for AQDS Reduction Assay/" \
+ "Photographs Just Prior to AQDS Introduction/"


destBaseDir = systemBaseDir \
+ "/BarstowLab Shared Folder/" \
+ "Analysis/2015-11-13 - Library Replication for AQDS Reduction Assay/" \
+ "Photographs Just Prior to AQDS Introduction/"


fileList = GenerateFileList(directory=sourceBaseDir, regex=".*\.jpg")

for file in fileList:
	RenameAndCopyImageFile(destBaseDir, sourceBaseDir + '/' + file, \
	diagnostics=True, preflight=False, organizeIntoDirectories=False)

