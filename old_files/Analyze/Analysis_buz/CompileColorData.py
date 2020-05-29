# ------------------------------------------------------------------------------------------------ #
# Code to analyze color data

import re
import pdb
from specutils11 import GenerateFileList
from macroscopeUtils7 import gWellList, ExportListForXML, indent, \
ImportUpdatedSudokuCatalogFromCSV, ExperimentWell, ImportWellColorsDictFromXML, \
WriteCondensedColorDataToXML
import xml.etree.ElementTree as ET

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# File information

updatedSudokuCatalogFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-11-14-19 - Shewanella Sudoku Library Replication for AQDS Reduction/Manually Curated Updated Catalog/Manually Curated Updated Catalog.csv'

imageBaseDir = "/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Image Data/"

wellColorsXMLFileName = "/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Condensed Data/wellColors.xml"

imageDirs = \
['AQDS01T', 'AQDS02T', 'AQDS03T', 'AQDS04T', 'AQDS05T', 'AQDS06T', 'AQDS07T', 'AQDS08T', 'AQDS09T', 'AQDS10T', 'AQDS11T', 'AQDS12T', 'AQDS13T', 'AQDS14T', 'AQDS15T', 'AQDS16T', 'AQDS17T', 'AQDS18T', 'AQDS19T', 'AQDS20T', 'AQDS21T', 'AQDS22T', 'AQDS23T', 'AQDS24T', 'AQDS25T', 'AQDS26T', 'AQDS27T', 'AQDS28T', 'AQDS29T', 'AQDS30T', 'AQDS31T', 'AQDS32T', 'AQDS33T', 'AQDS34T', 'AQDS35T', 'AQDS36T', 'AQDS37T', 'AQDS38T', 'AQDS39T', 'AQDS40T', 'AQDS41T', 'AQDS42T', 'AQDS43T', 'AQDS44T', 'AQDS45T', 'AQDS47BT', 'AQDS47T', 'AQDS48T', 'AQDS49T', 'AQDS50T', 'AQDS51T', 'AQDS52T', 'AQDS53BT', 'AQDS53T', 'AQDS54T', 'AQDS55T', 'AQDS56T', 'AQDS57T', 'AQDS58BT', 'AQDS58T', 'AQDS59T', 'AQDS60BT', 'AQDS60T', 'AQDS61BT', 'AQDS61T', 'AQDS62T', 'AQDS63BT', 'AQDS63T', 'AQDS64BT', 'AQDS64T', 'AQDS65T', 'AQDS66T', 'AQDS67T', 'AQDS68T', 'AQDS69T', 'AQDS70BT', 'AQDS70T', 'AQDS71T', 'AQDS72T', 'AQDS73T', 'AQDS74T', 'AQDS75T', 'AQDS76T', 'AQDS77T', 'AQDS78T']


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
wellColorsPostfix = r"_colors"
wellColorsFileExt = r".xml"
wellColorsFileRegex = r".*" + wellColorsPostfix + wellColorsFileExt
plateIDRe = re.compile(r"AQDS(\d+[v\d+]*[A-H]*)T")
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the updated, manually curated Sudoku catalog
sudokuCatalog = ImportUpdatedSudokuCatalogFromCSV(updatedSudokuCatalogFileName)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Go through the image directories and import the well color XML files

collectedWellColorsDict = {}

for imageDir in imageDirs:
	
	workingDir = imageBaseDir + '/' + imageDir
	wellColorsFileList = GenerateFileList(directory=workingDir, regex=wellColorsFileRegex)
	
	plateIDMatch = plateIDRe.match(imageDir)
					
	if plateIDMatch != None:
		plateIDFromImageDir = plateIDMatch.group(1)
	else:
		pdb.set_trace()
	
	# Make an entry in the collectedWellColorsDict for the plate
	# Working under the assumption that all images for a particular plate are organized into a
	# single directory.
	collectedWellColorsDict[plateIDFromImageDir] = {}
	currentEntry = collectedWellColorsDict[plateIDFromImageDir]
	
	for well in gWellList:
		currentEntry[well] = ExperimentWell(plateIDFromImageDir, well)
		
		currentEntry[well].dataDict['meanRed'] = []
		currentEntry[well].dataDict['meanGreen'] = []
		currentEntry[well].dataDict['meanBlue'] = []
		currentEntry[well].dataDict['medianRed'] = []
		currentEntry[well].dataDict['medianGreen'] = []
		currentEntry[well].dataDict['medianBlue'] = []
		currentEntry[well].dataDict['time'] = []
		
	i = 0
	while i < len(wellColorsFileList):
		
		wellColorsFileName = workingDir + '/' + wellColorsFileList[i]
		
		[wellPositionsDict, meanRedsDict, meanGreensDict, meanBluesDict, medianRedsDict, \
		medianGreensDict, medianBluesDict, wellSepX, wellSepY, originatingFile, timeStamp, \
		plateID] = \
		ImportWellColorsDictFromXML(wellColorsFileName)
		
		# This code normalizes the plate ID format and takes care of leading zeros
		# plateIDForCol = plate ID for collected data
		if plateID == plateIDFromImageDir:
			plateIDForCol = plateID
		elif str(int(plateID)).zfill(2) == str(int(plateIDFromImageDir)).zfill(2):
			plateIDForCol = str(int(plateID)).zfill(2)
		else:
			pdb.set_trace()
		
		print("Plate ID: " + str(plateIDForCol))
		
		for well in gWellList:
			currentEntry[well].dataDict['meanRed'].append(meanRedsDict[well])
			currentEntry[well].dataDict['meanGreen'].append(meanGreensDict[well])
			currentEntry[well].dataDict['meanBlue'].append(meanBluesDict[well])
			currentEntry[well].dataDict['medianRed'].append(medianRedsDict[well])
			currentEntry[well].dataDict['medianGreen'].append(medianGreensDict[well])
			currentEntry[well].dataDict['medianBlue'].append(medianBluesDict[well])
			currentEntry[well].dataDict['time'].append(float(timeStamp))
			
		i += 1


# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Update the collectedWellColorsDict with information from the manually curated Sudoku catalog

plateIDs = collectedWellColorsDict.keys()
for plate in plateIDs:
	currentPlate = collectedWellColorsDict[plate]
	
	if plate in sudokuCatalog.keys():
		plateNumberKey = plate
	elif str(int(plate)).zfill(1) in sudokuCatalog.keys():
		plateNumberKey = str(int(plate)).zfill(1)
	else:
		pdb.set_trace()
	
	currentSudokuCatalogEntry = sudokuCatalog[plateNumberKey]
	
	for well in gWellList:
		currentPlate[well].transposonCoord = currentSudokuCatalogEntry[well][0]
		currentPlate[well].featureName = currentSudokuCatalogEntry[well][1]
		currentPlate[well].occupied = currentSudokuCatalogEntry[well][2]
		currentPlate[well].occupancyTest = currentSudokuCatalogEntry[well][3]

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Write out the collectedWellColorsDict as an XML file

experimentalWellsLinearArray = []
for plate in collectedWellColorsDict.keys():
	for wellKey in collectedWellColorsDict[plate].keys():
		experimentalWellsLinearArray.append(collectedWellColorsDict[plate][wellKey])

WriteCondensedColorDataToXML(wellColorsXMLFileName, experimentalWellsLinearArray)
# ------------------------------------------------------------------------------------------------ #
