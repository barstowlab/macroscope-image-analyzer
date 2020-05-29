# ------------------------------------------------------------------------------------------------ #
# FindWellPositionsOnAssayPlate.py
# Author: Buz Barstow
# The purpose of this code is to identify wells in a series of images and write out their 
# coordinates to a text file
# Created: November 3rd 2015
# Last Updated: December 2nd 2015
# ------------------------------------------------------------------------------------------------ #

from specutils11 import GenerateFileList
from macroscopeUtils7 import ensure_dir, vecLength, \
FindWellPositionsIn96WellAssayPlateImage, MeasureAverageWellColorsIn96WellPlate, ProcessTimeStamp, \
ProcessFileNameExtension, ImportWellPositionsDictFromXML, WriteWellPositionsDictForImageAsXML, \
WriteAverageColorsDictForImageAsXML

import shutil
import re
import pdb


# ------------------------------------------------------------------------------------------------ #
# Image recognition parameter data
blueRegistrationColor=(100.0, 140.0, 200.)
pinkRegistrationColor=(212.0, 120.0, 205.)


# Plate geometry in pixels based upon template photograph
blueOriginInit = array([4526, 1459])
pinkOriginInit = array([4529, 2687])
wellRadius = 278/2.
wellSepX = 304
wellSepY = 304
originToA12Init = array([-222, -612])

A1CoordInit = array([918,732])

blueToPinkVector = pinkOriginInit - blueOriginInit
blueToPinkLength = vecLength(blueToPinkVector)

wellSepXInitRatio = wellSepX/blueToPinkLength
wellSepYInitRatio = wellSepY/blueToPinkLength
originToA12InitRatio = originToA12Init/blueToPinkLength
wellRadiusRatio = wellRadius/blueToPinkLength
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Some configuration information
writeWellPositions = True
moveImages = True
moveScaledMarkedImages = True

scaledMarkFilePostfix = r"_scaled_marked"
scaledMarkedFileExt = r".png"
wellCoordsPostfix = r"_well_coords"
wellCoordsFileExt = r".xml"
wellColorsPostfix = r"_colors"
wellColorsFileExt = r".xml"
wellCoordsFileRe = re.compile(r"(.*)"+wellCoordsPostfix)
plateIDRe = re.compile(r"AQDS(\d+[v\d+]*[A-H]*)T")
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# File information
imageBaseDir = "/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Image Data/"

# imageDirs = \
# ['AQDS01T', 'AQDS02T', 'AQDS03T', 'AQDS04T', 'AQDS05T', 'AQDS06T', 'AQDS07T', 'AQDS08T', 'AQDS09T', 'AQDS10T', 'AQDS11T', 'AQDS12T', 'AQDS13T', 'AQDS14T', 'AQDS15T', 'AQDS16T', 'AQDS17T', 'AQDS18T', 'AQDS19T', 'AQDS20T', 'AQDS21T', 'AQDS22T', 'AQDS23T', 'AQDS24T', 'AQDS25T', 'AQDS26T', 'AQDS27T', 'AQDS28T', 'AQDS29T', 'AQDS30T', 'AQDS31T', 'AQDS32T', 'AQDS33T', 'AQDS34T', 'AQDS35T', 'AQDS36T', 'AQDS37T', 'AQDS38T', 'AQDS39T', 'AQDS40T', 'AQDS41T', 'AQDS42T', 'AQDS43T', 'AQDS44T', 'AQDS45T', 'AQDS47BT', 'AQDS47T', 'AQDS48T', 'AQDS49T', 'AQDS50T', 'AQDS51T', 'AQDS52T', 'AQDS53BT', 'AQDS53T', 'AQDS54T', 'AQDS55T', 'AQDS56T', 'AQDS57T', 'AQDS58BT', 'AQDS58T', 'AQDS59T', 'AQDS60BT', 'AQDS60T', 'AQDS61BT', 'AQDS61T', 'AQDS62T', 'AQDS63BT', 'AQDS63T', 'AQDS64BT', 'AQDS64T', 'AQDS65T', 'AQDS66T', 'AQDS67T', 'AQDS68T', 'AQDS69T', 'AQDS70BT', 'AQDS70T', 'AQDS71T', 'AQDS72T', 'AQDS73T', 'AQDS74T', 'AQDS75T', 'AQDS76T', 'AQDS77T', 'AQDS78T']

imageDirs = \
['AQDS58BT', 'AQDS58T']


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Figure out the well positions and write out. 

fittedWellPositionsDict = {}

for imageDir in imageDirs:
	
	fittedWellPositionsDict[imageDir] = {}
	
	workingDir = imageBaseDir + '/' + imageDir
	fileList = GenerateFileList(directory=workingDir, regex=".*.jpg")
	fileList = sorted(fileList)

	for fileName in fileList:
		print(fileName)
		
		[fileNameNoExt, extension] = ProcessFileNameExtension(fileName)
		
		fullFileName = workingDir + '/' + fileName
		annotatedWellPositionsFigFileName = workingDir + '/' + fileNameNoExt \
		+ scaledMarkFilePostfix + scaledMarkedFileExt
		
		fittedWellPositionsFileName = workingDir + '/' + fileNameNoExt + wellCoordsPostfix \
		+ wellCoordsFileExt
		
		timeStamp = ProcessTimeStamp(fileName)
	
		[fittedWellPositionsDictForImage, wellSepXFitted, wellSepYFitted] = \
		FindWellPositionsIn96WellAssayPlateImage(fullFileName, blueRegistrationColor, \
		pinkRegistrationColor, wellSepXInitRatio, wellSepYInitRatio, originToA12InitRatio, \
		wellRadiusRatio, blueOriginInit, pinkOriginInit, diagnostics=False, scale=0.2, \
		useRingPositionsToGuessWellGrid=False, outputFigFileName=annotatedWellPositionsFigFileName)
	
		fittedWellPositionsDict[imageDir][fileNameNoExt] = fittedWellPositionsDictForImage
		
		WriteWellPositionsDictForImageAsXML(fittedWellPositionsFileName, \
		fittedWellPositionsDictForImage, timeStamp, wellSepXFitted, wellSepYFitted, fullFileName)
	
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Figure out the well colors and write out, along with data on well positions. 

wellColorsDict = {}

for imageDir in imageDirs:
	
	wellColorsDict[imageDir] = {}
	
	workingDir = imageBaseDir + '/' + imageDir
	imageFileList = GenerateFileList(directory=workingDir, regex=".*.jpg")
	imageFileList = sorted(imageFileList)
	
	wellPositionsFileList = GenerateFileList(directory=workingDir, regex=".*_well_coords.xml")
	wellPositionsFileList = sorted(wellPositionsFileList)
	
	if len(wellPositionsFileList) == len(imageFileList):
		i = 0
		while i < len(imageFileList):
			imageFileName = imageFileList[i]
			wellPositionsFileName = wellPositionsFileList[i]
			
			[imageFileNameNoExt, imageFileExtension] = ProcessFileNameExtension(imageFileName)
			[wellPositionsFileNameNoExt, wellPositionsFileNameExtension] = \
			ProcessFileNameExtension(wellPositionsFileName)
			
			wellCoordsFileMatch = wellCoordsFileRe.match(wellPositionsFileNameNoExt)
			
			if wellCoordsFileMatch != None:
				if wellCoordsFileMatch.group(1) == imageFileNameNoExt:
					[wellPositionsDict, wellSepX, wellSepY, originatingFile, timeStamp] = \
					ImportWellPositionsDictFromXML(workingDir + '/' \
					+  wellPositionsFileName)
					
					radiusForColorMeasurement = mean([wellSepX, wellSepY])*0.2
					
					wellAverageColorsDict = MeasureAverageWellColorsIn96WellPlate(\
					workingDir + '/' + imageFileName, \
					wellPositionsDict, radiusForColorMeasurement)
					
					plateIDMatch = plateIDRe.match(imageFileNameNoExt)
					
					if plateIDMatch != None:
						plateID = plateIDMatch.group(1)
					else:
						pdb.set_trace()
					
					wellAverageColorsFileName = workingDir + '/' + imageFileNameNoExt + \
					wellColorsPostfix + wellColorsFileExt
					
					# put in the code to write out the average colors here!
					WriteAverageColorsDictForImageAsXML(wellAverageColorsFileName, \
					wellAverageColorsDict, wellPositionsDict, wellSepX, wellSepY, originatingFile, \
					timeStamp, plateID)
					
					
				else:
					print('Image file name and well positions file name do not match!')
					pdb.set_trace()
			else:
				print('Format of well position file name does not match template!')
				pdb.set_trace()
			
			i += 1
# ------------------------------------------------------------------------------------------------ #



