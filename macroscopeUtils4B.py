# ------------------------------------------------------------------------------------------------ #
# A branched version of the macroscope utils
# Buz Barstow
# 17th November 2015
# Designed for compatibility with image recognition scripts requiring SimpleCV and python 2.
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #


import sympy
import re

# ------------------------------------------------------------------------------------------------ #
gOffset = sympy.N(-6.213577130399991E10)
gConstant = sympy.N(1.000002747119042E-7)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FileNameTimeStampToSecond(fileNameTimeStamp):
	from sympy import N
	
	seconds = sympy.N(fileNameTimeStamp)*gConstant + gOffset
	
	return seconds
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ConvertFileNameToSeconds(fileName):
	import re
	import os
	
	baseName = os.path.basename(fileName)
	
	timeStampRegex = re.compile("(\w*)\-(\d+)")
	timeStampSearch = timeStampRegex.search(baseName)

	if timeStampSearch != None:
		timeStamp = timeStampSearch.group(2)
	else:
		ex = TimeStampFailure('Time Stamp Failure')
		raise ex
	
	
	timeInSeconds = FileNameTimeStampToSecond(timeStamp)
	
	return timeInSeconds
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
class TimeStampFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class FileExtensionFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class SelfVectorFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class ExperimentCodeFailure(Exception):
	pass
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def ProcessTimeStamp(fileName):
	import re
	import os
	
	baseName = os.path.basename(fileName)
	
	timeStampRegex = re.compile("(\w*)\-(\d+)")
	timeStampSearch = timeStampRegex.search(baseName)

	if timeStampSearch != None:
		timeStamp = timeStampSearch.group(2)
	else:
		ex = TimeStampFailure('Time Stamp Failure')
		raise ex
		
	return timeStamp
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ProcessFileNameExtension(fileName):
	import re
	import os
	
	baseName = os.path.basename(fileName)
	
	extRegex = re.compile("\.(\w+)")
	extSearch = extRegex.search(baseName)

	if extSearch != None:
		extension = extSearch.group(1)
	else:
		ex = FileExtensionFailure('File Extension Failure')
		raise ex
		
	return extension
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)


def RenameAndCopyImageFile(destBaseDir, fileName, diagnostics=False, preflight=False, \
organizeIntoDirectories=True):
	import shutil
	import pdb
	
	[plateID, timeStamp, fileExtension] = ProcessPlateBarcode(fileName)
	
	try:
		if fileExtension[0] != '.':
			newFileName = plateID + '-' + timeStamp + '.' + fileExtension
		else:
			newFileName = plateID + '-' + timeStamp + fileExtension
	except:
		pdb.set_trace()
	
	
	if diagnostics == True:
		print fileName + ' --> ' + newFileName
	
	if preflight == False:
		if organizeIntoDirectories == True:
			ensure_dir(destBaseDir + '/' + plateID + '/' + newFileName)
			shutil.copy(fileName, destBaseDir + '/' + plateID + '/' + newFileName)
		else:
			ensure_dir(destBaseDir + '/' + newFileName)
			shutil.copy(fileName, destBaseDir + '/' + newFileName)
	
# ------------------------------------------------------------------------------------------------ #
	




# ------------------------------------------------------------------------------------------------ #
def IdentifyRowAndColumnNumber(xCoord, yCoord, xEdge1, yEdge1, dX, dY):
	xCoordR = xCoord - xEdge1
	yCoordR = yCoord - yEdge1
	
	yCoordScaled = yCoordR/dY
	yNumber = round(yCoordScaled, 0)
	
	xCoordScaled = xCoordR/dX
	xNumber = round(xCoordScaled, 0)
	
	xNumberAdj = int(abs(xNumber)+1)
	yNumberAdj = int(abs(yNumber)+1)
	
	return [xNumberAdj, yNumberAdj]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
RowNumberToLetterConversionDict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', \
'7':'G', '8':'H'}

def ConvertRowAndColumnNumberToWellID(queryRowNumber, queryColumnNumber):
	queryRowLetter = RowNumberToLetterConversionDict[str(queryRowNumber)]

	wellID = queryRowLetter + str(queryColumnNumber)
	return wellID
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GenerateCheckSum(text):
	import numpy
	
	textArray = []
	ordArray = []
	subtractionArray = []
	
	for letter in text:
		textArray.append(letter)
		ordArray.append(ord(letter))
		subtractionArray.append(ord(letter) - 55)
		
	j = 0
	checksum = 0
	subtractionArrayCopy = numpy.zeros(8, numpy.int)

	while j < min([len(subtractionArray), 8]):
		subtractionArrayCopy[j] = subtractionArray[j]
		j += 1
	
	j = 0
	while j < min([len(subtractionArray), 8]):
		if subtractionArrayCopy[j] < 10:
			subtractionArrayCopy[j] += 7
		checksum += subtractionArrayCopy[j]
		j += 1
	
	checksumMod = checksum%36
	if checksumMod < 10:
		checksumMod -= 7
	checkcar = chr(checksumMod+55)
	
	#return [textArray, ordArray, subtractionArray, subtractionArrayCopy, checksum, \
	#checksumMod, checkcar]
	
	return checkcar
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ProcessPlateBarcode(fileName):
	from SimpleCV import Image
	import pdb
	
	timeStamp = ProcessTimeStamp(fileName)
	[fileNameNoExt, fileExtension] = ProcessFileNameExtension(fileName)
	
# 	pdb.set_trace()
	
	img = Image(fileName)
	
	barcode = img.findBarcode()
	
	if barcode != None:
		if len(barcode[0].data) >= 2:	
			plateID = barcode[0].data[0:-1]
			checkSum = barcode[0].data[-1]	
			calculatedChecksum = GenerateCheckSum(plateID)
	
			if checkSum is not calculatedChecksum:
				print "Error in barcode check sum. File: " + fileName
		else:
			plateID = 'UNKNOWN'
	else:
		plateID = 'UNKNOWN'
	
	
	return [plateID, timeStamp, fileExtension]
		
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def PlotPlateGraphs(CollectedWellTimeCourseDict, plateODDataDict, dirList, baseDir, \
prefix='SOGRT|CTRL'):
	import matplotlib.pyplot as plt
	import pdb

	scale=1.4
	dpi=80
	xWidth = 11*scale
	yWidth = 8.5*scale


	for dir in dirList:
		WellTimeCourseDict = CollectedWellTimeCourseDict[dir]
		print dir
		
		pdb.set_trace()
	
		plateKeyForplateODDataDict = NormalizePlateIDFormat(dir, prefix=prefix,zfillNumber=3)
	
	
		# Generate a time course for each well:
		fig, axes = plt.subplots(nrows=8, ncols=12, figsize=(xWidth,yWidth),dpi=dpi)

		fig.subplots_adjust(hspace=0.15) # Vertical spacing between plots
		fig.subplots_adjust(wspace=0.25) # Horizontal spacing
		fig.subplots_adjust(left=0.05)
		fig.subplots_adjust(right=0.95)
		fig.subplots_adjust(top=0.95)
		fig.subplots_adjust(bottom=0.05)


		xCenter = xWidth*dpi/2
		yCenter = yWidth*0.8
		fig.text(0.5,0.975,plateKeyForplateODDataDict,style='oblique')

		yTicksArray = arange(0,300,50)
		yTickLabelsCol1 = arange(0,300,50)
		yTickLabelsCol2to12 = ['','','','','','','','']
		xTicksArray = [0,350,700]
		xTickLabelsRows1to7 = ['','','']

		for well in gWellList:
			[row, col] = ConvertWellIDToRowAndColumnNumber(well)
			plotNumber = ConvertWellIDToWellNumber(well)
			print str(row) + ', ' + str(col) + ', ' + str(plotNumber)
			print well
		
			subplot(8, 12, plotNumber+1)
	
			timesArray = array(WellTimeCourseDict[well].timeValues)/1E9
	
			if plateODDataDict[plateKeyForplateODDataDict][well] > 0.2:
				plot(timesArray, WellTimeCourseDict[well].meanRed, 'r')
				plot(timesArray, WellTimeCourseDict[well].meanGreen, 'g')
				plot(timesArray, WellTimeCourseDict[well].meanBlue, 'b')
			else:
				plot(timesArray, WellTimeCourseDict[well].meanRed, 'k')
				plot(timesArray, WellTimeCourseDict[well].meanGreen, 'k')
				plot(timesArray, WellTimeCourseDict[well].meanBlue, 'k')

	
			xlim([0,700])
			ylim([0,250])
			grid()
	
			if row == 1:
				title(col, fontsize=12)
				xticks(xTicksArray, xTickLabelsRows1to7, fontsize=10)
			elif row == 8:
				xticks(xTicksArray, fontsize=10)
			else:
				xticks(xTicksArray, xTickLabelsRows1to7, fontsize=10)
	
			if col == 1:
				ylabel(gRows[row-1], fontsize=12)
				yticks(yTickLabelsCol1, fontsize=10)
			else:
				yticks(yTicksArray, yTickLabelsCol2to12, fontsize=10)
	
		# 	plot(WellTimeCourseDict[well].timeValues, WellTimeCourseDict[well].medianRed, '--r')
		# 	plot(WellTimeCourseDict[well].timeValues, WellTimeCourseDict[well].medianGreen, '--g')
		# 	plot(WellTimeCourseDict[well].timeValues, WellTimeCourseDict[well].medianBlue, '--b')
		fig.savefig(baseDir+'/'+dir+'/'+dir+'.pdf',format='pdf')
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
def ExportPlateTimeCourses(CollectedWellTimeCourseDict, dirList, baseDir, prefix='SOGRT|CTRL'):
	import pdb
	from numpy import array, float, ones
	from vectorOutput import generateOutputMatrixWithHeaders, writeOutputMatrix

	for dir in dirList:
		WellTimeCourseDict = CollectedWellTimeCourseDict[dir]
		print dir
				
		headerList = []
		vectorList = []
		
		for well in gWellList:
			[row, col] = ConvertWellIDToRowAndColumnNumber(well)
			plotNumber = ConvertWellIDToWellNumber(well)
			
			if WellTimeCourseDict[well].hasODData == True:
				opticalDensities = \
				ones(len(WellTimeCourseDict[well].timeValues), float)*WellTimeCourseDict[well].OD
				
			if len(WellTimeCourseDict[well].genotype) > 0:
				genotypes = []
				i = 0
				while i < len(WellTimeCourseDict[well].timeValues):
					genotypes.append(WellTimeCourseDict[well].genotype)
					i += 1
			
			
			headerList.append(well+'_time')
			headerList.append(well+'_r')
			headerList.append(well+'_g')
			headerList.append(well+'_b')
			
			if WellTimeCourseDict[well].hasODData == True:
				headerList.append(well+'_od')
			
			#pdb.set_trace()
			if len(WellTimeCourseDict[well].genotype) > 0:
				headerList.append(well+'_gene')
			
			timesArray = list(array(WellTimeCourseDict[well].timeValues))
			
			vectorList.append(timesArray)
			vectorList.append(list(WellTimeCourseDict[well].meanRed))
			vectorList.append(list(WellTimeCourseDict[well].meanGreen))
			vectorList.append(list(WellTimeCourseDict[well].meanBlue))
			
			if WellTimeCourseDict[well].hasODData == True:
				vectorList.append(list(opticalDensities))
			
			if len(WellTimeCourseDict[well].genotype) > 0:
				vectorList.append(genotypes)
		
		#pdb.set_trace()
		oMatrix = generateOutputMatrixWithHeaders(vectorList, headerList, delimeter=',')
		
		writeOutputMatrix(baseDir+'/'+dir+'/'+dir+'.csv', oMatrix)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
####################################################################################################
# Functions and global variables
####################################################################################################
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
RowNumberToLetterConversionDict = {'1':'A', '2':'B', '3':'C', '4':'D', '5':'E', '6':'F', \
'7':'G', '8':'H'}

RowLetterToNumberConversionDict = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8}

gRows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
gColumns = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

gWellList = \
['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08', 'A09', 'A10', 'A11', 'A12',\
 'B01', 'B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B09', 'B10', 'B11', 'B12', \
 'C01', 'C02', 'C03', 'C04', 'C05', 'C06', 'C07', 'C08', 'C09', 'C10', 'C11', 'C12', \
 'D01', 'D02', 'D03', 'D04', 'D05', 'D06', 'D07', 'D08', 'D09', 'D10', 'D11', 'D12', \
 'E01', 'E02', 'E03', 'E04', 'E05', 'E06', 'E07', 'E08', 'E09', 'E10', 'E11', 'E12', \
 'F01', 'F02', 'F03', 'F04', 'F05', 'F06', 'F07', 'F08', 'F09', 'F10', 'F11', 'F12', \
 'G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10', 'G11', 'G12', \
 'H01', 'H02', 'H03', 'H04', 'H05', 'H06', 'H07', 'H08', 'H09', 'H10', 'H11', 'H12']


queryWellIDRe = re.compile('(\w)(\d*)')

def ConvertRowAndColumnNumberToWellID(queryRowNumber, queryColumnNumber):
	queryRowLetter = RowNumberToLetterConversionDict[str(queryRowNumber)]

	wellID = queryRowLetter + str(queryColumnNumber)
	return wellID

def ConvertWellIDToRowAndColumnNumber(wellID):
	test = queryWellIDRe.search(wellID)
	if test != None:
		queryRow = test.group(1)
		queryColumn = test.group(2)
	
		queryRowNumber = RowLetterToNumberConversionDict[queryRow]
		queryColumnNumber = int(queryColumn)
		
	return [queryRowNumber, queryColumnNumber]
	
def ConvertWellIDToWellNumber(wellID):
	wellNumber = gWellList.index(wellID)
	return wellNumber


def ConvertWellNumberToWellID(wellNumber):
	wellID = gWellList[wellNumber]
	return wellID
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ProcessTimeStamp(fileName):
	import re
	import os
	
	baseName = os.path.basename(fileName)
	
	timeStampRegex = re.compile("(\w*)\-(\d+)")
	timeStampSearch = timeStampRegex.search(baseName)

	if timeStampSearch != None:
		timeStamp = timeStampSearch.group(2)
	else:
		ex = TimeStampFailure('Time Stamp Failure')
		raise ex
		
	return timeStamp
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
####################################################################################################
# Vector processing functions
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def vecLength(x1):
	from numpy import sqrt
	length = sqrt(dotProduct(x1, x1))
	return length
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def dotProduct(vec1, vec2):
	i = 0
	dotProduct = 0.0
	while i < len(vec1):
		dotProduct = dotProduct + vec1[i]*vec2[i]
		i = i + 1
	return dotProduct
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def normalize(x1):
	import numpy
	normalized = numpy.zeros(len(x1))
	i = 0
	mag = numpy.sqrt(dotProduct(x1, x1))
	while i < len(x1):
		normalized[i] = x1[i] / mag
		i = i + 1
	
	return normalized
# ------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------ #

####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def MakeCircularMask(centerCoords, radius):
	from numpy import array, int, sqrt

	pixelCoordsSet = []
	
	# Using computer coordinate system in which origin of x,y system is at top left of window
	upperLeftPixel = array([centerCoords[0]-radius, centerCoords[1]-radius], int)
	upperRightPixel = array([centerCoords[0]+radius, centerCoords[1]-radius], int)
	lowerLeftPixel = array([centerCoords[0]-radius, centerCoords[1]+radius], int)
	lowerRightPixel = array([centerCoords[0]+radius, centerCoords[1]+radius], int)
	
	leftEdge = upperLeftPixel[0]
	rightEdge = upperRightPixel[0]
	topEdge = upperLeftPixel[1]
	bottomEdge = lowerRightPixel[1]
	
	# Scan left to right, top to bottom, picking up pixels
	
	x = leftEdge
	y = topEdge
	
	while x < rightEdge:
		y = topEdge
		while y < bottomEdge:
			radialDistance = vecLength(array([x,y])-(centerCoords))
			if radialDistance <= radius:
				pixelCoordsSet.append([x, y])
			y += 1
		x += 1
	
	return pixelCoordsSet
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def MoveMask(pixelCoordsSet, offsetVector):
	from copy import copy
	from numpy import array, int
	
	offsetPixelCoordsSet = copy(pixelCoordsSet)
	i = 0
	while i < len(offsetPixelCoordsSet):
		offsetPixelCoordsSet[i] = array(offsetPixelCoordsSet[i], int) + array(offsetVector, int)
		i += 1
	
	return offsetPixelCoordsSet
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GetColorValues(pixels, pixelCoordsSet, imageSize):
	pixelValueSet = []
	for pix in pixelCoordsSet:
		if (pix[0] >= 0 and pix[0] < imageSize[0]) and (pix[1] >= 0 and pix[1] < imageSize[1]):
			pixelValueSet.append(pixels[int(pix[0]), int(pix[1])])
	return pixelValueSet
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def GetMeanAndMedianColorValuesInCircle2(pixels, mask, imageSize, coord, radius):
	from numpy import array, mean, median, zeros
	
	pixelCircleSet = mask
	circleColorValueSet = GetColorValues(pixels, pixelCircleSet, imageSize)
	circleColorValueSet = array(circleColorValueSet)
	
	#print circleColorValueSet
	
	if len(circleColorValueSet) == 0:
		print "Error!"
		meanRed = 0
		meanGreen = 0
		meanBlue = 0
	
		medianRed = 0
		medianGreen = 0
		medianBlue = 0
	else:
		redValueSet = circleColorValueSet[:,0]
		greenValueSet = circleColorValueSet[:,1]
		blueValueSet = circleColorValueSet[:,2]

		meanRed = mean(redValueSet)
		meanGreen = mean(greenValueSet)
		meanBlue = mean(blueValueSet)
	
		medianRed = median(redValueSet)
		medianGreen = median(greenValueSet)
		medianBlue = median(blueValueSet)
	
	
	return [[meanRed, meanGreen, meanBlue], [medianRed, medianGreen, medianBlue]]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class Well:
	def __init__(self, plateID, wellID):
		self.colorValue = 0
		self.plateID = plateID
		self.wellID = wellID
		self.hasODData = False
		self.occupied = False
		self.OD = 0
		self.colorValues = []
		self.timeValues = []
		self.meanRed = []
		self.meanGreen = []
		self.meanBlue = []
		self.medianRed = []
		self.medianGreen = []
		self.medianBlue = []
		self.opticalDensities = []
		self.genotype = ''
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
####################################################################################################
# File utilities
####################################################################################################

# ------------------------------------------------------------------------------------------------ #
def ProcessFileNameExtension(fileName):
	import re
	import os
	import pdb
	
	baseName = os.path.basename(fileName)
	[fileName, fileExtension] = os.path.splitext(baseName)
	
# 	pdb.set_trace()
		
	return [fileName, fileExtension]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def NormalizeWellIDFormat(wellID):
	import re
	queryWellIDRe = re.compile('([A-Ha-h]+)(\d+)')
	
	queryResult = queryWellIDRe.search(wellID)
	
	wellIDNormalized = queryResult.group(1) + queryResult.group(2).zfill(2)
	
	return wellIDNormalized
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def NormalizePlateIDFormat(plateID, prefix='SOGRT',zfillNumber=3):
	import re
	import pdb
	
	#pdb.set_trace()
	
	queryPlateIDRe = re.compile('('+prefix+')_*(\d+)')
	
	queryResult = queryPlateIDRe.search(plateID)
	
	plateIDNormalized = queryResult.group(1) + queryResult.group(2).zfill(zfillNumber)
	
	return plateIDNormalized
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def WriteFittedWellPositions(fileName, dirName, FittedWellPositionsCollectionDict, saveDir):
# Write out the fitted well positions
	from numpy import array, float
	import pdb
	
	#pdb.set_trace()
	
	wellIDs = sorted(FittedWellPositionsCollectionDict[dirName][fileName][2].keys())
	[fileNameBase, fileExtension] = ProcessFileNameExtension(fileName)
	wellCoordsFileName = fileNameBase + '_well_coordinates.csv'
	
	fileHandle = open(saveDir + '/' + wellCoordsFileName, 'w')
	
	data = FittedWellPositionsCollectionDict[dirName][fileName]
	
	for well in wellIDs:
		scale = data[0]
		coord = array(data[2][well], float)/scale
		outputString = NormalizeWellIDFormat(well) + ',' + str(coord[0]) + ',' + str(coord[1]) + '\n'
		fileHandle.write(outputString)
	
	fileHandle.close()
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def ExtractPlateIDFromFileName(fileName):
	import re
	fileNameRe = re.compile('(\w+)\-(\d+)\.(\d+)')
	fileNameSearch = fileNameRe.search(fileName)
	plateID = ''
	if fileNameSearch != None:
		plateID = fileNameSearch.group(1)
	return plateID
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ImportPlateODData(plateDataFile, prefix='CTRL|SOGRT'):
	import re
	import pdb
	
	fileHandle = open(plateDataFile, 'r')
	data = fileHandle.readlines()
	fileHandle.close()

	currentPlate = ''
	plateDict = {}

	for line in data:
		print line.strip()
		dataLine = line.split(',')
		
		#pdb.set_trace()
		
		plate = NormalizePlateIDFormat(dataLine[0], prefix=prefix)
		wellID = dataLine[1]
		row = dataLine[2]
		column = dataLine[3]
		opticalDensity = dataLine[4]
	
		if plate != currentPlate:
			currentPlate = NormalizePlateIDFormat(plate, prefix=prefix)
			plateDict[currentPlate] = {}
	
		normalizedWellID = NormalizeWellIDFormat(wellID)
		#print currentPlate + ': ' + normalizedWellID + ': ' + opticalDensity
	
		plateDict[currentPlate][normalizedWellID] = float(opticalDensity)
	
	return plateDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AddODDataToTimeCourses(WellTimeCourseDict, plateODDataDict, plateODDataDictKey, \
prefix='SOGRT|CTRL', occupiedODThreshold=0.2):
	
	from macroscopeUtils4 import gWellList, NormalizePlateIDFormat
	from numpy import ones, float
	
	plateODDataDictKey = NormalizePlateIDFormat(plateODDataDictKey, prefix=prefix)
	
	plateODData = plateODDataDict[plateODDataDictKey]
	
	
	for well in gWellList:
		#opticalDensities = ones(len(WellTimeCourseDict[well].timeValues), float)*plateODData[well]
		#WellTimeCourseDict[well].opticalDensities = opticalDensities
		WellTimeCourseDict[well].OD = float(plateODData[well])
		WellTimeCourseDict[well].hasODData = True
		
		if plateODData[well] > occupiedODThreshold:
			WellTimeCourseDict[well].occupied = True
		else:
			WellTimeCourseDict[well].occupied = False
			
		

	return WellTimeCourseDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def AddODDataToCollectedTimeCourses(CollectedWellTimeCourseDict, plateODDataDict, dirList, \
prefix='SOGRT|CTRL'):
	
	for dir in dirList:
		WellTimeCourseDict = CollectedWellTimeCourseDict[dir]
		
		
		
		WellTimeCourseDict = AddODDataToTimeCourses(WellTimeCourseDict, plateODDataDict, \
		dir, prefix=prefix)

		CollectedWellTimeCourseDict[dir] = WellTimeCourseDict
		
		
	return CollectedWellTimeCourseDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def ImportPlateTimeCourse(fileName, plateID, occupiedODThreshold=0.2):
	
	import pdb
	from numpy import array, float
	from macroscopeUtils4 import Well
	
	dataHandle = open(fileName, 'r')
	data = dataHandle.readlines()
	dataHandle.close()
	
	
	WellTimeCourseDict = {}
	
	lines = []
	
	for line in data:
		line = line.strip()
		line = line.split(',')
		lines.append(line)
	
	keys = lines[0]
	lines.pop(0)
	
	data = array(lines, float)
	
	i = 0
	dataDict = {}
	while i < data.shape[1]:
		dataDict[keys[i]] = data[:,i]
		i += 1
	
	
	wells = []
	for key in keys:
		well = key.split('_')[0]
		if well not in wells:
			if well != '':
				wells.append(well)
	
	
	for well in wells:
		WellTimeCourseDict[well] = Well(plateID, well)
		WellTimeCourseDict[well].timeValues = dataDict[well+'_time']
		WellTimeCourseDict[well].meanRed = dataDict[well+'_r']
		WellTimeCourseDict[well].meanGreen = dataDict[well+'_g']
		WellTimeCourseDict[well].meanBlue = dataDict[well+'_b']
		
		if dataDict.has_key(well+'_od'):
			WellTimeCourseDict[well].OD = dataDict[well+'_od'][0]
			WellTimeCourseDict[well].hasODData = True
			if WellTimeCourseDict[well].OD > occupiedODThreshold:
				WellTimeCourseDict[well].occupied = True
			else:
				WellTimeCourseDict[well].occupied = False
		else:
			WellTimeCourseDict[well].hasODData = False
		
		#pdb.set_trace()
	
	return WellTimeCourseDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def FindMaxTimeValueInWellTimeCourseDict(WellTimeCourseDict):
	from macroscopeUtils4 import gWellList
	
	import pdb
	
	maxTimeValue = 0
	
	for well in gWellList:
		maxWellTime = max(WellTimeCourseDict[well].timeValues)
		if maxWellTime > maxTimeValue:
			maxTimeValue = maxWellTime
		
	return maxTimeValue
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def FindBlobSelfVectors(fblobs):
	from SimpleCV.Features.Features import FeatureSet
	
	
	selfVectors = []
	origins = []
	
	i = 0
	while i < len(fblobs):
		currentBlob = fblobs[i]
		currentOrigin = fblobs[i].coordinates()
		origins.append(currentBlob.coordinates())
		selfVectors.append([])
		
		j = 0
		while j < len(fblobs):
			if i != j:
				selfVectors[i].append(fblobs[j].coordinates()-fblobs[i].coordinates())
			j += 1
		i += 1
		
	return [origins, selfVectors]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def DrawWellPositions(WellPositionsDict, img, color, drawPoints=False):
	for key in WellPositionsDict.keys():
		coord = WellPositionsDict[key]
		img.drawText(key, x=coord[0], y=coord[1], color=color)
		if drawPoints:
			img.drawPoints([coord], color=color)

	img.show()
# ------------------------------------------------------------------------------------------------ #



