# ------------------------------------------------------------------------------------------------ #
# Interpolate the well color data and write out for use in plotting

from macroscopeUtils7 import ExperimentWell, ImportCondensedColorDataFromXML, \
CalculateRelativeDataDictArrayForExperimentalWell, ScaleDataDictArrayForExperimentalWell, \
GenerateWellDisplayCoordinatesForOriginalPlateLayout, DisplayWellColorsOnGrid, \
WriteCondensedColorDataToXML

from scipy import unique
from scipy.interpolate import interp1d

from numpy import argsort

import pdb

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the collected well color data from XML.
wellColorsXMLFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Condensed Data/' \
+ 'wellColors.xml'

collectedWellColorsDict = ImportCondensedColorDataFromXML(wellColorsXMLFileName)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Compute relative time arrays (time from earliest time point) for each experimental well and 
# convert
for plate in collectedWellColorsDict.keys():
	for wellKey in collectedWellColorsDict[plate].keys():
		
		well = collectedWellColorsDict[plate][wellKey]
		
		CalculateRelativeDataDictArrayForExperimentalWell(well, originalDataDictArray='time', \
		relativeDataDictArray='relativeTime')
	
		ScaleDataDictArrayForExperimentalWell(well, originalDataDictArray='relativeTime', \
		conversionFactor=1/3600, convertedDataDictArray='relativeTimeHours')
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# Find the maximum and minimum relative times, mean reds, mean greens and mean blues in the 
# experiment
minRelativeTime = 0
maxRelativeTime = 40
minColor = 150
maxColor = 151
lastTimesArray = []

for plate in collectedWellColorsDict.keys():
	for wellKey in collectedWellColorsDict[plate].keys():
		
		well = collectedWellColorsDict[plate][wellKey]
		
		relativeTime = well.dataDict['relativeTimeHours']
		maxRelativeTime = max(well.dataDict['relativeTimeHours'])
		lastTimesArray.append(maxRelativeTime)
		
		meanRed = array(well.dataDict['meanRed'], float)
		meanGreen = array(well.dataDict['meanGreen'], float)
		meanBlue = array(well.dataDict['meanBlue'], float)

		if max(relativeTime) > maxRelativeTime:
			maxRelativeTime = max(relativeTime)
		if min(relativeTime) < minRelativeTime:
			minRelativeTime = min(relativeTime)
		if max(meanRed) > maxColor:
			maxColor = max(meanRed)
		if min(meanRed) < minColor:
			minColor = min(meanRed)
		if max(meanGreen) > maxColor:
			maxColor = max(meanGreen)
		if min(meanGreen) < minColor:
			minColor = min(meanGreen)
		if max(meanBlue) > maxColor:
			maxColor = max(meanBlue)
		if min(meanBlue) < minColor:
			minColor = min(meanBlue)
	
xLim = [minRelativeTime, maxRelativeTime]
yLim = [minColor, maxColor]
lastTimes = unique(lastTimesArray)
earliestLastTime = min(lastTimes)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Set up an interpolation that figures out the colors in the wells at extremely fine time points
# between time 0 and the earliest last time.

relativeTimeHoursInterp = arange(0, earliestLastTime, 0.1)
for plate in collectedWellColorsDict.keys():
	for wellKey in collectedWellColorsDict[plate].keys():
		
		well = collectedWellColorsDict[plate][wellKey]
		
		relativeTime = well.dataDict['relativeTimeHours']
		meanRed = array(well.dataDict['meanRed'], float)
		meanGreen = array(well.dataDict['meanGreen'], float)
		meanBlue = array(well.dataDict['meanBlue'], float)
		
		interpFunctionMeanRed = interp1d(relativeTime, meanRed)
		interpFunctionMeanGreen = interp1d(relativeTime, meanGreen)
		interpFunctionMeanBlue = interp1d(relativeTime, meanBlue)
		
		interpMeanRed = interpFunctionMeanRed(relativeTimeHoursInterp)
		interpMeanGreen = interpFunctionMeanGreen(relativeTimeHoursInterp)
		interpMeanBlue = interpFunctionMeanBlue(relativeTimeHoursInterp)
		
		well.dataDict['relativeTimeHoursInterp'] = relativeTimeHoursInterp
		well.dataDict['interpMeanRed'] = interpMeanRed
		well.dataDict['interpMeanGreen'] = interpMeanGreen
		well.dataDict['interpMeanBlue'] = interpMeanBlue
	
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Write out the updated data

updatedWellColorsXMLFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Condensed Data/' \
+ 'wellColors_interpolated.xml'

experimentalWellsLinearArray = []
for plate in collectedWellColorsDict.keys():
	for wellKey in collectedWellColorsDict[plate].keys():
		experimentalWellsLinearArray.append(collectedWellColorsDict[plate][wellKey])

WriteCondensedColorDataToXML(updatedWellColorsXMLFileName, experimentalWellsLinearArray)
# ------------------------------------------------------------------------------------------------ #
