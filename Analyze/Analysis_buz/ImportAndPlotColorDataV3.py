# ------------------------------------------------------------------------------------------------ #
def FindMinAndMaxTimesAndColorValuesInExperimentalWellLinearArray(experimentWellsLinearArray):
# Find the maximum and minimum relative times, mean reds, mean greens and mean blues in the 
# experiment
	minRelativeTime = 0
	maxRelativeTime = 40
	minColor = 150
	maxColor = 151

	for well in experimentWellsLinearArray:
		relativeTime = array(well.dataDict['relativeTimeHours'], float)
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
	
	return [xLim,  yLim]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def OnPickExperimentalWell(event):
	
	from macroscopeUtils7 import ExperimentWellCircle
	
	import pdb
	
	artist = event.artist
	
	if isinstance(artist, ExperimentWellCircle):
		featureName = artist.experimentalWell.featureName
		plateID = artist.experimentalWell.plateID
		wellID = artist.experimentalWell.wellID

	print('Plate ' + str(plateID) + ', ' + str(wellID) + ': ' + str(featureName) )
# ------------------------------------------------------------------------------------------------ #








# ------------------------------------------------------------------------------------------------ #
# Import and plot compiled color data

from macroscopeUtils7 import ExperimentWell, ImportCondensedColorDataFromXML, \
CalculateRelativeDataDictArrayForExperimentalWell, ScaleDataDictArrayForExperimentalWell, \
GenerateWellDisplayCoordinatesForOriginalPlateLayout, DisplayWellColorsOnGrid, \
PlotExperimentWellTimeCourseGraphsOnGrid

from numpy import argsort

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the collected well color data from XML.
wellColorsXMLFileName = '../Condensed Data/wellColors_interpolated.xml'

timeCourseFigOutputDir = '../Time Courses/'

collectedWellColorsDict = ImportCondensedColorDataFromXML(wellColorsXMLFileName)
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# Generate a dictionary of well coords for the collectedWellColorsDict
[plateBoundariesX, plateBoundariesY, wellDisplayCoordsDict] = \
GenerateWellDisplayCoordinatesForOriginalPlateLayout(collectedWellColorsDict, platesPerRow=12, \
plotGridX=1440, origin=0, coordsDataDictName='wellDisplayCoordsArrangedByOriginalPlate')
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Arrange the experiment well objects into a single list ordered first by plate number and then 
# by well id.
experimentWellsLinearArray = []
plateKeys = sorted(collectedWellColorsDict.keys())
for plateKey in plateKeys:
	wellKeys = sorted(collectedWellColorsDict[plateKey].keys())
	for wellKey in wellKeys:
		experimentWellsLinearArray.append(collectedWellColorsDict[plateKey][wellKey])
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Go through the experimental wells and flag potential hits
hitsArray = []
for well in experimentWellsLinearArray:
	firstBlue = float(well.dataDict['meanBlue'][0])
	lastBlue = float(well.dataDict['meanBlue'][-1])
	
	status = well.occupancyTest
	occupancy = well.occupied
	
	if lastBlue > 0.9*firstBlue and status == 'OK' and occupancy == 'True':
		hitsArray.append(well)
		well.potentialHit = True

hitsFileHandle = open('../Condensed Data/hits.csv', 'w')
for hit in hitsArray:
	outputStr = hit.featureName + ',' + str(hit.transposonCoord) + ',' + hit.plateID + ',' \
	+ hit.wellID + '\n'
	hitsFileHandle.write(outputStr)
hitsFileHandle.close()
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Find the maximum and minimum relative times, mean reds, mean greens and mean blues in the 
# experiment
[xLim, yLim] = \
FindMinAndMaxTimesAndColorValuesInExperimentalWellLinearArray(experimentWellsLinearArray)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Display the well grids
wellColorGridFig = DisplayWellColorsOnGrid(experimentWellsLinearArray, \
coordsDataDictName='wellDisplayCoordsArrangedByOriginalPlate', diagnosticMarkUp=True, \
plateBoundariesX=plateBoundariesX, plateBoundariesY=plateBoundariesY, showGrid=True)

wellColorGridFig.canvas.mpl_connect('pick_event', OnPickExperimentalWell)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Plot out the time course graphs for each well
PlotExperimentWellTimeCourseGraphsOnGrid(experimentWellsLinearArray, 8, 12, xLim=xLim, \
yLim=yLim, saveDir=timeCourseFigOutputDir)

# ------------------------------------------------------------------------------------------------ #
