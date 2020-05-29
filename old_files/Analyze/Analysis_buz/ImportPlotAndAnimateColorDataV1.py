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
GenerateWellDisplayCoordinatesForOriginalPlateLayout, DisplayWellColorsOnGrid, ExperimentWellCircle

from matplotlib import animation
from matplotlib import pyplot as pyplot
from  matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

from numpy import argsort


# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the collected well color data from XML.
wellColorsXMLFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Condensed Data/' \
+ 'wellColors_interpolated.xml'

timeCourseFigOutputDir = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Time Courses/'

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
# Convert the interpolated colors to float arrays
i = 0
while i < len(experimentWellsLinearArray):
	workingWell = experimentWellsLinearArray[i]
	workingWell.dataDict['interpMeanRed'] = array(workingWell.dataDict['interpMeanRed'], float)
	workingWell.dataDict['interpMeanGreen'] = array(workingWell.dataDict['interpMeanGreen'], float)
	workingWell.dataDict['interpMeanBlue'] = array(workingWell.dataDict['interpMeanBlue'], float)
	
	workingWell.dataDict['interpMeanRedScaled'] = workingWell.dataDict['interpMeanRed']/255
	workingWell.dataDict['interpMeanGreenScaled'] = workingWell.dataDict['interpMeanGreen']/255
	workingWell.dataDict['interpMeanBlueScaled'] = workingWell.dataDict['interpMeanBlue']/255

	
	i += 1
	


# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Animate the color dot plot

animationFig = pyplot.figure()
coordsDataDictName = 'wellDisplayCoordsArrangedByOriginalPlate'
circlesArray = []
showGrid = True
maxFrames = len(experimentWellsLinearArray[0].dataDict['relativeTimeHoursInterp'])
timeStamps = experimentWellsLinearArray[0].dataDict['relativeTimeHoursInterp']


def init():
	wellCentersX = []
	wellCentersY = []
	i = 0
	while i < len(experimentWellsLinearArray):
	
		workingWell = experimentWellsLinearArray[i]
		firstRed = workingWell.dataDict['interpMeanRedScaled'][0]
		firstGreen = workingWell.dataDict['interpMeanGreenScaled'][0]
		firstBlue =workingWell.dataDict['interpMeanBlueScaled'][0]
	
		circleColor = (firstRed, firstGreen, firstBlue)
	
		wellCircleCenterX = workingWell.dataDict[coordsDataDictName].x
		wellCircleCenterY = workingWell.dataDict[coordsDataDictName].y
	
		circleRadius = workingWell.dataDict[coordsDataDictName].radius
		wellCircleCenter = (wellCircleCenterX, wellCircleCenterY)
	
		tempCircle = ExperimentWellCircle(wellCircleCenter, circleRadius, workingWell, \
		color=circleColor,\
		picker=1)
	
		animationFig.gca().add_artist(tempCircle)
	
		circlesArray.append(tempCircle)
		
		wellCentersX.append(wellCircleCenterX)
		wellCentersY.append(wellCircleCenterY)
		
		i += 1
	
	uniqueWellCentersX = sorted(unique(wellCentersX))
	wellCentersY = sorted(unique(wellCentersY))
	wellSepX = wellCentersX[1] -  wellCentersX[0]
	wellSepY = wellCentersY[1] -  wellCentersY[0]
	
	xlim = [min(wellCentersX)-wellSepX, max(wellCentersX)+wellSepX]
	ylim = [min(wellCentersY)-wellSepY, max(wellCentersY)+wellSepY]
	
	animationFig.gca().set_xlim(xlim)
	animationFig.gca().set_ylim(ylim)
	animationFig.gca().invert_yaxis()
	
	animationFig.gca().set_title('Frame Index: ' + str(0).zfill(3) + '. Time: ' \
	+ str(timeStamps[0]).zfill(1) + ' Hours.')
	
	if plateBoundariesX != None:
		animationFig.gca().set_xticks(plateBoundariesX)
	if plateBoundariesY != None: 
		animationFig.gca().set_yticks(plateBoundariesY)
	if showGrid == True:
		animationFig.gca().grid(True)	
	
	return


def updateExperimentalWellColorGridPlot(frameNumber):
	
	i = 0
	while i < len(circlesArray):
		newRed = circlesArray[i].experimentalWell.dataDict['interpMeanRedScaled'][frameNumber]
		newGreen = circlesArray[i].experimentalWell.dataDict['interpMeanGreenScaled'][frameNumber]
		newBlue = circlesArray[i].experimentalWell.dataDict['interpMeanBlueScaled'][frameNumber]
		newColor = (newRed, newGreen, newBlue)
		
		tempCircle = circlesArray[i]
		
		tempCircle.set_color(newColor)

		i += 1
		
	animationFig.gca().set_title('Frame Index: ' + str(frameNumber).zfill(3) + '. Time: ' \
	+ str(timeStamps[frameNumber]).zfill(1) + ' Hours.')



wellColorGridAnimation = \
FuncAnimation(animationFig, updateExperimentalWellColorGridPlot, init_func=init,\
interval=1, frames=maxFrames)
		

pyplot.show()

# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Buz Barstow'))

animationFileName = '/Users/buz/Dropbox (BarstowLab)/BarstowLab Shared Folder/Analysis/' \
+ '2015-11-14-19 - AQDS Reduction with Shewanella Sudoku Library Parts 1 and 2/Animation/' \
+ 'wellColorGrid.mp4'

wellColorGridAnimation.save(animationFileName)

# ------------------------------------------------------------------------------------------------ #

