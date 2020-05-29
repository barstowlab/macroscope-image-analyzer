# Script to grab data from condensed well array
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Circle

from macroscopeUtils7 import ImportCondensedColorDataFromXML, gWellList

from macroscopeUtils7 import ExperimentWell, ImportCondensedColorDataFromXML, \
CalculateRelativeDataDictArrayForExperimentalWell, ScaleDataDictArrayForExperimentalWell, \
GenerateWellDisplayCoordinatesForOriginalPlateLayout, DisplayWellColorsOnGrid, \
PlotExperimentWellTimeCourseGraphsOnGrid



fileName = '/Users/Don/Dropbox/group eight/2015-11-14-19 - AQDS Reduction with Shewanella\
 Sudoku Library Parts 1 and 2/Condensed Data/wellColors_interpolated.xml'


condensedData = ImportCondensedColorDataFromXML(fileName)

plateKeys = sorted(condensedData.keys())

transposonCoordData = []
transposonCoordArray = []

# ------------------------------------------------------------------------------------------------ #
# Import necessary data for 5 interpolated time points


timesToGet = ['0.0', '10.0', '20.0', '30.0', '39.5']


#n=0
for plateKey in plateKeys:
	
	plateData = condensedData[plateKey]
	
	for wellKey in gWellList:
		
		wellData = plateData[wellKey]
		wellCoord = wellData.transposonCoord
		
		if wellCoord != '-1':
		
			selectedMeanRedPoints = []
			selectedMeanGreenPoints = []
			selectedMeanBluePoints = []
			timeArray = wellData.dataDict['relativeTimeHoursInterp']
			
			for timeToGet in timesToGet:
				index = timeArray.index(timeToGet)
				selectedMeanRedPoints.append(wellData.dataDict['interpMeanRed'][index])
				selectedMeanGreenPoints.append(wellData.dataDict['interpMeanGreen'][index])
				selectedMeanBluePoints.append(wellData.dataDict['interpMeanBlue'][index])
				
		
		
			transposonCoordData.append([plateKey, wellKey, float(wellData.transposonCoord), \
			wellData.featureName, wellData.dataDict, \
			wellData.dataDict['relativeTimeHoursInterp'], \
			wellData.dataDict['interpMeanRed'], \
			wellData.dataDict['interpMeanGreen'], \
			wellData.dataDict['interpMeanBlue'], \
			selectedMeanRedPoints, \
			selectedMeanGreenPoints, \
			selectedMeanBluePoints, \
			timesToGet, \
			wellData])
			
			
dataLength = len(transposonCoordData)	

n = 0	
while n < dataLength:
	
	transposonCoordArray.append(float(transposonCoordData[n][2]))
	n += 1
	
	
# ------------------------------------------------------------------------------------------------ #
# Grab all relevant data for desired transposon range

rangeLower = 262500
rangeUpper = 287500
selectedTransposonCoordData = []
selectedTransposonArray = []

i = 0
while i < len(transposonCoordData):
	transposonCoord = transposonCoordData[i][2]
	
	if rangeLower <= transposonCoord <= rangeUpper and transposonCoord not in selectedTransposonArray:
		selectedTransposonCoordData.append(transposonCoordData[i])
		selectedTransposonArray.append(transposonCoord)
	i += 1

sortedSelectedTransposonArray = sorted(selectedTransposonArray)
sortedSelectedTransposonCoordData = sorted(selectedTransposonCoordData, key=lambda x:x[2])

		
# ------------------------------------------------------------------------------------------------ #
# Generate coordinates for data plot

# Find smallest gab between transposons

i=0
differences= []
while i < len(selectedTransposonCoordData)-1:
	
	differences.append((float(sortedSelectedTransposonCoordData[i+1][2]))\
	 -(float(sortedSelectedTransposonCoordData[i][2])))
	i+=1
gap = min(differences)


plotGridX = rangeUpper-rangeLower 
plotOriginX = rangeLower
plotOriginY = 0

# Concatenate array of relative X locations

colSepX = []
i = 1

while i < len(selectedTransposonCoordData):
	Xright = (float(sortedSelectedTransposonCoordData[i][2])-rangeLower)
	Xleft = (float(sortedSelectedTransposonCoordData[i-1][2])-rangeLower)
	colSepX.append(Xright-Xleft)
	i+=1

# Concatenate array of circle X locations (left to right)

circleXArray = [0]	
i=1

while i < len(selectedTransposonCoordData):
	circleXArray.append(circleXArray[i-1] + colSepX[i-1])
	i += 1
	
i=0

while i < len(selectedTransposonCoordData)-1:
	circleXArray[i]= circleXArray[i] + plotOriginX
	i += 1

# Concatenate array of circle Y locations (top to bottom)

colSepY = gap
circleYArray = []
i = 0

while i < len(timesToGet):
	circleYArray.append(plotOriginY + .2*colSepY + i*.5*colSepY)
	i += 1
	
circleYArray.reverse()


# ------------------------------------------------------------------------------------------------ #
# Plot out results

summaryFig = pyplot.figure(figsize=(11,2.5))
circleRadius = 0.4*gap
reds=[]
i = 0
while i < len(circleXArray):
	j = 0
	
	summaryFig.gca().text(circleXArray[i], circleYArray[0] + colSepY/4, \
	 sortedSelectedTransposonCoordData[i][3], fontsize=7, rotation=270, ha='right', va='bottom')
	
	# Draw wide ellipses to counteract stretch (circles calibrated to PlotGridX = 25,000)
	while j < len(circleYArray):

		currentRed = float(sortedSelectedTransposonCoordData[i][-5][j])
		currentGreen = float(sortedSelectedTransposonCoordData[i][-4][j])
		currentBlue = float(sortedSelectedTransposonCoordData[i][-3][j])
		circleColor = (currentRed/255., currentGreen/255., currentBlue/255.)
		wellCircleCenter = (circleXArray[i], circleYArray[j])
		reds[j]=reds.append(currentRed)
		tempCircle = Ellipse(wellCircleCenter, width=2*circleRadius, \
		 height=circleRadius/2, color=circleColor)
	
		summaryFig.gca().add_artist(tempCircle)
		
		j += 1
	i += 1

xlim(plotOriginX - gap, plotOriginX + plotGridX + gap)
ylim(plotOriginY, circleYArray[0] + colSepY)
show()

# Thank you for a great semester, Buz! 