# ------------------------------------------------------------------------------------------------ #
# Import and plot compiled color data

from macroscopeUtils7 import ImportCondensedColorDataFromXML, gWellList
from copy import deepcopy
import matplotlib.pyplot as pyplot
from matplotlib.pyplot import Circle
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Import the collected well color data from XML.
wellColorsXMLFileName = '../Condensed Data/' \
+ 'wellColors_interpolated.xml'

timeCourseFigOutputDir = '../Time Courses/'

collectedWellColorsDict = ImportCondensedColorDataFromXML(wellColorsXMLFileName)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Pick out the 
selectedGeneNames = [\
'SO_2083', 'lpxM', 'hypA', 'hypE', 'hypD', 'hypC', 'hypB', 'hypF', 'hyaE', 'hyaD', 'hyaC', 'hyaB', \
'hyaA', 'SO_2100', 'SO_2102', 'tufA', 'ccmE ', 'ccmD', 'ccmC', 'ccmB', 'ccmA ', 'ccmI ', 'ccmF ', \
'ccmG', 'ccmH', 'SO_0269', 'pncC', 'SO_1907', 'SO_1909', 'menA', 'SO_1911', 'tesB', 'SO_4571', \
'SO_4572', 'menD', 'menH', 'menC', 'menE', 'SO_4577', 'fimA', 'SO_4711', 'SO_4712', 'menF', \
'SO_4714', 'SO_4715', 'ccoG', 'SO_4738', 'menB', 'SO_4740', 'glmR', 'ackA', 'pta', 'SO_2917', \
'SO_2919', 'sdaR ', 'SO_1775', 'mtrB', 'mtrA', 'mtrC', 'omcA', 'mtrF', 'mtrE', 'mtrD', 'SO_1787', \
'folD', 'SO_4589', 'SO_4590', 'cymA', 'SO_4592', 'SO_4593']

plateKeys = collectedWellColorsDict.keys()

selectedGenesArray = []

for plateKey in plateKeys:
	for well in gWellList:
		currentWell = collectedWellColorsDict[plateKey][well]
		currentFeatureName = collectedWellColorsDict[plateKey][well].featureName
		if currentFeatureName in selectedGeneNames:
			selectedGenesArray.append(currentWell)

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
timesToGet = ['0.0', '10.0', '20.0', '30.0', '39.5']


# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Find the earliest last time
summarizedDataArray = []
for well in selectedGenesArray:
	
	tempDataDict = {}
	
	timeArray = well.dataDict['relativeTimeHoursInterp']
	selectedMeanRedPoints = []
	selectedMeanGreenPoints = []
	selectedMeanBluePoints = []
	
	for timeToGet in timesToGet:
		index = timeArray.index(timeToGet)
		selectedMeanRedPoints.append(well.dataDict['interpMeanRed'][index])
		selectedMeanGreenPoints.append(well.dataDict['interpMeanGreen'][index])
		selectedMeanBluePoints.append(well.dataDict['interpMeanBlue'][index])
	
	tempDataDict['interpMeanRed'] = selectedMeanRedPoints
	tempDataDict['interpMeanGreen'] = selectedMeanGreenPoints
	tempDataDict['interpMeanBlue'] = selectedMeanBluePoints
	tempDataDict['relativeTimeHoursInterp'] = timesToGet
	tempDataDict['featureName'] = well.featureName
	tempDataDict['plateID'] = well.plateID
	tempDataDict['wellID'] = well.wellID
	tempDataDict['transposonCoord'] = int(well.transposonCoord)
	
	summarizedDataArray.append(tempDataDict)

summarizedDataArray = sorted(summarizedDataArray, key=lambda k: k['transposonCoord']) 
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Generate coordinates for data plot
plotGridX=8500
circleXArray = []
i = 0
colSepX = plotGridX/(len(summarizedDataArray)-1)
plotOriginX = 0
plotOriginY = 0
colSepY = colSepX

while i < len(summarizedDataArray):
	circleXArray.append(plotOriginX + i*colSepX)
	i += 1

i = 0
circleYArray = []
while i < len(timesToGet):
	circleYArray.append(plotOriginY + i*colSepY)
	i += 1
	
circleYArray.reverse()
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Plot out results
summaryFig = pyplot.figure(figsize=(11,8.5))
circleRadius = colSepY*0.4

i = 0
while i < len(circleXArray):
	j = 0
	currentGene = summarizedDataArray[i]
	
	summaryFig.gca().text(circleXArray[i], circleYArray[0] + colSepY/2, currentGene['featureName'], \
	fontsize=5, rotation=270, ha='right', va='bottom')
	
	while j < len(circleYArray):
		currentRed = float(currentGene['interpMeanRed'][j])
		currentGreen = float(currentGene['interpMeanGreen'][j])
		currentBlue = float(currentGene['interpMeanBlue'][j])
		circleColor = (currentRed/255., currentGreen/255., currentBlue/255.)
		wellCircleCenter = (circleXArray[i], circleYArray[j])
		
		tempCircle = Circle(wellCircleCenter, radius=circleRadius, color=circleColor)
	
		summaryFig.gca().add_artist(tempCircle)
		
		j += 1
	i += 1

xlim(-colSepX, plotGridX + colSepY)
ylim(-colSepY, circleYArray[0] + colSepY)
show()
summaryFig.savefig('AQDSSummaryMap.eps')
# ------------------------------------------------------------------------------------------------ #
