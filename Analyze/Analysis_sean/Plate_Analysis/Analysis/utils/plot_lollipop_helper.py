# helper function for PlotLollipopGraph.py

# Functions for making lollipop graphs

# ------------------------------------------------------------------------------------------------ #
def GenerateCoordinatesForLollipopGraph(summarizedDataArray, timesToGet, plotGridX=8500):
    # Generate coordinates for data plot

    circleXArray = []
    i = 0
    colSepX = plotGridX / (len(summarizedDataArray) - 1)
    plotOriginX = 0
    plotOriginY = 0
    colSepY = colSepX

    while i < len(summarizedDataArray):
        circleXArray.append(plotOriginX + i * colSepX)
        i += 1

    i = 0
    circleYArray = []
    while i < len(timesToGet):
        circleYArray.append(plotOriginY + i * colSepY)
        i += 1

    circleYArray.reverse()

    return [circleXArray, circleYArray]


# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotLollipopGraph(summarizedDataArray, timesToGet, lollipopOutputDir, plotGridX=8500):
    import matplotlib.pyplot as pyplot
    from matplotlib.patches import Circle

    [circleXArray, circleYArray] = GenerateCoordinatesForLollipopGraph(summarizedDataArray, timesToGet, plotGridX=plotGridX)

    colSepX = plotGridX / (len(summarizedDataArray) - 1)
    summaryFig = pyplot.figure(figsize=(11, 8.5))
    plotOriginX = 0
    plotOriginY = 0
    colSepY = colSepX
    circleRadius = colSepY * 0.4

    i = 0
    while i < len(circleXArray):
        j = 0
        currentGene = summarizedDataArray[i]

        summaryFig.gca().text(circleXArray[i], circleYArray[0] + colSepY / 2,
                              currentGene['GeneName'],
                              fontsize=5, rotation=270, ha='right', va='bottom')

        while j < len(circleYArray):
            currentRed = float(currentGene['MeanRed'][j])
            currentGreen = float(currentGene['MeanGreen'][j])
            currentBlue = float(currentGene['MeanBlue'][j])
            circleColor = (currentRed / 255., currentGreen / 255., currentBlue / 255.)
            wellCircleCenter = (circleXArray[i], circleYArray[j])

            tempCircle = Circle(wellCircleCenter, radius=circleRadius, color=circleColor)

            summaryFig.gca().add_artist(tempCircle)

            j += 1
        i += 1

    pyplot.xlim(-colSepX, plotGridX + colSepY)
    pyplot.ylim(-colSepY, circleYArray[0] + colSepY)
    pyplot.show()
    summaryFig.savefig(lollipopOutputDir + '/' + 'AQDSSummaryMap.eps')

    return
# ------------------------------------------------------------------------------------------------ #