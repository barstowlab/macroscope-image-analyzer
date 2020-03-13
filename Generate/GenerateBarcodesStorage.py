import numpy as np


# ------------------------------------------------------------------------------------------------ #
def GenerateCheckSum(text):
	import numpy as np
	textArray = []
	ordArray = []
	subtractionArray = []
	
	for letter in text:
		textArray.append(letter)
		ordArray.append(ord(letter))
		subtractionArray.append(ord(letter) - 55)
		
	j = 0
	checksum = 0
	subtractionArrayCopy = np.zeros(8, np.int)

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






# outputFile = "REBarcodes4SIDE1000.csv"
outputFile = "REBarcodes4p2SIDE1000.csv"


delimeter = ","
lineBreak = "\n"

prefix = "RE-"
postfix = ""

# For the Diversified Biotech SIDE1000 labels, you have a grid of 4x39 labels.

MaxRows = 39
MaxColumns = 4
TotalLabels = 4*39

PaddingColumns = True
PaddingRows = False


# This is where I add in special barcodes
# sampleNumbers = arange(201, 357,1)
sampleNumbers = np.arange(757,804,1) 
	#TODO - Change this to accept input from user, and then do for number in range(...) in line 82, currently not inclusive



# ------------------------------------------------------------------------------------------------ #
# Generate the sample labels
sampleIDs = []
for number in sampleNumbers:
	sampleID = prefix + str(number).zfill(3) + 'T'
	sampleIDs.append(sampleID)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
labels = []
for sampleID in sampleIDs:
	barcode = sampleID + GenerateCheckSum(sampleID)
	labels.append("*" + barcode + "*")
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
row = 1
column = 1
labelCount = 1
TotalLabels = min([len(labels), TotalLabels])


rowData = []

while row <= MaxRows:
	rowText = ""
	
	column = 1
	while column <= MaxColumns and labelCount <= TotalLabels:
		label = labels[labelCount-1]
		labelCount += 1
		
		if column < MaxColumns:
			rowText += label + delimeter
			if PaddingColumns == True:
				rowText += delimeter
		else:
			rowText += label + lineBreak
			if PaddingRows == True:
				rowText += lineBreak
				
		column += 1
		
	rowData.append(rowText)
	row += 1
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Output the row data

fileHandle = open(outputFile, 'w')
for row in rowData:
	fileHandle.write(row)
	
fileHandle.close()
# ------------------------------------------------------------------------------------------------ #

