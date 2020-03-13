import numpy as np

# ------------------------------------------------------------------------------------------------ #
class PlatePosition:
	def __init__(self, number, rows, columns):
		self.number = number
		self.plate = None
		[self.row, self.column] = convertPositionNumberToRowAndColumn(number, rows, columns)
		
class Plate:
	def __init__(self, number):
		self.number = number


def shuffleArray(array):
	from random import shuffle
	shuffle(array)
	return array
	#TODO remove this function - unnecessary

def convertPositionNumberToRowAndColumn(positionNumber, rows, columns):
	from math import ceil
	
	row = ceil(positionNumber/columns)
	col = np.mod(positionNumber, columns)
	
	if col == 0:
		col = columns
	
	
	return [row, col]

# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def GenerateCheckSum(text, maxLength=8):	
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
	subtractionArrayCopy = np.zeros(maxLength, np.int)

	while j < min([len(subtractionArray), maxLength]):
		subtractionArrayCopy[j] = subtractionArray[j]
		j += 1
	
	j = 0
	while j < min([len(subtractionArray), maxLength]):
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

import numpy as np
import random


# NOTE What determines these numbers?
maxRows = 20
maxColumns = 21
maxPlates = 417

# i = 1 #TODO remove
platePositions = []
plates = []

for i in range(1, maxPlates+1):
# while i <= maxPlates: `	#TODO  remove
	platePositions.append(PlatePosition(i, maxRows, maxColumns))
	plates.append(Plate(i))	
	#i += 1 #TODO remove

# platesShuffled = shuffleArray(plates)
	#TODO why is this a thing
platesShuffled = plates

# i = 0 #TODO remove

for i in range(len(platesShuffled)):
# while i < len(platesShuffled): #TODO remove
	platePositions[i].plate = platesShuffled[i]
	# i += 1 #TODO remove


for platePosition in platePositions:
	outputString = str(platePosition.number) + '\t' + str(platePosition.plate.number) \
	+ '\t' + str(platePosition.row) + '\t' + str(platePosition.column)
	print(outputString)
	
outputString = ''
delimeter = '\t'

for platePosition in platePositions:
	outputString += str(platePosition.plate.number) + delimeter
	if (platePosition.number)%maxColumns == 0:
		outputString += '\n'
	
print(outputString)


sortedPlatePositions = sorted(platePositions, key=lambda platePosition: platePosition.plate.number)

for platePosition in sortedPlatePositions:
	outputString = 'P' + str(platePosition.plate.number).zfill(3) \
	+ 'R' + str(platePosition.row).zfill(2) \
	+ 'C' + str(platePosition.column).zfill(2) + 'T'
	barcode = outputString + GenerateCheckSum(outputString, maxLength=11)
	print("*" + barcode + "*")
	






