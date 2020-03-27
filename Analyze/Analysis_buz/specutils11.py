from scipy import *
from numpy import *
from pylab import *

# ---------------------------------------------------------------------------- #
def scaleres(p,y,x):
# A residuals function to calculate best scaling factor between two intensity 
# matrices
# Example of use:
# scalingCoeff = scipy.optimize.leastsq(scaleres, 1, args=(vector1, vector2))
# scalingCoeff will return a multiplication coefficient that will allow vector2
# to be multiplied to most closely match vector1
	err = y - x*p
	return err
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Moving average function
def MovingAverage(inputMatrix, width):
# Calculate the moving average of an input vector, in a window centered around 
# each element in the inputMatrix. The width of the window must be an odd number.
# For elements closer than width/2 to the edges of the input matrix, the moving
# average function takes a window of the same width, but shifts its center so as 
# not to step outside the array. For the 0th element, the average is calculated 
# from the first width elements. For the 1st element, the average is calculated 
# from the zeroth to the width element

	# Check that the moving average window size is an odd number
	if width%2 != 1:
		print("Moving average window must be an odd number.")
		return
	
	movingAverageMatrix = zeros(len(inputMatrix), float)
	i = int(width/2)
	
	# Calculate the moving average for elements in the center of the inputMatrix
	lowerLimit = -int(width/2)
	if width == 1:
		return inputMatrix
	else:
		upperLimit = width/2
	
	while i < (len(inputMatrix) - width/2):
		j = lowerLimit 
		while j <= upperLimit:
			movingAverageMatrix[i] = movingAverageMatrix[i] + inputMatrix[i+j]
			j = j+1
		movingAverageMatrix[i] = movingAverageMatrix[i] / width
		i = i+1
		
	# Handle the moving average of elements that are close the edges
	# The moving averages here will be the same as the moving average of the 
	# first element a distance of width/2 away from the edge
	
	# At the start of the output matrix
	j = 0
	while j < upperLimit:
		movingAverageMatrix[j] = movingAverageMatrix[upperLimit]
		j = j+1
	
	# At the end of the output matrix
	j = len(inputMatrix) - width/2
	referenceElement = len(inputMatrix) - width/2 - 1
	while j < len(inputMatrix):
		movingAverageMatrix[j] = movingAverageMatrix[referenceElement]
		j = j+1
	
	return movingAverageMatrix
# ---------------------------------------------------------------------------- #
	


# ---------------------------------------------------------------------------- #
# Function to apply correction factors to y data
def correct(datasets, scaleData):
	import numpy
	scaleDataY = scaleData[:,1]
	j = 0
	correctedDatasets = []
	while j < len(datasets):
		data = datasets[j]
		correctedData = float64(empty((len(data), 2)))
		k = 0
		while k < len(correctedData):
			cdat = float64(data[k,1])*float64(scaleDataY[k])
			correctedData[k,1] = cdat
			correctedData[k,0] = data[k,0]
			k = k+1
		correctedDatasets.append(correctedData)
		j = j+1
	return correctedDatasets
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Get filenames and collect data from files 
def importRawData(argv):
	import scipy
	import scipy.io.array_import
	import re
	
	tempMatch = re.compile("-*\d+C")
	numberMatch = re.compile("-*\d+")
		
	i = 0
	datalabels = []
	datasets = []
	temperatures = []
	while i < len(argv):
		filename = argv[i]
		datalabel = filename.split('.')
		datalabels.append(datalabel[0])
		data = scipy.io.array_import.read_array(filename,lines=(17,(18,2065)))
		datasets.append(data)
		
		tempMatchObj = tempMatch.search(datalabel[0])
		
		if tempMatchObj == None:
			temperature = "?"
		else:
			temperature = float(numberMatch.search(tempMatchObj.group()).group())
		
		temperatures.append(temperature)
		
		i = i+1
		
	return [datasets, datalabels, temperatures]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Get filenames and data from files 
def importRawData2(argv,startLine,endLine):
	import scipy
	from scipy.io.array_import import read_array 
	import re
	
	i = 0
	datalabels = []
	datasets = []
	while i < len(argv):
		filename = argv[i]
		datalabel = filename.split('.')
		datalabels.append(datalabel[0])
		data = read_array(filename,lines=(startLine,(startLine+1,endLine)))
		datasets.append(data)		
		i = i+1
		
	return [datasets, datalabels]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# Function to determine index of upper and lower wavelength limits
# The function will work with a multicolumn matrix, with the wavelength data in
# the first column or with a vector of wavelengths
def limits(lowerLimit, upperLimit, data):
	lower = lowerLimit
	upper = upperLimit
	
	if len(data.shape) == 1:
		x = data
	elif len(data.shape) > 1:
		x = data[:,0]
	
	i = 0
	lowerFound = False
	upperFound = False
	
	while i < len(x) and upperFound == False:
		if int(x[i]) == lower and lowerFound == False:
			lowerFound == True
			lowerIndex = i
		if int(x[i]) == upper and upperFound == False:
			upperFound == True
			upperIndex = i
		i = i+1

	return [lowerIndex, upperIndex]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to determine index of upper and lower wavelength limits
# The function will work with a multicolumn matrix, with the wavelength data in
# the first column or with a vector of wavelengths
# The function is a modified version of limits, that can handle floating point
# limits.
def limits2(lowerLimit, upperLimit, data, tol=0.2):
	lower = float(lowerLimit)
	upper = float(upperLimit)
	
	
	if len(data.shape) == 1:
		x = data
	elif len(data.shape) > 1:
		x = data[:,0]
	
	i = 0
	lowerFound = False
	upperFound = False
	
	while i < len(x) and upperFound == False:
		if (x[i] - tol) < lower < (x[i] + tol) and lowerFound == False:
			lowerFound == True
			lowerIndex = i
		if (x[i] - tol) < upper < (x[i] + tol) and upperFound == False:
			upperFound == True
			upperIndex = i
		i = i+1

	return [lowerIndex, upperIndex]
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
# Function to perform singular valule decomposition of a M x N numpy matrix, 
# and return the decomposed numpy matrices
# The first matrix to be returned is the MxN matrix of basis vectors
# The second is a NxN matrix with elements only on its diagonal
# The third matrix is a NxN matrix also
def SVD(inputMatrix):
	# Use the Linear algebra package from numpy to perform the singular value
	# decomposition
	svdResult = LinearAlgebra.singular_value_decomposition(inputMatrix)
	
	# Unfortunately, the numpy linear algebra routine returns the diagonal
	# matrix as a vector. The purpose of this function is to convert it to
	# a matrix
	size = svdResult[1].shape[0]
	diagonalMatrix = zeros([size, size],Float)
	i = 0
	while i < len(svdResult[1]):
		diagonalMatrix[i,i] = svdResult[1][i]
		i = i+1
	
	PMatrix = dot(diagonalMatrix,thirdMatrix)
	basisMatrix = svdResult[0]

	return [basisMatrix,PMatrix]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to return a truncated matrix multiplication
def TruncMatMult(m1, m2,N):
	if m1.shape[1] != m2.shape[0]:
		print("Matrix sizes do not match.")
		return None
	
	if N > m1.shape[1]:
		print("N must be less than the number of columns in matrix 1")
		return None
	
	a = zeros([m1.shape[0], m2.shape[1]],Float)
	
	q = 0
	k = 0
		
	while q < m1.shape[0]:
		k = 0
		while k < m2.shape[1]:
			j = 0
			while j < N:
				a[q,k] = a[q,k]+ m1[q,j]*m2[j,k]
				j = j+1
			k = k+1
		q = q+1
		
	return a
# ---------------------------------------------------------------------------- #

# -----------------------------------------------------------------------------#
# Functions for calculating a linear fit to a peak versus pressure plot
def linearFunction(x,p):
	fit = p[0] + p[1]*x
	return fit
	
def linearResiduals(p, y, x): 
	err = y - linearFunction(x,p) 
	return err

def linearFit(lowerLimit, upperLimit, x,y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p01 = (1.0/14500.0)
	
	p0 = array([525.0, p01])
	plsq = leastsq(linearResiduals, p0, args=(y, x), maxfev=2000)
	
# 	rplsq = []
# 	i=1
# 	while i <= len(plsq[0]):
# 		rplsq.append(plsq[0][-i])
# 		i = i+1
	
	return [plsq, linearFunction(x, plsq[0])]


# -----------------------------------------------------------------------------#
# Functions for calculating a polynomial fit to the peak of spectrum
def quadfunction(x,p): 
	fit = p[0] + p[1]*x + p[2]*(x**2)
	return fit 

def quadresiduals(p, y, x): 
	err = y - quadfunction(x,p) 
	return err

def quadfit(lowerLimit, upperLimit, x,y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025.,1051.,1.])
	plsq = leastsq(quadresiduals, p0, args=(y, x), maxfev=2000)
	
	rplsq = []
	i=1
	while i <= len(plsq[0]):
		rplsq.append(plsq[0][-i])
		i = i+1
		
	diff = scipy.polyder(rplsq,1)
	roots = scipy.roots(diff)
	
	peak = 0
	for root in roots:
		if root < x[-1] and root > x[0]:
			peak = root
	
	fit = []
	
	for i in  x:
		fit.append((i,quadfunction(i,plsq[0])))
	
	return float(real(peak))
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
def quadint(lowerLimit, upperLimit, x, y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025.,1051.,1.])
	plsq = leastsq(quadresiduals, p0, args=(y, x), maxfev=2000)
	
	rplsq = []
	i=1
	while i <= len(plsq[0]):
		rplsq.append(plsq[0][-i])
		i = i+1
		
	diff = scipy.polyder(rplsq,1)
	roots = scipy.roots(diff)
	
	peak = 0
	for root in roots:
		if root < x[-1] and root > x[0]:
			peak = root
	
	peakIntensity = quadfunction(peak,plsq[0])
		
	return float(real(peakIntensity))
# -----------------------------------------------------------------------------#


# -----------------------------------------------------------------------------#
def quadFitVector(lowerLimit, upperLimit, x,y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025,1051,-1,1,100])
	plsq = leastsq(quadresiduals, p0, args=(y, x), maxfev=2000)
	
	return quadfunction(x, plsq[0])
# -----------------------------------------------------------------------------#
	




# -----------------------------------------------------------------------------#
# Functions for calculating a polynomial fit to the peak of spectrum
def polyfunction(x,p): 
	fit = p[0] + p[1]*x + p[2]*(x**2) + p[3]*(x**3) + p[4]*(x**4)
	return fit 

def polyresiduals(p, y, x): 
	err = y - polyfunction(x,p) 
	return err

def polyfit(lowerLimit, upperLimit, x,y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025,1051,-1,1,100])
	plsq = leastsq(polyresiduals, p0, args=(y, x), maxfev=2000)
	
	rplsq = []
	i=1
	while i <= len(plsq[0]):
		rplsq.append(plsq[0][-i])
		i = i+1
		
	diff = scipy.polyder(rplsq,1)
	roots = scipy.roots(diff)
	
	peak = 0
	for root in roots:
		if root < x[-1] and root > x[0]:
			peak = root
	
	fit = []
	
	for i in  x:
		fit.append((i,polyfunction(i,plsq[0])))
	
	return float(real(peak))


def polyFitVector(lowerLimit, upperLimit, x,y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025,1051,-1,1,100])
	plsq = leastsq(polyresiduals, p0, args=(y, x), maxfev=2000)
	
	return polyfunction(x, plsq[0])
	

def polyint(lowerLimit, upperLimit, x, y):
	import scipy
	from scipy.optimize import leastsq
	
	y = y[lowerLimit:upperLimit]
	x = x[lowerLimit:upperLimit]
	
	maxY = max(y)
	
	p0 = array([-276025,1051,-1,1,100])
	plsq = leastsq(polyresiduals, p0, args=(y, x), maxfev=2000)
	
	rplsq = []
	i=1
	while i <= len(plsq[0]):
		rplsq.append(plsq[0][-i])
		i = i+1
		
	diff = scipy.polyder(rplsq,1)
	roots = scipy.roots(diff)
	
	peak = 0
	for root in roots:
		if root < x[-1] and root > x[0]:
			peak = root
	
	peakIntensity = polyfunction(peak,plsq[0])
		
	return float(real(peakIntensity))
# -----------------------------------------------------------------------------#


# ---------------------------------------------------------------------------- #
# Functions for performing center of spectral mass calculation
def cosm(lowerLimit, upperLimit, x, y):
	
	y = y[lowerLimit:upperLimit,1]
	x = x[lowerLimit:upperLimit,0]
	maxY = 1.1* max(y)
	
	# Integrate the data using a trapezoidal approximation
	integral = 0
	moment = 0
	i = 0
	while i < len(x)-1:
		dx = x[i+1] - x[i]
		dI = dx * (y[i] + 0.5*(y[i+1]-y[i]))
		integral = integral + dI
		moment = moment + 0.5*(x[i]+ x[i+1])*dI
		i = i+1
	cosm = moment / integral
		
	return cosm
# -----------------------------------------------------------------------------#

# -----------------------------------------------------------------------------#
# Functions for performing sort of datasets
# Sort order is pressure frozen lowest temperature to warmest, then cooling down 
# warmest temperature to lowest.
def temperatureArraySorter(list1, list2):
	import re
	temperature1 = list1[0]
	temperature2 = list2[0]
	
	dataLabel1 = list1[1]
	dataLabel2 = list2[1]
	
	PFMatch = re.compile("PF-*\d+C")
	CDMatch = re.compile("cooldown-*\d+C")
	
	dL1PFMatch = PFMatch.search(dataLabel1)
	dL1CDMatch = CDMatch.search(dataLabel1)
	dL2PFMatch = PFMatch.search(dataLabel2)
	dL2CDMatch = CDMatch.search(dataLabel2)
	
	if dL1PFMatch == None and dL1CDMatch != None:
		dL1State = 'CD'
	elif dL1PFMatch != None and dL1CDMatch == None:
		dL1State = 'PF'
	
	if dL2PFMatch == None and dL2CDMatch != None:
		dL2State = 'CD'
	elif dL2PFMatch != None and dL2CDMatch == None:
		dL2State = 'PF'
	
	if (dL1State == dL2State) and (temperature1 == temperature2):
		return 0
	elif (dL1State == dL2State):
		if dL1State == 'PF':
			return cmp(temperature1, temperature2)
		elif dL1State == 'CD':
			return -1*cmp(temperature1,temperature2)
	elif dL1State == 'PF' and dL2State == 'CD':
		return -1
	elif dL1State == 'CD' and dL2State == 'PF':
		return 1
# -----------------------------------------------------------------------------#

# ---------------------------------------------------------------------------- #
# Function to generate a filelist for import
def GenerateFileList(directory=".", regex=".*\.ProcSpec.dat", ignoreCase=True):
	import os
	import re
	fileList = os.listdir(directory)
	
	if ignoreCase==True:
		filePattern = re.compile(regex, re.IGNORECASE)
	else:
		filePattern = re.compile(regex)

	i = 0
	selectFiles = []

	while i < len(fileList):
		if filePattern.match(fileList[i]) != None:
			selectFiles.append(fileList[i])
		i = i+1
		
	return selectFiles
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# Function to import raw data, correct it and return matrix of corrected 
# intensities and wavelengths
def importData(fileList, \
scaleFactors="/Users/buz/bin/craneUSB2000XtalScaleFactors.dat"):
	
	import scipy
	
	# Load in correction factor data
	scaleData = scipy.io.array_import.read_array(scaleFactors,lines=(3,(4,2051)))
	
	# Get filenames and collect data from files 
	[datasets, datalabels, temperatures] = importRawData(fileList)

	# Apply correction factors
	correctedDatasets = correct(datasets, scaleData)
	
	# Sort datasets by temperature
	i = 0
	bigArray = []
	
	while i < len(temperatures):
		bigArray.append([temperatures[i], datalabels[i], correctedDatasets[i]])
		i = i+1
	
	bigArray.sort(temperatureArraySorter)
	
	i = 0
	while i < len(bigArray):
		temperatures[i] = bigArray[i][0]
		datalabels[i] = bigArray[i][1]
		correctedDatasets[i] = bigArray[i][2]
		i = i+1
	
	# Generate a matrix of wavelengths for use with the matrix of intensities
	wavelengthMatrix = zeros(len(correctedDatasets[0]),Float)
	i = 0
	while i < wavelengthMatrix.shape[0]:
		wavelengthMatrix[i] = correctedDatasets[0][i][0]
		i = i+1
	
	# Generate a matrix of intensities for use in singular value decomposition
	# Each set spectrum is given its own 2048 row column in the matrix
	intensitiesMatrix = zeros([len(correctedDatasets[0]),len(correctedDatasets)],Float)
	i = 0
	while i < intensitiesMatrix.shape[1]:
		j = 0
		while j < intensitiesMatrix.shape[0]:
			intensitiesMatrix[j,i] = correctedDatasets[i][j][1]
			j = j+1
		i = i+1
		
	return [wavelengthMatrix, intensitiesMatrix, temperatures]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to import raw data, leave it uncorrected and return matrix of 
# intensities and a matrix of wavelengths
def importData2(fileList,startLine=17,endLine=2065):
	
	import scipy
	
	# Get filenames and collect data from files 
	[datasets, datalabels] = importRawData2(fileList,startLine,endLine)

	# Generate a matrix of wavelengths for use with the matrix of intensities
	wavelengthMatrix = zeros(len(datasets[0]),Float)
	i = 0
	while i < wavelengthMatrix.shape[0]:
		wavelengthMatrix[i] = datasets[0][i][0]
		i = i+1
	
	# Generate a matrix of intensities for use in singular value decomposition
	# Each set spectrum is given its own 2048 row column in the matrix
	intensitiesMatrix = zeros([len(datasets[0]),len(datasets)],Float)
	i = 0
	while i < intensitiesMatrix.shape[1]:
		j = 0
		while j < intensitiesMatrix.shape[0]:
			intensitiesMatrix[j,i] = datasets[i][j][1]
			j = j+1
		i = i+1
		
	return [wavelengthMatrix, intensitiesMatrix, datalabels]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to automatically determine start and end lines of SpectraSuite data
# import raw data, leave it uncorrected and return matrix of intensities and a 
# matrix of wavelengths
def importData3(fileList, baseDir=""):
	
	#i = 0
	#while i < len(fileList):
	#	fileList[i] = baseDir + fileList[i]
	#	i += 1
	
	
	
	[startLines, endLines] = findStartEndLines(fileList,directory=baseDir)
	[datasets, datalabels] = importRawData3(fileList,startLines,endLines,\
	baseDir=baseDir)
	
	# Generate a matrix of wavelengths for use with the matrix of intensities
	wavelengthMatrix = zeros(len(datasets[0]),Float)
	i = 0
	while i < wavelengthMatrix.shape[0]:
		wavelengthMatrix[i] = datasets[0][i][0]
		i = i+1
	
	# Generate a matrix of intensities for use in singular value decomposition
	# Each set spectrum is given its own 2048 row column in the matrix
	intensitiesMatrix = zeros([len(datasets[0]),len(datasets)],Float)
	i = 0
	while i < intensitiesMatrix.shape[1]:
		j = 0
		while j < intensitiesMatrix.shape[0]:
			intensitiesMatrix[j,i] = datasets[i][j][1]
			j = j+1
		i = i+1
		
	return [wavelengthMatrix, intensitiesMatrix, datalabels]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# Function to automatically determine start and end lines of SpectraSuite data
# import raw data, correct it and return matrix of intensities and a 
# matrix of wavelengths
def importData4(fileList, baseDir="", \
scaleFactors="/Users/buz/bin/usb2000scalefactors.dat"):
	import scipy
	

	# Load in correction factor data
	scaleData = scipy.io.array_import.read_array(scaleFactors)
	
	
	[startLines, endLines] = findStartEndLines(fileList, directory=baseDir)
	[datasets, datalabels] = importRawData3(fileList,startLines,endLines,\
	baseDir=baseDir)
	
	# Apply correction factors
	correctedDatasets = correct(datasets, scaleData)
	
	
	# Generate a matrix of wavelengths for use with the matrix of intensities
	wavelengthMatrix = zeros(len(correctedDatasets[0]),Float)
	i = 0
	while i < wavelengthMatrix.shape[0]:
		wavelengthMatrix[i] = correctedDatasets[0][i][0]
		i = i+1
	
	# Generate a matrix of intensities for use in singular value decomposition
	# Each set spectrum is given its own 2048 row column in the matrix
	intensitiesMatrix = zeros([len(correctedDatasets[0]),\
	len(correctedDatasets)],Float)
	i = 0
	while i < intensitiesMatrix.shape[1]:
		j = 0
		while j < intensitiesMatrix.shape[0]:
			intensitiesMatrix[j,i] = correctedDatasets[i][j][1]
			j = j+1
		i = i+1
		
	return [wavelengthMatrix, intensitiesMatrix, datalabels]
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
# Get filenames and data from files with variable start and end lines
def importRawData3(argv,startLines,endLines,baseDir=""):
	import scipy
	from scipy.io.array_import import read_array 
	import re
	
	i = 0
	datalabels = []
	datasets = []
	while i < len(argv):
		filename = baseDir + argv[i]
		datalabel = filename.split('.')
		datalabels.append(datalabel[0])
		data = read_array(filename,\
		lines=(startLines[i]-1,(startLines[i],endLines[i])))
		datasets.append(data)		
		i = i+1
		
	return [datasets, datalabels]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def importChronosData(fileList, directory='.'):
	
	import re
		
	# Find the file starts
	[startLines, endLines] = findChronosDataStarts(fileList, directory=directory)
	
	# Make a pressure list
	pressureList = findPressureFromChronosFile(fileList, directory=directory)
	
	# Import the intensity data
	[wavelengthMatrix, iMatrix, stdErrorMatrix] = \
	importRawChronosData(fileList, startLines, endLines, directory=directory)
	
	
	
	return [pressureList, wavelengthMatrix, iMatrix, stdErrorMatrix]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to import raw Chronos data
def importRawChronosData(fileList, startLines, endLines, directory='.'):
	
	wavelengthMatrix = []
	iMatrix = []
	stdErrorsMatrix = []	
	
	i = 0
	while i < len(fileList):
		wavelengths = []
		intensities = []
		stdErrors = []
		
		fHandle = open(directory + '/' + fileList[i], 'r')
		data = fHandle.readlines()
		
		j = startLines[i] + 1
		
		while j < endLines[i]:
			dataCols = data[j].split()
			wavelengths.append(float(dataCols[0]))
			intensities.append(float(dataCols[1]))
			if len(dataCols) >= 3:
				stdErrors.append(float(dataCols[2]))
			
			j = j + 1
		
		wavelengthMatrix.append(wavelengths)
		iMatrix.append(intensities)
		stdErrorsMatrix.append(stdErrors)
		
		i = i + 1
	
	return [wavelengthMatrix, iMatrix, stdErrorsMatrix]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# Function to find pressure from header information in chronos data file
def findPressureFromChronosFile(fileList, directory='.'):
	import re
	
	pRegex = re.compile('\d+\.*\d*\s*psi')
	numberRe = re.compile('\d+\.*\d*')
	
	pressureList = []
	i = 0
	
	while i < len(fileList):
		fileHandle = open(directory + '/' + fileList[i], 'r')
		data = fileHandle.readlines()
		fileHandle.close()
		
		j = 0
		
		pressureFound = False
		
		while j < len(data):
			pressureMatch = pRegex.search(data[j])
			if pressureMatch != None and pressureFound == False:
				pressureData = pressureMatch.group()
				pressure = numberRe.search(pressureData).group()
				pressureList.append(pressure)
				pressureFound = True
			j = j + 1
			
		i = i + 1
	
	return pressureList
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Function to find start of data in chronos data file
def findChronosDataStarts(fileList, directory='.'):
	import re
	
	# Import chronos data
	dataArray = []
	
	for file in fileList:
		fileHandle = open(directory + '/' + file,'r')
		data = fileHandle.readlines()
		dataArray.append(data)
		fileHandle.close()

	
	# Find the start of the data section in the datafile
	dataStart = re.compile('[data]')
	
	dataStartLines = []
	dataEndLines = []
	
	j = 0
	while j < len(dataArray):
		data = dataArray[j]
		dataEndLine = len(data)
		
		i = 0
		while i < len(data):
			if dataStart.search(data[i]):
				dataStartLine = i
			i = i+1
		
		dataStartLines.append(dataStartLine)
		dataEndLines.append(dataEndLine)
		j = j + 1
	
	
	return [dataStartLines, dataEndLines]
# ---------------------------------------------------------------------------- #

	
	
# ---------------------------------------------------------------------------- #
# Function to automatically determine start and end lines of SpectraSuite data
def findStartEndLines(fileList, directory='.'):
	import re
	numberMatch = re.compile('\d')
	
	startLines = []
	endLines = []
	
	for file in fileList:
		fileHandle = open(directory + '/' + file,'r')
		lines = fileHandle.readlines()
		i = 0
		numberLines = []
		while i < len(lines):
			if numberMatch.match(lines[i]):
				numberLines.append(i+1)
			i = i+1
		fileHandle.close()
		startLines.append(numberLines[0])
		endLines.append(numberLines[-1])

	return [startLines, endLines]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Function to truncate wavelength and intensity matrices to cover reduced 
# wavelength range
def DataTruncate(wavelengthMatrix, intensitiesMatrix, lowerLimit, upperLimit):
	lims = limits(500,600,wavelengthMatrix)
	trWavelengthMatrix = wavelengthMatrix[lims[0]:lims[1]]
	trIntensitiesMatrix = intensitiesMatrix[lims[0]:lims[1],:]
	return [trWavelengthMatrix, trIntensitiesMatrix]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Function to return peaks of fluorescence spectra in an intensities matrix by 
# polynomial fitting
def QuadPeaks(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,upperLimit=540):
	i = 0
	length = len(intensitiesMatrix[0,:])
	peaks = zeros(length,Float)
	dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	while i < length:
		peak = quadfit(dataLimits[0], dataLimits[1], wavelengthMatrix,\
		intensitiesMatrix[:,i])
		peaks[i] = peak
		i = i+1
	return peaks
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to refine fit by polynomial fitting
# If function cannot find a peak inside the range of the lower to upper limits,
# it will not return an answer for that entry in the intensity matrix.
# It will return a list of all of the indices that could be fitted
def QuadPeakIntsRefine(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,\
upperLimit=540,fineFitWindow=5.0):
	
	i = 0
	length = intensitiesMatrix.shape
	if len(length) > 1:
		length = length[1]
	else:
		length = 1
	
	peaks = []
	intensities = []
	dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	trWMatrices = []
	peakFitMatrix = []
	indexList = []
	
	print("Length of iMatrix: " + str(length))
	
	while i < length:
		
		# Find the approximate position of the peak
		peakCoarse = quadfit(dataLimits[0], dataLimits[1], wavelengthMatrix, \
		intensitiesMatrix[:,i])
		
		if lowerLimit < peakCoarse < upperLimit:		
			# Calculate a fine window for the fit
			upperLimitFine = peakCoarse + fineFitWindow/2.0
			lowerLimitFine = peakCoarse - fineFitWindow/2.0
			
			dataLimitsFine = limits2(lowerLimitFine,upperLimitFine,wavelengthMatrix)
			
			peakFine = quadfit(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			intensity = quadint(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			
			peaks.append(peakFine)
			intensities.append(intensity)
			
			trWMatrix = wavelengthMatrix[dataLimitsFine[0]:dataLimitsFine[1]]
			trWMatrices.append(trWMatrix)
			
			quadfitvector = quadFitVector(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			peakFitMatrix.append(quadfitvector)
			indexList.append(i)
		
		i = i+1
	
	
	return [indexList, peaks, intensities, trWMatrices, peakFitMatrix]
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
# Function to return peaks of fluorescence spectra in an intensities matrix by 
# polynomial fitting
def PolyPeaks(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,upperLimit=540):
	i = 0
	length = len(intensitiesMatrix[0,:])
	peaks = zeros(length,Float)
	dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	while i < length:
		peak = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix,\
		intensitiesMatrix[:,i])
		peaks[i] = peak
		i = i+1
	return peaks

# Function to be used with single intensity vectors; useful for Chronos data
# where spectrum files often have different lengths
def PolyPeak(wavelengthVector,iVector,lowerLimit=520,upperLimit=540):
	i = 0
	length = len(iVector)
	dataLimits = limits2(lowerLimit,upperLimit,wavelengthVector)
	peak = polyfit(dataLimits[0], dataLimits[1], wavelengthVector,iVector)
	return peak


def PolyPeakInts(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,upperLimit=540):
	i = 0
	
	if len(intensitiesMatrix.shape) > 1:
		length = intensitiesMatrix.shape[1]
		peaks = zeros(length,Float)
		intensities = zeros(length, Float)
		dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	
		while i < length:
			peak = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix,\
			intensitiesMatrix[:,i])
			intensity = polyint(dataLimits[0], dataLimits[1], wavelengthMatrix,\
			intensitiesMatrix[:,i])
			peaks[i] = peak
			intensities[i] = intensity
			i = i+1
	else:
		dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
		peak = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix,\
		intensitiesMatrix)
		intensity = polyint(dataLimits[0], dataLimits[1], wavelengthMatrix,\
		intensitiesMatrix)
		peaks = peak
		intensities = intensity
		
	return [peaks, intensities]

def PolyPeakFits(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,upperLimit=540):
	i = 0
	
	if len(intensitiesMatrix.shape) > 1:	
		length = len(intensitiesMatrix[0,:])
		dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
		trWMatrix = wavelengthMatrix[dataLimits[0]:dataLimits[1]]
		peakFits = zeros((trWMatrix.shape[0], intensitiesMatrix.shape[1]),Float)
		
		while i < length:
			polyfitvector = polyFitVector(dataLimits[0], dataLimits[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			peakFits[:,i] = polyfitvector
			i = i+1
	else:
		dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
		trWMatrix = wavelengthMatrix[dataLimits[0]:dataLimits[1]]
		polyfitvector = polyFitVector(dataLimits[0], dataLimits[1], \
		wavelengthMatrix, intensitiesMatrix)
		peakFits = polyfitvector
		
	return [trWMatrix, peakFits]
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Function to refine fit by polynomial fitting
# If function cannot find a peak inside the range of the lower to upper limits,
# it will not return an answer for that entry in the intensity matrix.
# It will return a list of all of the indices that could be fitted
def PolyPeakIntsRefine(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,\
upperLimit=540,fineFitWindow=5.0, shapeByPass=False):
	
	i = 0
	length = intensitiesMatrix.shape
	if len(length) > 1:
		length = length[1]
	else:
		length = 1
	
	peaks = []
	intensities = []
	dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	trWMatrices = []
	peakFitMatrix = []
	indexList = []
	
	print("Length of iMatrix: " + str(length))
	
	while i < length:
		
		# Find the approximate position of the peak
		peakCoarse = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix, intensitiesMatrix[:,i])
		
		if lowerLimit < peakCoarse < upperLimit:		
			# Calculate a fine window for the fit
			upperLimitFine = peakCoarse + fineFitWindow/2.0
			lowerLimitFine = peakCoarse - fineFitWindow/2.0
			
			dataLimitsFine = limits2(lowerLimitFine,upperLimitFine,wavelengthMatrix)
			
			peakFine = polyfit(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			intensity = polyint(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			
			peaks.append(peakFine)
			intensities.append(intensity)
			
			trWMatrix = wavelengthMatrix[dataLimitsFine[0]:dataLimitsFine[1]]
			trWMatrices.append(trWMatrix)
			
			polyfitvector = polyFitVector(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, intensitiesMatrix[:,i])
			peakFitMatrix.append(polyfitvector)
			indexList.append(i)
		
		i = i+1
	
	
	return [indexList, peaks, intensities, trWMatrices, peakFitMatrix]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
# Function to refine fit by polynomial fitting
# If function cannot find a peak inside the range of the lower to upper limits,
# it will not return an answer for that entry in the intensity matrix.
# It will return a list of all of the indices that could be fitted

def PolyPeakIntsRefine2(wavelengthMatrix,intensitiesMatrix,lowerLimit=520,\
upperLimit=540,fineFitWindow=30.0):
	
	i = 0
	length = intensitiesMatrix.shape
	if len(length) > 1:
		length = length[1]
	else:
		length = 1
	
	peaks = []
	intensities = []
	dataLimits = limits(lowerLimit,upperLimit,wavelengthMatrix)
	trWMatrices = []
	peakFitMatrix = []
	indexList = []
	
	print("Length of iMatrix: " + str(length))
	
	while i < length:
		
		# Find the approximate position of the peak
		if length > 1:
			peakCoarse = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix, \
			intensitiesMatrix[:,i])
		else:
			peakCoarse = polyfit(dataLimits[0], dataLimits[1], wavelengthMatrix, intensitiesMatrix)
		
		print("Peak Coarse: " + str(peakCoarse))
		
		if lowerLimit < peakCoarse < upperLimit:		
			# Calculate a fine window for the fit
			upperLimitFine = peakCoarse + fineFitWindow/2.0
			lowerLimitFine = peakCoarse - fineFitWindow/2.0
			
			dataLimitsFine = limits2(lowerLimitFine,upperLimitFine,wavelengthMatrix,tol=1.0)
			
			print(str(dataLimitsFine))
			
			if length > 1:
				iMatrixTemp = intensitiesMatrix[:,i]
			else:
				iMatrixTemp = intensitiesMatrix
			
			peakFine = polyfit(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, iMatrixTemp)
			
			print("Peak Fine: " +  str(peakFine))
			
			intensity = polyint(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, iMatrixTemp)
			
			print("Peak Intensity: " + str(intensity))
			
			peaks.append(peakFine)
			intensities.append(intensity)
			
			trWMatrix = wavelengthMatrix[dataLimitsFine[0]:dataLimitsFine[1]]
			trWMatrices.append(trWMatrix)
			
			polyfitvector = polyFitVector(dataLimitsFine[0], dataLimitsFine[1], \
			wavelengthMatrix, iMatrixTemp)
			peakFitMatrix.append(polyfitvector)
			indexList.append(i)
		
		i = i+1
	
	
	return [indexList, peaks, intensities, trWMatrices, peakFitMatrix]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def BackgroundSubtraction(wMatrix, iMatrix, fitStartIndex=125, fitEndIndex=139, ):
	# Performs Lorentzian based background subtraction on an intensity matrix
	import numpy
	
	fitWavelengths = numpy.arange(483,650,0.25)
	fitIntensityMatrix = numpy.zeros([fitWavelengths.shape[0], iMatrix.shape[1]])
	sIMatrix = numpy.zeros([iMatrix.shape[0], iMatrix.shape[1]])
	plsqMatrix = []
	
	i = 0
	while i < iMatrix.shape[1]:
		iVector = iMatrix[:,i]
		plsq = scatteringFit(fitStartIndex, fitEndIndex, wMatrix, iVector)
		fitIntensityMatrix[:,i] = scatteringFunction(fitWavelengths, plsq[0])
		sIMatrix[:,i] = iVector - scatteringFunction(wMatrix, plsq[0])
		plsqMatrix.append(plsq)
		i = i+1
	
	return [sIMatrix, fitWavelengths, fitIntensityMatrix]
	
# -----------------------------------------------------------------------------#
# Functions for calculating a polynomial fit to the wings of the scattering peak
# in a spectrum	
def scatteringFunction(x,p): 
	fit = p[0]/( p[1] + p[2]*(x**2) )
	return fit 

def scatteringResiduals(p, y, x): 
	err = y - scatteringFunction(x,p) 
	return err

def scatteringFit(lowerIndex, upperIndex, x,y):
	import scipy
	
	y = y[lowerIndex:upperIndex]
	x = x[lowerIndex:upperIndex]
	
	p0 = array([57689.4, -111290.0, 0.461074])
	plsq = scipy.optimize.leastsq(scatteringResiduals, p0, args=(y, x), maxfev=20000)
	
	return plsq
# -----------------------------------------------------------------------------#



# ------------------------------------------------------------------------------------------------ #
# Background subtraction algorithm
def subtract_background2(dataDict, blankKey, dataKey, zeroOffset=0, UVvisCrossover=349, \
UVvisMisMatchCorrection=True):
# 2012-3-27
# Developed for background subtraction of data taken on Cary 300 UV/vis spectrophotometer
	from numpy import zeros,ones
	
	# Get the wavelength matrix
	
	blankMatrix = dataDict[blankKey]
	rawDataMatrix = dataDict[dataKey]
	
	wlBlank = blankMatrix[:,0]
	wlRawData = rawDataMatrix[:,0]
	
	blankData = blankMatrix[:,1]
	rawData = rawDataMatrix[:,1]
	
	# Test that the wavelengths of the blank and data are the same
	
	i = 0
	currentBlank = wlBlank[0]
	
	dataBGSVector = []
	wlBGSVector = []
	
	
	while i < len(wlRawData):
		j = 0
		while j < len(wlRawData):
			if wlRawData[j] == wlBlank[i]:
				wlBGSVector.append(wlRawData[j])
				dataBGSVector.append(rawData[j] - blankData[i]+zeroOffset)
			j += 1
		i += 1
		
	
	dataBGSArray = zeros(len(dataBGSVector), float)
	dataBGSArray = dataBGSVector
	
	dataBGS = zeros((len(wlBGSVector),2))
	dataBGS[:,0] = wlBGSVector
	dataBGS[:,1] = dataBGSVector
	
	if UVvisMisMatchCorrection == True:
		
		#print str(len(dataBGSArray))
		
		# Find the index of the UV vis crossover point
		crossoverIndex1 = findIndex(wlBGSVector, UVvisCrossover)
		#print crossoverIndex1
		crossoverIndex2 = findIndex(wlBGSVector, UVvisCrossover-1)
		crossoverIndex3 = findIndex(wlBGSVector, UVvisCrossover-2)
		#print crossoverIndex2
		uv200Index = findIndex(wlBGSVector, 200)
		
		# End = 800 nm
		
		# 200 nm is at the end of the vector
		if uv200Index > crossoverIndex2:
			uvStartIndex = crossoverIndex2
			uvEndIndex = uv200Index
			visibleStartIndex = 0
			visibleEndIndex = crossoverIndex2
		else:
		# 200 nm is at the start of the vector
			uvStartIndex = uv200Index
			uvEndIndex = crossoverIndex1
			visibleStartIndex = crossoverIndex1
			visibleEndIndex = -1
		
		#print uvStartIndex
		#print uvEndIndex
		#print visibleStartIndex
		#print visibleEndIndex
		
		
		
		difference = (dataBGSVector[crossoverIndex1] - dataBGSVector[crossoverIndex2]) \
		- (dataBGSVector[crossoverIndex2] - dataBGSVector[crossoverIndex3])
		
		dataBGSUVVMMVector = zeros(len(dataBGSArray), float)
		dataBGSUVVMMVector[uvStartIndex:uvEndIndex] = \
		dataBGSArray[uvStartIndex:uvEndIndex] + difference
		
		dataBGSUVVMMVector[visibleStartIndex:visibleEndIndex] = \
		dataBGSArray[visibleStartIndex:visibleEndIndex]
		
		dataBGS[:,1] = dataBGSUVVMMVector
		
	
	return dataBGS
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Background subtraction algorithm
def subtract_background(dataDict, blankKey, dataKey, zeroOffset=0):
	from numpy import zeros,ones
	
	# Get the wavelength matrix
	
	blankMatrix = dataDict[blankKey]
	rawDataMatrix = dataDict[dataKey]
	
	wlBlank = blankMatrix[:,0]
	wlRawData = rawDataMatrix[:,0]
	
	blankData = blankMatrix[:,1]
	rawData = rawDataMatrix[:,1]
	
	# Test that the wavelengths of the blank and data are the same
	
	i = 0
	currentBlank = wlBlank[0]
	
	dataBGSVector = []
	wlBGSVector = []
	
	while i < len(wlRawData):
		j = 0
		while j < len(wlRawData):
			if wlRawData[j] == wlBlank[i]:
				wlBGSVector.append(wlRawData[j])
				dataBGSVector.append(rawData[j] - blankData[i]+zeroOffset)
			j += 1
		i += 1
		
	
	dataBGS = zeros((len(wlBGSVector),2))
	dataBGS[:,0] = wlBGSVector
	dataBGS[:,1] = dataBGSVector
	
	return dataBGS
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def PlotAndCalculateBackgroundSubtractions(dataDict, blankKey, zeroOffsetsDict={},\
vLines=[215,280,570],plotOn=True):

	from matplotlib.pyplot import figure, axvline
	
	names = dataDict.keys()
	bgsDict = {}
	
	for name in names:
		
		
		offset = zeroOffsetsDict[name]
		
		bgsDict[name] = subtract_background(dataDict, blankKey, name, \
		zeroOffset=offset)
		
		if plotOn == True:
			figure()
			plot(dataDict[blankKey][:,0], dataDict[blankKey][:,1])
			plot(dataDict[name][:,0], dataDict[name][:,1])
			plot(bgsDict[name][:,0], bgsDict[name][:,1])
			
			for line in vLines:
				axvline(x=line)
			
			xlabel("Wavelength (nm)")
			ylabel("Absorbance")
			title(name)
			grid()
			
	return bgsDict
# ------------------------------------------------------------------------------------------------ #
	

# ------------------------------------------------------------------------------------------------ #
def findIndex(wlVector, number):
	i = 0
	length = len(wlVector)
	while i < length:
		if wlVector[i] == number:
			return i
		else:
			i += 1
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def CalculateSlopingBackground(data, x1, x2):
	xVector = data[:,0]
	yVector = data[:,1]
	
	x1Index = findIndex(xVector, x1)
	x2Index = findIndex(xVector, x2)
	
	y1 = yVector[x1Index]
	y2 = yVector[x2Index]
	
	a = (y2 - y1)/(x2 - x1)
	b = y2 - a*x2
	
	background = a*xVector+b
	
	return [a,b,background]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def SubtractSlopingBackground(data, background):
	
	from numpy import zeros, float
	
	xVector = data[:,0]
	yVector = data[:,1]
	
	yVectorNew = yVector-background
	
	dataNew = zeros(data.shape, float)
	dataNew[:,0] = xVector
	dataNew[:,1] = yVectorNew
	
	return dataNew
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def LoadCaryHeaders(fileName):
# Loads headers from CSV file from Cary
	dataHandle = open(fileName, 'r')
	
	data = dataHandle.readlines()
	
	dataHandle.close()
	
	headerLine = data[0]
	headerLine = headerLine.strip()
	headers = headerLine.split(',,')
	

	return headers
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PickHeaderPrefixes(headers):
	import re
	
	headerPrefixes = []
	
	# Look for the ending number '_digit' in the header
	
	postfixRe = re.compile('_\d+')
	
	
	
	for header in headers:
		objs = postfixRe.finditer(header)
		
		matchObjs = []
		matchCondition = True
		
		while matchCondition:
			try:
				matchObjs.append(objs.next())
			except StopIteration:
				matchCondition = False
				
		startIndex = matchObjs[-1].start()
		prefix = header[0:startIndex]
		
		headerPrefixes.append(prefix)
	
	return headerPrefixes
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def PickUniqueHeaderPrefixes(headerPrefixes, headers):
	import re
	
	uniquePrefixes = []
	uniquePrefixDict = {}
	
	prevPrefix = None
	currentPrefix = None
	
	i = 0
	
	if len(headers) != len(headerPrefixes):
		print("Something is wrong with lengths of headers and headerPrefixes!")
		return
	
	while i < len(headers):
		prefix = headerPrefixes[i]
		postfix = headers[i].split('_')[-1]
		
		if prefix != currentPrefix:
			uniquePrefixes.append(prefix)
			currentPrefix = prefix
			uniquePrefixDict[currentPrefix] = {}
			uniquePrefixDict[currentPrefix][postfix] = i
		else:
			uniquePrefixDict[currentPrefix][postfix] = i
		
		i +=1 
	
	
	return [uniquePrefixes, uniquePrefixDict]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def AverageSpectrumMatrices(spectrumList):
	spectrumVectorList = []
	
	i = 0
	while i < len(spectrumList):
		spectrumVectorList.append(spectrumList[i][:,1])
		i += 1
	
	averagedWavelengths = spectrumList[0][:,0]
	
	averagedSpectrumVector = zeros(len(averagedWavelengths))
	stdSpectrumVector = zeros(len(averagedWavelengths))
	
	
	i = 0
	while i < len(averagedSpectrumVector):
		j = 0
		spectralPoints = []
		while j < len(spectrumVectorList):
			spectralPoints.append(spectrumVectorList[j][i])
			j+=1
		averagedSpectrumVector[i] = mean(spectralPoints)
		stdSpectrumVector[i] = std(spectralPoints)
		i += 1
	
	averagedSpectrumMatrix = zeros([spectrumList[0].shape[0],3])
	averagedSpectrumMatrix[:,0] = averagedWavelengths
	averagedSpectrumMatrix[:,1] = averagedSpectrumVector
	averagedSpectrumMatrix[:,2] = stdSpectrumVector
	
	return averagedSpectrumMatrix
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def RetrieveSpectrum(data, uniquePrefixDict, key, cycleNumber):
	from numpy import zeros, float
	
	# Get the index in data associated with the key and cyclenumber
	wlIndex = int(uniquePrefixDict[key][cycleNumber])*2
	absIndex = wlIndex+1
	
	wavelengths = data[:,wlIndex]
	absorbances = data[:,absIndex]
	
	spectrumMatrix = zeros([len(wavelengths), 2], float)

	spectrumMatrix[:,0] = wavelengths
	spectrumMatrix[:,1] = absorbances

	return spectrumMatrix


def RetrieveSpectra(data, uniquePrefixDict, key, cycleNumbers):
	
	spectrumMatrices = []
	labels = []
	
	for number in cycleNumbers:
		specMatrix = RetrieveSpectrum(data, uniquePrefixDict, key, number)
		spectrumMatrices.append(specMatrix)
		labels.append(key + '_' + number)
		
	return [spectrumMatrices, labels]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def PlotSpectrumMatrices(spectrumMatrices, labels, titleText, newFig=True, legendOn=True):
	from matplotlib.pyplot import figure, plot, xlabel, ylabel, legend, title
	
	if newFig:
		figure()
	
	i = 0
	while i < len(spectrumMatrices):
		sMatrix = spectrumMatrices[i]
		label = labels[i]
		plot(sMatrix[:,0], sMatrix[:,1], label=label)
		i += 1
		
	xlabel("Wavelength (nm)")
	ylabel("Absorbance")
	title(titleText)
	
	if legendOn:
		legend()
	

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def RetrieveAndPlotSpectrumMatrices(data, uniquePrefixDict, key, cycleNumbers, title, newFig=True):
	[spectrumMatrices, labels] = RetrieveSpectra(data, uniquePrefixDict, key, cycleNumbers)
	PlotSpectrumMatrices(spectrumMatrices, labels, title, newFig=newFig)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def RetrieveAndAverageSpectrumMatrices(data, uniquePrefixDict, key, cycleNumbers):
	[spectrumMatrices, labels] = RetrieveSpectra(data, uniquePrefixDict, key, cycleNumbers)
	averagedSpectrumMatrix = AverageSpectrumMatrices(spectrumMatrices)
	return averagedSpectrumMatrix
# ------------------------------------------------------------------------------------------------ #



