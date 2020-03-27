# ---------------------------------------------------------------------------- #
# vectorUtils14.py
# A module for vector manipulation, originally written for protein structure
# analysis, but later used for lots of things
# Originally created: approximately 6th July 2007
# Last updated: 6th November 2015
# Updated to version 14 or python 3
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def dotProduct(vec1, vec2):
	i = 0
	dotProduct = 0.0
	while i < 3:
		dotProduct = dotProduct + vec1[i]*vec2[i]
		i = i + 1
	return dotProduct
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def nVecLink(x1, x2):
	import math
	i = 0
	nvec = [0.0, 0.0, 0.0]
	while i < 3:
		nvec[i] = x2[i] - x1[i]
		i = i + 1
	i = 0
	mag = math.sqrt(dotProduct(nvec, nvec))
	while i < 3:
		nvec[i] = nvec[i] / mag
		i = i + 1
	
	return nvec
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def vecLink(x1, x2):
	import math
	i = 0
	vec = [0.0, 0.0, 0.0]
	while i < 3:
		vec[i] = x2[i] - x1[i]
		i = i + 1
	i = 0
	
	return vec
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def normalize(x1):
	import math
	normalized = [0.0, 0.0, 0.0]
	i = 0
	mag = math.sqrt(dotProduct(x1, x1))
	while i < 3:
		normalized[i] = x1[i] / mag
		i = i + 1
	
	return normalized
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def vecLength(x1):
	from math import sqrt
	length = sqrt(dotProduct(x1, x1))
	return length
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def generateTransformMatrixForRotation(rotationMatrix):
# Function to generate input matrix for pymol transform_selection matrix
	outputMatrix = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
	0.0, 0.0, 0.0, 0.0]
	

	outputMatrix[0]  = rotationMatrix[0,0]
	outputMatrix[1]  = rotationMatrix[0,1]
	outputMatrix[2]  = rotationMatrix[0,2]
	outputMatrix[3]  = 0.0
	outputMatrix[4]  = rotationMatrix[1,0]
	outputMatrix[5]  = rotationMatrix[1,1]
	outputMatrix[6]  = rotationMatrix[1,2]
	outputMatrix[7]  = 0.0
	outputMatrix[8]  = rotationMatrix[2,0]
	outputMatrix[9]  = rotationMatrix[2,1]
	outputMatrix[10] = rotationMatrix[2,2]
	outputMatrix[11] = 0.0
	outputMatrix[12] = 0.0
	outputMatrix[13] = 0.0
	outputMatrix[14] = 0.0
	outputMatrix[15] = 0.0
	
	return outputMatrix
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def generateTransformMatrixForTranslation(translationVector):
# Function to generate input matrix for pymol transform_selection matrix
	outputMatrix = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, \
	0.0, 0.0, 0.0, 0.0]
	
	outputMatrix[0]  = 1.0
	outputMatrix[1]  = 0.0
	outputMatrix[2]  = 0.0
	outputMatrix[3]  = translationVector[0]
	outputMatrix[4]  = 0.0
	outputMatrix[5]  = 1.0
	outputMatrix[6]  = 0.0
	outputMatrix[7]  = translationVector[1]
	outputMatrix[8]  = 0.0
	outputMatrix[9]  = 0.0
	outputMatrix[10] = 1.0
	outputMatrix[11] = translationVector[2]
	outputMatrix[12] = 0.0
	outputMatrix[13] = 0.0
	outputMatrix[14] = 0.0
	outputMatrix[15] = 0.0
	
	return outputMatrix
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def calculateRotationMatrix(x1, x2, x3, y1, y2, y3):
# Calculates the rotation matrix, r, that maps the orthogonal coordinate system 
# defined by x1, x2 and x3 onto the orthogonal vectors y1, y2 and y3. 
# The coordinate vectors are normalized by this algorithm.

	from numpy import array, transpose, dot
	from numpy.oldnumeric.linear_algebra import inverse

	x1n = normalize(x1)
	x2n = normalize(x2)
	x3n = normalize(x3)
	
	y1n = normalize(y1)
	y2n = normalize(y2)
	y3n = normalize(y3)
	
	
	m = array([x1n, x2n, x3n])
	m = transpose(m)	

	
	n = array([y1n, y2n, y3n])
	n = transpose(n)	
	
	r = dot(m,inverse(n))
	
	return r


# ---------------------------------------------------------------------------- #








# ---------------------------------------------------------------------------- #
# A function for creating a rotation matrix that rotates a vector called
# "from" into another vector called "to".
# Input : fromVector[3], toVector[3] which both must be *normalized* 
# non-zero vectors
# Output: mtx[3][3] -- a 3x3 matrix in colum-major form
# Authors: Tomas Moller, John Hughes
#          "Efficiently Building a Matrix to Rotate One Vector to Another"
#          Journal of Graphics Tools, 4(4):1-4, 1999
# Adapted to Python by Buz Barstow

def fromToRotation(fromVector, toVector):
	epsilon = 0.000001
	from numpy import cross, dot, array, zeros, float
	v = cross(fromVector, toVector)
	e = dot(fromVector, toVector)
	u = zeros(3, float)
	mtx = zeros([3,3], float)
	
	if e < 0:
		f = -1*e
	else:
		f = e
	
	# from and to vectors are almost parallel
	if f > 1.0 - epsilon: 
		x = zeros(3)
		i = 0
		while i < 3:
			if fromVector[i] < 0.0:
				fromVector[i] = - fromVector[i]
			i = i + 1
		
		if x[0] < x[1]:
			if x[0] < x[2]:
				x[0] = 1.0
				x[1] = 0.0
				x[2] = 0.0
			else:
				x[2] = 1.0
				x[0] = 0.0
				x[1] = 0.0
		
		i = 0
		while i < 3:
			u[i] = x[i] - fromVector[i]
			v[i] = x[i] - toVector[i]
			i = i + 1
		
		c1 = 2.0/dot(u,u)
		c2 = 2.0/dot(v,v)
		c3 = c1*c2*dot(u,v)
		
		i = 0
		while i < 3:
			j = 0
			while j < 3:
				mtx[i,j] = -c1*u[i]*u[j] -c2*v[i]*v[j] + c3*v[i]*u[j]
				j = j + 1
			mtx[i,i] = mtx[i,i] + 1.0
			i = i + 1
	# the most common case, unless "from"="to", or "from"=-"to"
	else: 
		h = 1.0/(1.0 + e)
		
		hvx = h*v[0]
		hvz = h*v[2]
		hvxy = hvx * v[1]
		hvxz = hvx * v[2]
		hvyz = hvz * v[1]
		
		mtx[0,0] = e + (hvx*v[0])
		mtx[0,1] = hvxy - v[2]
		mtx[0,2] = hvxz + v[1]

		mtx[1,0] = hvxy + v[2]
		mtx[1,1] = e + h * v[1] * v[1]
		mtx[1,2] = hvyz - v[0]

		mtx[2,0] = hvxz - v[1]
		mtx[2,1] = hvyz + v[0]
		mtx[2,2] = e + hvz * v[2]
		
	return mtx
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def findCOMAndDirsTyr203(selection):
# Function to find the center of mass of the tyrosine 203 ring
# and the direction vectors, CE1 to CE2 and CZ to CG
# The input selection needs to be the tyrosine 203 ring. 
# The function will return the center of the 6 Carbon ring and 
# the direction vectors CZ to CG and CE1 to CE2.
	import re
	from pymol import cmd
	
	tyr203AtomsRe = re.compile('CG|CD1|CD2|CE1|CE2|CZ')
	centerOfMassTyr203 = [0.0, 0.0, 0.0]
	
	tyr203_model = cmd.get_model(selection)
	tyr203_atoms = tyr203_model.atom
	
	for atom in tyr203_atoms:
		search = tyr203AtomsRe.search(atom.name)
		if search != None:
			centerOfMassTyr203[0] = centerOfMassTyr203[0] + atom.coord[0] 
			centerOfMassTyr203[1] = centerOfMassTyr203[1] + atom.coord[1]
			centerOfMassTyr203[2] = centerOfMassTyr203[2] + atom.coord[2]
			if search.group() == 'CE1':
				CE1 = atom.coord
			if search.group() == 'CE2':
				CE2 = atom.coord
			if search.group() == 'CG':
				CG = atom.coord
			if search.group() == 'CZ':
				CZ = atom.coord
	
	centerOfMassTyr203[0] = centerOfMassTyr203[0]/6.0
	centerOfMassTyr203[1] = centerOfMassTyr203[1]/6.0
	centerOfMassTyr203[2] = centerOfMassTyr203[2]/6.0
	
	E1E2 = nVecLink(CE1, CE2)
	GZ = nVecLink(CG, CZ)
	
	return [centerOfMassTyr203, E1E2, GZ]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def findCOMAndDirsCroPhenol(selection):
# Function to find the center of mass of the Chromophore phenol
# and the direction vectors, CE1 to CE2 and CZ to CG2
# The input selection needs to contain the chromophore phenol ring. 
# The function will return the center of the 6 Carbon ring and 
# the direction vectors CZ to CG2 and CE1 to CE2.
	import re
	from pymol import cmd
	
	croAtomsRe = re.compile('CG2|CD1|CD2|CE1|CE2|CZ')
	centerOfMassCro = [0.0, 0.0, 0.0]
	
	cro_model = cmd.get_model(selection)
	cro_atoms = cro_model.atom
	
	for atom in cro_atoms:
		search = croAtomsRe.search(atom.name)
		if search != None:
			centerOfMassCro[0] = centerOfMassCro[0] + atom.coord[0] 
			centerOfMassCro[1] = centerOfMassCro[1] + atom.coord[1]
			centerOfMassCro[2] = centerOfMassCro[2] + atom.coord[2]
			if search.group() == 'CE1':
				CE1 = atom.coord
			if search.group() == 'CE2':
				CE2 = atom.coord
			if search.group() == 'CG2':
				CG2 = atom.coord
			if search.group() == 'CZ':
				CZ = atom.coord
	
	centerOfMassCro[0] = centerOfMassCro[0]/6.0
	centerOfMassCro[1] = centerOfMassCro[1]/6.0
	centerOfMassCro[2] = centerOfMassCro[2]/6.0
	
	E1E2 = nVecLink(CE1, CE2)
	G2Z = nVecLink(CG2, CZ)
	
	return [centerOfMassCro, E1E2, G2Z]
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
def moveCroToPositionAndOrientation(selection, phenolCOMPosition, \
targetE1E2Vector, targetGZVector, state=1):
	# Function to move the chromophore to a particular position and orientation
	
	import re
	from pymol import cmd
	from numpy import array, cos, sin, cross
	
	# Get the current position and orientation of the chromophore
	croCOM_Dir = findCOMAndDirsCroPhenol(selection)
	
	# Generate a transformation series that first takes the chromophore
	# to the origin, aligns it with the coordinate axes, 
	# rotates it to the desired orientation, and then translates the chromphore
	# to the desired location. 
	
	# Transformation to origin
	transVector = -1.0*array(croCOM_Dir[0])
	tTMatrix1 = generateTransformMatrixForTranslation(transVector)
	cmd.transform_selection(selection, tTMatrix1, state)
	
	# Start work on rotating the chromophore
	targetVector = normalize(cross(targetE1E2Vector, targetGZVector))
	croCOM_Dir = findCOMAndDirsCroPhenol(selection)
	origVector = normalize(cross(croCOM_Dir[1],croCOM_Dir[2]))
	
	rMatrix1 = fromToRotation(origVector, targetVector)
	tRMatrix1 = generateTransformMatrixForRotation(rMatrix1)
	cmd.transform_selection(selection, tRMatrix1, state)

	i = 0
	while i < 10:
		# Next, rotate the GZ vector
		croCOM_Dir = findCOMAndDirsCroPhenol(selection)
		rMatrix4 = fromToRotation(croCOM_Dir[2], targetGZVector)
		tRMatrix4 = generateTransformMatrixForRotation(rMatrix4)
		cmd.transform_selection(selection, tRMatrix4, state)
		
		# Rotate the E1E2 vector first
		croCOM_Dir = findCOMAndDirsCroPhenol(selection)
		rMatrix3 = fromToRotation(croCOM_Dir[1], targetE1E2Vector)
		tRMatrix3 = generateTransformMatrixForRotation(rMatrix3)
		cmd.transform_selection(selection, tRMatrix3, state)
		
		i = i + 1
 	
 	# Now, translate the center of mass
	tTMatrix2 = generateTransformMatrixForTranslation(phenolCOMPosition)
	cmd.transform_selection(selection, tTMatrix2, state)
# ---------------------------------------------------------------------------- #






# ---------------------------------------------------------------------------- #
def moveTyr203ToPositionAndOrientation(selection, phenolCOMPosition, \
targetE1E2Vector, targetGZVector, state=1):
	# Function to move the chromophore to a particular position and orientation
	
	import re
	from pymol import cmd
	from numpy import array, cos, sin
	
	# Get the current position and orientation of the chromophore
	tyrCOM_Dir = findCOMAndDirsTyr203(selection)
	
	# Generate a transformation series that first takes the chromophore
	# to the origin, aligns it with the coordinate axes, 
	# rotates it to the desired orientation, and then translates the chromphore
	# to the desired location. 
	
	# Transformation to origin
	transVector = -1.0*array(tyrCOM_Dir[0])
	tTMatrix1 = generateTransformMatrixForTranslation(transVector)
	cmd.transform_selection(selection, tTMatrix1, state)
	
 
	# Rotate to align with y axis with E1E2 vector of chromophore phenol
	tyrCOM_Dir = findCOMAndDirsTyr203(selection)
	rMatrix1 = fromToRotation(tyrCOM_Dir[1], [0,1,0])
	tRMatrix1 = generateTransformMatrixForRotation(rMatrix1)
	cmd.transform_selection(selection, tRMatrix1, state)

	# Recalculate direction vectors of chromophore and then rotate GZ vector
	# with the x axis
	tyrCOM_Dir = findCOMAndDirsTyr203(selection)
	rMatrix2 = fromToRotation(tyrCOM_Dir[2], [1,0,0])
	tRMatrix2 = generateTransformMatrixForRotation(rMatrix2)
	cmd.transform_selection(selection, tRMatrix2, state)
 	 	
 	# Rotate the E1E2 vector first
	tyrCOM_Dir = findCOMAndDirsTyr203(selection)
	rMatrix3 = fromToRotation(tyrCOM_Dir[1], targetE1E2Vector)
	tRMatrix3 = generateTransformMatrixForRotation(rMatrix3)
	cmd.transform_selection(selection, tRMatrix3, state)
 	
	# Next, rotate the GZ vector
	tyrCOM_Dir = findCOMAndDirsTyr203(selection)
	rMatrix4 = fromToRotation(tyrCOM_Dir[2], targetGZVector)
	tRMatrix4 = generateTransformMatrixForRotation(rMatrix4)
	cmd.transform_selection(selection, tRMatrix4, state)
 	
 	# Now, translate the center of mass
	tTMatrix2 = generateTransformMatrixForTranslation(phenolCOMPosition)
	cmd.transform_selection(selection, tTMatrix2, state)
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def moveRelative(selection, phenolCOMPosition, state=1):
	# Function to move the chromophore to a particular position and orientation
	from pymol import cmd
	from numpy import array, cos, sin

	# Now, translate the center of mass
	tTMatrix = generateTransformMatrixForTranslation(phenolCOMPosition)
	cmd.transform_selection(selection, tTMatrix, state)
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def moveRelative2(selection, phenolCOMPosition):
	# Function to move the chromophore to a particular position and orientation
	from pymol import cmd

	cmd.translate(vector=phenolCOMPosition, selection=selection)
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def generateRandomCoordinates(steps, standardError):
	from numpy import array, float
	from random import random
	from math import sqrt
	# Generate random coordinates in a sphere
	randomCoordinates = []
	
	# Make the first coordinate be at zero error
	randomCoordinates.append([0.0, 0.0, 0.0])
	
	i = 0
	while i <  steps:
		if random() > 0.5:
			signX = +1.0
		else:
			signX = -1.0
		
		if random() > 0.5:
			signY = +1.0
		else:
			signY = -1.0
		
		if random() > 0.5:
			signZ = +1.0
		else:
			signZ = -1.0
		

		randX = random()*signX
		randY = random()*signY
		randZ = random()*signZ
		
		if sqrt(randX**2 + randY**2 + randZ**2) <= 1:
			randCoords = array([randX, randY, randZ], float)
			randCoords = randCoords*standardError
			randomCoordinates.append(randCoords)
			i = i + 1
	
	return randomCoordinates
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def generateNormalRandCoordinates(steps, standardError, firstCoordNoError=True, sigma=1.0):
	# sigma is the standard deviation at which to accept random numbers
	from RandomArray import normal
	from numpy import array, float
	from random import random
	from math import sqrt
	# Generate random coordinates in a sphere
	randomCoordinates = []
	
	# Make the first coordinate be at zero error
	if firstCoordNoError == True:
		randomCoordinates.append([0.0, 0.0, 0.0])
	
	i = 0
	while i <  steps:

		randX = normal(0.0,standardError)
		randY = normal(0.0,standardError)
		randZ = normal(0.0,standardError)
		
		if randX*randX + randY*randY + randZ*randZ <= \
		standardError*standardError*sigma*sigma:
			randCoords = array([randX, randY, randZ], float)
			randomCoordinates.append(randCoords)
			i = i + 1
	
	return randomCoordinates
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def generateNormalRandomAngles(steps, standardError, length, \
firstAngleNoError=True, sigma=1.0):
	from RandomArray import normal
	from numpy import arctan, pi, sqrt
	
	randomLengths = []
	
	# Make the first coordinate be at zero error
	if firstAngleNoError == True:
		randomLengths.append(0.0)
	
	i = 0
	while i <  steps:
		randLength = normal(0.0,standardError)
		if sqrt(randLength*randLength) <= standardError*sigma:
			randomLengths.append(randLength)
			i = i + 1
		
	i = 0
	randomAngles = []
	while i < len(randomLengths):
		randomAngle = arctan(randomLengths[i]/length)
		randomAngles.append(randomAngle)
		i = i + 1
	
	return randomAngles
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def randSign():
	from random import random
	
	rand = random()
	
	if rand >= 0.5:
		sign = +1.0
	else:
		sign = -1.0

	return sign

# ---------------------------------------------------------------------------- #








# ---------------------------------------------------------------------------- #
def MonteCarloMove5(identifier, croSelection, tyr203Selection, \
standardError, steps, state=1):
# The first energies returned are with no positional error
	from numpy import array, float
	from yaehmop2 import tightbind, energyGap, makeYaehmopInput
	from pymol import cmd
	
	# Generate random coordinates in a sphere
	randomCoordinatesCro = generateNormalRandCoordinates(steps, standardError)
	
	# Calculate the energies of the unperturbed state
	i = 0
	homoEnergies = []
	lumoEnergies = []
	selectionForCalculation = croSelection + ' or ' + tyr203Selection

	# Send out the state for an energy calculation 
	filename = makeYaehmopInput(identifier +  '_' + str(i), \
	selectionForCalculation, charge=-1)
	
	tightbind(filename)
	[HOMOEnergy, LUMOEnergy] = energyGap(filename)
	homoEnergies.append(HOMOEnergy)
	lumoEnergies.append(LUMOEnergy)	
	
	# Use the random coordinates to move the selection around
	
	i = 1
	temp_cro_Name = 'temp_cro'
	temp_tyr_Name = 'temp_tyr'
	
	while i < len(randomCoordinatesCro):
		randCroCoord = randomCoordinatesCro[i]
			
		# Generate a copy of the chromophore selection
		cmd.create(temp_cro_Name, croSelection)
		cmd.create(temp_tyr_Name, tyr203Selection)
		cmd.select(temp_cro_Name, temp_cro_Name)
		cmd.select(temp_tyr_Name, temp_tyr_Name)
		
		# Move the temporary selection
		moveRelative(temp_cro_Name, randCroCoord, state=state)
		
		selectionForCalculation = temp_cro_Name + ' or ' +\
		temp_tyr_Name
		
		# Send out the state for an energy calculation 
		filename = makeYaehmopInput(identifier +  '_' + str(i), \
		selectionForCalculation, charge=-1)
		
		tightbind(filename)
		[HOMOEnergy, LUMOEnergy] = energyGap(filename)
		homoEnergies.append(HOMOEnergy)
		lumoEnergies.append(LUMOEnergy)	
		
		# Delete the copy of the chromophore selection
		cmd.delete(temp_cro_Name)
		cmd.delete(temp_tyr_Name)
				
		i = i + 1
	
	return [homoEnergies, lumoEnergies]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #

def atomicNumber(atomSymbol):
# Returns the atomic number of an atom with a given symbol
	atomicNumber = 0.0
	if atomSymbol == 'C':
		atomicNumber = 6
	if atomSymbol == 'O':
		atomicNumber = 8
	if atomSymbol == 'N':
		atomicNumber = 7
	if atomSymbol == 'H':
		atomicNumber = 1
	
	return atomicNumber
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def centerOfElectronDensity(selection):
	from pymol import cmd
	model = cmd.get_model(selection)
	x, y , z = 0.0, 0.0 , 0.0
	electronCount = 0.0
	for a in model.atom:
		aNumber = atomicNumber(a.symbol)
		x+= a.coord[0]*aNumber
		y+= a.coord[1]*aNumber
		z+= a.coord[2]*aNumber
		electronCount+= aNumber
	
	return [x/electronCount, y/electronCount, z/electronCount]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def centerOfMass(selection):
	from pymol import cmd
	
	x, y , z = 0.0, 0.0 , 0.0
	totalMass = 0.0
	
	model = cmd.get_model(selection)
	
	for a in model.atom:
		mass = atomicMass(a.symbol)
		x+= a.coord[0]*mass
		y+= a.coord[1]*mass
		z+= a.coord[2]*mass
		totalMass += mass
	
	return [x/totalMass, y/totalMass, z/totalMass]
# ---------------------------------------------------------------------------- #
		


# ---------------------------------------------------------------------------- #
def atomicMass(atomSymbol):
# Returns the atomic number of an atom with a given symbol
# Will return atomic mass of 0.0 if symbol is not recognized. 
	atomicMass = 0.0
	if atomSymbol == 'C':
		atomicMass = 12.0107
	if atomSymbol == 'O':
		atomicMass = 15.9994
	if atomSymbol == 'N':
		atomicMass = 14.0067
	if atomSymbol == 'H':
		atomicMass = 1.00794
	if atomSymbol == 'S':
		atomicMass = 32.065

	
	return atomicMass
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def calculateAxesAndDistancesCro(croSelection):
	from pymol import cmd
	from numpy import cross
	from math import sqrt
	
	eCenterCro=centerOfElectronDensity(croSelection)
	selection = "OH"
	cmd.select(selection, "name OH and " + croSelection)
	OHcoords = cmd.get_model(selection).atom[0].coord
	e1Axis = vecLink(eCenterCro, OHcoords)
	CCroOH = sqrt(dotProduct(e1Axis, e1Axis))
	e1Axis = normalize(e1Axis)
	
	selection = "CB2"
	cmd.select(selection, "name CB2 and " + croSelection)
	CB2coords = cmd.get_model(selection).atom[0].coord
	e2Axis = vecLink(eCenterCro, CB2coords)
	CCroCB = sqrt(dotProduct(e2Axis, e2Axis))
	e2Axis = normalize(e2Axis)
	
	e3Axis = normalize(cross(e1Axis, e2Axis))

	return [eCenterCro, e1Axis, e2Axis, e3Axis, CCroOH, CCroCB]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def calculateAxesAndDistancesTyr203(tyr203Selection):
	from pymol import cmd
	from numpy import cross
	from math import sqrt
	
	eCenterTyr203 = centerOfElectronDensity(tyr203Selection)
	selection = "OH"
	cmd.select(selection, "name OH and " + tyr203Selection)
	OHcoords = cmd.get_model(selection).atom[0].coord
	t1Axis = vecLink(eCenterTyr203, OHcoords)
	C203OH = sqrt(dotProduct(t1Axis, t1Axis))
	t1Axis = normalize(t1Axis)
	
	selection = "CE2"
	cmd.select(selection, "name CE2 and " + tyr203Selection)
	CE2coords = cmd.get_model(selection).atom[0].coord
	t2Axis = vecLink(eCenterTyr203, CE2coords)
	C203CE2 = sqrt(dotProduct(t2Axis, t2Axis))
	t2Axis = normalize(t2Axis)
	
	t3Axis = normalize(cross(t1Axis, t2Axis))

	return [eCenterTyr203, t1Axis, t2Axis, t3Axis, C203OH, C203CE2]
# ---------------------------------------------------------------------------- #






# ---------------------------------------------------------------------------- #
def MonteCarloMoveTilt(identifier, croSelection, tyr203Selection, \
standardError, steps, sigma=3.0, state=1):
	
	# Calculate the energy levels variation of the Citrine chromophore
	# under randomly generated positional and tilt errors
	from numpy import array, float
	from yaehmop2 import tightbind, energyGap, makeYaehmopInput
	from pymol import cmd
	from buz_rotate_cmd import buz_rotate_cmd
	
	# Calculate eCenterCro, e1, e2, e3 axes and CCroOH and CCroCB distances for chromophore
	[eCenterCro, e1Axis, e2Axis, e3Axis, CCroOH, CCroCB] = \
	calculateAxesAndDistancesCro(croSelection)
	
	# Calculate eCenter203, t1, t2, t3 axes and C203OH and C203CE2 distances for tyrosine 203
	[eCenter203, t1Axis, t2Axis, t3Axis, C203OH, C203CE2] = \
	calculateAxesAndDistancesTyr203(tyr203Selection)

	# Generate gaussian distributed random coordinates for the center of mass
	# of the chromophore and tyrosine 203
	randomCoordinatesCro = generateNormalRandCoordinates(steps, standardError, \
	sigma=sigma)
	randomCoordinatesTyr203 = generateNormalRandCoordinates(steps, \
	standardError, sigma=sigma)
	
	# Generate random tilt angles for chromophore; theta1 and theta2 
	# and second theta1; theta1_2
	randomAnglesTheta1_1 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	randomAnglesTheta2   = generateNormalRandomAngles(steps, standardError, \
	CCroCB, firstAngleNoError=True, sigma=sigma)
	randomAnglesTheta1_2 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	# Generate random tilt angles for tyrosine 203; phi1 and phi2 
	# and second phi1; phi1_2
	randomAnglesPhi1_1   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)
	randomAnglesPhi2     = generateNormalRandomAngles(steps, standardError, \
	C203CE2,firstAngleNoError=True, sigma=sigma)
	randomAnglesPhi1_2   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)	
	
	# Monte Carlo loop
	# Remember, that the first coordinate and angle set has no error
	i = 0
	homoEnergies = []
	lumoEnergies = []
	temp_cro_Name = 'temp_cro'
	temp_tyr_Name = 'temp_tyr'
	
	while i < steps+1:
		# Generate a copy of the chromophore selection
		cmd.create(temp_cro_Name, croSelection)
		cmd.create(temp_tyr_Name, tyr203Selection)
		cmd.select(temp_cro_Name, temp_cro_Name)
		cmd.select(temp_tyr_Name, temp_tyr_Name)

		theta1_1 = randomAnglesTheta1_1[i]
		theta2 = randomAnglesTheta2[i]
		theta1_2 = randomAnglesTheta1_2[i]
		phi1_1 = randomAnglesPhi1_1[i]
		phi2 = randomAnglesPhi2[i]
		phi1_2 = randomAnglesPhi1_2[i]
		
		randCroCoord = randomCoordinatesCro[i]
		randTyrCoord = randomCoordinatesTyr203[i]
		
		# Rotate chromophore by theta1 about e3
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e3Axis, theta1_1)
		# Rotate chromophore by theta2 about e1
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e1Axis, theta2)
		# Rotate chromophore by theta1_2 about e2
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e2Axis, theta1_2)
		# Rotate tyrosine 203 by phi1 about t3
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t3Axis, phi1_1)
		# Rotate tyrosine 203 by phi2 about t1
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t1Axis, phi2)
		# Rotate tyrosine 203 by phi1_2 about t2
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t2Axis, phi1_2)
		
		# Move center of mass of tyrosine 203 to randomly generated coordinates
		moveRelative(temp_tyr_Name, randTyrCoord, state=state)
		
		# Move center of mass of chromophore to randomly generated coordinates
		moveRelative(temp_cro_Name, randCroCoord, state=state)
		
		# Send geometry to yaehmop2 for energy level calculation
		selectionForCalculation = temp_cro_Name + ' or ' + temp_tyr_Name
		filename = makeYaehmopInput(identifier +  '_' + str(i), \
		selectionForCalculation, charge=-1)
		tightbind(filename)
	
		# Read in energy levels and save
		[HOMOEnergy, LUMOEnergy] = energyGap(filename)
		homoEnergies.append(HOMOEnergy)
		lumoEnergies.append(LUMOEnergy)	
	
		
		# Delete the temporary chromophore and tyrosine 203 selection
		cmd.delete(temp_cro_Name)
		cmd.delete(temp_tyr_Name)
		
		# Repeat
		i = i + 1
		
	
	return [homoEnergies, lumoEnergies]
	
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def generateRandomCoordsAndAngles(standardError, CCroOH, CCroCB, C203OH, \
C203CE2, CD1CD2):
	from RandomArray import normal
	from numpy import array, float, arctan, pi, sqrt
	from random import random
	from math import sqrt
	

	randCroX = normal(0.0,standardError)
	randCroY = normal(0.0,standardError)
	randCroZ = normal(0.0,standardError)
	randTyrX = normal(0.0,standardError)
	randTyrY = normal(0.0,standardError)
	randTyrZ = normal(0.0,standardError)
		
	randCroCoords = [randCroX, randCroY, randCroZ]
	randTyrCoords = [randTyrX, randTyrY, randTyrZ]
	
	
	randLengthTheta1_1 = normal(0.0,standardError)
	randLengthTheta2_1 = normal(0.0,standardError)
	randLengthTheta1_2 = normal(0.0,standardError)
	randLengthPhi1_1 = normal(0.0,standardError)
	randLengthPhi2_1 = normal(0.0,standardError)
	randLengthPhi1_2 = normal(0.0,standardError)
	
	croPhenolLength = normal(0.0,standardError)
	
	theta1_1 = arctan(randLengthTheta1_1/CCroOH)
	theta2_1 = arctan(randLengthTheta2_1/CCroCB)
	theta1_2 = arctan(randLengthTheta1_2/CCroOH)
	phi1_1 = arctan(randLengthPhi1_1/C203OH)
	phi2_1 = arctan(randLengthPhi2_1/C203CE2)
	phi1_2 = arctan(randLengthPhi1_2/C203OH)
	
	croPhenolAngle = arctan(croPhenolLength/(0.5*CD1CD2))
	
	thetaAngles = [theta1_1, theta2_1, theta1_2]
	phiAngles = [phi1_1, phi2_1, phi1_2]
	
	return [randCroCoords, randTyrCoords, thetaAngles, phiAngles, \
	croPhenolAngle]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def perturbChromophoreAndTyr203(temp_cro_Name, temp_tyr_Name, randCroCoord, \
randTyrCoord, thetaAngles, phiAngles, eCenterCro, eCenter203, e1Axis, e2Axis,\
e3Axis, t1Axis, t2Axis, t3Axis, CG2Coord, nCG2CB2Vec, croPhenolAngle):
	from buz_rotate_cmd import buz_rotate_cmd
	from pymol import cmd
	
	# Rotate the chromophore phenol by croPhenolAngle
	name = selectChromophorePhenol(temp_cro_Name)
	buz_rotate_cmd(name, CG2Coord, nCG2CB2Vec, croPhenolAngle)

	[theta1_1, theta2, theta1_2] = thetaAngles
	[phi1_1, phi2, phi1_2] = phiAngles

	# Rotate chromophore by theta1 about e3
	buz_rotate_cmd(temp_cro_Name, eCenterCro, e3Axis, theta1_1)
	# Rotate chromophore by theta2 about e1
	buz_rotate_cmd(temp_cro_Name, eCenterCro, e1Axis, theta2)
	# Rotate chromophore by theta1_2 about e2
	buz_rotate_cmd(temp_cro_Name, eCenterCro, e2Axis, theta1_2)
	# Rotate tyrosine 203 by phi1 about t3
	buz_rotate_cmd(temp_tyr_Name, eCenter203, t3Axis, phi1_1)
	# Rotate tyrosine 203 by phi2 about t1
	buz_rotate_cmd(temp_tyr_Name, eCenter203, t1Axis, phi2)
	# Rotate tyrosine 203 by phi1_2 about t2
	buz_rotate_cmd(temp_tyr_Name, eCenter203, t2Axis, phi1_2)
	
	# Move center of mass of tyrosine 203 to randomly generated coordinates
	cmd.translate(selection=temp_tyr_Name, \
	vector=[randTyrCoord[0], randTyrCoord[1], randTyrCoord[2]])
		
	# Move center of mass of chromophore to randomly generated coordinates
	cmd.translate(selection=temp_cro_Name, \
	vector=[randCroCoord[0], randCroCoord[1], randCroCoord[2]])
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def insideErrorVolume(temp_cro_Name, temp_tyr_Name, croSelection, \
tyr203Selection, standardError, sigmaCoE, sigmaAtom):
	from pymol import cmd
	from math import sqrt
	
	structureOK = True
		
	# Check that the center of electron density of the perturbed chromophore
	# is no more than than sigmaCoE*standard error away from the recorded
	# position
	
	eCenterCroOrig = calculateAxesAndDistancesCro(croSelection)[0]
	eCenterTyr203Orig = calculateAxesAndDistancesTyr203(tyr203Selection)[0]
	eCenterCroPerturbed = calculateAxesAndDistancesCro(temp_cro_Name)[0]
	eCenterTyr203Perturbed = calculateAxesAndDistancesTyr203(temp_tyr_Name)[0]
	
	eCenterCroPerturbation = vecLength(vecLink(eCenterCroOrig, eCenterCroPerturbed))
	eCenterTyrPerturbation = vecLength(vecLink(eCenterTyr203Orig, eCenterTyr203Perturbed))
	
	if eCenterCroPerturbation > sigmaCoE*standardError or \
	eCenterTyrPerturbation > sigmaCoE*standardError:
		structureOK = False
	
	# Check that positions of atoms in perturbed chromophore are not more than
	# sigma numbers of standard errors away from the original structure

	perturbedCro = cmd.get_model(temp_cro_Name)
	perturbedTyr = cmd.get_model(temp_tyr_Name)
	originalCro = cmd.get_model(croSelection)
	originalTyr = cmd.get_model(tyr203Selection)
	
	i = 0
	while i < len(perturbedCro.atom) and structureOK == True:
		if perturbedCro.atom[i].symbol != 'H':
			length = vecLength(vecLink(perturbedCro.atom[i].coord, originalCro.atom[i].coord))
			if perturbedCro.atom[i].name != originalCro.atom[i].name:
				print("Error, atom names are not the same!")
			if length > standardError*sigmaAtom:
				structureOK = False
		i = i + 1
	
	# Check that positions of atoms in perturbed tyrosine are not more than
	# sigma numbers of standard errors away from the original structure
	i = 0
	while i < len(perturbedTyr.atom) and structureOK == True:
		if perturbedTyr.atom[i].symbol != 'H':
			length = vecLength(vecLink(perturbedTyr.atom[i].coord, originalTyr.atom[i].coord))
			if perturbedTyr.atom[i].name != perturbedTyr.atom[i].name:
				print("Error, atom names are not the same!")
			if length > standardError*sigmaAtom:
				structureOK = False
		i = i + 1

	return structureOK
# ---------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
def selectChromophorePhenol(chromophoreSelection):
	from pymol import cmd
	
	name = 'phenol_' + chromophoreSelection
	
	cmd.select(name, \
	'(name CG2 or name CD1 or name CD2 or name CE1 or name CE2 or name CZ ' + \
	'or name OH or name H01 or name H02 or name H03 or name H05) and ' + \
	chromophoreSelection)
	
	return name
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def calculateAxesAndDistancesCroPhenol(croSelection):
	# Calculate the positions of the CG2 atom, CG2 to CB2 vector and CD1 to CD2
	# distance
	
	from pymol import cmd
	from math import sqrt
	
	name = selectChromophorePhenol(croSelection)
	model = cmd.get_model(name)
	
	print(str(croSelection))
	
	for atom in model.atom:
		if atom.name == 'CG2':
			CG2Coord = atom.coord
			print(str(CG2Coord))
		if atom.name == 'CD1':
			CD1Coord = atom.coord
			print(str(CD1Coord))
		if atom.name == 'CD2':
			CD2Coord = atom.coord
			print(str(CD2Coord))
			
	# Select the CB2 atom
	cmd.select('CB2', 'name CB2 and ' + croSelection)
	model2 = cmd.get_model('CB2')
	CB2Coord = model2.atom[0].coord
	print(str(CB2Coord))
	
	nCG2CB2Vec = nVecLink(CG2Coord, CB2Coord)
	CD1CD2 = vecLink(CD1Coord, CD2Coord)
	CD1CD2 = sqrt(dotProduct(CD1CD2, CD1CD2))
	
	return [CG2Coord, nCG2CB2Vec, CD1CD2]
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def MonteCarloMoveTilt2(identifier, croSelection, tyr203Selection, \
standardError, steps, sigma=3.0, state=1):
# As well as moving and tilting the chromphore as a rigid body, this Monte Carlo
# algorithm also rotates the chromophore phenol and imidazolinone relative to
# one another about the bond that connects the CG2 atom to the CB2 atom 
	
	# Calculate the energy levels variation of the Citrine chromophore
	# under randomly generated positional and tilt errors
	from numpy import array, float
	from yaehmop2 import tightbind, energyGap, makeYaehmopInput
	from pymol import cmd
	from buz_rotate_cmd import buz_rotate_cmd
	
	# Calculate eCenterCro, e1, e2, e3 axes and CCroOH and CCroCB distances 
	# for chromophore
	[eCenterCro, e1Axis, e2Axis, e3Axis, CCroOH, CCroCB] = \
	calculateAxesAndDistancesCro(croSelection)
	
	# Calculate the positions of the CG2 atom, CG2 to CB2 vector and CD1 to CD2
	# distance
	[CG2Coord, nCG2CB2Vec, CD1CD2] = \
	calculateAxesAndDistancesCroPhenol(croSelection)
	
	# Calculate eCenter203, t1, t2, t3 axes and C203OH and C203CE2 distances 
	# for tyrosine 203
	[eCenter203, t1Axis, t2Axis, t3Axis, C203OH, C203CE2] = \
	calculateAxesAndDistancesTyr203(tyr203Selection)

	# Generate gaussian distributed random coordinates for the center of mass
	# of the chromophore and tyrosine 203
	randomCoordinatesCro = generateNormalRandCoordinates(steps, standardError, \
	sigma=sigma)
	randomCoordinatesTyr203 = generateNormalRandCoordinates(steps, \
	standardError, sigma=sigma)
	
	# Generate random tilt angles for chromophore; theta1 and theta2 
	# and second theta1; theta1_2
	randomAnglesTheta1_1 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesTheta2   = generateNormalRandomAngles(steps, standardError, \
	CCroCB, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesTheta1_2 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	
	# Generate random tilt angles for tyrosine 203; phi1 and phi2 
	# and second phi1; phi1_2
	
	randomAnglesPhi1_1   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesPhi2     = generateNormalRandomAngles(steps, standardError, \
	C203CE2,firstAngleNoError=True, sigma=sigma)
	
	randomAnglesPhi1_2   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)
	
	# Generate a random tilt angle for the chromophore phenol
	randomAnglesCroPhenol = generateNormalRandomAngles(steps, standardError, \
	CD1CD2/2.0, firstAngleNoError=True, sigma=sigma)
	
	# Monte Carlo loop
	# Remember, that the first coordinate and angle set has no error
	i = 0
	homoEnergies = []
	lumoEnergies = []
	temp_cro_Name = 'temp_cro'
	temp_tyr_Name = 'temp_tyr'
	
	while i < steps+1:
		# Generate a copy of the chromophore selection
		cmd.create(temp_cro_Name, croSelection)
		cmd.create(temp_tyr_Name, tyr203Selection)
		cmd.select(temp_cro_Name, temp_cro_Name)
		cmd.select(temp_tyr_Name, temp_tyr_Name)

		theta1_1 = randomAnglesTheta1_1[i]
		theta2 = randomAnglesTheta2[i]
		theta1_2 = randomAnglesTheta1_2[i]
		phi1_1 = randomAnglesPhi1_1[i]
		phi2 = randomAnglesPhi2[i]
		phi1_2 = randomAnglesPhi1_2[i]
		
		croPhenolAngle = randomAnglesCroPhenol[i]
		
		randCroCoord = randomCoordinatesCro[i]
		randTyrCoord = randomCoordinatesTyr203[i]
		
		# Rotate the chromophore phenol by croPhenolAngle
		name = selectChromophorePhenol(temp_cro_Name)
		buz_rotate_cmd(name, CG2Coord, nCG2CB2Vec, croPhenolAngle)
		
		
		# Rotate chromophore by theta1 about e3
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e3Axis, theta1_1)
		# Rotate chromophore by theta2 about e1
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e1Axis, theta2)
		# Rotate chromophore by theta1_2 about e2
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e2Axis, theta1_2)
		# Rotate tyrosine 203 by phi1 about t3
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t3Axis, phi1_1)
		# Rotate tyrosine 203 by phi2 about t1
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t1Axis, phi2)
		# Rotate tyrosine 203 by phi1_2 about t2
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t2Axis, phi1_2)
		
		
		# Move center of mass of tyrosine 203 to randomly generated coordinates
		moveRelative(temp_tyr_Name, randTyrCoord, state=state)
		
		# Move center of mass of chromophore to randomly generated coordinates
		moveRelative(temp_cro_Name, randCroCoord, state=state)
		
		# Send geometry to yaehmop2 for energy level calculation
		selectionForCalculation = temp_cro_Name + ' or ' + temp_tyr_Name
		filename = makeYaehmopInput(identifier +  '_' + str(i), \
		selectionForCalculation, charge=-1)
		tightbind(filename)
	
		# Read in energy levels and save
		[HOMOEnergy, LUMOEnergy] = energyGap(filename)
		homoEnergies.append(HOMOEnergy)
		lumoEnergies.append(LUMOEnergy)	
	
		
		# Delete the temporary chromophore and tyrosine 203 selection
		cmd.delete(temp_cro_Name)
		cmd.delete(temp_tyr_Name)
		
		# Repeat
		i = i + 1
		
	
	return [homoEnergies, lumoEnergies]
	
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def MonteCarloMoveTilt3(identifier, croSelection, tyr203Selection, \
standardError, steps, sigma=3.0, state=1):
# As well as moving and tilting the chromphore as a rigid body, this Monte Carlo
# algorithm also rotates the chromophore phenol and imidazolinone relative to
# one another about the bond that connects the CG2 atom to the CB2 atom 
	
	# Calculate the energy levels variation of the Citrine chromophore
	# under randomly generated positional and tilt errors
	from numpy import array, float
	from yaehmop2 import tightbind, energyGap, makeYaehmopInput
	from pymol import cmd
	from buz_rotate_cmd import buz_rotate_cmd
	
	# Calculate eCenterCro, e1, e2, e3 axes and CCroOH and CCroCB distances 
	# for chromophore
	[eCenterCro, e1Axis, e2Axis, e3Axis, CCroOH, CCroCB] = \
	calculateAxesAndDistancesCro(croSelection)
	
	# Calculate the positions of the CG2 atom, CG2 to CB2 vector and CD1 to CD2
	# distance
	[CG2Coord, nCG2CB2Vec, CD1CD2] = \
	calculateAxesAndDistancesCroPhenol(croSelection)
	
	# Calculate eCenter203, t1, t2, t3 axes and C203OH and C203CE2 distances 
	# for tyrosine 203
	[eCenter203, t1Axis, t2Axis, t3Axis, C203OH, C203CE2] = \
	calculateAxesAndDistancesTyr203(tyr203Selection)

	# Generate gaussian distributed random coordinates for the center of mass
	# of the chromophore and tyrosine 203
	randomCoordinatesCro = generateNormalRandCoordinates(steps, standardError, \
	sigma=sigma)
	randomCoordinatesTyr203 = generateNormalRandCoordinates(steps, \
	standardError, sigma=sigma)
	
	# Generate random tilt angles for chromophore; theta1 and theta2 
	# and second theta1; theta1_2
	randomAnglesTheta1_1 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesTheta2   = generateNormalRandomAngles(steps, standardError, \
	CCroCB, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesTheta1_2 = generateNormalRandomAngles(steps, standardError, \
	CCroOH, firstAngleNoError=True, sigma=sigma)
	
	# Generate random tilt angles for tyrosine 203; phi1 and phi2 
	# and second phi1; phi1_2
	
	randomAnglesPhi1_1   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)
	
	randomAnglesPhi2     = generateNormalRandomAngles(steps, standardError, \
	C203CE2,firstAngleNoError=True, sigma=sigma)
	
	randomAnglesPhi1_2   = generateNormalRandomAngles(steps, standardError, \
	C203OH, firstAngleNoError=True, sigma=sigma)
	
	# Generate a random tilt angle for the chromophore phenol
	randomAnglesCroPhenol = generateNormalRandomAngles(steps, standardError, \
	CD1CD2/2.0, firstAngleNoError=True, sigma=sigma)
	
	# Monte Carlo loop
	# Remember, that the first coordinate and angle set has no error
	i = 0
	homoEnergies = []
	lumoEnergies = []
	temp_cro_Name = 'temp_cro'
	temp_tyr_Name = 'temp_tyr'
	
	while i < steps+1:
		# Generate a copy of the chromophore selection
		cmd.create(temp_cro_Name, croSelection)
		cmd.create(temp_tyr_Name, tyr203Selection)
		cmd.select(temp_cro_Name, temp_cro_Name)
		cmd.select(temp_tyr_Name, temp_tyr_Name)

		theta1_1 = randomAnglesTheta1_1[i]
		theta2 = randomAnglesTheta2[i]
		theta1_2 = randomAnglesTheta1_2[i]
		phi1_1 = randomAnglesPhi1_1[i]
		phi2 = randomAnglesPhi2[i]
		phi1_2 = randomAnglesPhi1_2[i]
		
		croPhenolAngle = randomAnglesCroPhenol[i]
		
		randCroCoord = randomCoordinatesCro[i]
		randTyrCoord = randomCoordinatesTyr203[i]
		
		# Rotate the chromophore phenol by croPhenolAngle
		name = selectChromophorePhenol(temp_cro_Name)
		buz_rotate_cmd(name, CG2Coord, nCG2CB2Vec, croPhenolAngle)
		
		
		# Rotate chromophore by theta1 about e3
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e3Axis, theta1_1)
		# Rotate chromophore by theta2 about e1
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e1Axis, theta2)
		# Rotate chromophore by theta1_2 about e2
		buz_rotate_cmd(temp_cro_Name, eCenterCro, e2Axis, theta1_2)
		# Rotate tyrosine 203 by phi1 about t3
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t3Axis, phi1_1)
		# Rotate tyrosine 203 by phi2 about t1
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t1Axis, phi2)
		# Rotate tyrosine 203 by phi1_2 about t2
		buz_rotate_cmd(temp_tyr_Name, eCenter203, t2Axis, phi1_2)
		
		
		# Move center of mass of tyrosine 203 to randomly generated coordinates
		print(str(randTyrCoord))
		cmd.translate(vector=[randTyrCoord[0], randTyrCoord[1],randTyrCoord[2]],\
		selection=temp_tyr_Name)
		#moveRelative2(temp_tyr_Name, randTyrCoord)
		
		# Move center of mass of chromophore to randomly generated coordinates
		print(str(randCroCoord))
		cmd.translate(vector=[randCroCoord[0], randCroCoord[1],randCroCoord[2]],\
		selection=temp_cro_Name)
		#moveRelative2(temp_cro_Name, randCroCoord)
		
		# Send geometry to yaehmop2 for energy level calculation
		selectionForCalculation = temp_cro_Name + ' or ' + temp_tyr_Name
		filename = makeYaehmopInput(identifier +  '_' + str(i), \
		selectionForCalculation, charge=-1)
		tightbind(filename)
	
		# Read in energy levels and save
		[HOMOEnergy, LUMOEnergy] = energyGap(filename)
		homoEnergies.append(HOMOEnergy)
		lumoEnergies.append(LUMOEnergy)	
	
		
		# Delete the temporary chromophore and tyrosine 203 selection
		cmd.delete(temp_cro_Name)
		cmd.delete(temp_tyr_Name)
		
		# Repeat
		i = i + 1
		
	
	return [homoEnergies, lumoEnergies]
	
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def performEHTCalculation(identifier, selectionForCalculation):
	from yaehmop2 import tightbind, energyGap, makeYaehmopInput
	
	filename = makeYaehmopInput(identifier, selectionForCalculation, charge=-1)
	
	tightbind(filename)
		
	# Read in energy levels and save
	[HOMOEnergy, LUMOEnergy] = energyGap(filename)
	
	return [HOMOEnergy, LUMOEnergy]
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
def MonteCarloMoveVerify(identifier, croSelection, tyr203Selection, \
standardError, steps, sigmaCoE=1.0, sigmaAtom=1.0):

# Calculate the energy levels variation of the Citrine chromophore
# under randomly generated positional and tilt errors
# Verify that following tilting, deformation and translation of the chromophore
# that the chrompohore still lies within the error volume defined by
# sigma*standard error.

	from numpy import array, float
	from pymol import cmd
	
	# Calculate eCenterCro, e1, e2, e3 axes and CCroOH and CCroCB distances 
	# for chromophore
	[eCenterCro, e1Axis, e2Axis, e3Axis, CCroOH, CCroCB] = \
	calculateAxesAndDistancesCro(croSelection)
	
	# Calculate the positions of the CG2 atom, CG2 to CB2 vector and CD1 to CD2
	# distance
	[CG2Coord, nCG2CB2Vec, CD1CD2] = \
	calculateAxesAndDistancesCroPhenol(croSelection)
	
	# Calculate eCenter203, t1, t2, t3 axes and C203OH and C203CE2 distances 
	# for tyrosine 203
	[eCenter203, t1Axis, t2Axis, t3Axis, C203OH, C203CE2] = \
	calculateAxesAndDistancesTyr203(tyr203Selection)
	
	homoEnergies = []
	lumoEnergies = []
	temp_cro_Name = 'temp_cro'
	temp_tyr_Name = 'temp_tyr'
	
	# Perform calculation of energy of unperturbed state
	selection = croSelection + ' or ' + tyr203Selection
	id = identifier +  '_' + str(0)
	[HOMOEnergy, LUMOEnergy] = performEHTCalculation(id, selection)
	homoEnergies.append(HOMOEnergy)
	lumoEnergies.append(LUMOEnergy)
	
	# Monte Carlo loop
	
	i = 1
	while i <= steps:
		# Generate a copy of the chromophore selection
		cmd.create(temp_cro_Name, croSelection)
		cmd.create(temp_tyr_Name, tyr203Selection)

		# Generate random angles and random coordinates for chromophore and 
		# tyrosine 203
		[randCroCoord, randTyrCoord, thetaAngles, phiAngles, croPhenolAngle] = \
		generateRandomCoordsAndAngles(standardError, CCroOH, CCroCB, C203OH, \
		C203CE2, CD1CD2)
		
		# Perturb chromophore and tyrosine 203 by random angles
		perturbChromophoreAndTyr203(temp_cro_Name, temp_tyr_Name, randCroCoord, \
		randTyrCoord, thetaAngles, phiAngles, eCenterCro, eCenter203, e1Axis, \
		e2Axis, e3Axis, t1Axis, t2Axis, t3Axis, CG2Coord, nCG2CB2Vec, croPhenolAngle)
		
		# Check that perturbated structure still lies within error volume
		structureOK = \
		insideErrorVolume(temp_cro_Name, temp_tyr_Name, \
		croSelection, tyr203Selection, standardError, sigmaCoE=sigmaCoE, \
		sigmaAtom=sigmaAtom)
		
		print(str(structureOK))
		
		if structureOK == True:	
			# If it does, perform energy level calculation
			# Send geometry to yaehmop2 for energy level calculation
			selection = temp_cro_Name + ' or ' + temp_tyr_Name
			id = identifier +  '_' + str(i)
			
			[HOMOEnergy, LUMOEnergy] = performEHTCalculation(id, selection)
			
			homoEnergies.append(HOMOEnergy)
			lumoEnergies.append(LUMOEnergy)
			
			# Increment step count
			i = i + 1
		
		# Delete perturbed chromophore and tyrosine 203
		cmd.delete(temp_cro_Name)
		cmd.delete(temp_tyr_Name)
		
	return [homoEnergies, lumoEnergies]
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
def momentOfInertiaMatrix(coords, masses, center):
# Program to calculate the moment of inertia matrix of a group of atoms about
# a point
	from numpy import zeros, float, array
	
	tCoords = []
	
	# Transform all of the coordinates into the center of mass frame
	for coord in coords:
		tCoords.append(array(coord) - array(center))
	
	# Allocate the moment of inertia matrix
	inertiaMatrix = zeros([3,3], float)
	
	# Allocate the matrix elements
	j = 0 
	while j < len(coords):
		inertiaMatrix[0,0] += masses[j]*(tCoords[j][1]**2 + tCoords[j][2]**2)
		inertiaMatrix[1,1] += masses[j]*(tCoords[j][0]**2 + tCoords[j][2]**2)
		inertiaMatrix[2,2] += masses[j]*(tCoords[j][0]**2 + tCoords[j][1]**2)
		
		
		inertiaMatrix[0,1] += -1.0*masses[j]*tCoords[j][0]*tCoords[j][1]
		inertiaMatrix[0,2] += -1.0*masses[j]*tCoords[j][0]*tCoords[j][2]
		inertiaMatrix[1,2] += -1.0*masses[j]*tCoords[j][1]*tCoords[j][2]
		
		inertiaMatrix[1,0] = inertiaMatrix[0,1]
		inertiaMatrix[2,0] = inertiaMatrix[0,2]
		inertiaMatrix[2,1] = inertiaMatrix[1,2]
		
		j += 1
		
	return inertiaMatrix
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
def momentOfInertiaMatrixPymolSelection(selection):
# Returns the moment of inertia matrix of the selection about its center of mass
	from pymol import cmd
	
	# Calculate the center of mass of the object
	center = centerOfMass(selection)

	# Get the coordinates of the selection from pymol
	model = cmd.get_model(selection)
	
	masses = []
	coords = []
	
	for atom in model.atom:
		coords.append(atom.coord)
		masses.append(atomicMass(atom.symbol))
	
	
	inertiaMatrix = momentOfInertiaMatrix(coords, masses, center)

	return inertiaMatrix
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def VectorToOutputString(vector):
	i = 0
	outputString = ''
	while i < len(vector):
		outputString += str(vector[i]) + '\t'
		i += 1
	
	return outputString
# ---------------------------------------------------------------------------- #
