import numpy as np
from tkinter import filedialog
import tkinter as tk

# ------------------------------------------------------------------------------------------------ #

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
       
    return checkcar
# ------------------------------------------------------------------------------------------------ #

def generate_barcode(plate_type, output_file, prefix, postfix, MaxRows, MaxColumns, PaddingColumns, \
    PaddingRows, sampleNumbersStart, sampleNumbersEnd):

    delimeter = ","
    lineBreak = "\n"

    TotalLabels = MaxRows * MaxColumns

    sampleNumbers = np.arange(sampleNumbersStart, sampleNumbersEnd, 1)

    # Generate the sample labels
    sampleIDs = []
    for number in sampleNumbers:
        sampleID = prefix + str(number).zfill(3) + 'T'
        sampleIDs.append(sampleID)

    labels = []
    for sampleID in sampleIDs:
        barcode = sampleID + GenerateCheckSum(sampleID)
        labels.append("*" + barcode + "*")

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

    fileHandle = open(output_file, 'w')
    for row in rowData:
        fileHandle.write(row)
        
    fileHandle.close()
    print("done")
    return

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #

def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

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

# ------------------------------------------------------------------------------------------------ #
def ProcessTimeStamp(fileName):
    import re
    import os
	
    baseName = os.path.basename(fileName)

    timeStampRegex = re.compile("(\w*)\-(\d+)")
    timeStampSearch = timeStampRegex.search(baseName)

    if timeStampSearch != None:
        timeStamp = timeStampSearch.group(2)
    else:
        raise Exception('Time Stamp Failure')
		
    return timeStamp
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ProcessFileNameExtension(fileName):
    import re
    import os
    import pdb
	
    baseName = os.path.basename(fileName)
    [fileName, fileExtension] = os.path.splitext(baseName)
	
# 	pdb.set_trace()
		
    return [fileName, fileExtension]
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
def ProcessPlateBarcode(fileName):
    import pyzbar.pyzbar as pyzbar
    import numpy as np
    import cv2

    timeStamp = ProcessTimeStamp(fileName)
    processed_result = ProcessFileNameExtension(fileName)
    fileNameNoExt = processed_result[0]
    fileExtension = processed_result[1]
    image = cv2.imread(fileName)
    decodedObjects = pyzbar.decode(image)
	
    if decodedObjects != None:
        barcode = decodedObjects[0][0][:-1].decode("utf-8")
        if len(barcode) > 1:
            print(barcode)
            plateID = barcode[:-1]
            checkSum = barcode[-1]
            print(plateID)
            calculatedChecksum = GenerateCheckSum(plateID)
		
            if checkSum != calculatedChecksum:
                print(checkSum, calculatedChecksum)
                print("Error in barcode check sum. File: " + fileName)
        else:
            plateID = 'UNKNOWN'
    else:
        plateID = 'UNKNOWN'

    return [plateID, timeStamp, fileExtension]



# ------------------------------------------------------------------------------------------------ #

def RenameAndCopyImageFile(destBaseDir, fileName, diagnostics=False, preflight=False, \
organizeIntoDirectories=True):
    import shutil
    import pdb
	
    [plateID, timeStamp, fileExtension] = ProcessPlateBarcode(fileName)
	
    try:
        if fileExtension[0] != '.':
            newFileName = plateID + '-' + timeStamp + '.' + fileExtension
        else:
            newFileName = plateID + '-' + timeStamp + fileExtension
    except:
        pdb.set_trace()
	
	
    if diagnostics == True:
        print(fileName + ' --> ' + newFileName)
	
    if preflight == False:
        if organizeIntoDirectories == True:
            ensure_dir(destBaseDir + '/' + plateID + '/' + newFileName)
            shutil.copy(fileName, destBaseDir + '/' + plateID + '/' + newFileName)
        else:
            ensure_dir(destBaseDir + '/' + newFileName)
        shutil.copy(fileName, destBaseDir + '/' + newFileName)
	
# ------------------------------------------------------------------------------------------------ #


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


def organize_by_barcode(plate_type, source, dest):
    fileList = GenerateFileList(directory=source, regex=".*\.jpg")
    print(source)
    print(dest)
    print(plate_type)


    if plate_type == "Assay":
        for file in fileList:
            print(file)
            RenameAndCopyImageFile(dest, source + '/' + file, \
            diagnostics=True, preflight=False, organizeIntoDirectories=True)
    elif plate_type == "Storage":
        for file in fileList:
            RenameAndCopyImageFile(dest, source + '/' + file, \
            diagnostics=True, preflight=False, organizeIntoDirectories=False)
    else:
        raise Exception("Bad plate type")

    print("done")
    return