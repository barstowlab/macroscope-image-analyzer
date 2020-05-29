from tkinter import filedialog

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

source = filedialog.askdirectory()
dest = filedialog.askdirectory()
fileList = GenerateFileList(directory=source, regex=".*\.jpg")

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

def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

def ProcessPlateBarcode(fileName):
    import pyzbar.pyzbar as pyzbar
    import numpy as np
    import cv2

    timeStamp = ProcessTimeStamp(fileName)
    # print(fileName)
    processed_result = ProcessFileNameExtension(fileName)
    fileNameNoExt = processed_result[0]
    fileExtension = processed_result[1]
    image = cv2.imread(fileName)
    # print(image)
    decodedObjects = pyzbar.decode(image)
	
    if decodedObjects == []:
        print("No Barcode found in: " + fileName )
        plateID = 'UNKNOWN'
    elif decodedObjects != None:
        # print(decodedObjects)
        barcode = decodedObjects[0][0][:-1].decode("utf-8")
        if len(barcode) > 1:
            # print(barcode)
            plateID = barcode[:-1]
            checkSum = barcode[-1]
            # print(plateID)
            calculatedChecksum = GenerateCheckSum(plateID)
		
            # if checkSum != calculatedChecksum:
                # print(checkSum, calculatedChecksum)
                # print("Error in barcode check sum. File: " + fileName)
        else:
            plateID = 'UNKNOWN'
    else:
        plateID = 'UNKNOWN'

    return [plateID, timeStamp, fileExtension]

def RenameAndCopyImageFile(destBaseDir, fileName, diagnostics=False, preflight=False, \
organizeIntoDirectories=True):
    import shutil
    import pdb
    import os
	
    [plateID, timeStamp, fileExtension] = ProcessPlateBarcode(fileName)

    print(plateID)
    if plateID == "AQDS25":
        try:
            if fileExtension[0] != '.':
                newFileName = plateID + '-' + timeStamp + '.' + fileExtension
            else:
                newFileName = plateID + '-' + timeStamp + fileExtension
        except:
            pdb.set_trace()
        
        
        # if diagnostics == True:
            # print(fileName + ' --> ' + newFileName)
        
        if preflight == False:
            # if organizeIntoDirectories == True:
            #     ensure_dir(destBaseDir + '/' + plateID + '/' + newFileName)
            #     shutil.copy(fileName, destBaseDir + '/' + plateID + '/' + newFileName)
            # else:
            #     ensure_dir(destBaseDir + '/' + newFileName)
            print('fileName', fileName)
            print('destBaseDir', destBaseDir)
            shutil.copy(fileName, destBaseDir + '/' + os.path.basename(fileName))

for file in fileList:
    # print(file)
    RenameAndCopyImageFile(dest, source + '/' + file, diagnostics=True, preflight=False, organizeIntoDirectories=True)