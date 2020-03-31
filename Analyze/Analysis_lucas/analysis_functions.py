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