
def ProcessPlateBarcode(fileName):
    import pyzbar.pyzbar as pyzbar
    import numpy as numpy
    import cv2

    timeStamp = ProcessTimeStamp(fileName)
	[fileNameNoExt, fileExtension] = ProcessFileNameExtension(fileName)

    decodedObjects = pyzbar.decode(image)
    if decodedObjects != None:
        barcode = decodedObjects[0][0][:-1].decode("utf-8")
        if len(barcode > 1):
            plateID = barcode[:-1]
            checkSum = barcode[-1]
            calculatedChecksum = GenerateCheckSum(plateID)
        
            if checkSum != calculatedChecksum:
                print("Error in barcode check sum. File: " + fileName)\
        else:
            plateID = 'UNKNOWN'
    else:
        plateID = 'UNKNOWN'


if __name__ == '__main__':
    ProcessPlateBarcode('test.jpg')