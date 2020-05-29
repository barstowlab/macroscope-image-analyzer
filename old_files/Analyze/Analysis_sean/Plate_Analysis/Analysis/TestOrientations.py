#DO NOT USE

import tkinter as Tk
from PIL import Image, ImageTk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imshow, annotate, savefig, close

# generates list of files from the given "directory" that meet the given regular expression "regex"
def GenerateFileList(directory=".", regex=".*\.ProcSpec.dat", ignoreCase=True):
    import os
    import re
    fileList = os.listdir(directory)

    if ignoreCase == True:
        filePattern = re.compile(regex, re.IGNORECASE)
    else:
        filePattern = re.compile(regex)

    i = 0
    selectFiles = []

    while i < len(fileList):
        if filePattern.match(fileList[i]) != None:
            selectFiles.append(fileList[i])
        i = i + 1

    return selectFiles

workingDir = '../../../../../Photographic Data/Farshid/ElectronUptake'
fileList = GenerateFileList(directory=workingDir, regex=".*" + '.jpg')
fileList = sorted(fileList)

scale = 10
img = Image.open(workingDir + '/' + fileList[0])
width0,height0 = img.size
img = img.resize((width0//scale, height0//scale))
width0,height0 = img.size

off = 0

# for file in fileList:
#     img = Image.open(workingDir + '/' + file)
#     width,height = img.size
#     if width != width0 or height != height0:
#         off += 1
#         print(file)

img2 = Image.open(workingDir + '/IMG-1573273023.jpg')
img2 = img2.rotate(90, expand=True)
#img2.save(workingDir + '/IMG-1573273023_1.jpg', quality=95)
width,height = img2.size
img2 = img2.resize((width//scale, height//scale))

width,height = img2.size

root = Tk.Tk()
canvas = Tk.Canvas(root, width=width0, height=height0)
plate = ImageTk.PhotoImage(img, master=canvas)
canvas.create_image(width0 / 2, height0 / 2, image=plate)
canvas.pack()
Tk.mainloop()

root = Tk.Tk()
canvas = Tk.Canvas(root, width=width, height=height)
plate = ImageTk.PhotoImage(img2, master=canvas)
canvas.create_image(width / 2, height / 2, image=plate)
canvas.pack()
Tk.mainloop()
