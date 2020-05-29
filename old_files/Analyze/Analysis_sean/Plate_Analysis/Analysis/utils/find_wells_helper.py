# helper functions for FindWellPositionsImproved


# takes a csv file name containing information about the wells as an input
def parse_well_info(plate_info_file):

    import csv

    wells_info = {}
    with open(plate_info_file, newline='') as csvfile:
        well_reader = csv.reader(csvfile, delimiter=',')
        headers = next(well_reader)
        for well in well_reader:
            well_info = {}
            for i in range(1, len(headers)):
                well_info[headers[i].replace(" ", "")] = well[i]
            wells_info[well[0]] = well_info
    return wells_info


# generates xml file for each plate
def find_wells(plate_names, img_types, folder_location, well_info, new_pixel_width, save_location, scaled_img_postfix,
               well_coords_postfix):
    fittedWellPositionDict = {}
    well_colors_dict = {}
    print(plate_names)
    for plate in plate_names:

        output_img_file_name = save_location + '/' + plate + scaled_img_postfix
        fittedWellPositionDict[plate], well_radius = locate_wells(plate, img_types, folder_location, new_pixel_width,
                                                                  output_img_file_name)
        well_colors_dict[plate] = get_well_colors(plate, img_types, folder_location, fittedWellPositionDict[plate], int(well_radius / 2))
        well_coords_file_name = save_location + '/' + plate + well_coords_postfix
        WriteWellPositionsDictForImageAsXML(well_coords_file_name, plate, fittedWellPositionDict[plate], folder_location + '/' + plate,
                                            well_info, well_radius, well_colors_dict[plate])
        WriteWellPositionsDictForImageAsCSV(well_coords_file_name, plate, well_info[plate], well_colors_dict[plate])
    return


# creates GUI for finding well locations
def locate_wells(plate, img_types, folder_location, new_pixel_width, output_file_name):
    import tkinter as Tk
    from PIL import Image, ImageTk
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import imshow, annotate, savefig, close

    # gets all images
    workingDir = folder_location + '/' + plate
    fileList = GenerateFileList(directory=workingDir, regex=".*" + img_types[plate])
    fileList = sorted(fileList)

    root = Tk.Tk()
    well_locs = {}  # dictionary of pixel locations for each well

    # imports image and resizes it
    image = Image.open(workingDir + '/' + fileList[0])
    width,height = image.size
    control_panel_width = 400
    wall_offset = 100
    scaling = new_pixel_width / width
    new_height = int(scaling * height)
    image = image.resize((new_pixel_width, new_height))

    mode = 0  # mode 0 means no locations defined, 1 means A1 defined, 2 means both are defined
    first_well_x = 0
    first_well_y = 0
    last_well_x = 0
    last_well_y = 0
    crs = 4  # radius in pixels for circles drawn on the plate
    row_num = 8
    col_num = 12
    radius = 1
    spacing = 0
    shapes = []  # keeps track of all the shapes drawn on canvas

    canvas = Tk.Canvas(root, width=new_pixel_width + control_panel_width, height=new_height)

    # Undoes last well location defined
    def undo_well_click():
        nonlocal mode
        if mode == 1:  # gets rid of first well definition
            canvas.delete(shapes.pop())
            mode -= 1
        elif mode == 2:  # gets rid of second well definition
            while len(shapes) > 1:
                canvas.delete(shapes.pop())
            mode -= 1

    # If well locations have been defined, plots where every well is
    def test_wells_on_click():
        nonlocal row_num
        nonlocal col_num
        nonlocal radius
        nonlocal spacing
        if mode == 2:
            row_num = int(row_entry.get())
            col_num = int(column_entry.get())
            w0 = abs(first_well_x - last_well_x)
            h0 = abs(first_well_y - last_well_y)
            spacing = h0 * col_num / (row_num - 1) - w0
            radius = (w0 - spacing * (col_num - 1)) / 2 / col_num
            x_dir = 1
            y_dir = 1
            if first_well_x > last_well_x:
                x_dir = -1
            if first_well_y > last_well_y:
                y_dir = -1
            for r in range(1, row_num + 1):
                for c in range(1, col_num + 1):
                    dr = radius / 2
                    x = first_well_x + x_dir * ((c-1)*spacing + (2*c - 1)*radius)
                    y = first_well_y + y_dir * ((r-1)*spacing + (2*r - 2)*radius)
                    shapes.append(canvas.create_oval(x-dr, y-dr, x+dr, y+dr))
                    well_locs[get_well_id(r,c)] = np.array([x,y]) / scaling

    # saves the scaled and annotated image, exits window and moves on in the function
    def approve_wells_on_click():
        for file in fileList[0:1]: # for now will only do this for one
            new_image = Image.open(workingDir + '/' + file)
            new_image = new_image.resize((new_pixel_width, new_height))
            if mode == 2:
                imshow(new_image)
                for r in range(1, row_num + 1):
                    for c in range(1, col_num + 1):
                        annotate(s=get_well_id(r,c), xy=well_locs[get_well_id(r,c)] * scaling)
                print(output_file_name, file)
                save_spot = output_file_name[:len(output_file_name) - 4] + file[:len(file)-4] + output_file_name[len(output_file_name)-4:]
                savefig(save_spot)
                close()
        root.destroy()

    # defines known locations on the plate
    def click_plate(event):
        nonlocal mode
        nonlocal first_well_x
        nonlocal first_well_y
        nonlocal last_well_x
        nonlocal last_well_y
        if mode == 0:
            shapes.append(canvas.create_oval(event.x-crs,event.y-crs,event.x+crs,event.y+crs,fill='blue'))
            first_well_x = event.x
            first_well_y = event.y
            mode += 1

        elif mode == 1:
            shapes.append(canvas.create_oval(event.x-crs,event.y-crs,event.x+crs,event.y+crs, fill='green'))
            last_well_x = event.x
            last_well_y = event.y
            mode += 1
        elif mode == 2 and len(shapes) > 2:
            # gets row and column number of well pressed
            column = int(np.floor(abs(event.x - first_well_x) / (2*radius + spacing))) + 1
            row = int(np.floor(abs(event.y - first_well_y) / (2*radius + spacing))) + 1
            if column <= col_num and row <= row_num and column > 0 and row > 0:
                loc = well_locs[get_well_id(row, column)] * scaling
                red_pixels, green_pixels, blue_pixels = get_pixels(image.load(), int(loc[0]), int(loc[1]), int(radius/2))
                plt.figure()
                plt.hist(red_pixels)
                plt.hist(green_pixels)
                plt.hist(blue_pixels)
                plt.show()

    plate = ImageTk.PhotoImage(image, master=canvas)
    canvas.create_image(new_pixel_width / 2, new_height / 2, image=plate)
    canvas.bind("<Button-1>", click_plate)
    canvas.pack()

    leftx = (new_pixel_width + wall_offset) / (new_pixel_width + control_panel_width)

    Tk.Label(canvas, text="Enter Row #").place(relx=leftx,rely=0.1)
    row_entry = Tk.Entry(canvas)
    row_entry.place(relx=leftx,rely=0.2)
    row_entry.insert(0,'8')
    Tk.Label(canvas, text="Enter Col #").place(relx=leftx,rely=0.3)
    column_entry = Tk.Entry(canvas)
    column_entry.place(relx=leftx,rely=0.4)
    column_entry.insert(0,'12')
    undo_wells = Tk.Button(canvas, text="Undo Well Define", command=undo_well_click)
    undo_wells.place(relx=leftx, rely=0.5)
    test_wells = Tk.Button(canvas, text="Test Well Locations", command=test_wells_on_click)
    test_wells.place(relx=leftx,rely=0.6)
    approve_wells = Tk.Button(canvas, text="Approve Well Locations", command=approve_wells_on_click)
    approve_wells.place(relx=leftx,rely=0.7)

    Tk.mainloop()

    return well_locs, radius / scaling


# gets all pixels within given radius and given location
def get_pixels(pixels, x, y, radius):
    import numpy as np
    red_pixels = []
    green_pixels = []
    blue_pixels = []
    # print("i", y - radius, y + radius + 1)
    for i in range(y - radius, y + radius + 1):
        delta_x = max(0, int(np.sqrt(radius ** 2 - (y - i) ** 2)) - 1)
        for j in range(x - delta_x, x + delta_x + 1):
            red_pixels.append(pixels[j, i][0])
            green_pixels.append(pixels[j, i][1])
            blue_pixels.append(pixels[j, i][2])
    return red_pixels, green_pixels, blue_pixels


# get colors from well
def get_well_colors(plate, img_types, folder_location, fittedWellPositionDict, radius):
    from PIL import Image
    import numpy as np

    # defines yellow and white
    yellow = np.array([225, 153, 0])
    white = np.array([255, 255, 255])

    # gets all images
    workingDir = folder_location + '/' + plate
    fileList = GenerateFileList(directory=workingDir, regex=".*" + img_types[plate])
    fileList = sorted(fileList)
    well_colors_dict = {}
    for well in fittedWellPositionDict:
        well_colors_dict[well] = {'mean_red': np.array([]), 'median_red': np.array([]),
                                  'mean_green': np.array([]), 'median_green': np.array([]),
                                  'mean_blue': np.array([]), 'median_blue': np.array([]),
                                  'mean_yellow': np.array([]), 'median_yellow': np.array([]), 'times': np.array([])}
    for file in fileList:
        img = Image.open(workingDir + '/' + file)
        time = ProcessTimeStamp(workingDir + '/' + file)
        pixels = img.load()
        for well in fittedWellPositionDict:
            # gets all pixels in well within a certain radius
            x = int(fittedWellPositionDict[well][0])
            y = int(fittedWellPositionDict[well][1])
            red_pixels, green_pixels, blue_pixels = get_pixels(pixels, x, y, -radius)

            # finds the mean and median pixel values for these colors
            # print(red_pixels)
            r_mean = np.mean(red_pixels)
            r_median = np.median(red_pixels)
            g_mean = np.mean(green_pixels)
            g_median = np.median(green_pixels)
            b_mean = np.mean(blue_pixels)
            b_median = np.median(blue_pixels)
            mean_color = np.array([r_mean, g_mean, b_mean]) - white
            median_color = np.array([r_median, g_median, b_median]) - white
            mean_yellow = mean_color.dot(yellow.T-white.T) / np.linalg.norm(yellow - white)
            median_yellow = median_color.dot(yellow.T-white.T) / np.linalg.norm(yellow - white)
            well_colors_dict[well]['mean_red'] = np.append(well_colors_dict[well]['mean_red'], r_mean)
            well_colors_dict[well]['mean_green'] = np.append(well_colors_dict[well]['mean_green'], g_mean)
            well_colors_dict[well]['mean_blue'] = np.append(well_colors_dict[well]['mean_blue'], b_mean)
            well_colors_dict[well]['mean_yellow'] = np.append(well_colors_dict[well]['mean_yellow'], mean_yellow)
            well_colors_dict[well]['median_red'] = np.append(well_colors_dict[well]['median_red'], r_median)
            well_colors_dict[well]['median_green'] = np.append(well_colors_dict[well]['median_green'], g_median)
            well_colors_dict[well]['median_blue'] = np.append(well_colors_dict[well]['median_blue'], b_median)
            well_colors_dict[well]['median_yellow'] = np.append(well_colors_dict[well]['median_yellow'], median_yellow)

            well_colors_dict[well]['times'] = np.append(well_colors_dict[well]['times'], time)

    return well_colors_dict


# translates row and column numbers to well id
def get_well_id(r, c):
    return chr(ord('A') + r - 1) + str(c).zfill(2)


# writes out well information in xml form
def WriteWellPositionsDictForImageAsXML(fileName, plate, fittedWellPositionsDict, originatingFolder, well_info, well_radius,
                                        well_colors_dict):

    import xml.etree.ElementTree as ET

    plateImageTreeRoot = ET.Element('plateImage')
    plateImageTreeRoot.set('originatingFolder', originatingFolder)
    plateImageTreeRoot.set('wellRadius', str(well_radius))

    for well in fittedWellPositionsDict:
        wellSubElement = ET.SubElement(plateImageTreeRoot, 'well')
        wellSubElement.set('wellID', well)
        wellSubElement.set('wellPosition', ExportListForXML(fittedWellPositionsDict[well]))
        if well in well_info[plate]:
            for key in well_info[plate][well]:
                wellSubElement.set(key, well_info[plate][well][key])
        if well in well_colors_dict:
            for key in well_colors_dict[well]:
                wellSubElement.set(key, ExportListForXML(well_colors_dict[well][key]))

    indent(plateImageTreeRoot)
    plateImageTree = ET.ElementTree(plateImageTreeRoot)
    plateImageTree.write(fileName)

    return


# writes out well information in csv form
def WriteWellPositionsDictForImageAsCSV(fileName, plate, well_info, well_colors_dict):
    import numpy as np
    import csv

    # first converts well tracking of colors to gene tracking of colors
    relevant_colors = ['median_yellow']
    gene_colors_dict = {}
    times = np.array([])
    for well in well_info:
        gene = well_info[well]['GeneName']
        times = well_colors_dict[well]['times']
        if gene not in gene_colors_dict:
            gene_colors_dict[gene] = {}
        for color in relevant_colors:
            new_arr = np.reshape(well_colors_dict[well][color], (1, len(well_colors_dict[well][color])))
            if color in gene_colors_dict[gene]:
                gene_colors_dict[gene][color] = np.vstack((gene_colors_dict[gene][color], new_arr))
            else:
                gene_colors_dict[gene][color] = new_arr

    # now converts to the csv and takes average and standard deviation of values
    csv_file_name = fileName[:len(fileName)-3] + 'csv'
    with open(csv_file_name, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        csvwriter.writerow(['GeneName'] + ExportListForCSV(times))
        for gene in gene_colors_dict:
            for color in gene_colors_dict[gene]:
                mean_colors = np.mean(gene_colors_dict[gene][color], axis=0)
                std_colors = np.std(gene_colors_dict[gene][color], axis=0)
                row = [gene + ' ' + color] + ExportListForCSV(mean_colors)
                csvwriter.writerow(row)
                new_row = [gene + ' ' + color + ' std'] + ExportListForCSV(std_colors)
                csvwriter.writerow(new_row)


# method for indenting some text when writing to xml (from Buz's old code)
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


# translates array for csv
def ExportListForCSV(original_list):
    new_arr = []
    for val in original_list:
        new_arr.append(str(val))
    return new_arr


# translates array for xml (also from Buz)
def ExportListForXML(listToExport, delimeter=','):

    outputStr = ''

    i = 0
    while i < len(listToExport):
        outputStr += str(listToExport[i])
        if i < len(listToExport) - 1:
            outputStr += delimeter
        i += 1

    return outputStr


# gets absolute time image was taken (it's in the fileName)
def ProcessTimeStamp(fileName):
    import re
    import os

    baseName = os.path.basename(fileName)

    timeStampRegex = re.compile("(\w*)\-(\d+)")
    timeStampSearch = timeStampRegex.search(baseName)

    if timeStampSearch != None:
        timeStamp = timeStampSearch.group(2)
    else:
        ex = 0
        raise 0

    return timeStamp
# -----------------------


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


# reads xml file created in FindWellPositionsImproved.py
def import_gene_color_change_dict(color_info_file):
    import xml.etree.ElementTree as ET
    import numpy as np

    tree = ET.parse(color_info_file)
    root = tree.getroot()
    gene_series = {}

    wells = root.findall('well')
    color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'median_yellow',
                  'mean_yellow']

    for well in wells:
        gene = well.attrib['GeneName']
        times = np.array(well.attrib['times'].split(','),float)
        times = times - times[0]
        ordered_indices = np.argsort(times)
        if gene in gene_series:
            for color in color_keys:
                color_series = np.array(well.attrib[color].split(','), float)
                gene_series[gene][color] = np.vstack((gene_series[gene][color],
                                                     np.reshape(color_series[ordered_indices], (1, len(times)))))
        else:
            gene_series[gene] = {}
            gene_series[gene]["times"] = times
            for color in color_keys:
                color_series = np.array(well.attrib[color].split(','), float)
                gene_series[gene][color] = np.reshape(color_series[ordered_indices], (1, len(times)))

    return gene_series


# reorients image to have the width be the longest dimension and maximizes similarities in red to old image
# to make sure they are all facing the same direction
def orient_images(plate_names, img_types, folder_location, new_folder_location):
    from PIL import Image, ImageTk
    import os

    for plate in plate_names:
        # gets all images
        workingDir = folder_location + plate
        fileList = GenerateFileList(directory=workingDir, regex=".*" + img_types[plate])
        fileList = sorted(fileList)

        # makes sure image has a larger width than height
        image = Image.open(workingDir + '/' + fileList[0])
        width, height = image.size
        if height > width:
            image = image.rotate(90)
        new_width = 200
        new_height = int(new_width / width * height)
        image_small = image.resize((new_width, new_height))

        # creates directory if neccessary
        if not os.path.exists(new_folder_location + plate):
            os.makedirs(new_folder_location + plate)

        image.save(new_folder_location + plate + '/' + fileList[0])

        # uses the first image as a baseline to compare to all the other images
        for file in fileList[1:]:
            new_image = Image.open(workingDir + '/' + file)
            width, height = new_image.size
            if height > width:
                im1 = new_image.rotate(90)
                temp = width
                width = height
                height = temp
            else:
                im1 = new_image
            im2 = im1.rotate(180)
            im1_small = im1.resize((new_width,new_height))
            im2_small = im2.resize((new_width, new_height))

            cor1 = calc_red_cor(image_small, im1_small, new_width, new_height)
            cor2 = calc_red_cor(image_small, im2_small, new_width, new_height)
            if cor1 > cor2:
                im1.save(new_folder_location + plate + '/' + file)
            else:
                im2.save(new_folder_location + plate + '/' + file)


# calculates cross correlation of im1 and im2 in the red color channel
def calc_red_cor(im1, im2, width, height):
    import math
    pix1 = im1.load()
    pix2 = im2.load()
    pix1_norm = 0
    pix2_norm = 0
    red_cor = 0
    for i in range(width):
        for j in range(height):
            red_cor += pix1[i,j][0] * pix2[i,j][0]
            pix1_norm += pix1[i,j][0]**2
            pix2_norm += pix2[i,j][0]**2
    pix1_norm = math.sqrt(pix1_norm / (width * height))
    pix2_norm = math.sqrt(pix2_norm / (width * height))
    return red_cor / (width * height) / pix1_norm / pix2_norm
