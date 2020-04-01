# Author: Sean Medin, Zhengxing Xue
# This code is used to plot the lollipop graph of the color change for various mutants in a plate
# Created: 12-06-19 by Zhengxing Xue
# Last Updated: 11-07-19 by Zhengxing Xue
# ------------------------------------------------------------------------------------------------ #

from utils.find_wells_helper import *
from utils.plot_lollipop_helper import *
import numpy as np
import sys
from utils.input import get_input
import matplotlib.pyplot as plt
import tkinter as Tk

# reads input file
argv = sys.argv
inputParameters = ["color_info_file", "times", "save_location"]
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../inputs/' + proj_name + '/'
file_name = file_intro + 'PlotLollipopGraph.inp'
inputParameterValues = get_input(file_name, inputParameters)

color_info_file = inputParameterValues["color_info_file"]
timesToGet = [float(i) for i in inputParameterValues["times"].split(',')]
timesToGet.reverse()
save_loc = inputParameterValues["save_location"]

# creates dictionary of genes with time series of color values
gene_series = import_gene_color_change_dict(color_info_file)
all_genes = list(gene_series.keys())

# generates GUI for plotting
root = Tk.Tk()

# creates gene plotting options
gene_options = Tk.Listbox(root, selectmode=Tk.MULTIPLE, exportselection=0)
gene_options.pack()
for gene in all_genes:
    gene_options.insert(Tk.END, gene)


def plot_stuff():

    genes_chosen = []
    for gen in list(gene_options.curselection()):
        genes_chosen.append(all_genes[int(gen)])

    summarizedDataArray = []

    for gene_name in genes_chosen:

        tempDataDict = {}

        timeArray = np.array(gene_series[gene_name]['times'])
        timeArray = (timeArray - timeArray[0]) / 60 # converts to minutes
        print(timeArray)

        selectedMeanRedPoints = []
        selectedMeanGreenPoints = []
        selectedMeanBluePoints = []

        for timeToGet in timesToGet:
            index = np.argmin(np.abs(timeArray - float(timeToGet)))
            print(index)
            print(timeToGet)

            selectedMeanRedPoints.append(np.mean(gene_series[gene_name]['mean_red'], axis=0)[index])
            selectedMeanGreenPoints.append(np.mean(gene_series[gene_name]['mean_green'], axis=0)[index])
            selectedMeanBluePoints.append(np.mean(gene_series[gene_name]['mean_blue'], axis=0)[index])

        tempDataDict['MeanRed'] = selectedMeanRedPoints
        tempDataDict['MeanGreen'] = selectedMeanGreenPoints
        tempDataDict['MeanBlue'] = selectedMeanBluePoints
        tempDataDict['RelativeTimeInterp'] = timesToGet
        tempDataDict['GeneName'] = gene_name

        summarizedDataArray.append(tempDataDict)

    PlotLollipopGraph(summarizedDataArray, timesToGet, save_loc, plotGridX=8500)
    root.destroy()


# button for creating the plots
plot_stuff_button = Tk.Button(root, text="Make Lollipop Graphs", command=plot_stuff)
plot_stuff_button.pack()

Tk.mainloop()