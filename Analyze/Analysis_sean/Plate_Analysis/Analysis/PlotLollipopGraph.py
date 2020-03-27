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
inputParameters = ["color_info_file"]
if len(argv) != 2:
    print("Error: IndexSummarize takes 1 argument, the input file name")
    sys.exit(-1)
proj_name = sys.argv[1]
file_intro = '../inputs/' + proj_name + '/'
file_name = file_intro + 'PlotColorTimeSeries.inp'
inputParameterValues = get_input(file_name, inputParameters)

color_info_file = inputParameterValues["color_info_file"]

# creates dictionary of genes with time series of color values
gene_series = import_gene_color_change_dict(color_info_file)
all_genes = list(gene_series.keys())
color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'mean_yellow',
              'median_yellow']

selectedGeneNames = ['Q_WT', 'Blank']

selectedGenesArray = []

for gene in selectedGeneNames:

    selectedGenesArray.append(gene_series[gene])

times = gene_series['Blank']['times']

timesToGet = [times[0], times[10], times[20], times[30],times[40]]

list.reverse(timesToGet)

summarizedDataArray = []

for gene, gene_name in zip(selectedGenesArray, selectedGeneNames):

    tempDataDict = {}

    timeArray = gene['times']

    selectedMeanRedPoints = []
    selectedMeanGreenPoints = []
    selectedMeanBluePoints = []

    for timeToGet in timesToGet:
        index = np.where (timeArray == float(timeToGet))

        selectedMeanRedPoints.append(np.mean(gene_series[gene_name]['mean_red'], axis=0)[index])
        selectedMeanGreenPoints.append(np.mean(gene_series[gene_name]['mean_green'], axis=0)[index])
        selectedMeanBluePoints.append(np.mean(gene_series[gene_name]['mean_blue'], axis=0)[index])

    tempDataDict['MeanRed'] = selectedMeanRedPoints
    tempDataDict['MeanGreen'] = selectedMeanGreenPoints
    tempDataDict['MeanBlue'] = selectedMeanBluePoints
    tempDataDict['relativeTimeInterp'] = timesToGet
    tempDataDict['GeneName'] = gene_name

    summarizedDataArray.append(tempDataDict)

PlotLollipopGraph(summarizedDataArray, timesToGet, '/Users/barstowlab/Downloads/test', plotGridX=8500)