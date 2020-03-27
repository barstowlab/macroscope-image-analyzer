# Author: Sean Medin
# This code is used to plot the time series of the color change for various mutants in a plate
# Created: 11-07-19 by Sean Medin
# Last Updated: 11-07-19 by Sean Medin
# ------------------------------------------------------------------------------------------------ #

from utils.find_wells_helper import *
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

# generates GUI for plotting
root = Tk.Tk()

# creates gene plotting options
gene_options = Tk.Listbox(root, selectmode=Tk.MULTIPLE, exportselection=0)
gene_options.pack()
for gene in all_genes:
    gene_options.insert(Tk.END, gene)

# creates color plotting options
color_options = Tk.Listbox(root, selectmode=Tk.MULTIPLE, exportselection=0)
color_options.pack()
for color in color_keys:
    color_options.insert(Tk.END, color)

# generates time interval plotting options
Tk.Label(root, text="Enter Start Time (minutes)").pack()
lower_bound = Tk.Entry(root)
lower_bound.pack()
Tk.Label(root, text="Enter End Time (minutes)").pack()
higher_bound = Tk.Entry(root)
higher_bound.pack()

# function for plotting stuff
def plot_stuff():
    genes_chosen = []
    colors_chosen = []
    for gen in list(gene_options.curselection()):
        genes_chosen.append(all_genes[int(gen)])
    for col in list(color_options.curselection()):
        colors_chosen.append(color_keys[int(col)])
        
    times = gene_series[genes_chosen[0]]['times']
    min_time_str = lower_bound.get()
    max_time_str = higher_bound.get()
    min_time = min(times)
    max_time = max(times)
    if min_time_str.isdigit():
        min_time = int(min_time_str) * 60
    if max_time_str.isdigit():
        max_time = int(max_time_str) * 60

    plt.figure()
    idxes = np.logical_and(min_time <= times, max_time >= times)
    items_plotted = []
    for gene in genes_chosen:
        for color in colors_chosen:
            mean_colors = np.mean(gene_series[gene][color], axis=0)
            std_colors = np.std(gene_series[gene][color], axis=0)
#           plt.errorbar(times[idxes], mean_colors[idxes], std_colors[idxes])
            plt.plot(times[idxes], mean_colors[idxes])
            items_plotted.append(gene + " " + color)
    plt.legend(items_plotted)
    plt.title("Comparing colors for various genes")
    plt.xlabel('time (seconds)')
    plt.ylabel('degree of color')
    plt.show()

# button for creating the plots
plot_stuff_button = Tk.Button(root, text="Plot Stuff", command=plot_stuff)
plot_stuff_button.pack()

Tk.mainloop()
