def num_from_well_ID(wellID):
    return (ord(wellID[0])-65) * 12 + int(wellID[1:3])


def import_gene_color_change_dict(color_info_file):
    import xml.etree.ElementTree as ET
    import numpy as np
    print(color_info_file)

    tree = ET.parse(color_info_file)
    root = tree.getroot()
    gene_series = {}

    wells = root.findall('well')
    color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'median_yellow',
                  'mean_yellow']

    for well in wells:
        wellID = well.attrib['wellID']
        gene = well.attrib['GeneName']
        well_row = ord(wellID[0]) - 65
        well_column = int(wellID[1:3])
        times = np.array(well.attrib['times'].split(','),float)
        times = times - times[0]
        ordered_indices = np.argsort(times)
        # if gene in gene_series:
        #     for color in color_keys:
        #         color_series = np.array(well.attrib[color].split(','), float)
        #         gene_series[gene][color] = np.vstack((gene_series[gene][color],
        #                                              np.reshape(color_series[ordered_indices], (1, len(times)))))
        # else:
        #     gene_series[gene] = {}
        #     gene_series[gene]["times"] = times
        #     for color in color_keys:
        #         color_series = np.array(well.attrib[color].split(','), float)
        #         gene_series[gene][color] = np.reshape(color_series[ordered_indices], (1, len(times)))
        
        gene_series[wellID] = {}
        gene_series[wellID]["gene"] = gene
        gene_series[wellID]["times"] = times
        for color in color_keys:
            color_series = np.array(well.attrib[color].split(','), float)
            gene_series[wellID][color] = np.reshape(color_series[ordered_indices], (1, len(times)))
        gene_series[wellID]["well_row"] = well_row
        gene_series[wellID]["well_column"] = well_column

    return gene_series

def get_input(file_name, inputParameters):
    import re
    import pdb

    file = open(file_name,'r')
    inputData = file.readlines()

    # Go through all lines and search for input parameters

    inputParameterValues = {}

    for item in inputParameters:
        result = inputParameterSearch(inputData, item)
        if result != None:
            inputParameterValues[item] = result

    return inputParameterValues
