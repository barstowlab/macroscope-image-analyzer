import numpy as np
from tkinter import filedialog
import tkinter as tk
import analysis_functions
import csv

tk.TK_SILENCE_DEPRECATION=1

###############################################

def create_notification(message):
    notification = tk.Tk()
    notification.title("Notification")
    notification.geometry("400x50")

    tk.Label(notification, text = message).pack()

    tk.Button(notification, text = "Close", command = notification.destroy).pack()
    notification.mainloop()

###############################################

# Generate Screen
def open_generate():
    master.withdraw()
    master.destroy()
    
    # Generate Variables
    generate = tk.Tk()
    generate.title("Generate Barcodes")
    generate.geometry("500x650")

    plate_type = tk.StringVar(generate)
    output_file = tk.StringVar(generate)
    barcode_prefix = tk.StringVar(generate)
    barcode_postfix = tk.StringVar(generate)
    max_rows = tk.IntVar(generate)
    max_columns = tk.IntVar(generate)
    padding_columns = tk.BooleanVar(generate)
    padding_rows = tk.BooleanVar(generate)
    barcode_number_start = tk.IntVar(generate)
    barcode_number_end = tk.IntVar(generate)
    dest_directory = tk.StringVar(generate)

    def run_generate():
        if output_file.get() == "":
            create_notification("ERROR: Must have created an Output File Name")
            return
        if dest_directory.get() == "":
            create_notification("ERROR: Must select a destination")
            return

        generate.destroy()

        create_notification("Saved file to Local Directory")

        output_file_and_directory = dest_directory.get() + "/" + output_file.get()

        analysis_functions.generate_barcode(plate_type.get(), output_file_and_directory, \
            barcode_prefix.get(), barcode_postfix.get(), max_rows.get(), max_columns.get(), \
                padding_columns.get(), padding_rows.get(), barcode_number_start.get(), \
                    barcode_number_end.get()+1)
        return




    # Plate Type
    tk.Label(generate, text = "Type of Plate:").pack()
    plate_type.set("Assay")
    tk.Radiobutton(generate, text = "Assay", variable = plate_type, value = "Assay").pack()
    tk.Radiobutton(generate, text = "Storage", variable = plate_type, value = "Storage").pack()

    # Prefix
    tk.Label(generate, text = "Barcode Prefix").pack()
    tk.Entry(generate, textvariable = barcode_prefix).pack()

    # Postfix
    tk.Label(generate, text = "Barcode Postfix").pack()
    tk.Entry(generate, textvariable = barcode_postfix).pack()

    # Output File Name
    tk.Label(generate, text = "Output File name (end with '.csv'):").pack()
    tk.Entry(generate, textvariable = output_file).pack()

    # Max Rows:
    tk.Label(generate, text = "Max Rows").pack()
    max_rows.set(39)
    tk.Entry(generate, textvariable = max_rows).pack()
    
    # Max Columns:
    tk.Label(generate, text = "Max Columns").pack()
    max_columns.set(4)
    tk.Entry(generate, textvariable = max_columns).pack()

    # Padding Columns
    tk.Label(generate, text = "Padding Colums?").pack()
    padding_columns.set(True)
    tk.Radiobutton(generate, text = "True", variable = padding_columns, value = True).pack()
    tk.Radiobutton(generate, text = "False", variable = padding_columns, value = False).pack()
    
    # Padding Rows
    tk.Label(generate, text = "Padding Rows?").pack()
    padding_rows.set(False)
    tk.Radiobutton(generate, text = "True", variable = padding_rows, value = True).pack()
    tk.Radiobutton(generate, text = "False", variable = padding_rows, value = False).pack()

    # Barcode Number Start
    tk.Label(generate, text = "First Barcode Number").pack()
    tk.Entry(generate, textvariable = barcode_number_start).pack()

    # Barcode Number End
    tk.Label(generate, text = "Final Barcode Number").pack()
    tk.Entry(generate, textvariable = barcode_number_end).pack()




    # Save Location, destination directory
    tk.Label(generate, text = dest_directory.get()).pack()
    def select_dest_folder():
        dest_directory.set(filedialog.askdirectory())

    tk.Button(generate, text = "Select Destination", command = select_dest_folder).pack()
    
    # Generate Barcodes - finished with inputs
    tk.Button(generate, text = "Generate Barcodes", command = run_generate).pack()
    

    generate.mainloop()

###############################################

# Organize Screen

def open_organize():
    master.withdraw()
    master.destroy()
    
    # Generate Variables
    organize = tk.Tk()
    organize.title("Organize by Barcodes")
    organize.geometry("500x300")
    UNSET_DESTINATION = "Select Destination"

    plate_type = tk.StringVar(organize)
    source_directory = tk.StringVar(organize)
    source_directory.set(UNSET_DESTINATION)
    dest_directory = tk.StringVar(organize)
    dest_directory.set(UNSET_DESTINATION)

    def run_organize():
        if source_directory.get() ==  "" or source_directory.get() == UNSET_DESTINATION:
            create_notification("ERROR: Must designate a source directory")
            return
        if dest_directory.get() ==  "" or dest_directory.get() == UNSET_DESTINATION:
            create_notification("ERROR: Must designate a destination directory")
            return
        analysis_functions.organize_by_barcode(plate_type.get(), source_directory.get(),\
            dest_directory.get())
        return


    # Plate Type
    tk.Label(organize, text = "Type of Plate:").pack()
    plate_type.set("Assay")
    tk.Radiobutton(organize, text = "Assay", variable = plate_type, value = "Assay").pack()
    tk.Radiobutton(organize, text = "Storage", variable = plate_type, value = "Storage").pack()

    # Destination Directory directory
    dest_text = tk.Label(organize, textvariable = dest_directory)
    dest_text.pack()
    def select_dest_folder():
        dest_directory.set(filedialog.askdirectory())

    tk.Button(organize, text = "Select Destination", command = select_dest_folder).pack()

    # Source directory
    souce_text = tk.Label(organize , textvariable = source_directory)
    souce_text.pack()
    def select_source_folder():
        source_directory.set(filedialog.askdirectory())

    tk.Button(organize, text = "Select Source", command = select_source_folder).pack()

    # Organize by Barcodes - finished with inputs
    tk.Button(organize, text = "Organize by Barcodes", command = run_organize).pack()

    
     
    organize.mainloop()

###############################################

# Analyze Screen

def open_analyze():
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.backend_bases import key_press_handler
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    import temp_functions
    import os

    master.withdraw()
    master.destroy()
    
    # Generate Variables
    analyze = tk.Tk()
    analyze.title("Analyze Images")
    analyze.geometry("1600x800")

    source_directory = tk.StringVar(analyze)
    single_or_multiple = tk.StringVar(analyze)
    single_or_multiple.set("multiple")

    # Source directory
    souce_text = tk.Label(analyze , textvariable=source_directory)
    souce_text.pack()

    lbox_files = tk.Listbox(analyze, exportselection = 0)
    lbox_files.pack(side=tk.LEFT, padx =20)

    def select_source_folder():
        source_directory.set(filedialog.askdirectory())

        lbox_files.delete(0, tk.END)
        file_list = os.listdir(source_directory.get())
        for item in file_list:
            lbox_files.insert(tk.END, item)
        lbox_files.update()
        return

    tk.Button(analyze, text = "Select Source", command=select_source_folder).pack(side = tk.LEFT, padx=10)
    # tk.Button(analyze, text = "fill", command = fill_file_list_box).pack(side = tk.LEFT, padx=10)


    color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'mean_yellow',
              'median_yellow']

    lbox_colors = tk.Listbox(analyze, exportselection = 0)
    lbox_colors.pack(side=tk.LEFT, padx =20)
    for color in color_keys:
        lbox_colors.insert(tk.END, color)
    lbox_colors.update()

    

##############
    def onclick(event):
        ix = event.xdata

        axss = []
        for a in range(8):
            _a = a
            for b in range(12):
                _b = b
                axss.append(axs[_a,_b])

        for i, ax in enumerate(axss):
            if ax == event.inaxes:
                i += 1
                a = chr(int((i-1)/12)  + 65)
                b = (i-1) % 12 + 1

                if b < 10:
                    b = "0" + str(b)
                else:
                    b = str(b)

                print(a,b)
                wellID = a + b

                abx = plt.figure()
                abx.add_subplot(111)

                file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
                color = lbox_colors.get(lbox_colors.curselection()[0])
                color_info_file = file
                gene_series = temp_functions.import_gene_color_change_dict(color_info_file)

                times = gene_series["A01"]['times']
                min_time = min(times)
                max_time = max(times)
                idxes = np.logical_and(min_time <= times, max_time >= times)

                print("COLORS for " + wellID)
                mean_colors = np.mean(gene_series[wellID][color], axis=0)

                plt.title(gene_series[wellID]["gene"])
                plt.plot(times[idxes], mean_colors[idxes])
                abx.show()
        return

##############

    # Add in graphs

    figure_frame = tk.Frame(analyze, background='black', padx=2, pady=2)

    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)

    fig, axs = plt.subplots(8,12, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
    fig.figsize=[16,10]
    fig.dpi = 120
    
    canvas = FigureCanvasTkAgg(fig, master=figure_frame)
    fig.tight_layout()
    canvas.draw()



    def update_graph():
        file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
        print(file)

        color = lbox_colors.get(lbox_colors.curselection()[0])
        color_info_file = file
        gene_series = temp_functions.import_gene_color_change_dict(color_info_file)

        all_genes = list(gene_series.keys())
        times = gene_series["A01"]['times']
        min_time = min(times)
        max_time = max(times)
        idxes = np.logical_and(min_time <= times, max_time >= times)

        for wellID in gene_series:
            a = gene_series[wellID]["well_row"] - 0
            b = gene_series[wellID]["well_column"] - 1
            x = a*12+b
            print(a,b,x)
            mean_colors = np.mean(gene_series[wellID][color], axis=0)
            axs[a ,b].cla()

            axs[a,b].tick_params(
                axis='x',
                which='both',
                bottom=False,    
                top=False,      
                labelbottom=False)
            axs[a,b].tick_params(
                axis='y',     
                which='both',     
                left=False,    
                right=False,      
                labelbottom=False)

            if gene_series[wellID]["gene"] == "Blank":
                axs[a, b].plot(times[idxes], mean_colors[idxes], color = "gray")
                print(type(axs[a,b]))
            else:
                axs[a, b].plot(times[idxes], mean_colors[idxes], color = "blue")

        canvas.draw()


    canvas.get_tk_widget().pack(padx=4, pady=4)
    figure_frame.pack()

    tk.Button(analyze, text = "Update Graph", command=update_graph).pack()

    clicker = fig.canvas.mpl_connect('button_press_event', onclick)

    analyze.mainloop()


###############################################

# Set up Main Selection Screen
master = tk.Tk()
master.title("Jurgensen Well Analyzer Demo")
master.geometry("300x150")

tk.Label(master, text = "Select which function you want to perform").pack()
tk.Button(master, text = "Generate Barcodes", command = open_generate).pack()
tk.Button(master, text = "Organize by Barcodes", command = open_organize).pack()
tk.Button(master, text = "Analyze Plates", command = open_analyze).pack()

master.mainloop()

