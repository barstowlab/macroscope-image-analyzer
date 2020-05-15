import numpy as np
from tkinter import filedialog
import tkinter as tk
import analysis_functions
import threading
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
    from tkinter import ttk, HORIZONTAL
    import matplotlib
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.backend_bases import key_press_handler
    from matplotlib.figure import Figure
    import matplotlib.pyplot as plt
    import temp_functions
    import os
    matplotlib.use("TkAgg")

    master.withdraw()
    master.destroy()
    
    # Generate Variables
    analyze = tk.Tk()
    analyze.title("Analyze Images")
    analyze.geometry("1600x800")

    source_directory = tk.StringVar(analyze)
    current_data = tk.StringVar(analyze)
    single_or_multiple = tk.StringVar(analyze)
    single_or_multiple.set("multiple")

    left_frame = tk.Frame(analyze, relief=tk.RAISED, bd=2)
    figure_frame = tk.Frame(analyze, background='black', padx=2, pady=2)
    right_frame = tk.Frame(analyze, relief=tk.RAISED, bd=2)


############### LEFT FRAME   ###########################
    # Source directory
    souce_text = tk.Label(left_frame, textvariable=source_directory)
    souce_text.grid(row=0,column=0, sticky="ew", padx=5, pady=5)

    lbox_files = tk.Listbox(left_frame, exportselection = 0)
    lbox_files.grid(row=1,column=0, sticky="ew", padx=5, pady=5)

    def select_source_folder():
        source_directory.set(filedialog.askdirectory())

        lbox_files.delete(0, tk.END)
        file_list = os.listdir(source_directory.get())
        for item in file_list:
            lbox_files.insert(tk.END, item)
        lbox_files.update()
        return

    tk.Button(left_frame, text = "Select Source", command=select_source_folder).grid(\
        row=2,column=0, sticky="ew", padx=5, pady=5)
    # tk.Button(analyze, text = "fill", command = fill_file_list_box).pack(side = tk.LEFT, padx=10)


    color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'mean_yellow',
              'median_yellow']

    lbox_colors = tk.Listbox(left_frame, exportselection = 0)
    lbox_colors.grid(row=3,column=0, sticky="ew", padx=5, pady=5)
    for color in color_keys:
        lbox_colors.insert(tk.END, color)
    lbox_colors.update()


    # Namespace and execute functions
    text_entry = tk.Text(left_frame, height=10, width=30).grid(row=4,column=0, sticky="ew", padx=5, pady=5)


############### RIGHT FRAME  ###########################
    operation = tk.StringVar(analyze)
    operation.set("view")
    tk.Label(right_frame, text = "Click Operation:").grid(row=0,column=0, sticky="ew", padx=5, pady=5)
    tk.Radiobutton(right_frame, text = "Select Hits", variable=operation, value = "select").grid(row=1,column=0, sticky="ew", padx=5, pady=5)
    tk.Radiobutton(right_frame, text = "View Graph", variable=operation, value = "view").grid(row=2,column=0, sticky="ew", padx=5, pady=5)
    lbox_hits = tk.Listbox(right_frame, exportselection = 0)
    lbox_hits.grid(row=3,column=0, sticky="ew", padx=5, pady=5)


    hit_list = []

    def update_hit_list():
        lbox_hits.delete(0,tk.END)
        for item in hit_list:
            lbox_hits.insert(tk.END,item)
        lbox_hits.update()
        return

    def add_hit(gene):
        if gene not in hit_list:
            hit_list.append(gene)
            update_hit_list()
        return
    
    def remove_hit():
        gene = lbox_hits.get(lbox_hits.curselection()[0])
        if gene in hit_list:
            hit_list.remove(gene)
            update_hit_list()
        return

    def save_hits():
        save_location = filedialog.asksaveasfilename()
        with open(save_location, 'w') as f:
            for hit in hit_list:
                f.writelines(hit + "\n")


    tk.Button(right_frame, text="Remove Hit", command=remove_hit).grid(row=4,column=0,sticky="ew",padx=5,pady=5)
    tk.Button(right_frame, text="Save Hits", command=save_hits).grid(row=5,column=0,sticky="ew",padx=5,pady=5)

    
############### MIDDLE FRAME ###########################
    def onclick(event):
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

                if operation.get() == "select":
                    file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
                    # color = lbox_colors.get(lbox_colors.curselection()[0])
                    color_info_file = file
                    gene_series = temp_functions.import_gene_color_change_dict(color_info_file)
                    add_hit(gene_series[wellID]["gene"])
                    

                if operation.get() == "view":
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
                    plt.ylim([0,300])
                    plt.plot(times[idxes], mean_colors[idxes])
                    abx.show()
                else:
                    print("invalid value for variable 'operation'")

        return


    fig, axs = plt.subplots(8,12, sharex='all', sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})

    for a in range(len(axs)):
        for b in range(len(axs[a])):
            axs[a,b].axes.get_yaxis().set_visible(False)
            axs[a,b].axes.get_xaxis().set_visible(False)

    fig.figsize=[16,10]
    fig.dpi = 120
    
    canvas = FigureCanvasTkAgg(fig, master=figure_frame)
    fig.tight_layout()
    canvas.draw()


    def update_graph():
        try:
            file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
            current_data.set(lbox_files.get(lbox_files.curselection()[0]))
            print(file)

            color = lbox_colors.get(lbox_colors.curselection()[0])
            color_info_file = file
            gene_series = temp_functions.import_gene_color_change_dict(color_info_file)
            all_genes = list(gene_series.keys())
            times = gene_series["A01"]['times']
            min_time = min(times)
            max_time = max(times)
            idxes = np.logical_and(min_time <= times, max_time >= times)

            current_data.set("Processing '" + lbox_files.get(lbox_files.curselection()[0]) + "' ..." )
            figure_frame.update()
            for wellID in gene_series:
                a = gene_series[wellID]["well_row"] - 0
                b = gene_series[wellID]["well_column"] - 1
                x = a*12+b
                print(a,b,x)
                mean_colors = np.mean(gene_series[wellID][color], axis=0)
                axs[a ,b].cla()

                axs[a,b].tick_params(
                    axis='x',
                    bottom=False,    
                    top=False)
                axs[a,b].tick_params(
                    axis='y',     
                    left=False,    
                    right=False)
                # axs[a,b].set_title("demo")
                if gene_series[wellID]["gene"] != "Blank":
                    axs[a,b].text(0,1,gene_series[wellID]["gene"], fontsize=5, style='italic')
                axs[a,b].set_ylim([0,300])
                axs[a,b].axes.get_yaxis().set_visible(False)
                axs[a,b].axes.get_xaxis().set_visible(False)
                if gene_series[wellID]["gene"] == "Blank":
                    axs[a, b].plot(times[idxes], mean_colors[idxes], color = "gray")
                else:
                    axs[a, b].plot(times[idxes], mean_colors[idxes], color = "blue")
            current_data.set(lbox_files.get(lbox_files.curselection()[0]) + "  by  " + color)
            canvas.draw()
        except IndexError:
            current_data.set("Error: No File or Color Selected")
        except:
            current_data.set("Error: Issue with file: '" + lbox_files.get(lbox_files.curselection()[0]) + "'")

    tk.Label(figure_frame, textvariable=current_data).grid(row=1,column=0, sticky="ew", padx=5, pady=5)
    tk.Button(figure_frame, text = "Update Graph", command=update_graph).grid(row=2,column=0, sticky="ew", padx=5, pady=5)

    clicker = fig.canvas.mpl_connect('button_press_event', onclick)
    canvas.get_tk_widget().grid(row=0,column=0, sticky="ew", padx=5, pady=5)


############### COMBINE      ###########################


    def on_closing():
        plt.close("all")
        analyze.destroy()

    analyze.protocol("WM_DELETE_WINDOW", on_closing)

    left_frame.grid(row=0,column=0,sticky="ew",padx=5,pady=5)
    figure_frame.grid(row=0,column=1, sticky="ew", padx=5, pady=5)
    right_frame.grid(row=0, column=2, sticky="ew", padx=5, pady=5)

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

