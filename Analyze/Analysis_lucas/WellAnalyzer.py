import numpy as np
from tkinter import filedialog
import tkinter as tk
import analysis_functions
import csv

TK_SILENCE_DEPRECATION=1

def create_notification(message):
    notification = tk.Tk()
    notification.title("Notification")
    notification.geometry("400x50")

    tk.Label(notification, text = message).pack()

    tk.Button(notification, text = "Close", command = notification.destroy).pack()
    notification.mainloop()



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
        print("HELLO")
        generate.destroy()
        create_notification("Saved file to Local Directory")

        print("MODDLE")
        output_file_and_directory = dest_directory.get() + "/" + output_file.get()
        print(output_file_and_directory)
        print("GOODBYE")

        analysis_functions.generate_barcode(plate_type.get(), output_file_and_directory, \
            barcode_prefix.get(), barcode_postfix.get(), max_rows.get(), max_columns.get(), \
                padding_columns.get(), padding_rows.get(), barcode_number_start.get(), \
                    barcode_number_end.get()+1)




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
    def select_dest_folder():
        dest_directory.set(filedialog.askdirectory())

    tk.Button(generate, text = "Select Destination", command = select_dest_folder).pack()
    
    # Generate Barcodes - finished with inputs
    tk.Button(generate, text = "Generate Barcodes", command = run_generate).pack()
    

    generate.mainloop()

# Organize Screen


# Analyze Screen




# Set up Main Selection Screen
master = tk.Tk()
master.title("Jurgensen Well Analyzer Demo")
master.geometry("500x300")

choice_label = tk.Label(master, text = "Select which function you want to perform")
generate_button = tk.Button(master, text = "Generate Barcodes", command = open_generate)
choice_label.pack()
generate_button.pack()


organize_button = tk.Button(master, text = "Organize Images by Barcodes")

analyze_button = tk.Button(master, text = "Analyze Images")

master.mainloop()