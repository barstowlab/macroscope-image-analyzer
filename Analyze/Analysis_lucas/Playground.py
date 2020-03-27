import numpy as np
from tkinter import filedialog
import tkinter as tk
import csv




# Set up the Tkinter Window
master = tk.Tk()
master.title("Jurgensen Well Analyzer Demo")
master.geometry("500x300")

# Select Source Folder
source_dir = "Select Source Directory"
def select_source_folder():
    source_dir = filedialog.askdirectory()
    return

# Select Destination Folder
dest_dir = tk.StringVar(master)
dest_dir.set("Select Destination Directory")
def select_dest_folder():
    dest_dir.set(filedialog.askdirectory())
    print(dest_dir)
    return

source_text = tk.Label(master, text = source_dir)
source_text.pack()
source_button = tk.Button(master, text = "Source Folder", command = select_source_folder)
source_button.pack()

dest_text = tk.Label(master, textvariable = dest_dir)
dest_text.pack()
dest_button = tk.Button(master, text = "Dest Folder", command = select_dest_folder)
dest_button.pack()


# Run main loop
master.mainloop()