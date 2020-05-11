# # import numpy as np
# # from tkinter import filedialog
# # import tkinter as tk
# # import csv




# # # Set up the Tkinter Window
# # master = tk.Tk()
# # master.title("Jurgensen Well Analyzer Demo")
# # master.geometry("500x300")

# # # Select Source Folder
# # source_dir = "Select Source Directory"
# # def select_source_folder():
# #     source_dir = filedialog.askdirectory()
# #     return

# # # Select Destination Folder
# # dest_dir = tk.StringVar(master)
# # dest_dir.set("Select Destination Directory")
# # def select_dest_folder():
# #     dest_dir.set(filedialog.askdirectory())
# #     print(dest_dir)
# #     return

# # source_text = tk.Label(master, text = source_dir)
# # source_text.pack()
# # source_button = tk.Button(master, text = "Source Folder", command = select_source_folder)
# # source_button.pack()

# # dest_text = tk.Label(master, textvariable = dest_dir)
# # dest_text.pack()
# # dest_button = tk.Button(master, text = "Dest Folder", command = select_dest_folder)
# # dest_button.pack()


# # # Run main loop
# # master.mainloop()

# import tkinter

# from matplotlib.backends.backend_tkagg import (
#     FigureCanvasTkAgg, NavigationToolbar2Tk)
# # Implement the default Matplotlib key bindings.
# from matplotlib.backend_bases import key_press_handler
# from matplotlib.figure import Figure

# import numpy as np


# root = tkinter.Tk()
# root.wm_title("Embedding in Tk")

# fig = Figure(figsize=(5, 4), dpi=100)
# t = np.arange(0, 3, .01)
# fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))

# canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
# canvas.draw()

# toolbar = NavigationToolbar2Tk(canvas, root)
# toolbar.update()


# def on_key_press(event):
#     print("you pressed {}".format(event.key))
#     key_press_handler(event, canvas, toolbar)


# canvas.mpl_connect("key_press_event", on_key_press)

# button = tkinter.Button(master=root, text="Quit", command=root.quit)

# # Packing order is important. Widgets are processed sequentially and if there
# # is no space left, because the window is too small, they are not displayed.
# # The canvas is rather flexible in its size, so we pack it last which makes
# # sure the UI controls are displayed as long as possible.
# button.pack(side=tkinter.BOTTOM)
# canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

# tkinter.mainloop()

import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np


root = tk.Tk()
# frames to create the black border:
fig_frame1 = tk.Frame(root, background='black', padx=2, pady=2)  

fig1, ax1 = plt.subplots(figsize=(2, 2))
t = np.arange(0, 2*np.pi, 0.1)

ax1.plot(t, np.cos(t))
fig1.tight_layout()
fig1_canvas = FigureCanvasTkAgg(fig1, master=fig_frame1)  # set master of fig1_canvas to the border frame
fig1_canvas.get_tk_widget().pack(padx=10, pady=10)  # change padx and pady to choose the thickness of the border


frame1 = tk.Frame()

tk.Label(frame1, text='hello').pack()

frame1.grid(row=0, column=0, rowspan=2)
# put the border frames in the root window
fig_frame1.grid(row=0, column=1)

root.mainloop()