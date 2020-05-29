from tkinter import *
from tkinter import ttk
import time

MAX = 30

root = Tk()
root.geometry('{}x{}'.format(400, 100))
progress_var = DoubleVar() #here you have ints but when calc. %'s usually floats
theLabel = Label(root, text="Sample text to show")
theLabel.pack()
progressbar = ttk.Progressbar(root, variable=progress_var, maximum=MAX)
progressbar.pack(fill=X, expand=1)


def loop_function():

    k = 0
    while k <= MAX:
    ### some work to be done
        progress_var.set(k)
        k += 1
        time.sleep(0.02)
        root.update_idletasks()
    root.after(100, loop_function)

loop_function()
root.mainloop()
# from tkinter import *
# from tkinter.messagebox import *
# from tkinter.filedialog import *
  
# class Notepad: 
  
#     __root = Tk() 
  
#     # default window width and height 
#     __thisWidth = 300
#     __thisHeight = 300
#     __thisTextArea = Text(__root) 
#     __thisMenuBar = Menu(__root) 
#     __thisFileMenu = Menu(__thisMenuBar, tearoff=0) 
#     __thisEditMenu = Menu(__thisMenuBar, tearoff=0) 
#     __thisHelpMenu = Menu(__thisMenuBar, tearoff=0) 
      
#     # To add scrollbar 
#     __thisScrollBar = Scrollbar(__thisTextArea)      
#     __file = None
  
#     def __init__(self,**kwargs): 
  
#         # Set icon 
#         try: 
#                 self.__root.wm_iconbitmap("Notepad.ico")  
#         except: 
#                 pass
  
#         # Set window size (the default is 300x300) 
  
#         try: 
#             self.__thisWidth = kwargs['width'] 
#         except KeyError: 
#             pass
  
#         try: 
#             self.__thisHeight = kwargs['height'] 
#         except KeyError: 
#             pass
  
#         # Set the window text 
#         self.__root.title("Untitled - Notepad") 
  
#         # Center the window 
#         screenWidth = self.__root.winfo_screenwidth() 
#         screenHeight = self.__root.winfo_screenheight() 
      
#         # For left-alling 
#         left = (screenWidth / 2) - (self.__thisWidth / 2)  
          
#         # For right-allign 
#         top = (screenHeight / 2) - (self.__thisHeight /2)  
          
#         # For top and bottom 
#         self.__root.geometry('%dx%d+%d+%d' % (self.__thisWidth, 
#                                               self.__thisHeight, 
#                                               left, top))  
  
#         # To make the textarea auto resizable 
#         self.__root.grid_rowconfigure(0, weight=1) 
#         self.__root.grid_columnconfigure(0, weight=1) 
  
#         # Add controls (widget) 
#         self.__thisTextArea.grid(sticky = N + E + S + W) 
          
#         # To open new file 
#         self.__thisFileMenu.add_command(label="New", 
#                                         command=self.__newFile)     
          
#         # To open a already existing file 
#         self.__thisFileMenu.add_command(label="Open", 
#                                         command=self.__openFile) 
          
#         # To save current file 
#         self.__thisFileMenu.add_command(label="Save", 
#                                         command=self.__saveFile)     
  
#         # To create a line in the dialog         
#         self.__thisFileMenu.add_separator()                                          
#         self.__thisFileMenu.add_command(label="Exit", 
#                                         command=self.__quitApplication) 
#         self.__thisMenuBar.add_cascade(label="File", 
#                                        menu=self.__thisFileMenu)      
          
#         # To give a feature of cut  
#         self.__thisEditMenu.add_command(label="Cut", 
#                                         command=self.__cut)              
      
#         # to give a feature of copy     
#         self.__thisEditMenu.add_command(label="Copy", 
#                                         command=self.__copy)          
          
#         # To give a feature of paste 
#         self.__thisEditMenu.add_command(label="Paste", 
#                                         command=self.__paste)          
          
#         # To give a feature of editing 
#         self.__thisMenuBar.add_cascade(label="Edit", 
#                                        menu=self.__thisEditMenu)      
          
#         # To create a feature of description of the notepad 
#         self.__thisHelpMenu.add_command(label="About Notepad", 
#                                         command=self.__showAbout)  
#         self.__thisMenuBar.add_cascade(label="Help", 
#                                        menu=self.__thisHelpMenu) 
  
#         self.__root.config(menu=self.__thisMenuBar) 
  
#         self.__thisScrollBar.pack(side=RIGHT,fill=Y)                     
          
#         # Scrollbar will adjust automatically according to the content         
#         self.__thisScrollBar.config(command=self.__thisTextArea.yview)      
#         self.__thisTextArea.config(yscrollcommand=self.__thisScrollBar.set) 
      
          
#     def __quitApplication(self): 
#         self.__root.destroy() 
#         # exit() 
  
#     def __showAbout(self): 
#         showinfo("Notepad","Mrinal Verma") 
  
#     def __openFile(self): 
          
#         self.__file = askopenfilename(defaultextension=".txt", 
#                                       filetypes=[("All Files","*.*"), 
#                                         ("Text Documents","*.txt")]) 
  
#         if self.__file == "": 
              
#             # no file to open 
#             self.__file = None
#         else: 
              
#             # Try to open the file 
#             # set the window title 
#             self.__root.title(os.path.basename(self.__file) + " - Notepad") 
#             self.__thisTextArea.delete(1.0,END) 
  
#             file = open(self.__file,"r") 
  
#             self.__thisTextArea.insert(1.0,file.read()) 
  
#             file.close() 
  
          
#     def __newFile(self): 
#         self.__root.title("Untitled - Notepad") 
#         self.__file = None
#         self.__thisTextArea.delete(1.0,END) 
  
#     def __saveFile(self): 
  
#         if self.__file == None: 
#             # Save as new file 
#             self.__file = asksaveasfilename(initialfile='Untitled.txt', 
#                                             defaultextension=".txt", 
#                                             filetypes=[("All Files","*.*"), 
#                                                 ("Text Documents","*.txt")]) 
  
#             if self.__file == "": 
#                 self.__file = None
#             else: 
                  
#                 # Try to save the file 
#                 file = open(self.__file,"w") 
#                 file.write(self.__thisTextArea.get(1.0,END)) 
#                 file.close() 
                  
#                 # Change the window title 
#                 self.__root.title(os.path.basename(self.__file) + " - Notepad") 
                  
              
#         else: 
#             file = open(self.__file,"w") 
#             file.write(self.__thisTextArea.get(1.0,END)) 
#             file.close() 
  
#     def __cut(self): 
#         self.__thisTextArea.event_generate("<<Cut>>") 
  
#     def __copy(self): 
#         self.__thisTextArea.event_generate("<<Copy>>") 
  
#     def __paste(self): 
#         self.__thisTextArea.event_generate("<<Paste>>") 
  
#     def run(self): 
  
#         # Run main application 
#         self.__root.mainloop() 
  
  
  
  
# # Run main application 
# notepad = Notepad(width=600,height=400) 
# notepad.run() 

# # from tkinter import Button, Tk, HORIZONTAL

# # from tkinter.ttk import Progressbar
# # import time
# # import threading

# # class MonApp(Tk):
# #     def __init__(self):
# #         super().__init__()


# #         self.btn = Button(self, text='Traitement', command=self.traitement)
# #         self.btn.grid(row=0,column=0)
# #         self.progress = Progressbar(self, orient=HORIZONTAL,length=100,  mode='indeterminate')


# #     def traitement(self):
# #         def real_traitement():
# #             self.progress.grid(row=1,column=0)
# #             self.progress.start()
# #             time.sleep(5)
# #             self.progress.stop()
# #             self.progress.grid_forget()

# #             self.btn['state']='normal'

# #         self.btn['state']='disabled'
# #         threading.Thread(target=real_traitement).start()

# # if __name__ == '__main__':

# #     app = MonApp()
# #     app.mainloop()





# # # import numpy as np
# # # from tkinter import filedialog
# # # import tkinter as tk
# # # import analysis_functions
# # # import csv

# # # tk.TK_SILENCE_DEPRECATION=1

# # # ###############################################

# # # def create_notification(message):
# # #     notification = tk.Tk()
# # #     notification.title("Notification")
# # #     notification.geometry("400x50")

# # #     tk.Label(notification, text = message).pack()

# # #     tk.Button(notification, text = "Close", command = notification.destroy).pack()
# # #     notification.mainloop()

# # # ###############################################

# # # # Generate Screen
# # # def open_generate():
# # #     master.withdraw()
# # #     master.destroy()
    
# # #     # Generate Variables
# # #     generate = tk.Tk()
# # #     generate.title("Generate Barcodes")
# # #     generate.geometry("500x650")

# # #     plate_type = tk.StringVar(generate)
# # #     output_file = tk.StringVar(generate)
# # #     barcode_prefix = tk.StringVar(generate)
# # #     barcode_postfix = tk.StringVar(generate)
# # #     max_rows = tk.IntVar(generate)
# # #     max_columns = tk.IntVar(generate)
# # #     padding_columns = tk.BooleanVar(generate)
# # #     padding_rows = tk.BooleanVar(generate)
# # #     barcode_number_start = tk.IntVar(generate)
# # #     barcode_number_end = tk.IntVar(generate)
# # #     dest_directory = tk.StringVar(generate)

# # #     def run_generate():
# # #         if output_file.get() == "":
# # #             create_notification("ERROR: Must have created an Output File Name")
# # #             return
# # #         if dest_directory.get() == "":
# # #             create_notification("ERROR: Must select a destination")
# # #             return

# # #         generate.destroy()

# # #         create_notification("Saved file to Local Directory")

# # #         output_file_and_directory = dest_directory.get() + "/" + output_file.get()

# # #         analysis_functions.generate_barcode(plate_type.get(), output_file_and_directory, \
# # #             barcode_prefix.get(), barcode_postfix.get(), max_rows.get(), max_columns.get(), \
# # #                 padding_columns.get(), padding_rows.get(), barcode_number_start.get(), \
# # #                     barcode_number_end.get()+1)
# # #         return




# # #     # Plate Type
# # #     tk.Label(generate, text = "Type of Plate:").pack()
# # #     plate_type.set("Assay")
# # #     tk.Radiobutton(generate, text = "Assay", variable = plate_type, value = "Assay").pack()
# # #     tk.Radiobutton(generate, text = "Storage", variable = plate_type, value = "Storage").pack()

# # #     # Prefix
# # #     tk.Label(generate, text = "Barcode Prefix").pack()
# # #     tk.Entry(generate, textvariable = barcode_prefix).pack()

# # #     # Postfix
# # #     tk.Label(generate, text = "Barcode Postfix").pack()
# # #     tk.Entry(generate, textvariable = barcode_postfix).pack()

# # #     # Output File Name
# # #     tk.Label(generate, text = "Output File name (end with '.csv'):").pack()
# # #     tk.Entry(generate, textvariable = output_file).pack()

# # #     # Max Rows:
# # #     tk.Label(generate, text = "Max Rows").pack()
# # #     max_rows.set(39)
# # #     tk.Entry(generate, textvariable = max_rows).pack()
    
# # #     # Max Columns:
# # #     tk.Label(generate, text = "Max Columns").pack()
# # #     max_columns.set(4)
# # #     tk.Entry(generate, textvariable = max_columns).pack()

# # #     # Padding Columns
# # #     tk.Label(generate, text = "Padding Colums?").pack()
# # #     padding_columns.set(True)
# # #     tk.Radiobutton(generate, text = "True", variable = padding_columns, value = True).pack()
# # #     tk.Radiobutton(generate, text = "False", variable = padding_columns, value = False).pack()
    
# # #     # Padding Rows
# # #     tk.Label(generate, text = "Padding Rows?").pack()
# # #     padding_rows.set(False)
# # #     tk.Radiobutton(generate, text = "True", variable = padding_rows, value = True).pack()
# # #     tk.Radiobutton(generate, text = "False", variable = padding_rows, value = False).pack()

# # #     # Barcode Number Start
# # #     tk.Label(generate, text = "First Barcode Number").pack()
# # #     tk.Entry(generate, textvariable = barcode_number_start).pack()

# # #     # Barcode Number End
# # #     tk.Label(generate, text = "Final Barcode Number").pack()
# # #     tk.Entry(generate, textvariable = barcode_number_end).pack()




# # #     # Save Location, destination directory
# # #     tk.Label(generate, text = dest_directory.get()).pack()
# # #     def select_dest_folder():
# # #         dest_directory.set(filedialog.askdirectory())

# # #     tk.Button(generate, text = "Select Destination", command = select_dest_folder).pack()
    
# # #     # Generate Barcodes - finished with inputs
# # #     tk.Button(generate, text = "Generate Barcodes", command = run_generate).pack()
    

# # #     generate.mainloop()

# # # ###############################################

# # # # Organize Screen

# # # def open_organize():
# # #     master.withdraw()
# # #     master.destroy()
    
# # #     # Generate Variables
# # #     organize = tk.Tk()
# # #     organize.title("Organize by Barcodes")
# # #     organize.geometry("500x300")
# # #     UNSET_DESTINATION = "Select Destination"

# # #     plate_type = tk.StringVar(organize)
# # #     source_directory = tk.StringVar(organize)
# # #     source_directory.set(UNSET_DESTINATION)
# # #     dest_directory = tk.StringVar(organize)
# # #     dest_directory.set(UNSET_DESTINATION)

# # #     def run_organize():
# # #         if source_directory.get() ==  "" or source_directory.get() == UNSET_DESTINATION:
# # #             create_notification("ERROR: Must designate a source directory")
# # #             return
# # #         if dest_directory.get() ==  "" or dest_directory.get() == UNSET_DESTINATION:
# # #             create_notification("ERROR: Must designate a destination directory")
# # #             return
# # #         analysis_functions.organize_by_barcode(plate_type.get(), source_directory.get(),\
# # #             dest_directory.get())
# # #         return


# # #     # Plate Type
# # #     tk.Label(organize, text = "Type of Plate:").pack()
# # #     plate_type.set("Assay")
# # #     tk.Radiobutton(organize, text = "Assay", variable = plate_type, value = "Assay").pack()
# # #     tk.Radiobutton(organize, text = "Storage", variable = plate_type, value = "Storage").pack()

# # #     # Destination Directory directory
# # #     dest_text = tk.Label(organize, textvariable = dest_directory)
# # #     dest_text.pack()
# # #     def select_dest_folder():
# # #         dest_directory.set(filedialog.askdirectory())

# # #     tk.Button(organize, text = "Select Destination", command = select_dest_folder).pack()

# # #     # Source directory
# # #     souce_text = tk.Label(organize , textvariable = source_directory)
# # #     souce_text.pack()
# # #     def select_source_folder():
# # #         source_directory.set(filedialog.askdirectory())

# # #     tk.Button(organize, text = "Select Source", command = select_source_folder).pack()

# # #     # Organize by Barcodes - finished with inputs
# # #     tk.Button(organize, text = "Organize by Barcodes", command = run_organize).pack()

    
     
# # #     organize.mainloop()

# # # ###############################################

# # # # Analyze Screen

# # # def open_analyze():
# # #     import matplotlib
# # #     from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# # #     from matplotlib.backend_bases import key_press_handler
# # #     from matplotlib.figure import Figure
# # #     import matplotlib.pyplot as plt
# # #     import temp_functions
# # #     import os

# # #     matplotlib.use("Agg")

# # #     master.withdraw()
# # #     master.destroy()
    
# # #     # Generate Variables
# # #     analyze = tk.Tk()
# # #     analyze.title("Analyze Images")
# # #     analyze.geometry("1600x800")

# # #     source_directory = tk.StringVar(analyze)
# # #     single_or_multiple = tk.StringVar(analyze)
# # #     single_or_multiple.set("multiple")

# # #     left_frame = tk.Frame(analyze, relief=tk.RAISED, bd=2)
# # #     figure_frame = tk.Frame(analyze, background='black', padx=2, pady=2)
# # #     right_frame = tk.Frame(analyze, relief=tk.RAISED, bd=2)


# # # ########################## LEFT FRAME ################
# # #     # Source directory
# # #     souce_text = tk.Label(left_frame, textvariable=source_directory)
# # #     souce_text.grid(row=0,column=0, sticky="ew", padx=5, pady=5)

# # #     lbox_files = tk.Listbox(left_frame, exportselection = 0)
# # #     lbox_files.grid(row=1,column=0, sticky="ew", padx=5, pady=5)

# # #     def select_source_folder():
# # #         source_directory.set(filedialog.askdirectory())

# # #         lbox_files.delete(0, tk.END)
# # #         file_list = os.listdir(source_directory.get())
# # #         for item in file_list:
# # #             lbox_files.insert(tk.END, item)
# # #         lbox_files.update()
# # #         return

# # #     tk.Button(left_frame, text = "Select Source", command=select_source_folder).grid(\
# # #         row=2,column=0, sticky="ew", padx=5, pady=5)
# # #     # tk.Button(analyze, text = "fill", command = fill_file_list_box).pack(side = tk.LEFT, padx=10)


# # #     color_keys = ['mean_green', 'mean_red', 'mean_blue', 'median_green', 'median_red', 'median_blue', 'mean_yellow',
# # #               'median_yellow']

# # #     lbox_colors = tk.Listbox(left_frame, exportselection = 0)
# # #     lbox_colors.grid(row=3,column=0, sticky="ew", padx=5, pady=5)
# # #     for color in color_keys:
# # #         lbox_colors.insert(tk.END, color)
# # #     lbox_colors.update()


# # # ############### RIGHT FRAME ###########################
# # #     operation = tk.StringVar(analyze)
# # #     operation.set("view")
# # #     tk.Label(right_frame, text = "Click Operation:").grid(row=0,column=0, sticky="ew", padx=5, pady=5)
# # #     tk.Radiobutton(right_frame, text = "Select Hits", variable=operation, value = "select").grid(row=1,column=0, sticky="ew", padx=5, pady=5)
# # #     tk.Radiobutton(right_frame, text = "View Graph", variable=operation, value = "view").grid(row=2,column=0, sticky="ew", padx=5, pady=5)
# # #     lbox_hits = tk.Listbox(right_frame, exportselection = 0)
# # #     lbox_hits.grid(row=3,column=0, sticky="ew", padx=5, pady=5)


# # #     hit_list = []

# # #     def update_hit_list():
# # #         lbox_hits.delete(0,tk.END)
# # #         for item in hit_list:
# # #             lbox_hits.insert(tk.END,item)
# # #         lbox_hits.update()
# # #         return

# # #     def add_hit(gene):
# # #         if gene not in hit_list:
# # #             hit_list.append(gene)
# # #             update_hit_list()
# # #         return
    
# # #     def remove_hit():
# # #         gene = lbox_hits.get(lbox_hits.curselection()[0])
# # #         if gene in hit_list:
# # #             hit_list.remove(gene)
# # #             update_hit_list()
# # #         return

# # #     tk.Button(right_frame, text="Remove Hit", command=remove_hit).grid(row=4,column=0,sticky="ew", padx=5, pady=5)
    
    






    
# # # ############### MIDDLE FRAME ##########################
# # #     def onclick(event):
# # #         axss = []
# # #         for a in range(8):
# # #             _a = a
# # #             for b in range(12):
# # #                 _b = b
# # #                 axss.append(axs[_a,_b])

# # #         for i, ax in enumerate(axss):
# # #             if ax == event.inaxes:
# # #                 i += 1
# # #                 a = chr(int((i-1)/12)  + 65)
# # #                 b = (i-1) % 12 + 1

# # #                 if b < 10:
# # #                     b = "0" + str(b)
# # #                 else:
# # #                     b = str(b)

# # #                 print(a,b)
# # #                 wellID = a + b

# # #                 if operation.get() == "select":
# # #                     file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
# # #                     # color = lbox_colors.get(lbox_colors.curselection()[0])
# # #                     color_info_file = file
# # #                     gene_series = temp_functions.import_gene_color_change_dict(color_info_file)
# # #                     add_hit(gene_series[wellID]["gene"])
                    

# # #                 if operation.get() == "view":
# # #                     abx = plt.figure()
# # #                     abx.add_subplot(111)

# # #                     file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
# # #                     color = lbox_colors.get(lbox_colors.curselection()[0])
# # #                     color_info_file = file
# # #                     gene_series = temp_functions.import_gene_color_change_dict(color_info_file)

# # #                     times = gene_series["A01"]['times']
# # #                     min_time = min(times)
# # #                     max_time = max(times)
# # #                     idxes = np.logical_and(min_time <= times, max_time >= times)

# # #                     print("COLORS for " + wellID)
# # #                     mean_colors = np.mean(gene_series[wellID][color], axis=0)

# # #                     plt.title(gene_series[wellID]["gene"])
# # #                     plt.ylim([0,300])
# # #                     plt.plot(times[idxes], mean_colors[idxes])
# # #                     abx.show()
# # #                 else:
# # #                     print("invalid value for variable 'operation'")

# # #         return


# # #     fig, axs = plt.subplots(8,12, sharex='all', sharey='all', gridspec_kw={'hspace': 0, 'wspace': 0})

# # #     for a in range(len(axs)):
# # #         for b in range(len(axs[a])):
# # #             axs[a,b].axes.get_yaxis().set_visible(False)
# # #             axs[a,b].axes.get_xaxis().set_visible(False)

# # #     fig.figsize=[16,10]
# # #     fig.dpi = 120
    
# # #     canvas = FigureCanvasTkAgg(fig, master=figure_frame)
# # #     fig.tight_layout()
# # #     canvas.draw()



# # #     def update_graph():
# # #         file = source_directory.get() + "/" + lbox_files.get(lbox_files.curselection()[0])
# # #         print(file)

# # #         color = lbox_colors.get(lbox_colors.curselection()[0])
# # #         color_info_file = file
# # #         gene_series = temp_functions.import_gene_color_change_dict(color_info_file)

# # #         all_genes = list(gene_series.keys())
# # #         times = gene_series["A01"]['times']
# # #         min_time = min(times)
# # #         max_time = max(times)
# # #         idxes = np.logical_and(min_time <= times, max_time >= times)

# # #         for wellID in gene_series:
# # #             a = gene_series[wellID]["well_row"] - 0
# # #             b = gene_series[wellID]["well_column"] - 1
# # #             x = a*12+b
# # #             print(a,b,x)
# # #             mean_colors = np.mean(gene_series[wellID][color], axis=0)
# # #             axs[a ,b].cla()

# # #             axs[a,b].tick_params(
# # #                 axis='x',
# # #                 bottom=False,    
# # #                 top=False)
# # #             axs[a,b].tick_params(
# # #                 axis='y',     
# # #                 left=False,    
# # #                 right=False)
# # #             # axs[a,b].set_title("demo")
# # #             if gene_series[wellID]["gene"] != "Blank":
# # #                 axs[a,b].text(0,1,gene_series[wellID]["gene"], fontsize=5, style='italic')
# # #             axs[a,b].set_ylim([0,300])
# # #             axs[a,b].axes.get_yaxis().set_visible(False)
# # #             axs[a,b].axes.get_xaxis().set_visible(False)
# # #             if gene_series[wellID]["gene"] == "Blank":
# # #                 axs[a, b].plot(times[idxes], mean_colors[idxes], color = "gray")
# # #             else:
# # #                 axs[a, b].plot(times[idxes], mean_colors[idxes], color = "blue")
# # #         canvas.draw()


# # #     canvas.get_tk_widget().grid(row=0,column=1, sticky="ew", padx=5, pady=5)
# # #     left_frame.grid(row=0,column=0,sticky="ew",padx=5,pady=5)
# # #     figure_frame.grid(row=0,column=1, sticky="ew", padx=5, pady=5)
# # #     right_frame.grid(row=0, column=2, sticky="ew", padx=5, pady=5)

# # #     tk.Text(ana)
# # #     tk.Button(analyze, text = "Update Graph", command=update_graph).grid(row=1,column=1, sticky="ew", padx=5, pady=5)

# # #     clicker = fig.canvas.mpl_connect('button_press_event', onclick)

# # #     analyze.mainloop()


# # # ###############################################

# # # # Set up Main Selection Screen
# # # master = tk.Tk()
# # # master.title("Jurgensen Well Analyzer Demo")
# # # master.geometry("300x150")

# # # tk.Label(master, text = "Select which function you want to perform").pack()
# # # tk.Button(master, text = "Generate Barcodes", command = open_generate).pack()
# # # tk.Button(master, text = "Organize by Barcodes", command = open_organize).pack()
# # # tk.Button(master, text = "Analyze Plates", command = open_analyze).pack()

# # # master.mainloop()













# # # # import tkinter as tk
# # # # from tkinter.filedialog import askopenfilename, asksaveasfilename

# # # # def open_file():
# # # #     """Open a file for editing."""
# # # #     filepath = askopenfilename(
# # # #         filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
# # # #     )
# # # #     if not filepath:
# # # #         return
# # # #     txt_edit.delete(1.0, tk.END)
# # # #     with open(filepath, "r") as input_file:
# # # #         text = input_file.read()
# # # #         txt_edit.insert(tk.END, text)
# # # #     window.title(f"Simple Text Editor - {filepath}")

# # # # def save_file():
# # # #     """Save the current file as a new file."""
# # # #     filepath = asksaveasfilename(
# # # #         defaultextension="txt",
# # # #         filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
# # # #     )
# # # #     if not filepath:
# # # #         return
# # # #     with open(filepath, "w") as output_file:
# # # #         text = txt_edit.get(1.0, tk.END)
# # # #         output_file.write(text)
# # # #     window.title(f"Simple Text Editor - {filepath}")

# # # # window = tk.Tk()
# # # # window.title("Simple Text Editor")
# # # # window.rowconfigure(0, minsize=800, weight=1)
# # # # window.columnconfigure(1, minsize=800, weight=1)

# # # # txt_edit = tk.Text(window)
# # # # fr_buttons = tk.Frame(window, relief=tk.RAISED, bd=2)
# # # # btn_open = tk.Button(fr_buttons, text="Open", command=open_file)
# # # # btn_save = tk.Button(fr_buttons, text="Save As...", command=save_file)

# # # # btn_open.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
# # # # btn_save.grid(row=1, column=0, sticky="ew", padx=5)

# # # # fr_buttons.grid(row=0, column=0, sticky="ns")
# # # # txt_edit.grid(row=0, column=1, sticky="nsew")

# # # # window.mainloop()
# # # # # # import numpy as np
# # # # # # from tkinter import filedialog
# # # # # # import tkinter as tk
# # # # # # import csv




# # # # # # # Set up the Tkinter Window
# # # # # # master = tk.Tk()
# # # # # # master.title("Jurgensen Well Analyzer Demo")
# # # # # # master.geometry("500x300")

# # # # # # # Select Source Folder
# # # # # # source_dir = "Select Source Directory"
# # # # # # def select_source_folder():
# # # # # #     source_dir = filedialog.askdirectory()
# # # # # #     return

# # # # # # # Select Destination Folder
# # # # # # dest_dir = tk.StringVar(master)
# # # # # # dest_dir.set("Select Destination Directory")
# # # # # # def select_dest_folder():
# # # # # #     dest_dir.set(filedialog.askdirectory())
# # # # # #     print(dest_dir)
# # # # # #     return

# # # # # # source_text = tk.Label(master, text = source_dir)
# # # # # # source_text.pack()
# # # # # # source_button = tk.Button(master, text = "Source Folder", command = select_source_folder)
# # # # # # source_button.pack()

# # # # # # dest_text = tk.Label(master, textvariable = dest_dir)
# # # # # # dest_text.pack()
# # # # # # dest_button = tk.Button(master, text = "Dest Folder", command = select_dest_folder)
# # # # # # dest_button.pack()


# # # # # # # Run main loop
# # # # # # master.mainloop()

# # # # # import tkinter

# # # # # from matplotlib.backends.backend_tkagg import (
# # # # #     FigureCanvasTkAgg, NavigationToolbar2Tk)
# # # # # # Implement the default Matplotlib key bindings.
# # # # # from matplotlib.backend_bases import key_press_handler
# # # # # from matplotlib.figure import Figure

# # # # # import numpy as np


# # # # # root = tkinter.Tk()
# # # # # root.wm_title("Embedding in Tk")

# # # # # fig = Figure(figsize=(5, 4), dpi=100)
# # # # # t = np.arange(0, 3, .01)
# # # # # fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))

# # # # # canvas = FigureCanvasTkAgg(fig, master=root)  # A tk.DrawingArea.
# # # # # canvas.draw()

# # # # # toolbar = NavigationToolbar2Tk(canvas, root)
# # # # # toolbar.update()


# # # # # def on_key_press(event):
# # # # #     print("you pressed {}".format(event.key))
# # # # #     key_press_handler(event, canvas, toolbar)


# # # # # canvas.mpl_connect("key_press_event", on_key_press)

# # # # # button = tkinter.Button(master=root, text="Quit", command=root.quit)

# # # # # # Packing order is important. Widgets are processed sequentially and if there
# # # # # # is no space left, because the window is too small, they are not displayed.
# # # # # # The canvas is rather flexible in its size, so we pack it last which makes
# # # # # # sure the UI controls are displayed as long as possible.
# # # # # button.pack(side=tkinter.BOTTOM)
# # # # # canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

# # # # # tkinter.mainloop()

# # # # # import tkinter as tk
# # # # # import matplotlib.pyplot as plt
# # # # # from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# # # # # import numpy as np


# # # # # root = tk.Tk()
# # # # # # frames to create the black border:
# # # # # fig_frame1 = tk.Frame(root, background='black', padx=2, pady=2)  

# # # # # fig1, ax1 = plt.subplots(figsize=(2, 2))
# # # # # t = np.arange(0, 2*np.pi, 0.1)

# # # # # ax1.plot(t, np.cos(t))
# # # # # fig1.tight_layout()
# # # # # fig1_canvas = FigureCanvasTkAgg(fig1, master=fig_frame1)  # set master of fig1_canvas to the border frame
# # # # # fig1_canvas.get_tk_widget().pack(padx=10, pady=10)  # change padx and pady to choose the thickness of the border


# # # # # frame1 = tk.Frame()

# # # # # tk.Label(frame1, text='hello').pack()

# # # # # frame1.grid(row=0, column=0, rowspan=2)
# # # # # # put the border frames in the root window
# # # # # fig_frame1.grid(row=0, column=1)

# # # # # root.mainloop()