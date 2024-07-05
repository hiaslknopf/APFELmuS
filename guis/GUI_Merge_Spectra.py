__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import tkinter as tk
from tkinter import ttk, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter.messagebox as messagebox

import matplotlib.pyplot as plt
import os
import sys

import ctypes

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output, FileTranslator
from guis.welcome_message import welcome_message

""" A simple GUI to merge 2 or 3 spectra of different gains.
    The input data is already analyzed .csv files. Its too complicated to include the whole analysis in this little GUI.

    The merged spectrum can also be saved as an analyzed .csv file.
 """

class GUI_Merge:

    myappid = 'APFELmuS_Merge' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        # Main window setup
        self.root = root
        self.root.title("APFELmuS - Merge Spectra")

        # Load the APFELmuS logo
        root.iconbitmap("ressources/logo.ico")

        # Create the widgets
        self.create_load_data_widgets()
        self.create_set_parameters_widgets()
        self.create_info_widgets()
        self.create_plot_widgets()
        self.create_buttons()
    
    ################################################################
    #----------------------GUI SETUP-------------------------------#
    ################################################################

    def create_load_data_widgets(self):
        pass

    def create_set_parameters_widgets(self):
        pass

    def create_info_widgets(self):
        pass

    def create_plot_widgets(self):
        pass

    def create_buttons(self):
        pass

    ################################################################
    #----------------------GUI FUNCTIONALITY-----------------------#
    ################################################################

    def welcome_message(self):

        tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
        fig, ax, canvas = self.setup_plot('Welcome', tab_text)
        canvas.draw()

        self.plot_objects['Welcome'] = (fig, ax, canvas)

        welcome_message(ax, "Merge Spectra Tool")
    
    def reset_GUI(self):
        pass

    def browse_input_files(self):
        pass

    def print_input_info(self):
        pass

    def browse_output_file(self):
        pass

    def save_output_file(self):
        pass

    def slider_change(self, event):
        pass

    def run_merge_spectra(self):
        pass

    def print_output_info(self):
        pass

    ################################################################
    #----------------------PLOTTING-------------------------------#
    ################################################################

    def setup_plot(self, tab_name, tab_text):
        pass

    def show_plots(self):
        pass

    def close_plots(self):
        pass

root = tk.Tk()
app = GUI_Merge(root)

root.mainloop()