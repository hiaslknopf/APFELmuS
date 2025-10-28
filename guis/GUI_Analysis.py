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

""" A simple GUI for showing and saving plots on the fly from Spe, MCA and root files + Calibration files. """

class GUI_Analysis:

    myappid = 'APFELMUS_Analysis' # Some arbitrary string

    if os.name == 'nt':
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        # Main window setup
        self.root = root
        self.root.title("APFELmuS - Spectrum Viewer")

        # Load the APFELmuS logo
        if os.name == 'nt':
            root.iconbitmap("ressources/logo.ico")

        # Create the widgets
        self.create_cal_widgets()
        self.create_plot_selection_widgets()
        self.create_output_widgets()
        self.create_plot_window()

        # Initialize the GUI plotting window
        self.welcome_message()

        self.last_browsed_path = os.path.dirname(os.path.abspath(__file__)) 
    
    def create_cal_widgets(self):

        ################################################################
        #-----------------File and Calibration section-----------------#
        ################################################################

        self.frame_file_cal = ttk.Frame(root, borderwidth=2, relief="solid")
        self.frame_file_cal.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        #Browse for file
        ttk.Label(self.frame_file_cal, text="File:").grid(row=0, column=0, sticky=tk.W)
        self.entry_file_path = ttk.Entry(self.frame_file_cal)
        self.entry_file_path.grid(row=0, column=1, sticky=tk.W)
        self.browse_file_path = ttk.Button(self.frame_file_cal, text="Browse", command=self.browse_file)
        self.browse_file_path.grid(row=0, column=2, sticky=tk.W)

        #Browse for calibration file
        ttk.Label(self.frame_file_cal, text="Linearization:").grid(row=1, column=0, sticky=tk.W)
        self.entry_linearization_path = ttk.Entry(self.frame_file_cal)
        self.entry_linearization_path.grid(row=1, column=1, sticky=tk.W)
        self.browse_linearization_path = ttk.Button(self.frame_file_cal, text="Browse", command=self.browse_linearization)
        self.browse_linearization_path.grid(row=1, column=2, sticky=tk.W)

        # Empty line
        ttk.Label(self.frame_file_cal, text="").grid(row=2, column=0, columnspan=3, sticky=tk.W)

        #Enter calibration factor
        ttk.Label(self.frame_file_cal, text="Calibration factor:").grid(row=3, column=0, sticky=tk.W)
        self.entry_cal_factor = ttk.Entry(self.frame_file_cal)
        self.entry_cal_factor.grid(row=3, column=1, sticky=tk.W)
        ttk.Label(self.frame_file_cal, text="keVum-1mV-1").grid(row=3, column=2, sticky=tk.W)

        # Checkbox for calibration method
        self.calibration_method_var = tk.BooleanVar(value=False)
        self.calibration_method_check = ttk.Checkbutton(self.frame_file_cal, text="Calibrate with edge + chord length",
                                                        variable=self.calibration_method_var, command=self.toggle_calib_method)
        self.calibration_method_check.grid(row=4, column=0, columnspan=3, sticky=tk.W)

        #Enter edge position (mV)
        ttk.Label(self.frame_file_cal, text="Edge position:").grid(row=5, column=0, sticky=tk.W)
        self.entry_edge_pos = ttk.Entry(self.frame_file_cal, state=tk.DISABLED)
        self.entry_edge_pos.grid(row=5, column=1, sticky=tk.W)
        ttk.Label(self.frame_file_cal, text="mV").grid(row=5, column=2, sticky=tk.W)

        #Enter the mean chord length (µm)
        ttk.Label(self.frame_file_cal, text="Mean chord length:").grid(row=6, column=0, sticky=tk.W)
        self.entry_mean_chord_length = ttk.Entry(self.frame_file_cal, state=tk.DISABLED)
        self.entry_mean_chord_length.grid(row=6, column=1, sticky=tk.W)
        ttk.Label(self.frame_file_cal, text="µm").grid(row=6, column=2, sticky=tk.W)

        # Variable to trace changes in the file entry
        self.file_var = tk.StringVar()
        self.file_var.trace_add("write", self.check_file_entry)
        self.entry_file_path["textvariable"] = self.file_var

        # Two dropdown menus for info_dict detector and material
        ttk.Label(self.frame_file_cal, text=" Detector:").grid(row=7, column=0, sticky=tk.W)
        self.detector_var = tk.StringVar()
        self.detector_var.set("silicon")
        self.detector_dropdown = ttk.OptionMenu(self.frame_file_cal, self.detector_var, "silicon", "silicon", "diamond", "sic")
        self.detector_dropdown.configure(state="disabled")
        self.detector_dropdown.grid(row=7, column=1, sticky=tk.W)

        ttk.Label(self.frame_file_cal, text=" Particle:").grid(row=8, column=0, sticky=tk.W)
        self.particle_var = tk.StringVar()
        self.particle_var.set("proton")
        self.material_dropdown = ttk.OptionMenu(self.frame_file_cal, self.particle_var, "proton", "proton", "carbon", "helium")
        self.material_dropdown.configure(state="disabled")
        self.material_dropdown.grid(row=8, column=1, sticky=tk.W)

        # Another dropdown for the database selection
        ttk.Label(self.frame_file_cal, text=" Database:").grid(row=9, column=0, sticky=tk.W)
        self.database_var = tk.StringVar()
        self.database_var.set("ICRU")
        self.database_dropdown = ttk.OptionMenu(self.frame_file_cal, self.database_var, "ICRU", "ICRU", "NIST", "SRIM")
        self.database_dropdown.configure(state="disabled")
        self.database_dropdown.grid(row=9, column=1, sticky=tk.W)

        # Empty line
        ttk.Label(self.frame_file_cal, text="").grid(row=9, column=0, columnspan=3, sticky=tk.W)

        # Option to cut off spectrum
        self.cutoff_var = tk.BooleanVar(value=False)
        self.cutoff_check = ttk.Checkbutton(self.frame_file_cal, text="Cut off spectrum:", variable=self.cutoff_var, command=self.toggle_cutoff)
        self.cutoff_check.grid(row=10, column=0, sticky=tk.W)

        # Entry for cutoff value
        self.entry_cutoff = ttk.Entry(self.frame_file_cal, state=tk.DISABLED)
        self.entry_cutoff.grid(row=10, column=1, sticky=tk.W)
        ttk.Label(self.frame_file_cal, text="channel").grid(row=10, column=2, sticky=tk.W)     

    def create_plot_selection_widgets(self):

        ################################################################
        #-----------------------Plot section---------------------------#
        ################################################################

        self.frame_plot = ttk.Frame(root, borderwidth=2, relief="solid")
        self.frame_plot.grid(row=2, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        ttk.Label(self.frame_plot, text="Select plots").grid(row=0, column=0, sticky=tk.W)
        self.plot_options = ["Original data (no log binning)", "f(y)", "F(y)", "d(y)", "yf(y)", "yd(y)"]
        self.check_vars = [tk.BooleanVar(value=False) for _ in range(len(self.plot_options))]
        self.check_vars[0].set(True)  # Default plot original spectrum

        for i, option in enumerate(self.plot_options):
            check_button = ttk.Checkbutton(self.frame_plot, text=option, variable=self.check_vars[i])
            check_button.grid(row=i + 1, column=0, padx=5, pady=5, columnspan=3, sticky=tk.W)

        ttk.Label(self.frame_plot, text="Select plots").grid(row=0, column=0, sticky=tk.W)

        # Empty line
        ttk.Label(self.frame_plot, text="").grid(row=8, column=0, columnspan=3, sticky=tk.W)

        # Log binning checkmark
        self.log_binning_var = tk.BooleanVar(value=True)
        self.log_binning_check = ttk.Checkbutton(self.frame_plot, text="Log binning", variable=self.log_binning_var, command=self.toggle_entry_num_bins)
        self.log_binning_check.grid(row=9, column=0, sticky=tk.W)

        # Number of bins entry
        self.num_bins_placeholder = 'Enter num bins'
        self.entry_num_bins = ttk.Entry(self.frame_plot, textvariable=tk.StringVar(), state=tk.NORMAL)
        self.entry_num_bins.grid(row=9, column=1, sticky=tk.W)
        self.entry_num_bins.insert(0, self.num_bins_placeholder)  # Insert placeholder text
        self.entry_num_bins.bind("<FocusIn>", self.clear_num_bins_placeholder)  # Clear placeholder on focus
        self.entry_num_bins.bind("<FocusOut>", self.restore_num_bins_placeholder)  # Restore placeholder if entry is empty

        # Set xlim
        self.xlim_var = tk.BooleanVar(value=False)
        self.xlim_check = ttk.Checkbutton(self.frame_plot, text="Set x-axis limits (keV/um)", variable=self.xlim_var, command=self.toggle_xlim)
        self.xlim_check.grid(row=10, column=0, sticky=tk.W)

        # xlim entry two parts: from to
        self.entry_xlim_from = ttk.Entry(self.frame_plot, textvariable=tk.StringVar(), state=tk.NORMAL)
        self.entry_xlim_from.grid(row=10, column=1, sticky=tk.W)
        self.entry_xlim_to = ttk.Entry(self.frame_plot, textvariable=tk.StringVar(), state=tk.NORMAL)
        self.entry_xlim_to.grid(row=10, column=2, sticky=tk.W)

        # Placeholder text for xlim entry
        self.xlim_from_placeholder = 'From'
        self.xlim_to_placeholder = 'To'

        self.entry_xlim_from.insert(0, self.xlim_from_placeholder)
        self.entry_xlim_from.bind("<FocusIn>", self.clear_xlim_placeholder)
        self.entry_xlim_from.bind("<FocusOut>", self.restore_xlim_placeholder)

        self.entry_xlim_to.insert(0, self.xlim_to_placeholder)
        self.entry_xlim_to.bind("<FocusIn>", self.clear_xlim_placeholder)
        self.entry_xlim_to.bind("<FocusOut>", self.restore_xlim_placeholder)

        # Disable xlim entry
        self.entry_xlim_from.config(state=tk.DISABLED)
        self.entry_xlim_to.config(state=tk.DISABLED)
    
        # Change to step plot
        self.step_var = tk.BooleanVar(value=False)
        self.step_check = ttk.Checkbutton(self.frame_plot, text="Step plot", variable=self.step_var)
        self.step_check.grid(row=11, column=0, sticky=tk.W)

        # Plot mean values
        self.plot_means_var = tk.BooleanVar(value=False)
        self.plot_means_check = ttk.Checkbutton(self.frame_plot, text="Plot mean values (y_F, y_D)", variable=self.plot_means_var)
        self.plot_means_check.grid(row=12, column=0, sticky=tk.W)
        
    def create_output_widgets(self):

        ################################################################
        #-----------------------Output section-------------------------#
        ################################################################

        self.frame_output = ttk.Frame(root, borderwidth=2, relief="solid")
        self.frame_output.grid(row=3, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)
        
        ttk.Label(self.frame_output, text="Save plots:").grid(row=0, column=0, sticky=tk.W)
        self.save_plots_var = tk.BooleanVar(value=False)
        self.check_save_plots = ttk.Checkbutton(self.frame_output, variable=self.save_plots_var, command=self.check_plot_out_entry).grid(row=0, column=1, sticky=tk.W)

        ttk.Label(self.frame_output, text="Output:").grid(row=1, column=0, sticky=tk.W)
        self.entry_plot_out = ttk.Entry(self.frame_output,state=tk.DISABLED)
        self.entry_plot_out.grid(row=1, column=1, sticky=tk.W)
        self.browse_out_button = ttk.Button(self.frame_output, text="Browse", command=self.browse_plot_out,state=tk.DISABLED)
        self.browse_out_button.grid(row=1, column=2, sticky=tk.W)

        ttk.Label(self.frame_output, text="Save csv:").grid(row=2, column=0, sticky=tk.W)
        self.save_csv_var = tk.BooleanVar(value=False)
        self.check_save_csv = ttk.Checkbutton(self.frame_output, variable=self.save_csv_var, command=self.check_csv_out_entry).grid(row=2, column=1, sticky=tk.W)

        ttk.Label(self.frame_output, text="Output:").grid(row=3, column=0, sticky=tk.W)
        self.entry_csv_out = ttk.Entry(self.frame_output,state=tk.DISABLED)
        self.entry_csv_out.grid(row=3, column=1, sticky=tk.W)
        self.browse_csv_out_button = ttk.Button(self.frame_output, text="Browse", command=self.browse_csv_out,state=tk.DISABLED)
        self.browse_csv_out_button.grid(row=3, column=2, sticky=tk.W)

        ################################################################
        #-----------------------Button section-------------------------#
        ################################################################

        self.button_frame = ttk.Frame(root)
        self.button_frame.grid(row=11, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.run_button = ttk.Button(self.button_frame, text="RUN", command=self.run)
        self.run_button.grid(row=0, column=0, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.close_plots_button = ttk.Button(self.button_frame, text="CLOSE PLOTS", command=self.close_plots)
        self.close_plots_button.grid(row=0, column=1, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.reset_button = ttk.Button(self.button_frame, text="RESET", command=self.reset_GUI)
        self.reset_button.grid(row=0, column=2, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

    def create_plot_window(self):

        ################################################################
        #-----------------------Plot window----------------------------#
        ################################################################

        # Frame for right side
        self.frame_plot_window = ttk.Frame(root)
        self.frame_plot_window.grid(row=0, column=1, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.frame_plot_window)
        self.notebook.grid(row=0, column=0, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a Figure and Axes for the plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot_window)
        self.canvas_widget = self.canvas.get_tk_widget()
 
        # Initialize the check_vars and plot_tabs attributes
        self.plot_objects = {}
        self.plot_tabs = {}

        # Empty row
        ttk.Label(self.frame_plot_window, text="").grid(row=1, column=0, sticky=tk.W)

        # Create textbox widget below the notebook with scrollbar for longer texts
        self.info_textbox = tk.Text(self.frame_plot_window, width=50, height=10, wrap=tk.WORD, bg="white")
        scrollbar = ttk.Scrollbar(self.frame_plot_window, orient=tk.VERTICAL, command=self.info_textbox.yview)
        self.info_textbox.config(yscrollcommand=scrollbar.set)
        scrollbar.grid(row=20, column=1, padx=0, pady=10, sticky=tk.W+tk.E+tk.N+tk.S, rowspan=12)
        self.info_textbox.grid(row=20, column=0, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.initial_text = "Welcome to the APFELmuS Spectrum Viewer\nPlease select your options on the left ..."
        self.info_textbox.insert(1.0, self.initial_text)

################################################################
#----------------------GUI FUNCTIONALITY-----------------------#
################################################################

    def welcome_message(self):
        """ Greeting message when the GUI is started or reset """

        tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
        fig, ax, canvas = self.setup_plot('Welcome', tab_text)
        canvas.draw()

        self.plot_objects['Welcome'] = (fig, ax, canvas)

        welcome_message(ax, "Analysis Toolkit")
    
    def plot_csv_data(self, campaign):

        # Set the plot window
        self.close_plots()
        tab_text_str = campaign.measurements[self.meas_name].y_axis
        tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
        fig, ax, canvas = self.setup_plot(tab_text_str, tab_text)

        # Set xlims if specified
        if self.xlim_var.get():
            xlim_from = float(self.entry_xlim_from.get())
            xlim_to = float(self.entry_xlim_to.get())
            ax.set_xlim(xlim_from, xlim_to)
        else:
            xlim_from = None
            xlim_to = None

        # Plot the data
        Output.plot_single(campaign.measurements[self.meas_name], name=f'{self.meas_name}',
                                   step=self.step_var.get(), y_F=False, y_D=False, xlim=[xlim_from, xlim_to],
                                   show_plot=False, gui=True, fig=fig, ax=ax)
             
        toolbar = NavigationToolbar2Tk(canvas, tab_text, pack_toolbar=False)
        toolbar.grid(row=1, column=0, padx=10, pady=10, sticky=tk.S)
        canvas.draw()

        # Store the plot objects in the dictionary (for closing the plots later)
        self.plot_objects[tab_text_str] = (fig, ax, canvas)
            
        # Remove the welcome message tab
        if 'Welcome' in self.plot_tabs:
            tab_text = self.plot_tabs.pop('Welcome')
            index = self.notebook.index(tab_text)
            self.notebook.forget(index)
            self.plot_objects.pop('Welcome')

    def reset_GUI(self):
        """ Reset the GUI to an initial state """

        # Reset all checkmarks
        for var in self.check_vars:
            var.set(False)
        self.log_binning_var.set(True)
        self.save_plots_var.set(False)
        self.cutoff_var.set(False)

        # Reset file and calibration entries and button
        self.calibration_method_check.state(['!selected'])
        self.entry_file_path.delete(0, tk.END)
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_edge_pos.delete(0, tk.END)
        self.entry_mean_chord_length.delete(0, tk.END)
        self.entry_cal_factor.delete(0, tk.END)
        self.entry_num_bins.delete(0, tk.END)
        self.entry_cutoff.delete(0, tk.END)

        # Reset dropdown
        self.detector_var.set("silicon")
        self.particle_var.set("proton")
        self.database_var.set("ICRU")
        
        # Lock entries
        self.entry_cal_factor.config(state=tk.NORMAL)
        self.entry_linearization_path.config(state=tk.NORMAL)
        self.browse_linearization_path.config(state=tk.NORMAL)
        self.entry_edge_pos.config(state=tk.DISABLED)
        self.entry_mean_chord_length.config(state=tk.DISABLED)
        self.entry_num_bins.config(state=tk.NORMAL)
        self.entry_plot_out.config(state=tk.DISABLED)
        self.entry_csv_out.config(state=tk.DISABLED)
        self.browse_out_button.config(state=tk.DISABLED)
        self.entry_cutoff.config(state=tk.DISABLED)

        # Select log binning
        self.log_binning_check.state(['!selected'])
        self.log_binning_check.config(state=tk.NORMAL)
        self.entry_num_bins.config(state=tk.NORMAL)
        self.log_binning_var.set(True)

        # Deselect step plot
        self.step_var.set(False)

        # Reset text box
        self.info_textbox.delete(1.0, tk.END)
        self.info_textbox.insert(tk.END, self.initial_text)

        # Deselect plot and csv output
        self.entry_plot_out.delete(0, tk.END)
        self.entry_csv_out.delete(0, tk.END)
        self.save_plots_var.set(False)
        self.check_plot_out_entry()
        self.save_csv_var.set(False)
        self.check_csv_out_entry()

        # Reset the plot tabs
        for tab_name in list(self.plot_tabs.keys()):
            tab_text = self.plot_tabs.pop(tab_name)
            index = self.notebook.index(tab_text)
            self.notebook.forget(index)
        
        # Activate step plot
        self.step_check.config(state=tk.NORMAL)
        self.step_var.set(False)

        # xlim values
        self.xlim_check.config(state=tk.NORMAL)
        self.xlim_var.set(False)
        self.entry_xlim_from.config(state=tk.DISABLED)
        self.entry_xlim_to.config(state=tk.DISABLED)
        self.entry_xlim_from.delete(0, tk.END)
        self.entry_xlim_to.delete(0, tk.END)
        self.entry_xlim_from.insert(0, self.xlim_from_placeholder)
        self.entry_xlim_to.insert(0, self.xlim_to_placeholder)

        # Mean values
        self.plot_means_check.config(state=tk.NORMAL)
        self.plot_means_var.set(False)

        # Reset the plot objects
        self.plot_objects = {}
        self.plot_tabs = {}

        # Reset the  run button
        self.run_button.config(state=tk.NORMAL)

        # Reset the plot canvas
        plt.close('all')
        self.welcome_message()

    def check_csv_out_entry(self, *args):
        """ Lock csv output fields if no output is selected """

        if self.save_csv_var.get():
            self.entry_csv_out.config(state=tk.NORMAL)
            self.browse_csv_out_button.config(state=tk.NORMAL)
        else:
            self.entry_csv_out.config(state=tk.DISABLED)
            self.browse_csv_out_button.config(state=tk.DISABLED)

    def check_plot_out_entry(self, *args):
        """ Lock plot output fields if no output is selected """

        if self.save_plots_var.get():
            self.entry_plot_out.config(state=tk.NORMAL)
            self.browse_out_button.config(state=tk.NORMAL)
        else:
            self.entry_plot_out.config(state=tk.DISABLED)
            self.browse_out_button.config(state=tk.DISABLED)

    def check_file_entry(self, *args):
        """ Lock input options for different types of files
        Currently supports .MCA, .Spe, .root and analyzed .csv files """

        file_content = self.file_var.get()

        # Option 1: Linearized .MCA file (APFELmuS)
        if ".MCA" in file_content.upper():
            # Lock linearization entry
            self.entry_linearization_path.config(state=tk.DISABLED)
            self.browse_linearization_path.config(state=tk.DISABLED)

            # Lock detector and particle dropdowns
            self.detector_dropdown.configure(state="disabled")
            self.material_dropdown.configure(state="disabled")

        # Option 2: Raw .Spe file (MAESTRO)
        elif ".Spe" in file_content.lower():
            # Unlock linearization entry
            self.entry_linearization_path.config(state=tk.NORMAL)
            self.browse_linearization_path.config(state=tk.NORMAL)

        # Option 3: .root file (Geant4 or GATE simulation)
        elif ".root" in file_content.lower():
            # Disable linearization and calibration entries
            self.entry_linearization_path.config(state=tk.DISABLED)
            self.browse_linearization_path.config(state=tk.DISABLED)
            self.entry_cal_factor.config(state=tk.DISABLED)
            # Activate chord_length entry
            self.entry_mean_chord_length.config(state=tk.NORMAL)
        
        # Option 4: Analyzed .csv file (APFELmuS)
        elif ".csv" in file_content.lower():
            # Check if this is really microdosimetric data
            with open(file_content, 'r') as f:
                lines = f.readlines()
                if not 'APFELmuS' in lines[0]:
                    messagebox.showinfo("Error", "This is not an analyzed APFELmuS microdosimetric data file")
                    raise ValueError("This is not an analyzed APFELmuS microdosimetric data file")
                else:
                    # Create Measurement object
                    campaign = MicroDosimetry()
                    campaign.read_file(file_content)
                    self.meas_name = file_content.split(os.path.sep)[-1].split('/')[-1].split('.csv')[0]

                    # Lock all fields
                    self.setup_csv_plot(campaign)
                    
                    # Plot the data right away
                    self.plot_csv_data(campaign)

    def setup_csv_plot(self, campaign):
        """ Prepare the GUI for plotting already analyzed .csv data"""

        # Lock all linearization and calibration fields
        self.entry_linearization_path.config(state=tk.DISABLED)
        self.browse_linearization_path.config(state=tk.DISABLED)
        self.entry_cal_factor.config(state=tk.DISABLED)
        self.calibration_method_check.state(['!selected'])
        self.entry_edge_pos.config(state=tk.DISABLED)
        self.entry_mean_chord_length.config(state=tk.DISABLED)
        self.cutoff_check.state(['!selected'])
        self.entry_cutoff.config(state=tk.DISABLED)

        # Set the detector and particle dropdowns to the values in the file
        self.detector_var.set(campaign.measurements[self.meas_name].detector)
        self.particle_var.set(campaign.measurements[self.meas_name].particle)

        # Lock all plot options
        self.log_binning_check.state(['!selected'])
        self.log_binning_check.config(state=tk.DISABLED)
        self.entry_num_bins.config(state=tk.DISABLED)
        self.step_check.state(['!selected'])
        self.step_check.config(state=tk.NORMAL)
        self.plot_means_check.state(['!selected'])
        self.plot_means_check.config(state=tk.DISABLED)

        # Lock the run button -> Plots the data automatically
        self.run_button.config(state=tk.DISABLED)

        # Lock the csv and plot output fields -> Unnecessary
        self.save_csv_var.set(False)
        self.check_csv_out_entry()
        self.save_plots_var.set(False)
        self.check_plot_out_entry()        

        # Select the plot option specified in the y_axis
        # Disable all other options
        for i, option in enumerate(self.plot_options):
            if campaign.measurements[self.meas_name].x_axis in ['ENERGY', 'CHANNEL']:
                self.check_vars[0].set(True)
            elif option in campaign.measurements[self.meas_name].y_axis:
                self.check_vars[i].set(True)
            else:
                self.check_vars[i].set(False)
    def toggle_cutoff(self):
        """ Check if a channel cutoff should be applied """

        if self.cutoff_var.get():
            self.entry_cutoff.config(state=tk.NORMAL)
        else:
            self.entry_cutoff.config(state=tk.DISABLED)

    def toggle_entry_num_bins(self):
        """ Check if the number of bins entry should be enabled """

        if self.log_binning_var.get():
            self.entry_num_bins.config(state=tk.NORMAL)
        else:
            self.entry_num_bins.config(state=tk.DISABLED)
    
    def toggle_calib_method(self):
        """ Check if the calibration method should be applied """

        if self.calibration_method_var.get():
            self.entry_edge_pos.config(state=tk.NORMAL)
            self.entry_mean_chord_length.config(state=tk.NORMAL)
            self.entry_cal_factor.config(state=tk.DISABLED)

            # Check if file is .MCA
            if ".MCA" in self.entry_file_path.get().upper():
                self.detector_dropdown.configure(state="disabled")
                self.material_dropdown.configure(state="disabled")
                self.database_dropdown.configure(state="enabled")
            else:
                self.detector_dropdown.configure(state="enabled")
                self.material_dropdown.configure(state="enabled")
                self.database_dropdown.configure(state="enabled")

        else:
            self.entry_edge_pos.config(state=tk.DISABLED)
            self.entry_mean_chord_length.config(state=tk.DISABLED)
            self.entry_cal_factor.config(state=tk.NORMAL)
            self.detector_dropdown.configure(state="disabled")
            self.material_dropdown.configure(state="disabled")
            self.database_dropdown.configure(state="disabled")

    def toggle_xlim(self):
        """ Check if the xlim entry should be enabled """

        if self.xlim_var.get():
            self.entry_xlim_from.config(state=tk.NORMAL)
            self.entry_xlim_to.config(state=tk.NORMAL)
        else:
            self.entry_xlim_from.config(state=tk.DISABLED)
            self.entry_xlim_to.config(state=tk.DISABLED)
    
    def clear_xlim_placeholder(self, event):
        """ Clear the placeholder text when the field is clicked """

        if self.entry_xlim_from.get() == self.xlim_from_placeholder:
            self.entry_xlim_from.delete(0, tk.END)
        if self.entry_xlim_to.get() == self.xlim_to_placeholder:
            self.entry_xlim_to.delete(0, tk.END)
    
    def restore_xlim_placeholder(self, event):
        """ Restore the placeholder text if the entry is empty """

        if not self.entry_xlim_from.get():
            self.entry_xlim_from.insert(0, self.xlim_from_placeholder)
        if not self.entry_xlim_to.get():
            self.entry_xlim_to.insert(0, self.xlim_to_placeholder)

    def clear_num_bins_placeholder(self, event):
        """ Clear the placeholder text when the field is clicked """

        if self.entry_num_bins.get() == self.num_bins_placeholder:
            self.entry_num_bins.delete(0, tk.END)

    def restore_num_bins_placeholder(self, event):
        """ Restore the placeholder text if the entry is empty """

        if not self.entry_num_bins.get():
            self.entry_num_bins.insert(0, self.num_bins_placeholder)

    def browse_file(self):
        """ Browse for the file to be analyzed """    

        file_path = filedialog.askopenfilename(initialdir=self.last_browsed_path,
                                               filetypes=[("All files", "*.*"),
                                                          ("MCA files", "*.MCA"),
                                                          ("MAESTRO files", "*.Spe"),
                                                          ("ROOT simulation data", "*.root"),
                                                          ("Analyzed Microdosimetry data", "*.csv")])
        
        self.last_browsed_path = os.path.dirname(file_path)
        
        # Empty previous entries
        self.entry_file_path.delete(0, tk.END)
        self.entry_file_path.insert(0, file_path)

    def browse_linearization(self):
        """ Browse for the linearization file """

        cal_path = filedialog.askopenfilename(initialdir=self.last_browsed_path,
                                              filetypes=[("Linearization csv files", "*.csv"), ("All files", "*.*")])
        
        self.last_browsed_path = os.path.dirname(cal_path)

        # Check if this is an APFELmuS linearization file
        with open(cal_path, 'r') as f:
            lines = f.readlines()
            if not 'CHANNEL' in lines[0]:
                messagebox.showinfo("Error", "This is not an APFELmuS linearization file")
                raise ValueError("This is not an APFELmuS linearization file")

        # Empty previous entries
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_linearization_path.insert(0, cal_path)

    def browse_plot_out(self):
        """ Browse for the output directory of the plot """

        output_path = filedialog.askdirectory(initialdir=self.last_browsed_path)
        self.last_browsed_path = output_path

        # Empty previous entries
        self.entry_plot_out.delete(0, tk.END)
        self.entry_plot_out.insert(0, output_path)
    
    def browse_csv_out(self):
        """ Browse for the output directory of the csv file """

        output_path = filedialog.askdirectory(initialdir=self.last_browsed_path)
        self.last_browsed_path = output_path

        # Empty previous entries
        self.entry_csv_out.delete(0, tk.END)
        self.entry_csv_out.insert(0, output_path)

    def close_plots(self):
        """ Close all open plots and open the welcome message """

        for tab_name in list(self.plot_tabs.keys()):
            tab_text = self.plot_tabs.pop(tab_name)
            index = self.notebook.index(tab_text)
            self.notebook.forget(index)
        
        # Reset RUN button
        self.run_button.config(state=tk.NORMAL)
        
        self.plot_objects = {}
        self.plot_tabs = {}

        plt.close('all')

        self.welcome_message()

        # Change textbox
        message = "Closed plots. Ready to start again ..."
        self.info_textbox.delete(1.0, tk.END)
        self.info_textbox.insert(1.0, message)

#################################################################
#----------------------PLOTTING FUNCTIONALITY-------------------#
#################################################################

    def calibrate_specs(self, campaign):
        """ Calibrate spectra according to the selected method and parameters """

        # Option ROOT: Attach a chord length and transform to a lineal energy axis
        if self.extension == '.root' and campaign.measurements[self.meas_name].x_axis == 'ENERGY':
            Calibrate.get_chord_length(campaign.measurements[self.meas_name], 'slab', float(self.entry_mean_chord_length.get()), plot=False)
            Calibrate.lineal_energy_axis(campaign.measurements[self.meas_name], chord_length='mean')
        
        # Option MCA: Calibrate with a calibration factor or edge + chord length
        if self.extension == '.MCA':

            if self.calibration_method_var.get():
                if self.entry_mean_chord_length.get() == "":
                    raise ValueError("Mean chord length has to be specified")
                if self.entry_edge_pos.get() == "":
                    raise ValueError("Calibration edge has to be specified")
                
                database = self.database_var.get()
                # Check if this combination is possible and set value to SRIM (exists for all detectors)
                if database == 'NIST' and self.detector_var.get() == 'sic':
                    messagebox.showinfo("Selected atabase does not exists", "NIST database not available for SiC detectors - Set to SRIM calculation")
                    database = 'SRIM'
                    self.database_var.set('SRIM')
                if database == 'ICRU' and self.detector_var.get() == 'sic':
                    messagebox.showinfo("Selected atabase does not exists", "ICRU database not available for SiC detectors - Set to SRIM calculation")
                    database = 'SRIM'
                    self.database_var.set('SRIM')
                if database == 'NIST' and self.particle_var.get() == 'carbon':
                    messagebox.showinfo("Selected atabase does not exists", "NIST database is only available for protons and helium particles - Set to SRIM calculation")
                    database = 'SRIM'
                    self.database_var.set('SRIM')
                
                # Calibration with edge + chord length
                Calibrate.get_chord_length(campaign.measurements[self.meas_name], 'slab', float(self.entry_mean_chord_length.get()), plot=False) 
                ymax, _, _ = Calibrate.get_stopping_power(campaign.measurements[self.meas_name], chord_length='mean', database=database, precision=0.05, plot=False)
                Calibrate.scale_energy_axis(campaign.measurements[self.meas_name], float(self.entry_edge_pos.get()), ymax, energy_axis='lineal', chord_length='mean')
            
            else:
                # Calibration not specified -> Default calibration factor
                if self.entry_cal_factor.get() == "":
                    cal_factor = 3
                    messagebox.showinfo("Calibration factor", f"Assumed calibration factor: {cal_factor}")
                else:
                    cal_factor = float(self.entry_cal_factor.get())
                
                Calibrate.scale_energy_axis_with_factor(campaign.measurements[self.meas_name], cal_factor)
    
    def calculate_means(self, campaign):
        """ Calculate the mean values for the selected representations """

        # Spectra are allready calibrated at this point
        if self.plot_means_var.get():
            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
            Spectrum.probability_density(campaign.measurements[self.meas_name])
            Spectrum.normalize_linear_spectrum(campaign.measurements[self.meas_name])

            self.y_F, self.y_D = Output.means_from_fy(campaign.measurements[self.meas_name])

            # Return the proper spectrum for further plotting
            Spectrum.retrieve_original_spectrum(campaign.measurements[self.meas_name])

            # Shit solution but you need to cutoff and calibrate again to calculate future spectra
            if self.cutoff_var.get():
                    Spectrum.cutoff(campaign.measurements[self.meas_name], channels=int(self.entry_cutoff.get()))
  
            self.calibrate_specs(campaign)

    def setup_plot(self, tab_name, tab_text):
        """ Create a new tab for the plot and set up the canvas """

        # Just a white screen for a smoother transition between plots
        self.canvas.draw()

        # Add the new tab to the notebook
        self.notebook.add(tab_text, text=tab_name)
        self.notebook.select(tab_text)
        self.plot_tabs[tab_name] = tab_text

        #Open up a new canvas
        fig, ax = plt.subplots()
        canvas = FigureCanvasTkAgg(fig, master=tab_text)
        canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Include space for the toolbar
        #fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        fig.subplots_adjust(bottom=0.1)     

        return fig, ax, canvas

    def show_plot(self, index, option, campaign):
        """ Plot one selected option and show it in a new tab """

        selected_option = self.check_vars[index].get()

        if selected_option:
            # Create a new tab for each plot
            tab_name = f"{option}"
            if tab_name not in self.plot_tabs:
                tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
                tab_text.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

                print('\n---------------------------------------------------')
                print('Plotting index ' + str(index) + ' ' + str(option) + '\n')
                
                # Start fresh every time with the original spectrum - Shit performance but less errors
                if self.extension == '.csv':
                    print('Just print csv')
                else:
                    Spectrum.retrieve_original_spectrum(campaign.measurements[self.meas_name])

                    # Cut off spectrum (if selected)
                    if self.cutoff_var.get():
                        if self.entry_cutoff.get() == "Enter cutoff value":
                            raise ValueError("Cutoff value has to be specified")
                        Spectrum.cutoff(campaign.measurements[self.meas_name], channels=int(self.entry_cutoff.get()))

                    # Calibrate spectra (But the original ADC channels)
                    if index != 0:
                        self.calibrate_specs(campaign)

                    # Calculate mean values if not deselected
                    if self.plot_means_var.get() and index !=0:
                        self.calculate_means(campaign)

                    # Calculate selected representations
                    #                     
                    # Logarithmic binning (if selected)
                    if self.log_binning_var.get() and index != 0:
                        if self.entry_num_bins.get() == self.num_bins_placeholder: # Default number of bins=50
                            messagebox.showinfo("Number of bins", "Assumed 50 log bins")
                            num_bins = 50
                        else:
                            num_bins = int(self.entry_num_bins.get())
                        
                        if index == 1: #f(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                            Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                        if index == 2: #F(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                            Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                        if index == 3: #d(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                            Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                        if index == 4: #yf(y) 
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                            Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                            Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                        if index == 5: #yd(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                            Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                            Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                    
                        # Normalize spectra to area=1 (except for the original and the F(y) spectrum)
                        if not index in [0,2]:
                            Spectrum.normalize_log_spectrum(campaign.measurements[self.meas_name])
                    
                    else:

                        if index == 1: #f(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                        if index == 2: #F(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                        if index == 3: #d(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                        if index == 4: #yf(y) 
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                            Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                        if index == 5: #yd(y)
                            Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                            Spectrum.probability_density(campaign.measurements[self.meas_name])
                            Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                    
                        # Normalize spectra to area=1 (except for the original and the F(y) spectrum)
                        if not index in [0,2]:
                            Spectrum.normalize_linear_spectrum(campaign.measurements[self.meas_name])
                    
                    # Setup plot canvas
                    fig, ax, canvas = self.setup_plot(tab_name, tab_text)

                    # Save the plot and csv output (if selected)
                    if self.save_plots_var.get():
                        output_path = self.entry_plot_out.get()
                        if output_path == "":
                            raise ValueError("Output path has to be specified")
                    else:
                        output_path = False
                    
                    if self.save_csv_var.get():
                        output_csv = self.entry_csv_out.get()
                        Output.csv_output(campaign.measurements[self.meas_name], output_csv, name=f'{self.meas_name}_{self.plot_options[index]}')
                        if output_csv == "":
                            raise ValueError("Output path has to be specified")
                    
                    # Set the x-axis limits (if selected)
                    if index == 0:
                        scale = 'lin'
                        xlim = [0, campaign.measurements[self.meas_name].num_channels]
                    else:
                        scale = 'log'
                        if self.xlim_var.get():
                            xlim = [float(self.entry_xlim_from.get()), float(self.entry_xlim_to.get())]
                        else:
                            xlim = [0.1, 1000]

                    # Display the mean values in the plot (if selected)
                    if self.plot_means_var.get() and index != 0:
                        # Plot the mean values
                        y_F = self.y_F
                        y_D = self.y_D
                    
                    else:
                        y_F = False
                        y_D = False

                    # Finally, plot the data
                    Output.plot_single(campaign.measurements[self.meas_name], name=f'{self.meas_name}_{self.plot_options[index]}', output_path=output_path,
                                    step=self.step_var.get(), y_F=y_F, y_D=y_D, scale=scale, xlim=xlim, 
                                    show_plot=False, gui=True, fig=fig, ax=ax)
                
                    toolbar = NavigationToolbar2Tk(canvas, tab_text, pack_toolbar=False)
                    toolbar.grid(row=1, column=0, padx=10, pady=10, sticky=tk.S)
                    canvas.draw()

                    # Store the new plot objects in the dictionary
                    self.plot_objects[tab_name] = (fig, ax, canvas)
        else:
            # If a checkmark is unselected between pressing RUN, remove the tab from the notebook
            tab_name = f"{option}"
            if tab_name in self.plot_tabs:
                tab_text = self.plot_tabs.pop(tab_name)
                index = self.notebook.index(tab_text)
                self.notebook.forget(index)
                self.plot_objects.pop(tab_name)
            
            # Remove the welcome message (if it is still there)
            if 'Welcome' in self.plot_tabs:
                tab_text = self.plot_tabs.pop('Welcome')
                index = self.notebook.index(tab_text)
                self.notebook.forget(index)
                self.plot_objects.pop('Welcome')

    def run(self):
        """ Run the analysis with the selected options and plot the data """

        # Close all open plots
        self.close_plots()

        #Create Microdosimetry object and read in data
        campaign = MicroDosimetry()
        self.datafile = self.entry_file_path.get()
        self.filename, self.extension = os.path.splitext(os.path.normpath(self.datafile))
        try:
            self.meas_name = self.filename.split(os.path.sep)[-1].split(self.extension)[0]
        except:
            messagebox.showinfo("Error", "No file selected or format corrupted")
            raise ValueError("No file selected or format corrupted")

        # Write message to the textbox
        self.info_textbox.delete(1.0, tk.END)
        message_run = "Calculating the selected plots. Please wait ...\n "
        self.info_textbox.insert(1.0, message_run)
        self.info_textbox.update()

        # Option .Spe: Translate MAESTRO file to MCA file and add metainformation
        if self.extension == '.Spe':
            if self.entry_linearization_path.get() == "":
                raise ValueError("Linearization file has to be specified")
                
            info_dict = {'DETECTOR': self.detector_var.get(), 'PARTICLE': self.particle_var.get(), 'GAIN': 'None'}

            # Translate the MAESTRO file to MCA -> Stored in a folder 'temp' in the APFELmuS repo folder
            FileTranslator.translate_MAESTRO_file(input_MAESTRO=self.datafile, input_linearization=self.entry_linearization_path.get(), output_path='temp', name=self.meas_name, info_dict=info_dict)
            self.extension = '.MCA' #Change the extension to .MCA (for GUI purposes)
            campaign.read_file(f'temp/{self.meas_name}.MCA') #Read in the translated MCA file

        elif self.extension == '.root':
            if self.entry_mean_chord_length.get() == "":
                raise ValueError("Mean chord length has to be specified")
            campaign.read_file(self.datafile)
        
        elif self.extension == '.MCA':
            campaign.read_file(self.datafile)

            # get particle and material and set dropdowns
            self.detector_var.set(campaign.measurements[self.meas_name].detector)
            self.particle_var.set(campaign.measurements[self.meas_name].particle)
        
        elif self.extension == '.csv':
            campaign.read_file(self.datafile)
            self.setup_csv_plot(campaign)
            self.plot_csv_data(campaign)

        if not self.extension == '.csv':
            # Check selected plot options
            selected_plots = []
            for i in range(len(self.plot_options)):
                selected_plots.append(self.check_vars[i].get())
            if not any(selected_plots):
                messagebox.showinfo("No plot selected", "No plot selected")
            else:
                #Iterate through the selected plot options
                for i, option in enumerate(self.plot_options):
                    self.show_plot(i, option, campaign)
            
            # Write the output of the analysis to the textbox
            self.write_output()
        
    def write_output(self):
        """ Write the output of the analysis to the textbox """

        l1 = f"Analysis completed: {self.meas_name}"
        if self.plot_means_var.get():
            l2 = f"Mean values y_F = {self.y_F:.3f} keV/um, y_D = {self.y_D:.3f} keV/um"
        else:
            l2 = ""
        l3 = f"Selected plots: {', '.join([option for i, option in enumerate(self.plot_options) if self.check_vars[i].get()])}\n"
        l4 = "Analysis choices:"
        l5 = f"Log binning: {self.log_binning_var.get()}"
        l6 = f"Step plot: {self.step_var.get()}"
        if self.cutoff_var.get():
            l7 = f"Cutoff value: {self.entry_cutoff.get()}"
        else:
            l7 = "No cutoff"
        if self.save_plots_var.get():
            l8 = f"Save plots: {self.entry_plot_out.get()}"
        else:  
            l8 = "No plot output"
        if self.save_csv_var.get():
            l9 = f"Save csv: {self.entry_csv_out.get()}"
        else:
            l9 = "No csv output"
        
        self.info_textbox.delete(1.0, tk.END)
        self.info_textbox.insert(1.0, f"{l1}\n{l2}\n{l3}\n{l4}\n{l5}\n{l6}\n{l7}\n{l8}\n{l9}")
        
root = tk.Tk()
app = GUI_Analysis(root)

# Fixed Window size
#root.resizable(False, False)
root.mainloop()