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
from lib import Spectrum, Output
from guis.welcome_message import welcome_message

""" A simple GUI to merge 2 or 3 spectra of different gains.
    The input data is already analyzed .csv files. Its too complicated to include the whole analysis in this little GUI.

    The merged spectrum can also be saved as an analyzed .csv file.
 """
### Im großen und ganzen funktionierts, muss aber noch bissl poliert werden
# Deppensicher machen: Ausnahmen etc
# Gscheide guesses for initial Parameter
# Plots schauen Kacke aus
# Iwie control panel für binning und normalization
# Gscheider Text output während rechnen und nachher usw

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
        self.create_plot_widgets()
        self.create_buttons()

        # Initialize the GUI
        self.welcome_message()

        self.last_browsed_path = os.path.dirname(__file__)
    
    ################################################################
    #----------------------GUI SETUP-------------------------------#
    ################################################################

    def create_load_data_widgets(self):

        self.frame_left_side = ttk.Frame(self.root)
        self.frame_left_side.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')
        
        self.frame_load_data = ttk.Frame(self.frame_left_side, borderwidth=2, relief='solid')
        self.frame_load_data.grid(row=0, column=0, padx=10, pady=10, sticky='new')

        # Browse for files
        self.frame_input_files = ttk.Frame(self.frame_load_data)
        self.frame_input_files.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        self.label_file_path = ttk.Label(self.frame_input_files, text="Select the input files:")
        self.label_file_path.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        self.entry_file_path = ttk.Entry(self.frame_input_files, width=30)  
        self.entry_file_path.grid(row=0, column=1, padx=10, pady=10, sticky='w')
        self.browse_file_path = ttk.Button(self.frame_input_files, text="Browse", command=self.browse_input_files)
        self.browse_file_path.grid(row=0, column=2, padx=10, pady=10, sticky='w')

        # Textbox showing information about the loaded spectra
        self.frame_textbox = ttk.Frame(self.frame_load_data)
        self.frame_textbox.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.text_input_info = tk.Text(self.frame_textbox, width=60, height=10, wrap=tk.WORD, bg="white", state=tk.NORMAL)
        self.text_input_info.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        self.initial_text = "Load input files to merge spectra."
        self.text_input_info.insert(tk.END, self.initial_text)

    def create_set_parameters_widgets(self):

        self.frame_set_parameters = ttk.Frame(self.frame_left_side, borderwidth=2, relief='solid')
        self.frame_set_parameters.grid(row=1, column=0, padx=10, pady=10, sticky='ew')
        
        # Set overlap region 1
        self.label_overlap1 = ttk.Label(self.frame_set_parameters, text="Overlap Region 1:", state=tk.DISABLED)
        self.label_overlap1.grid(row=1, column=0, padx=10, pady=10, sticky='w')
        self.entry_overlap1_from = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_overlap1_from.grid(row=1, column=1, padx=10, pady=10, sticky='w')
        self.label_overlap1_to = ttk.Label(self.frame_set_parameters, text="to", state=tk.DISABLED)
        self.label_overlap1_to.grid(row=1, column=2, padx=10, pady=10, sticky='w')
        self.entry_overlap1_to = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_overlap1_to.grid(row=1, column=3, padx=10, pady=10, sticky='w')
        self.overlap_1_unit = ttk.Label(self.frame_set_parameters, text="keV/µm", state=tk.DISABLED)
        self.overlap_1_unit.grid(row=1, column=4, padx=10, pady=10, sticky='w')

        # Set overlap region 2
        self.label_overlap2 = ttk.Label(self.frame_set_parameters, text="Overlap Region 2:", state=tk.DISABLED)
        self.label_overlap2.grid(row=3, column=0, padx=10, pady=10, sticky='w')
        self.entry_overlap2_from = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_overlap2_from.grid(row=3, column=1, padx=10, pady=10, sticky='w')
        self.label_overlap2_to = ttk.Label(self.frame_set_parameters, text="to", state=tk.DISABLED)
        self.label_overlap2_to.grid(row=3, column=2, padx=10, pady=10, sticky='w')
        self.entry_overlap2_to = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_overlap2_to.grid(row=3, column=3, padx=10, pady=10, sticky='w')
        self.overlap_2_unit = ttk.Label(self.frame_set_parameters, text="keV/µm", state=tk.DISABLED)
        self.overlap_2_unit.grid(row=3, column=4, padx=10, pady=10, sticky='w')

        # Empty row
        ttk.Label(self.frame_set_parameters, text="").grid(row=4, column=0, sticky=tk.W)

        # Scaling Factor 1
        self.label_scaling1 = ttk.Label(self.frame_set_parameters, text="Scaling Factor 1:", state=tk.DISABLED)
        self.label_scaling1.grid(row=5, column=0, padx=10, pady=10, sticky='w')
        self.entry_scaling1 = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_scaling1.grid(row=5, column=1, padx=10, pady=10, sticky='w')
        
        # Scaling Factor 2
        self.label_scaling2 = ttk.Label(self.frame_set_parameters, text="Scaling Factor 2:", state=tk.DISABLED)
        self.label_scaling2.grid(row=6, column=0, padx=10, pady=10, sticky='w')
        self.entry_scaling2 = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_scaling2.grid(row=6, column=1, padx=10, pady=10, sticky='w')

        # Empty row
        ttk.Label(self.frame_set_parameters, text="").grid(row=7, column=0, sticky=tk.W)

        # Stitchting point 1
        self.label_stitch1 = ttk.Label(self.frame_set_parameters, text="Stitching Point 1:", state=tk.DISABLED)
        self.label_stitch1.grid(row=8, column=0, padx=10, pady=10, sticky='w')
        self.entry_stitch1 = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_stitch1.grid(row=8, column=1, padx=10, pady=10, sticky='w')
        self.stitch_1_unit = ttk.Label(self.frame_set_parameters, text="keV/µm", state=tk.DISABLED)
        self.stitch_1_unit.grid(row=8, column=2, padx=10, pady=10, sticky='w')

        # Stitchting point 2
        self.label_stitch2 = ttk.Label(self.frame_set_parameters, text="Stitching Point 2:", state=tk.DISABLED)
        self.label_stitch2.grid(row=9, column=0, padx=10, pady=10, sticky='w')
        self.entry_stitch2 = ttk.Entry(self.frame_set_parameters, width=5, state=tk.DISABLED)
        self.entry_stitch2.grid(row=9, column=1, padx=10, pady=10, sticky='w')
        self.stitch_2_unit = ttk.Label(self.frame_set_parameters, text="keV/µm", state=tk.DISABLED)
        self.stitch_2_unit.grid(row=9, column=2, padx=10, pady=10, sticky='w')

    def create_plot_widgets(self):

        self.frame_right_side = ttk.Frame(self.root)
        self.frame_right_side.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')
        
        self.frame_plot = ttk.Frame(self.frame_right_side, borderwidth=2, relief='solid')
        self.frame_plot.grid(row=0, column=1, padx=10, pady=10, sticky='nsew')

        # Create a notebook to switch between different plots
        self.notebook = ttk.Notebook(self.frame_plot)
        self.notebook.grid(row=0, column=0, rowspan=12, padx=10, pady=10, sticky='nsew')

        # Create a Figure and Axes for the plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot)
        self.canvas_widget = self.canvas.get_tk_widget()
 
        # Initialize the check_vars and plot_tabs attributes
        self.plot_objects = {}

        # Empty row
        ttk.Label(self.frame_plot, text="").grid(row=1, column=0, sticky=tk.W)

        # Create textbox widget below the notebook with scrollbar for longer texts
        self.plot_info_textbox = tk.Text(self.frame_plot, width=50, height=10, wrap=tk.WORD, bg="white")
        scrollbar = ttk.Scrollbar(self.frame_plot, orient=tk.VERTICAL, command=self.plot_info_textbox.yview)
        self.plot_info_textbox.config(yscrollcommand=scrollbar.set)
        scrollbar.grid(row=20, column=1, padx=0, pady=10, sticky=tk.W+tk.E+tk.N+tk.S, rowspan=12)
        self.plot_info_textbox.grid(row=20, column=0, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

    def create_buttons(self):
        
        self.frame_buttons = ttk.Frame(self.frame_left_side, borderwidth=2, relief='solid')
        self.frame_buttons.grid(row=2, column=0, columnspan=2, padx=10, pady=10, sticky='nsew')

        # Create the buttons
        self.run_button = ttk.Button(self.frame_buttons, text="RUN", command=self.run, state=tk.DISABLED)
        self.run_button.grid(row=0, column=0, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.reset_button = ttk.Button(self.frame_buttons, text="RESET", command=self.reset_GUI)
        self.reset_button.grid(row=0, column=1, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.save_button = ttk.Button(self.frame_buttons, text="SAVE AS .csv", command=self.save_output_file, state=tk.DISABLED)
        self.save_button.grid(row=0, column=2, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.set_default_values_button = ttk.Button(self.frame_buttons, text="Set Default\nStart Values", command=self.set_default_values, state=tk.DISABLED)
        self.set_default_values_button.grid(row=0, column=3, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

    ################################################################
    #----------------------GUI FUNCTIONALITY-----------------------#
    ################################################################

    def welcome_message(self):

        fig, ax, canvas, _ = self.setup_plot('Welcome')
        canvas.draw()

        self.plot_objects['Welcome'] = (fig, ax, canvas)

        welcome_message(ax, "Merge Spectra Tool")
    
    def set_default_values(self):
        """ Put in some starting values for a first look """

        # Load the measurements if it is not already done
        if not hasattr(self, 'campaign'):
            self.load_measurements()
            print("Loaded measurements")

        # Some kind of function guessing parameters from the input spectra
        #TODO: Implement this shit

        # Momentane Notlösung
        if self.num_files == 3:
            nums = [2, 5, 10, 20, 3.5, 15]
        if self.num_files == 2:
            nums = [10,20,0,0,35,0]

        self.entry_overlap1_from.delete(0, tk.END)
        self.entry_overlap1_from.insert(0, nums[0])
        self.entry_overlap1_to.delete(0, tk.END)
        self.entry_overlap1_to.insert(0, nums[1])

        if self.num_files == 3:
            self.entry_overlap2_from.delete(0, tk.END)
            self.entry_overlap2_from.insert(0, nums[2])
            self.entry_overlap2_to.delete(0, tk.END)
            self.entry_overlap2_to.insert(0, nums[3])

        self.entry_scaling1.delete(0, tk.END)
        self.entry_scaling1.insert(0, 1.0)

        if self.num_files == 3:
            self.entry_scaling2.delete(0, tk.END)
            self.entry_scaling2.insert(0, 1.0)

        self.entry_stitch1.delete(0, tk.END)
        self.entry_stitch1.insert(0, nums[4])

        if self.num_files == 3:
            self.entry_stitch2.delete(0, tk.END)
            self.entry_stitch2.insert(0, nums[5])

    def load_measurements(self):
        # Load the measurements
        self.campaign = MicroDosimetry()
        self.meas_list_plots = []

        for file in self.files:
            filename, extension = os.path.splitext(os.path.normpath(file))
            try:
                meas_name = filename.split(os.path.sep)[-1].split(extension)[0]
            except:
                messagebox.showinfo("Error", "No file selected or format corrupted")
                raise ValueError("No file selected or format corrupted")
        
            # Check that these are actually analyzed .csv files
            #TODO: Check this

            # Check gains and stuff?

            # Files should not be log binned !!! Have to have a linear scale
            #TODO: Check this

            # Read in
            self.campaign.read_file(file)
            self.meas_list_plots.append(self.campaign.measurements[meas_name])

    def run(self):
        """ Run the merging of the spectra starting from the input files """
        
        if not self.files:
            messagebox.showerror("Error", "No input files selected.")
            return

        self.save_button.config(state=tk.DISABLED)
        
        # Close all plots
        self.close_plots(welcome=False)

        # Load the measurements if it is not already done
        if not hasattr(self, 'campaign'):
            self.load_measurements()

        # Run the merge command with the parameters set in GUI
        if len(self.files) == 2:
            merged_spec, info_dict = Spectrum.merge_spectra(self.meas_list_plots,
                                                            overlap_regions = [[float(self.entry_overlap1_from.get()),float(self.entry_overlap1_to.get())]],
                                                            scaling_factors = [float(self.entry_scaling1.get())],
                                                            stitching_points = [float(self.entry_stitch1.get())],
                                                            testplots=False)
            
        elif len(self.files) == 3:
            merged_spec, info_dict = Spectrum.merge_spectra(self.meas_list_plots,
                                                            overlap_regions = [[float(self.entry_overlap1_from.get()), float(self.entry_overlap1_to.get())],
                                                                               [float(self.entry_overlap2_from.get()),float(self.entry_overlap2_to.get())]],
                                                            scaling_factors = [float(self.entry_scaling1.get()), float(self.entry_scaling2.get())],
                                                            stitching_points = [float(self.entry_stitch1.get()), float(self.entry_stitch2.get())],
                                                            testplots=False)
            
        # Update scaling overlap regions, scaling factors and stitching points
        self.print_output_info()

        self.entry_overlap1_from.delete(0, tk.END)
        self.entry_overlap1_to.delete(0, tk.END)
        self.entry_overlap2_from.delete(0, tk.END)
        self.entry_overlap2_to.delete(0, tk.END)

        self.entry_overlap1_from.insert(0, round(info_dict['lower_overlap']['lower_energy_low'],2))
        self.entry_overlap1_to.insert(0, round(info_dict['lower_overlap']['upper_energy_high'],2))

        if len(self.files) == 3:
            self.entry_overlap2_from.insert(0, round(info_dict['upper_overlap']['lower_energy_low'], 2))
            self.entry_overlap2_to.insert(0, round(info_dict['upper_overlap']['upper_energy_high'], 2))
        
        self.entry_scaling1.delete(0, tk.END)
        self.entry_scaling1.insert(0, round(info_dict['scaling_factor_low'], 2))

        if len(self.files) == 3:
            self.entry_scaling2.delete(0, tk.END)
            self.entry_scaling2.insert(0, round(info_dict['scaling_factor_high'], 2))
        
        self.entry_stitch1.delete(0, tk.END)
        self.entry_stitch1.insert(0, round(info_dict['stitching_point_low'], 2))

        if len(self.files) == 3:
            self.entry_stitch2.delete(0, tk.END)
            self.entry_stitch2.insert(0, round(info_dict['stitching_point_high'], 2))
        
        # Add the new spectrum to the campaign
        self.campaign.add_measurement('merged', merged_spec)
        Spectrum.logarithmic_binning(self.campaign.measurements['merged'], 60) #TODO: Add these as options somewhere
        Spectrum.normalize_spectrum(self.campaign.measurements['merged'])

        self.add_plot_tab(None, info_dict=info_dict,
                          tab_name="merged", tab_text='Merged',
                          single=True, meas_name="merged")

        # Final stitched spectrum
        self.add_plot_tab(Spectrum._merge_plot_stitched, info_dict=info_dict, tab_name="stitched",
                          tab_text='Stitched', single=False, meas_name=None)
        
        # Merged scaled spectrum
        self.add_plot_tab(Spectrum._merge_plot_scaled, info_dict=info_dict, tab_name="merged_scaled",
                          tab_text='Merged Scaled', single=False, meas_name=None)
        
        # Merged unscaled spectrum
        self.add_plot_tab(Spectrum._merge_plot_unscaled, info_dict=info_dict, tab_name="merged_unscaled",
                          tab_text='Merged Unscaled', single=False, meas_name=None)
        
        # Original Spectra
        for i, meas_name in enumerate(self.campaign.measurements.keys()):

            if not meas_name == 'merged':
                gain = self.campaign.measurements[meas_name].gain
                self.add_plot_tab(None, info_dict=info_dict, tab_name=f"original_{gain}",
                                  tab_text=f'Gain {gain}', single=True, meas_name=meas_name)

        # Enable the save button
        self.save_button.config(state=tk.NORMAL)

    def reset_GUI(self):
        """ Reset the GUI to its initial state """
        
        # Clear the input fields
        self.entry_file_path.delete(0, tk.END)

        self.text_input_info.delete(1.0, tk.END)
        self.text_input_info.insert(tk.END, self.initial_text)

        # Clear the parameters
        self.entry_overlap1_from.delete(0, tk.END)
        self.entry_overlap1_to.delete(0, tk.END)
        self.entry_overlap2_from.delete(0, tk.END)
        self.entry_overlap2_to.delete(0, tk.END)
        self.entry_scaling1.delete(0, tk.END)
        self.entry_scaling2.delete(0, tk.END)
        self.entry_stitch1.delete(0, tk.END)
        self.entry_stitch2.delete(0, tk.END)

        # Lock all entries
        self.label_overlap1.config(state=tk.DISABLED)
        self.entry_overlap1_from.config(state=tk.DISABLED)
        self.label_overlap1_to.config(state=tk.DISABLED)
        self.entry_overlap1_to.config(state=tk.DISABLED)
        self.overlap_1_unit.config(state=tk.DISABLED)

        self.label_overlap2.config(state=tk.DISABLED)
        self.entry_overlap2_from.config(state=tk.DISABLED)
        self.label_overlap2_to.config(state=tk.DISABLED)
        self.entry_overlap2_to.config(state=tk.DISABLED)
        self.overlap_2_unit.config(state=tk.DISABLED)

        self.label_scaling1.config(state=tk.DISABLED)
        self.entry_scaling1.config(state=tk.DISABLED)
        self.label_scaling2.config(state=tk.DISABLED)
        self.entry_scaling2.config(state=tk.DISABLED)

        self.label_stitch1.config(state=tk.DISABLED)
        self.entry_stitch1.config(state=tk.DISABLED)
        self.stitch_1_unit.config(state=tk.DISABLED)
        self.label_stitch2.config(state=tk.DISABLED)
        self.entry_stitch2.config(state=tk.DISABLED)
        self.stitch_2_unit.config(state=tk.DISABLED)

        # Disable buttons
        self.run_button.config(state=tk.DISABLED)
        self.save_button.config(state=tk.DISABLED)
        self.set_default_values_button.config(state=tk.DISABLED)

        # Close all plots
        self.close_plots(welcome=True)
        
    def browse_input_files(self):
        """ Browse for two or three analyzed .csv files to merge """

        file_path = filedialog.askopenfilenames(initialdir=self.last_browsed_path,
                                                title="Select the input files",
                                                filetypes=[("Analyzed CSV files", "*.csv")],
                                                multiple=True)
        
        normpath = [os.path.normpath(file) for file in file_path]
        self.last_browsed_path = os.path.dirname(normpath[0])

        self.num_files = len(file_path)

        # Enable the parameter entries
        self.label_overlap1.config(state=tk.NORMAL)
        self.entry_overlap1_from.config(state=tk.NORMAL)
        self.label_overlap1_to.config(state=tk.NORMAL)
        self.entry_overlap1_to.config(state=tk.NORMAL)
        self.overlap_1_unit.config(state=tk.NORMAL)

        if self.num_files == 3:
            self.label_overlap2.config(state=tk.NORMAL)
            self.entry_overlap2_from.config(state=tk.NORMAL)
            self.label_overlap2_to.config(state=tk.NORMAL)
            self.entry_overlap2_to.config(state=tk.NORMAL)
            self.overlap_2_unit.config(state=tk.NORMAL)
        
        self.label_scaling1.config(state=tk.NORMAL)
        self.entry_scaling1.config(state=tk.NORMAL)

        if self.num_files == 3:
            self.label_scaling2.config(state=tk.NORMAL)
            self.entry_scaling2.config(state=tk.NORMAL)

        self.label_stitch1.config(state=tk.NORMAL)
        self.entry_stitch1.config(state=tk.NORMAL)
        self.stitch_1_unit.config(state=tk.NORMAL)

        if self.num_files == 3:
            self.label_stitch2.config(state=tk.NORMAL)
            self.entry_stitch2.config(state=tk.NORMAL)
            self.stitch_2_unit.config(state=tk.NORMAL)

        if self.num_files > 3:
            messagebox.showerror("Sorry this doesnt work :(", "WTF dude, how many amps do you have. Please select only 2 or 3 files.")
            self.entry_file_path.delete(0, tk.END)
            self.entry_file_path.insert(0, "Please select 2 or 3 files.")
            return
        else:
            self.entry_file_path.delete(0, tk.END)
            self.entry_file_path.insert(0, file_path)
            self.print_input_info()
        
        # Enable buttons
        self.run_button.config(state=tk.NORMAL)
        self.set_default_values_button.config(state=tk.NORMAL)

    def print_input_info(self):
        """ Print out names of the files + metainformation """

        file_path = self.entry_file_path.get()
        self.files = file_path.split(" ")

        # Clear the textbox
        self.text_input_info.delete(1.0, tk.END)

        # Put in file names
        for file in self.files:
            # Strip off the path and the extension
            file = os.path.basename(file)
            file = file.split(".")[0]
            self.text_input_info.insert(tk.END, file + "\n")

        self.text_input_info.insert(tk.END, "\n")

        #TODO: Actually retrieve the information and print it out

        # x-axis
        self.text_input_info.insert(tk.END, "x axis: Energy [keV/µm]\n")

        # y-axis
        self.text_input_info.insert(tk.END, "y axis: Counts [a.u.]\n")

        # Gains
        self.text_input_info.insert(tk.END, "Gains: [blablabla]\n")
        
    def browse_output_file(self):
        """ Browse for a location to save the merged spectrum """

        file_path = filedialog.asksaveasfilename(initialdir=self.last_browsed_path,
                                                title="Select the output file",
                                                filetypes=[("Analyzed CSV files", "*.csv")],
                                                defaultextension=".csv",
                                                confirmoverwrite=True)
        
        self.last_browsed_path = os.path.dirname(file_path)
        
        return file_path

    def save_output_file(self):
        """ Save the merged spectrum as a .csv file """

        self.output_path = self.browse_output_file()

        # Save the file
        for meas_name in self.campaign.measurements.keys():
            if meas_name == 'merged':
                Output.csv_output(self.campaign.measurements[meas_name], self.output_path)

    def print_output_info(self):
        """ Information about the current stats of stitching and scaling """

        # TODO: Print something interesting

        pass

    ################################################################
    #----------------------PLOTTING--------------------------------#
    ################################################################

    def setup_plot(self, tab_text):
        """ Create a new tab for the plot and set up the canvas """

        # Just a white screen for a smoother transition between plots
        self.canvas.draw()

        # Add the new tab to the notebook
        new_tab = ttk.Frame(self.notebook)
        self.notebook.add(new_tab, text=tab_text)
        #self.notebook.select(tab_text)

        #Open up a new canvas
        fig, ax = plt.subplots(figsize=(8, 6))  
        canvas = FigureCanvasTkAgg(fig, master=new_tab)
        canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Include space for the toolbar
        #fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        fig.subplots_adjust(bottom=0.1)     

        return fig, ax, canvas, new_tab

    def add_plot_tab(self, plot_function, info_dict, tab_name, tab_text, single=False, meas_name=None):
        
        # Possibly setting plot options and parameters before
        # TODO: Add these options to the GUI and implement them

        # Setup plot canvas
        fig, ax, canvas, new_tab = self.setup_plot(tab_text)

        # For single spectra use standard plot function
        if single:
            Output.plot_single(self.campaign.measurements[meas_name],
                               show_plot=False, gui=True, fig=fig, ax=ax) #TODO: Add all options needed
        else:
            # Run the plot function
            plot_function(self.meas_list_plots, info_dict, ax)

        toolbar = NavigationToolbar2Tk(canvas, new_tab, pack_toolbar=False)
        toolbar.grid(row=1, column=0, padx=10, pady=10, sticky=tk.S)
        canvas.draw()
 
        # Store the new plot objects in the dictionary
        self.plot_objects[tab_name] = (fig, ax, canvas)

    def close_plots(self, welcome=False):
        """ Close all open plots and open the welcome message """

        # Close all tabs
        for tab in self.notebook.tabs():
            self.notebook.forget(tab)
        
        # Reset RUN button
        self.run_button.config(state=tk.NORMAL)
        
        self.plot_objects = {}

        plt.close('all')

        if welcome:
            self.welcome_message()

            # Change textbox
            message = "Closed plots. Ready to start again ..."
            self.plot_info_textbox.delete(1.0, tk.END)
            self.plot_info_textbox.insert(1.0, message)

root = tk.Tk()
app = GUI_Merge(root)

root.mainloop()