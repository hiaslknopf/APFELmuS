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

    myappid = 'APFELMUS_Analysis' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        # Main window setup
        self.root = root
        self.root.title("APFELmuS - Spectrum Viewer")

        # Load the APFELmuS logo
        root.iconbitmap("ressources/logo.ico")

        # Create the widgets
        self.create_cal_widgets()
        self.create_plot_selection_widgets()
        self.create_output_widgets()
        self.create_plot_window()

        self.welcome_message() 
    
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
        self.browse_linearization_path = ttk.Button(self.frame_file_cal, text="Browse", command=self.browse_cal)
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
        self.calibration_method_check = ttk.Checkbutton(self.frame_file_cal, text="Calibrate with edge + chord length", variable=self.calibration_method_var, command=self.toggle_calib_method)
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

        # Empty line
        ttk.Label(self.frame_file_cal, text="").grid(row=7, column=0, columnspan=3, sticky=tk.W)

        # Variable to trace changes in the file entry
        self.file_var = tk.StringVar()
        self.file_var.trace_add("write", self.check_file_entry)
        self.entry_file_path["textvariable"] = self.file_var

        # Option to cut off spectrum
        self.cutoff_var = tk.BooleanVar(value=False)
        self.cutoff_check = ttk.Checkbutton(self.frame_file_cal, text="Cut off spectrum:", variable=self.cutoff_var, command=self.toggle_cutoff)
        self.cutoff_check.grid(row=9, column=0, sticky=tk.W)

        # Entry for cutoff value
        self.entry_cutoff = ttk.Entry(self.frame_file_cal, state=tk.DISABLED)
        self.entry_cutoff.grid(row=9, column=1, sticky=tk.W)
        ttk.Label(self.frame_file_cal, text="channel").grid(row=9, column=2, sticky=tk.W)

        # Empty line
        ttk.Label(self.frame_file_cal, text="").grid(row=10, column=0, columnspan=3, sticky=tk.W)

        # Two dropdown menus for info_dict detector and material
        ttk.Label(self.frame_file_cal, text="Detector:").grid(row=11, column=0, sticky=tk.W)
        self.detector_var = tk.StringVar()
        self.detector_var.set("silicon")
        self.detector_dropdown = ttk.OptionMenu(self.frame_file_cal, self.detector_var, "silicon", "silicon", "diamond", "sic")
        self.detector_dropdown.grid(row=11, column=1, sticky=tk.W)

        ttk.Label(self.frame_file_cal, text="Particle:").grid(row=12, column=0, sticky=tk.W)
        self.particle_var = tk.StringVar()
        self.particle_var.set("proton")
        self.material_dropdown = ttk.OptionMenu(self.frame_file_cal, self.particle_var, "proton", "proton", "carbon", "alpha")
        self.material_dropdown.grid(row=12, column=1, sticky=tk.W)      

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

        # Change to step plot
        self.step_var = tk.BooleanVar(value=False)
        self.step_check = ttk.Checkbutton(self.frame_plot, text="Step plot", variable=self.step_var)
        self.step_check.grid(row=10, column=0, sticky=tk.W)
        
    def create_output_widgets(self):

        ################################################################
        #-----------------------Output section-------------------------#
        ################################################################

        self.frame_output = ttk.Frame(root, borderwidth=2, relief="solid")
        self.frame_output.grid(row=3, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)
        
        ttk.Label(self.frame_output, text="Save plots:").grid(row=0, column=0, sticky=tk.W)
        self.save_plots_var = tk.BooleanVar(value=False)
        self.check_save_plots = ttk.Checkbutton(self.frame_output, variable=self.save_plots_var, command=self.check_out_entry).grid(row=0, column=1, sticky=tk.W)

        ttk.Label(self.frame_output, text="Output:").grid(row=1, column=0, sticky=tk.W)
        self.entry_out = ttk.Entry(self.frame_output,state=tk.DISABLED)
        self.entry_out.grid(row=1, column=1, sticky=tk.W)
        self.browse_out = ttk.Button(self.frame_output, text="Browse", command=self.browse_out,state=tk.DISABLED)
        self.browse_out.grid(row=1, column=2, sticky=tk.W)

        ################################################################
        #-----------------------Button section-------------------------#
        ################################################################

        self.button_frame = ttk.Frame(root)
        self.button_frame.grid(row=11, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.run_button = ttk.Button(self.button_frame, text="RUN", command=self.run)
        self.run_button.grid(row=0, column=0, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.run_button = ttk.Button(self.button_frame, text="CLOSE PLOTS", command=self.close_plots)
        self.run_button.grid(row=0, column=1, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        self.reset_button = ttk.Button(self.button_frame, text="RESET", command=self.reset_GUI)
        self.reset_button.grid(row=0, column=2, padx=5, pady=5, ipadx=10, ipady=10, sticky=tk.W+tk.E+tk.N+tk.S)

    def create_plot_window(self):

        ################################################################
        #-----------------------Plot window----------------------------#
        ################################################################

        # Create notebook for tabs
        self.notebook = ttk.Notebook(root)
        self.notebook.grid(row=0, column=1, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)
        
        self.frame_plot_window = ttk.Frame(root)
        self.frame_plot_window.grid(row=0, column=3, rowspan=12, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Create a Figure and Axes for the plot
        self.fig, self.ax = plt.subplots(figsize=(6, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame_plot_window)
        self.canvas_widget = self.canvas.get_tk_widget()
 
        # Initialize the check_vars and plot_tabs attributes
        self.plot_objects = {}
        self.plot_tabs = {}       

################################################################
#----------------------GUI FUNCTIONALITY-----------------------#
################################################################

    def welcome_message(self):

        tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
        fig, ax, canvas = self.setup_plot('Welcome', tab_text)
        canvas.draw()

        self.plot_objects['Welcome'] = (fig, ax, canvas)

        welcome_message(ax, "Analysis Toolkit")

    def reset_GUI(self):
        """ Reset the GUI to its initial state """

        # Reset all checkmarks
        for var in self.check_vars:
            var.set(False)
        self.log_binning_var.set(True)
        self.save_plots_var.set(False)
        self.cutoff_var.set(False)

        # Reset file and calibration entries
        self.entry_file_path.delete(0, tk.END)
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_edge_pos.delete(0, tk.END)
        self.entry_mean_chord_length.delete(0, tk.END)
        self.entry_cal_factor.delete(0, tk.END)

        # Reset dropdown
        self.detector_var.set("silicon")
        self.particle_var.set("carbon")

        # Reset the output entry
        self.entry_out.delete(0, tk.END)
        
        # Lock entries
        self.entry_linearization_path.config(state=tk.NORMAL)
        self.browse_linearization_path.config(state=tk.NORMAL)
        self.entry_edge_pos.config(state=tk.DISABLED)
        self.entry_mean_chord_length.config(state=tk.DISABLED)
        self.entry_num_bins.config(state=tk.NORMAL)
        self.entry_out.config(state=tk.DISABLED)
        self.browse_out.config(state=tk.DISABLED)
        self.entry_cutoff.config(state=tk.DISABLED)

        # Reset the plot tabs
        for tab_name in list(self.plot_tabs.keys()):
            tab_text = self.plot_tabs.pop(tab_name)
            index = self.notebook.index(tab_text)
            self.notebook.forget(index)
        
        # Reset the plot objects
        self.plot_objects = {}
        self.plot_tabs = {}

        # Reset the canvas
        self.welcome_message()

    def check_out_entry(self, *args):
        """ Lock output browse if no output is selected """

        self.close_plots()

        if self.save_plots_var.get():
            self.entry_out.config(state=tk.NORMAL)
            self.browse_out.config(state=tk.NORMAL)
        else:
            self.entry_out.config(state=tk.DISABLED)
            self.browse_out.config(state=tk.DISABLED)

    def check_file_entry(self, *args):
        """ Lock input options for different file types """
        file_content = self.file_var.get()

        # Check if the content contains the ".MCA" extension
        if ".MCA" in file_content.upper():
            # Disable the calibration entry if the condition is met
            self.entry_linearization_path.config(state=tk.DISABLED)
            self.browse_linearization_path.config(state=tk.DISABLED)
        elif ".Spe" in file_content.lower():
            # Enable the calibration entry if the condition is met
            self.entry_linearization_path.config(state=tk.NORMAL)
            self.browse_linearization_path.config(state=tk.NORMAL)
        elif ".root" in file_content.lower():
            # Disable the calibration entry if the condition is met
            self.entry_linearization_path.config(state=tk.DISABLED)
            self.browse_linearization_path.config(state=tk.DISABLED)
    
    def toggle_cutoff(self):
        if self.cutoff_var.get():
            self.entry_cutoff.config(state=tk.NORMAL)
        else:
            self.entry_cutoff.config(state=tk.DISABLED)

    def toggle_entry_num_bins(self):
        if self.log_binning_var.get():
            self.entry_num_bins.config(state=tk.NORMAL)
        else:
            self.entry_num_bins.config(state=tk.DISABLED)
    
    def toggle_calib_method(self):
        if self.calibration_method_var.get():
            self.entry_edge_pos.config(state=tk.NORMAL)
            self.entry_mean_chord_length.config(state=tk.NORMAL)
            self.entry_cal_factor.config(state=tk.DISABLED)
        else:
            self.entry_edge_pos.config(state=tk.DISABLED)
            self.entry_mean_chord_length.config(state=tk.DISABLED)
            self.entry_cal_factor.config(state=tk.NORMAL)

    def clear_num_bins_placeholder(self, event):
        # Clear the placeholder text when the entry gets focus
        if self.entry_num_bins.get() == self.num_bins_placeholder:
            self.entry_num_bins.delete(0, tk.END)

    def restore_num_bins_placeholder(self, event):
        # Restore the placeholder text if the entry is empty
        if not self.entry_num_bins.get():
            self.entry_num_bins.insert(0, self.num_bins_placeholder)

    def browse_file(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("All files", "*.*"),("MCA files", "*.MCA"),("MAESTRO files", "*.Spe"),("ROOT files", "*.root")])
        self.entry_file_path.delete(0, tk.END)
        self.entry_file_path.insert(0, file_path)

    def browse_cal(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        cal_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("csv files", "*.csv"), ("All files", "*.*")])
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_linearization_path.insert(0, cal_path)

    def browse_out(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        output_path = filedialog.askdirectory(initialdir=script_directory)
        self.entry_out.delete(0, tk.END)
        self.entry_out.insert(0, output_path)

    def close_plots(self):
        for tab_name in list(self.plot_tabs.keys()):
            tab_text = self.plot_tabs.pop(tab_name)
            index = self.notebook.index(tab_text)
            self.notebook.forget(index)
        
        self.plot_objects = {}
        self.plot_tabs = {}

        self.welcome_message()

#################################################################
#----------------------PLOTTING FUNCTIONALITY-------------------#
#################################################################

    def calibrate_specs(self, campaign):
        #Calibrate
        if self.extension == '.root':
            Calibrate.lineal_energy_axis(campaign.measurements[self.meas_name], chord_length='mean')
        if self.extension == '.MCA':

            if self.calibration_method_var.get():
                if self.entry_mean_chord_length.get() == "":
                    raise ValueError("Mean chord length has to be specified")
                if self.entry_edge_pos.get() == "":
                    raise ValueError("Calibration edge has to be specified")
                Calibrate.get_chord_length(campaign.measurements[self.meas_name], 'slab', float(self.entry_mean_chord_length.get()), plot=False) 
                ymax, _ = Calibrate.get_stopping_power(campaign.measurements[self.meas_name], chord_length='mean', database='SRIM', precision=0.01, plot=False)
                Calibrate.scale_energy_axis(campaign.measurements[self.meas_name], float(self.entry_edge_pos.get()), ymax, energy_axis='lineal', chord_length='mean')
            
            else:
                if self.entry_cal_factor.get() == "":
                    cal_factor = 3
                    messagebox.showinfo("Calibration factor", f"Assumed calibration factor: {cal_factor}")
                else:
                    cal_factor = float(self.entry_cal_factor.get())
                
                Calibrate.scale_energy_axis_with_factor(campaign.measurements[self.meas_name], cal_factor)
    
    def setup_plot(self, tab_name, tab_text):

        # Just a white screen for smoother transition
        self.canvas.draw()

        # Add the tab to the notebook
        self.notebook.add(tab_text, text=tab_name)
        self.notebook.select(tab_text)
        self.plot_tabs[tab_name] = tab_text

        #Open up a new canvas
        fig, ax = plt.subplots()
        canvas = FigureCanvasTkAgg(fig, master=tab_text)
        canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        return fig, ax, canvas

    def show_plot(self, index, option, campaign):

        selected_option = self.check_vars[index].get()

        if selected_option:
            # Create a new tab for each plot
            tab_name = f"{option}"
            if tab_name not in self.plot_tabs:
                tab_text = tk.Text(self.notebook, width=40, height=20, wrap=tk.WORD, bg="white", state=tk.DISABLED)
                tab_text.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

                print('\n-----------ENTERED PLOTTING---------------')
                # Start fresh with the original spectrum - Shit performance but it works
                Spectrum.retrieve_original_spectrum(campaign.measurements[self.meas_name])
                
                if self.cutoff_var.get():
                    if self.entry_cutoff.get() == "Enter cutoff value":
                        raise ValueError("Cutoff value has to be specified")
                    Spectrum.cutoff(campaign.measurements[self.meas_name], int(self.entry_cutoff.get()))

                if index != 0:
                    self.calibrate_specs(campaign)

                print('Plotting index ' + str(index) + ' ' + str(option) + '\n')

                if index == 1:
                    #f(y)
                    Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                    Spectrum.probability_density(campaign.measurements[self.meas_name])
                if index == 2:
                    #F(y)
                    Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                if index == 3:
                    #d(y)
                    Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                    Spectrum.probability_density(campaign.measurements[self.meas_name])
                if index == 4:
                    #yf(y)
                    Spectrum.probability_function(campaign.measurements[self.meas_name], 'F')
                    Spectrum.probability_density(campaign.measurements[self.meas_name])
                    Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                if index == 5:
                    #yd(y)
                    Spectrum.probability_function(campaign.measurements[self.meas_name], 'D')
                    Spectrum.probability_density(campaign.measurements[self.meas_name])
                    Spectrum.weighted_probability_density(campaign.measurements[self.meas_name])
                
                if self.log_binning_var.get() and index != 0:
                    # Log bin the spectrum
                    if self.entry_num_bins.get() == self.num_bins_placeholder:
                        messagebox.showinfo("Number of bins", "Assumed 50 log bins")
                        num_bins = 50
                    else:
                        num_bins = int(self.entry_num_bins.get())

                    Spectrum.logarithmic_binning(campaign.measurements[self.meas_name], num_bins)
                
                if not index in [0,2]:
                    Spectrum.normalize_spectrum(campaign.measurements[self.meas_name])
                
                fig, ax, canvas = self.setup_plot(tab_name, tab_text)

                if self.save_plots_var.get():
                    output_path = self.entry_out.get()
                    if output_path == "":
                        raise ValueError("Output path has to be specified")
                else:
                    output_path = False

                if index == 0:
                    scale = 'lin'
                else:
                    scale = 'log'

                Output.plot_single(campaign.measurements[self.meas_name], name=f'{self.meas_name}_{self.plot_options[index]}', output_path=output_path,
                                   step=self.step_var.get(), mean=False, scale=scale,
                                   show_plot=False, gui=True, fig=fig, ax=ax)
             
                toolbar = NavigationToolbar2Tk(canvas, tab_text, pack_toolbar=False)
                toolbar.grid(row=1, column=0, padx=10, pady=10, sticky=tk.S)
                canvas.draw()

                # Store the plot objects in the dictionary
                self.plot_objects[tab_name] = (fig, ax, canvas)
        else:
            # If the checkmark is unselected, remove the tab from the notebook
            tab_name = f"{option}"
            if tab_name in self.plot_tabs:
                tab_text = self.plot_tabs.pop(tab_name)
                index = self.notebook.index(tab_text)
                self.notebook.forget(index)
                self.plot_objects.pop(tab_name)
            
            # Remove the welcome message
            if 'Welcome' in self.plot_tabs:
                tab_text = self.plot_tabs.pop('Welcome')
                index = self.notebook.index(tab_text)
                self.notebook.forget(index)
                self.plot_objects.pop('Welcome')

    def run(self):

        # Close the open plots
        self.close_plots()

        #Create object and read in
        campaign = MicroDosimetry()
        self.datafile = self.entry_file_path.get()
        self.filename, self.extension = os.path.splitext(os.path.normpath(self.datafile))
        self.meas_name = self.filename.split(os.path.sep)[-1].split(self.extension)[0]

        if self.extension == '.Spe':
            if self.entry_linearization_path.get() == "":
                raise ValueError("Calibration file has to be specified")
                
            info_dict = {'DETECTOR': self.detector_var.get(), 'PARTICLE': self.particle_var.get(), 'GAIN': 'None'}

            FileTranslator.translate_MAESTRO_file(input_MAESTRO=self.datafile, input_linearization=self.entry_linearization_path.get(), output_path='temp', name=self.meas_name, info_dict=info_dict)
            self.extension = '.MCA' #The new extension
            campaign.read_file(f'temp/{self.meas_name}.MCA') #The new file
        else:
            campaign.read_file(self.datafile)

        test = []
        for i in range(len(self.plot_options)):
            test.append(self.check_vars[i].get())
        if not any(test):
            messagebox.showinfo("No plot selected", "No plot selected")
        else:
            #Iterate through the selected plot options
            for i, option in enumerate(self.plot_options):
                self.show_plot(i, option, campaign)
        
root = tk.Tk()
app = GUI_Analysis(root)

# Fixed Window size
#root.geometry("1035x715")
#root.resizable(False, False)

root.mainloop()