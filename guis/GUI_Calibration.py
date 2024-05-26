__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import os
import sys

import ctypes

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Spectrum, Calibrate, FileTranslator
from guis.welcome_message import welcome_message

class GUI_Calibration:

    myappid = 'APFELMUS_Calibration' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        self.root = root
        self.root.title("APFElmuS - Edge Calibration")
        self.create_widgets()

        # Load the APFELmuS logo
        root.iconbitmap("ressources/logo.ico")
    
    def create_sliders(self):
        # Sliders
        self.sliders = {
            "Fit bound LOW": (1, self.num_channels/4),
            "Fit bound HIGH": (self.num_channels/4, self.num_channels-1),
            "Spectrum cutoff": (0, self.num_channels-500)
        }

        self.slider_vars = {name: tk.IntVar() for name in self.sliders}

        #Set initial values for sliders
        self.slider_vars["Spectrum cutoff"].set(0)
        self.slider_vars["Fit bound LOW"].set(1)
        self.slider_vars["Fit bound HIGH"].set(self.num_channels-1)

        i = 0
        for label, (min_val, max_val) in self.sliders.items():
            ttk.Label(self.slider_frame, text=label + " [CH]:").grid(row=i, column=0, sticky=tk.W)
            ttk.Label(self.slider_frame, text='\t').grid(row=i, column=1, sticky=tk.W)

            slider = ttk.Scale(self.slider_frame, from_=min_val, to=max_val, orient=tk.HORIZONTAL, variable=self.slider_vars[label], command=lambda val, name=label: self.on_slider_change(val, name))
            slider.grid(row=i, column=2, sticky=tk.W+tk.E)
        
            i += 1

        # Text entry and label for the slider value
        self.low_bound_entry = ttk.Entry(self.slider_frame, textvariable=self.slider_vars["Fit bound LOW"], width=5)
        self.low_bound_entry.grid(row=0, column=3, sticky=tk.W+tk.E)

        self.high_bound_entry = ttk.Entry(self.slider_frame, textvariable=self.slider_vars["Fit bound HIGH"], width=5)
        self.high_bound_entry.grid(row=1, column=3, sticky=tk.W+tk.E)

        self.cutoff_entry = ttk.Entry(self.slider_frame, textvariable=self.slider_vars["Spectrum cutoff"], width=5)
        self.cutoff_entry.grid(row=2, column=3, sticky=tk.W+tk.E)

    def create_widgets(self):
        self.left_panel = ttk.Frame(root)
        self.left_panel.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        #Subframe for file and calibration file
        self.file_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.file_frame.grid(row=0, column=0, columnspan=3, sticky=tk.W+tk.E)

        #Browse for file
        ttk.Label(self.file_frame, text="File:").grid(row=0, column=0, sticky=tk.W)
        self.entry_file_path = ttk.Entry(self.file_frame)
        self.entry_file_path.grid(row=0, column=1, sticky=tk.W)
        self.browse_file_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_file)
        self.browse_file_path.grid(row=0, column=2, sticky=tk.W)

        #Browse for calibration file
        ttk.Label(self.file_frame, text="Linearization:").grid(row=1, column=0, sticky=tk.W)
        self.entry_linearization_path = ttk.Entry(self.file_frame)
        self.entry_linearization_path.grid(row=1, column=1, sticky=tk.W)
        self.browse_linearization_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_file)
        self.browse_linearization_path.grid(row=1, column=2, sticky=tk.W)

        # Empty line
        ttk.Label(self.file_frame, text="").grid(row=2, column=0, sticky=tk.W)
        
        # Two dropdown menus for info_dict detector and material
        ttk.Label(self.file_frame, text="Detector:").grid(row=3, column=0, sticky=tk.W)
        self.detector_var = tk.StringVar()
        self.detector_var.set("silicon")
        self.detector_dropdown = ttk.OptionMenu(self.file_frame, self.detector_var, "silicon", "silicon", "diamond", "sic")
        self.detector_dropdown.grid(row=3, column=1, sticky=tk.W)

        ttk.Label(self.file_frame, text="Particle:").grid(row=4, column=0, sticky=tk.W)
        self.particle_var = tk.StringVar()
        self.particle_var.set("proton")
        self.material_dropdown = ttk.OptionMenu(self.file_frame, self.particle_var, "proton", "proton", "carbon", "alpha")
        self.material_dropdown.grid(row=4, column=1, sticky=tk.W)

        # Entry for chord length
        ttk.Label(self.file_frame, text="Chord length:").grid(row=5, column=0, sticky=tk.W)
        self.entry_chord_length = ttk.Entry(self.file_frame)
        self.entry_chord_length.grid(row=5, column=1, sticky=tk.W)
        ttk.Label(self.file_frame, text="um").grid(row=5, column=2, sticky=tk.W)

        # Empty line between frames
        ttk.Label(self.left_panel, text="").grid(row=1, column=0, sticky=tk.W)

        # New frame for sliders
        self.slider_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.slider_frame.grid(row=2, column=0, columnspan=3, sticky=tk.W+tk.E)

        self.num_channels = 500 #Initial value

        # Sliders
        self.create_sliders()
        
        # Empty line between frames
        ttk.Label(self.left_panel, text="").grid(row=3, column=0, sticky=tk.W)

        # Frame for buttons
        self.button_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.button_frame.grid(row=4, column=0, columnspan=3, sticky=tk.W+tk.E)

        # Bottom Buttons
        self.run_button = ttk.Button(self.button_frame, text="RUN", command=self.run)
        self.run_button.grid(row=0, column=0, sticky=tk.W+tk.E, ipadx=50, ipady=10)
        self.rerun_button = ttk.Button(self.button_frame, text="Rerun Edge Calibration", command=self.rerun_edge, state='disabled')
        self.rerun_button.grid(row=1, column=0, sticky=tk.W+tk.E, ipadx=20, ipady=10)
        self.save_button = ttk.Button(self.button_frame, text="Save", command=self.save_values, state='disabled')
        self.save_button.grid(row=2, column=0, sticky=tk.W+tk.E, ipadx=20, ipady=10)
        self.reset_button = ttk.Button(self.button_frame, text="Reset GUI", command=self.reset_GUI)
        self.reset_button.grid(row=3, column=0, sticky=tk.W+tk.E, ipadx=20, ipady=10)

        # Right panel for plots and results
        self.right_panel = ttk.Frame(root, borderwidth=2, relief="solid")
        self.right_panel.grid(row=0, column=1, padx=20, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Notebook for plots
        self.notebook = ttk.Notebook(self.right_panel)
        self.notebook.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)
        self.add_initial_tab()

        # Results textbox below the notebook
        ttk.Label(self.right_panel, text="Results:").grid(row=1, column=0, sticky=tk.W)
        self.textbox = tk.Text(self.right_panel, height=12.5, width=30)
        self.textbox.grid(row=1, column=0, sticky=tk.W+tk.E)

        # Write sth in the textbox
        self.textbox.insert(tk.END, "Welcome to the edge calibration tool!\nPlease specify a file and the linearization file.\n\n")
        self.textbox.insert(tk.END, "After running the calibration, the results will be displayed here...")

        # Variable to trace changes in the file entry
        self.file_var = tk.StringVar()
        self.file_var.trace_add("write", self.check_file_entry)
        self.entry_file_path.config(textvariable=self.file_var)

    def browse_file(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("All files", "*.*"),("MCA files", "*.MCA"),("MAESTRO files", "*.Spe")])
        self.entry_file_path.delete(0, tk.END)
        self.entry_file_path.insert(0, file_path)

    def rerun_edge(self):

        # Cut the spectrum
        Spectrum.cutoff(self.campaign.measurements[self.meas_name], int(self.slider_vars["Spectrum cutoff"].get()))

        # Calculate the edge position and stopping power
        _, self.edge_dict = Calibrate.get_edge_pos(self.campaign.measurements[self.meas_name], 'hTC',
                               fit_bounds=[int(self.slider_vars["Fit bound LOW"].get()), int(self.slider_vars["Fit bound HIGH"].get())],
                               check_plot=False)
        
        print(f'Fit bounds: {int(self.slider_vars["Fit bound LOW"].get())} - {int(self.slider_vars["Fit bound HIGH"].get())}')
            
        # Update the edge plot (Close the old one)
        for i in range(len(self.notebook.tabs())-1, -1, -1):
            self.notebook.forget(i)

        self.add_plot_tab(Calibrate._plot_edge_calibration, self.edge_dict, title="Edge Spectrum")
        self.add_plot_tab(Calibrate._plot_stopping_power, self.stop_pow_dict, title="Stopping Power")

        # Update the text box
        self.update_textbox(self.edge_dict, self.stop_pow_dict)

    def reset_GUI(self):
        # Reset all entries and sliders
        self.entry_file_path.delete(0, tk.END)
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_chord_length.delete(0, tk.END)
        self.slider_vars["Fit bound LOW"].set(1)
        self.slider_vars["Fit bound HIGH"].set(4095)
        self.slider_vars["Spectrum cutoff"].set(0)

        # Reset the notebook
        for i in range(len(self.notebook.tabs())-1, -1, -1):
            self.notebook.forget(i)	
        self.add_initial_tab()

        # Unlock the linearization entry
        self.entry_linearization_path.config(state='normal')
        self.browse_linearization_path.config(state='normal')

        # Delete results from textbox
        self.textbox.delete(1.0, tk.END)

        # Reset the buttons
        self.run_button.config(state='normal')
        self.rerun_button.config(state='disabled')
        self.save_button.config(state='disabled')

        # Write sth in the textbox
        self.textbox.insert(tk.END, "Starting fresh ... \n")

    def run(self):

        # Check if file is specified
        if self.entry_file_path.get() == '':
            messagebox.showwarning("Warning", "Please specify a file")
            self.reset_GUI()

        # Read in file and create the MicroDosimetry object
        self.campaign = MicroDosimetry()
        self.datafile = self.entry_file_path.get()
        self.filename, self.extension = os.path.splitext(os.path.normpath(self.datafile))
        self.meas_name = self.filename.split(os.path.sep)[-1].split(self.extension)[0]

        if self.extension == '.Spe':
            if self.entry_linearization_path.get() == '':
                raise ValueError("Calibration file has to be specified")
                
            info_dict = {'DETECTOR': self.detector_var.get(), 'PARTICLE': self.particle_var.get(), 'GAIN': 'None'}

            FileTranslator.translate_MAESTRO_file(input_MAESTRO=self.datafile, input_linearization=self.entry_file_path.get(), output_path='temp', name=self.meas_name, info_dict=info_dict)
            self.extension = '.MCA' #The new extension
            self.campaign.read_file(f'temp/{self.meas_name}.MCA') #The new file
        else:
            self.campaign.read_file(self.datafile)
        
        # Update sliders
        self.num_channels = self.campaign.measurements[self.meas_name].num_channels
        self.create_sliders()

        # Attach chord length to the measurement
        try:
            chord_length = float(self.entry_chord_length.get())
            self.campaign.measurements[self.meas_name].mean_chord_length = chord_length
        except:
            messagebox.showwarning("Warning", "Please specify a chord length")
            self.reset_GUI()
            return

        # Info for tetxbox
        self.particle = self.campaign.measurements[self.meas_name].particle
        self.detector = self.campaign.measurements[self.meas_name].detector

        # Get slider values
        try:
            self.slider_vars["Fit bound LOW"].set(int(self.low_bound_entry.get()))
            self.slider_vars["Fit bound HIGH"].set(int(self.high_bound_entry.get()))
            self.slider_vars["Spectrum cutoff"].set(int(self.cutoff_entry.get()))
        except:
            pass

        # Calculate the edge position and stopping power
        _, self.edge_dict = Calibrate.get_edge_pos(self.campaign.measurements[self.meas_name], 'hTC',
                               fit_bounds=[int(self.slider_vars["Fit bound LOW"].get()), int(self.slider_vars["Fit bound HIGH"].get())],
                               check_plot=False)
        
        _, _, self.stop_pow_dict = Calibrate.get_stopping_power(self.campaign.measurements[self.meas_name], plot=False)                                                           

        # Close the initial tab
        self.notebook.forget(0)

        # Create the edge plot
        self.add_plot_tab(Calibrate._plot_edge_calibration, self.edge_dict, title="Edge Spectrum")

        # Create the stopping power plot
        self.add_plot_tab(Calibrate._plot_stopping_power, self.stop_pow_dict, title="Stopping Power")

        # Update the text box
        self.update_textbox(self.edge_dict, self.stop_pow_dict)

        # Enable buttons
        self.rerun_button.config(state='normal')
        self.save_button.config(state='normal')

        # Lock the run button
        self.run_button.config(state='disabled')

    def check_file_entry(self, *args):
        # Lock entries for already linearized files
         
        if self.entry_file_path.get().endswith('.Spe'):
            self.entry_linearization_path.config(state='normal')
            self.browse_linearization_path.config(state='normal')

            self.detector_dropdown.config(state='normal')
            self.material_dropdown.config(state='normal')
        else:
            self.entry_linearization_path.config(state='disabled')
            self.browse_linearization_path.config(state='disabled')
            self.detector_dropdown.config(state='disabled')
            self.material_dropdown.config(state='disabled')

    def save_values(self):
        # Save the values from the textbox to a file
        save_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])

        #Give it a default name if the user cancels
        if save_path == '':
            save_path = f"{self.meas_name}_edge_calibration_{self.particle}_on_{self.detector}.txt"

        with open(save_path, 'w') as file:
            file.write(self.textbox.get(1.0, tk.END))

    def update_textbox(self, edge_dict, stop_pow_dict):
        # Update the text box with the results of the edge calibration
        self.textbox.delete(1.0, tk.END)

        self.textbox.insert(tk.END, f"Calibration Results: {self.particle} on {self.detector}\n")
        self.textbox.insert(tk.END, "------------------------\n")
        self.textbox.insert(tk.END, f"Scaling factor lineal energy: {stop_pow_dict['ymax']/edge_dict['hTC']:.3f} keV/um/mV\n")
        self.textbox.insert(tk.END, f"Maximum energy loss: y_max = {stop_pow_dict['ymax']:.3f} keV/um and LET_max = {stop_pow_dict['Lmax']:.3f} keV/um\n")
        self.textbox.insert(tk.END, "hTC taken as reference, Scaling from mV to keV/um\n")
        self.textbox.insert(tk.END, "------------------------\n")

        self.textbox.insert(tk.END, "Edge Results\n")
        self.textbox.insert(tk.END, "------------------------\n")
        
        self.textbox.insert(tk.END, f"Edge position hTC: {edge_dict['hTC']:.3f} mV\n")
        self.textbox.insert(tk.END, f"Edge position hFlex: {edge_dict['hFlex']:.3f} mV\n")
        self.textbox.insert(tk.END, f"Edge position hD: {edge_dict['hDD']:.3f} mV\n")
        self.textbox.insert(tk.END, f"Fit bounds (channel): {edge_dict['fit_bounds']}\n")
        self.textbox.insert(tk.END, "------------------------\n")

    def add_initial_tab(self):
        initial_tab = ttk.Frame(self.notebook)
        self.notebook.add(initial_tab, text="Welcome")

        # Einen Beispiel-Plot erstellen
        fig = Figure(figsize=(7, 4.5), dpi=100)
        ax = fig.add_subplot(111)
        welcome_message(ax, "Edge Calibration Tool")
        canvas = FigureCanvasTkAgg(fig, master=initial_tab)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)

        # NavigationToolbar hinzufügen (optional)
        toolbar = NavigationToolbar2Tk(canvas, initial_tab, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky=tk.W+tk.E)
        toolbar.update()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)

    def add_plot_tab(self, plot_function, info_dict, title="New Tab"):
        # Erstellt einen neuen Tab im Notizbuch
        new_tab = ttk.Frame(self.notebook)
        self.notebook.add(new_tab, text=title)
        
        # Erstellt eine Figure und einen Axes-Subplot für den Plot
        fig = Figure(figsize=(7, 4.5), dpi=100)
        ax = fig.add_subplot(111)
        
        # Ruft die externe Plot-Funktion auf
        plot_function(info_dict, ax)
        
        # Integriert den Plot in den Tab
        canvas = FigureCanvasTkAgg(fig, master=new_tab)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)
        
        # Fügt optional eine NavigationToolbar hinzu
        toolbar = NavigationToolbar2Tk(canvas, new_tab, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky=tk.W+tk.E)
        toolbar.update()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)

    def on_slider_change(self, value, name):
        # Check if left slider is smaller than right slider and left slider bigger than cutoff
        # If not set to initial value and throw popup messagebox

        int_value = int(float(value))
        self.slider_vars[name].set(int_value)

        if name == "Fit bound LOW":
            if int_value >= self.slider_vars["Fit bound HIGH"].get():
                self.slider_vars[name].set(1)
                self.slider_vars["Fit bound HIGH"].set(4095)
                messagebox.showwarning("Warning", "Left slider has to be smaller than Fit bound HIGH")
        elif name == "Fit bound HIGH":
            if int_value <= self.slider_vars["Fit bound LOW"].get():
                self.slider_vars[name].set(4095)
                self.slider_vars["Fit bound LOW"].set(1)
                messagebox.showwarning("Warning", "Right slider has to be bigger than Fit bound LOW")
        elif name == "Spectrum cutoff":
            if int_value >= self.slider_vars["Fit bound LOW"].get():
                self.slider_vars[name].set(0)
                messagebox.showwarning("Warning", "Cutoff channel has to be smaller than Fit bound LOW")

root = tk.Tk()
app = GUI_Calibration(root)

# Fixed window size
root.resizable(False, False)
root.mainloop()