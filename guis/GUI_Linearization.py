__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import os
import sys
import pandas as pd
import numpy as np

import ctypes

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import Linearization, FileTranslator
from guis.welcome_message import welcome_message

#TODO: Theres some weird bug with deleting the old plots, Just hit buttons twice
#TODO: Its a pain in the ass to enter the mV values manually, maybe you find another way to do this

class GUI_Linearization:

    myappid = 'APFELMUS_Linearization' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        self.root = root
        self.root.title("APFELmuS - Calibration File Generator")
        self.create_widgets()

        # Load the APFELmuS logo
        root.iconbitmap("ressources/logo.ico")

    def create_widgets(self):
        self.left_panel = ttk.Frame(root)
        self.left_panel.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        #Subframe for pulser spectrum and pulser cal
        self.file_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.file_frame.grid(row=0, column=0, columnspan=3, sticky=tk.W+tk.E)

        #Browse for file
        ttk.Label(self.file_frame, text="Pulser File (*.Spe):").grid(row=0, column=0, sticky=tk.W)
        self.entry_pulser_file_path = ttk.Entry(self.file_frame)
        self.entry_pulser_file_path.grid(row=0, column=1, sticky=tk.W)
        self.browse_pulser_file_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_pulser_spec)
        self.browse_pulser_file_path.grid(row=0, column=2, sticky=tk.W)

        #Checkbutton for using calibration file
        self.use_calibration = tk.BooleanVar()
        self.use_calibration.set(False)
        self.checkbutton_use_calibration = ttk.Checkbutton(self.file_frame, text="Use calibration file\n(otherwise 1:1 response assumed)", variable=self.use_calibration, command=self.toggle_calibration)
        self.checkbutton_use_calibration.grid(row=1, column=0, sticky=tk.W)

        #Browse for calibration file
        self.entry_cal_label = ttk.Label(self.file_frame, text="Pulser Cal (*.csv):", state='disabled').grid(row=2, column=0, sticky=tk.W)
        self.entry_cal_path = ttk.Entry(self.file_frame, state='disabled')
        self.entry_cal_path.grid(row=2, column=1, sticky=tk.W)
        self.browse_cal_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_pulser_cal, state='disabled')
        self.browse_cal_path.grid(row=2, column=2, sticky=tk.W)

        # Empty line between frames
        ttk.Label(self.left_panel, text="").grid(row=1, column=0, sticky=tk.W)

        # New frame for mV list, cutoff and fit method
        self.specs_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.specs_frame.grid(row=2, column=0, columnspan=3, sticky=tk.W+tk.E)

        # Entry for mV list (textbox)
        ttk.Label(self.specs_frame, text="mV list (comma separated):").grid(row=0, column=0, sticky=tk.W)
        self.entry_mV_list = tk.Text(self.specs_frame, height=5, width=30)
        self.entry_mV_list.grid(row=0, column=1, sticky=tk.W)

        # Button to load mV list from file
        self.browse_mV_list_button = ttk.Button(self.specs_frame, text="Load from file", command=self.browse_mV_list)
        self.browse_mV_list_button.grid(row=2, column=0, sticky=tk.W)
        
        # Empty row
        ttk.Label(self.specs_frame, text="").grid(row=3, column=0, sticky=tk.W)

        # Slider and entry for cutoff
        ttk.Label(self.specs_frame, text="Spectrum cutoff:").grid(row=4, column=0, sticky=tk.W)

        self.slider_var = tk.IntVar()
        self.slider_cutoff = ttk.Scale(self.specs_frame, from_=0, to=250, orient=tk.HORIZONTAL, variable=self.slider_var)
        self.slider_cutoff.grid(row=4, column=0, sticky=tk.W+tk.E)
        
        self.entry_cutoff = ttk.Entry(self.specs_frame, width=5, textvariable=self.slider_var)
        self.entry_cutoff.grid(row=4, column=1, sticky=tk.W)

        self.slider_cutoff.bind("<ButtonRelease-1>", self.on_slider_cutoff_change)
        self.entry_cutoff.bind("<Return>", self.on_entry_cutoff_change)
        
        # Empty row
        ttk.Label(root, text="").grid(row=3, column=0, sticky=tk.W)

        # Frame for buttons
        self.button_frame = ttk.Frame(self.left_panel, borderwidth=2, relief="solid")
        self.button_frame.grid(row=4, column=0, columnspan=3, sticky=tk.W+tk.E)

        # Bottom Buttons
        self.run_button = ttk.Button(self.button_frame, text="RUN", command=self.run)
        self.run_button.grid(row=0, column=0, sticky=tk.W+tk.E, ipadx=50, ipady=10)
        self.rerun_button = ttk.Button(self.button_frame, text="Rerun Peak Finding", command=self.rerun_peaks, state='disabled')
        self.rerun_button.grid(row=0, column=1, sticky=tk.W+tk.E, ipadx=20, ipady=10)
        self.save_button = ttk.Button(self.button_frame, text="Save", command=self.save_values, state='disabled')
        self.save_button.grid(row=1, column=0, sticky=tk.W+tk.E, ipadx=20, ipady=10)
        self.reset_button = ttk.Button(self.button_frame, text="Reset GUI", command=self.reset_GUI)
        self.reset_button.grid(row=1, column=1, sticky=tk.W+tk.E, ipadx=20, ipady=10)

        # Right panel for plots and buttons
        self.right_panel = ttk.Frame(root, borderwidth=2, relief="solid")
        self.right_panel.grid(row=0, column=1, padx=20, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)

        # Notebook for plots
        self.notebook = ttk.Notebook(self.right_panel)
        self.notebook.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W+tk.E+tk.N+tk.S)
        self.add_initial_tab()
    
    def add_initial_tab(self):
        initial_tab = ttk.Frame(self.notebook)
        self.notebook.add(initial_tab, text="Welcome")

        # Einen Beispiel-Plot erstellen
        fig = Figure(figsize=(8, 3.5), dpi=100)
        ax = fig.add_subplot(111)
        welcome_message(ax, "Linearization File Generator")
        canvas = FigureCanvasTkAgg(fig, master=initial_tab)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)

    def browse_pulser_spec(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("MAESTRO files", "*.Spe")], multiple=True)
        self.entry_pulser_file_path.delete(0, tk.END)
        self.entry_pulser_file_path.insert(0, file_path)

    def browse_pulser_cal(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("Pulser Calibration", "*.csv")], multiple=True)
        self.entry_cal_path.delete(0, tk.END)
        self.entry_cal_path.insert(0, file_path)

    def toggle_calibration(self):
        if self.use_calibration.get():
            self.entry_cal_path.config(state='normal')
            self.browse_cal_path.config(state='normal')
        else:
            self.entry_pulser_file_path.config(state='disabled')
            self.browse_pulser_file_path.config(state='disabled')

    def add_plot_tab(self, plot_function, title="New Tab", arg_dict={}):
        # Erstellt einen neuen Tab im Notizbuch
        new_tab = ttk.Frame(self.notebook)
        self.notebook.add(new_tab, text=title)
        
        # Erstellt eine Figure und einen Axes-Subplot für den Plot
        fig = Figure(figsize=(9, 4), dpi=100)
        ax = fig.add_subplot(111)
        
        # Ruft die externe Plot-Funktion auf (mit argumenten in arg_dict)
        plot_function(ax=ax, **arg_dict)

        # Integriert den Plot in den Tab
        canvas = FigureCanvasTkAgg(fig, master=new_tab)
        canvas.draw()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)
        
        # Fügt optional eine NavigationToolbar hinzu
        toolbar = NavigationToolbar2Tk(canvas, new_tab, pack_toolbar=False)
        toolbar.grid(row=1, column=0, sticky=tk.W+tk.E)
        toolbar.update()
        canvas.get_tk_widget().grid(row=0, column=0, sticky=tk.W+tk.E+tk.N+tk.S)

    def on_slider_cutoff_change(self, event):

        # Get the slider value
        self.slider_var.set(int(self.slider_cutoff.get()))
        self.rerun_peaks()

    def on_entry_cutoff_change(self, event):
        try:
            cutoff = int(self.entry_cutoff.get())
            if cutoff < 0 or cutoff > 4095:
                raise ValueError
            self.slider_var.set(cutoff)
            self.rerun_peaks()
        except ValueError:
            messagebox.showwarning("Warning", "Please enter a valid cutoff value (0-4095)")

    def run(self):
        
        if self.entry_pulser_file_path.get() == '':
            messagebox.showwarning("Warning", "Please specify a file")
            self.reset_GUI()

        # If calibration file is used, check if it exists
        if self.use_calibration.get() and self.entry_cal_path.get() == '':
            messagebox.showwarning("Warning", "Please specify a calibration file")
            self.reset_GUI()

        # Use pulser file or one to one response
        if self.use_calibration.get():
            self.calibration_file = self.entry_cal_path.get()
        else:
            # Jump to project directory
            os.chdir(projekt_dir)
            self.calibration_file = 'ressources/1to1_response.csv'

        # Read in the pulser spectrum
        self.spectrum_df, _, self.num_channels = FileTranslator._read_MAESTRO_file(self.entry_pulser_file_path.get())

        # Run the linearization
        self.rerun_peaks()

        # Enable buttons
        self.rerun_button.config(state='normal')
        self.save_button.config(state='normal')
        self.run_button.config(state='disabled')

    def rerun_peaks(self):

        # Get the mV list, message if empty
        mV_list = self.entry_mV_list.get("1.0", tk.END).split(',')
        if mV_list == ['']:
            messagebox.showwarning("Warning", "Please specify a list of mV values")
            self.reset_GUI()

        # Get the slider value (if typed in)
        try:
            self.slider_var.set(int(self.entry_cutoff.get()))
        except:
            pass
        
        # Convert the mV list to floats and a sorted numpy array
        mV_list = np.array([float(mV) for mV in mV_list])
        self.mV_list = np.sort(mV_list)

        # Read in the pulser calibration
        self.pulser_df = Linearization._read_pulser(self.calibration_file, self.mV_list, testplot=False)

        # Fit the peaks
        self.cutoff = self.slider_var.get()
        self.peaks, self.popt_list = Linearization._fit_peaks(self.spectrum_df, testplot=False, cutoff=self.cutoff)

        # Check if all peaks were found
        if len(self.peaks) != len(self.mV_list):
            messagebox.showwarning("Warning",
                                   f"The number of peaks found ({len(self.peaks)}) does not match the number of mV values ({len(self.mV_list)})")
        
        # Get noise figure from average sigma of the gaussian fit
        self.sigma = np.mean([popt[2] for popt in self.popt_list])

        # Fit the linearization curve
        columns = ['CHANNEL', 'INPUT [mV]']
        self.linearization_df = pd.DataFrame(index=np.arange(self.num_channels), columns=columns)
        
        mV = self.pulser_df['mV_measured'].values
        lin_curve_channel, lin_curve_mV = Linearization._fit_linearization(self.peaks, mV, self.num_channels, method='interpol')

        self.linearization_df['CHANNEL'] = lin_curve_channel.astype(int)
        self.linearization_df['INPUT [mV]'] = np.round(lin_curve_mV, 3)

        # Close all open tabs
        for i in range(0, len(self.notebook.tabs())):
            self.notebook.forget(i)

        # Show the result plots
        self.add_plot_tab(Linearization._plot_peaks, title="Peaks", arg_dict={'spectrum': self.spectrum_df, 'positions': self.peaks, 'popt_list': self.popt_list})
        self.add_plot_tab(Linearization._plot_linearization, title="Linearization",
                  arg_dict={'channel': self.peaks, 'mV': mV, 'lin_curve_channel': lin_curve_channel, 'lin_curve_mV': lin_curve_mV,
                         'method': 'interpol', 'num_channels': self.num_channels})
        self.add_plot_tab(Linearization._plot_pulser_cal, title="Pulser Calibration", arg_dict={'pulser_df': self.pulser_df, 'label': 'Values'})

        messagebox.showinfo("Average noise figure", f"Average noise figure from fit (FWHM): {self.sigma:.2f} channels")

    def save_values(self):
        
        # Browse for output path
        script_directory = os.path.dirname(os.path.abspath(__file__))
        output_file_path = filedialog.asksaveasfilename(initialdir=script_directory, title="Select output file",
                                                        defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if output_file_path == '':
            print("No output path selected")
        
        #Write the calibration file
        self.linearization_df.to_csv(f'{output_file_path}', sep=',', index=False)
        print(f'Wrote linearization file: {output_file_path}')

    def reset_GUI(self):
        # Reset buttons
        self.run_button.config(state='normal')
        self.rerun_button.config(state='disabled')
        self.save_button.config(state='disabled')
        self.use_calibration.set(False)

        # Reset entry fields
        self.entry_pulser_file_path.delete(0, tk.END)
        self.entry_pulser_file_path.config(state='normal')
        self.browse_pulser_file_path.config(state='normal')
        self.entry_cal_path.delete(0, tk.END)
        self.entry_cal_path.config(state='disabled')
        self.browse_cal_path.config(state='disabled')
        self.entry_mV_list.delete("1.0", tk.END)
        self.entry_cutoff.delete(0, tk.END)
        self.slider_var.set(0)

        # Close all open tabs
        for i in range(0, len(self.notebook.tabs())):
            self.notebook.forget(i)
        
        # Add initial tab
        self.add_initial_tab()

    def browse_mV_list(self):
        file_path = filedialog.askopenfilename(initialdir=os.path.dirname(os.path.abspath(__file__)), filetypes=[("Text files", "*.txt")])
        self.get_voltages_from_file(file_path)

    def get_voltages_from_file(self, file_path):
        """ Make sure the format matches the Pulser GUI output """

        with open(file_path, 'r') as file:
            lines = file.readlines()

            for i, line in enumerate(lines):
                if 'Voltage' in line:
                    voltage_line = lines[i+1]
                    print(voltage_line)
                    break
                else:
                    voltage_line = ''
            
        if voltage_line == '':
            messagebox.showwarning("Warning", "No valid voltage sequence found")
            return

        # Remove the brackets and the newline character
        voltage_line = voltage_line.replace('[', '').replace(']', '').replace('\n', '')
        print(voltage_line)
        
        # Write the voltages into the entry field
        if self.entry_mV_list.get("1.0", tk.END) != '':
            self.entry_mV_list.delete("1.0", tk.END)
        
        self.entry_mV_list.insert("1.0", voltage_line)

root = tk.Tk()
app = GUI_Linearization(root)

# Fixed window size
root.resizable(False, False)
root.mainloop()