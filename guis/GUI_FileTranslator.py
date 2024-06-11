__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import os
import sys

import ctypes

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import FileTranslator

class GUI_FileTranslator:

    myappid = 'APFELMUS_FileTranslator' # arbitrary string
    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

    def __init__(self, root):
        self.root = root
        self.root.title("APFELmuS - File Translator")
        self.create_widgets()

        # Load the APFELmuS logo
        root.iconbitmap("ressources/logo.ico")

    def create_widgets(self):

        # Frame for file browsing
        self.file_frame = ttk.Frame(root, borderwidth=2, relief="solid")
        self.file_frame.grid(row=0, column=0, columnspan=3, sticky=tk.W+tk.E)

        # Enter file to translate
        ttk.Label(self.file_frame, text="File(s) to translate").grid(row=0, column=0, sticky=tk.E)
        ttk.Label(self.file_frame, text=" (*.Spe):").grid(row=0, column=1, sticky=tk.E)
        self.entry_file_path = ttk.Entry(self.file_frame)
        self.entry_file_path.grid(row=0, column=2, sticky=tk.W)
        self.browse_file_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_spe_file)
        self.browse_file_path.grid(row=0, column=3, sticky=tk.W)

        #Browse for calibration file
        ttk.Label(self.file_frame, text="Linearization").grid(row=1, column=0, sticky=tk.E)
        ttk.Label(self.file_frame, text=" (*.csv):").grid(row=1, column=1, sticky=tk.E)
        self.entry_linearization_path = ttk.Entry(self.file_frame)
        self.entry_linearization_path.grid(row=1, column=2, sticky=tk.W)
        self.browse_linearization_path = ttk.Button(self.file_frame, text="Browse", command=self.browse_linearization_file)
        self.browse_linearization_path.grid(row=1, column=3, sticky=tk.W)

        # Empty row
        ttk.Label(self.file_frame, text="").grid(row=2, column=0, sticky=tk.W)

        # Dropdowns for particle and detector
        ttk.Label(self.file_frame, text="Detector:").grid(row=3, column=0, sticky=tk.E)
        self.detector_var = tk.StringVar()
        self.detector_var.set("silicon")
        self.detector_dropdown = ttk.OptionMenu(self.file_frame, self.detector_var, "silicon", "silicon", "diamond", "sic")
        self.detector_dropdown.grid(row=3, column=1, sticky=tk.W)

        ttk.Label(self.file_frame, text="Particle:").grid(row=3, column=2, sticky=tk.E)
        self.particle_var = tk.StringVar()
        self.particle_var.set("proton")
        self.material_dropdown = ttk.OptionMenu(self.file_frame, self.particle_var, "proton", "proton", "carbon", "helium")
        self.material_dropdown.grid(row=3, column=3, sticky=tk.W)

        # Empty row
        ttk.Label(root, text="").grid(row=1, column=0, sticky=tk.W)

        # Frame for the buttons
        self.button_frame = ttk.Frame(root, borderwidth=2, relief="solid")
        self.button_frame.grid(row=1, column=0, columnspan=2, sticky=tk.W+tk.E)

        # Translate button
        self.translate_button = ttk.Button(self.button_frame, text="Translate", command=self.translate_file)
        self.translate_button.grid(row=0, column=0, sticky=tk.W, ipadx=50, ipady=15)

        # Reset button
        self.reset_button = ttk.Button(self.button_frame, text="Reset", command=self.reset)
        self.reset_button.grid(row=0, column=1, sticky=tk.W, ipadx=50, ipady=15)

    def translate_file(self):
        
        info_dict = {'DETECTOR': self.detector_var.get(), 'GAIN': 'None', 'PARTICLE': self.particle_var.get()}

        # Popup window for output path
        if self.entry_file_path.get() == "" or self.entry_linearization_path.get() == "":
            messagebox.showerror("Error", "Please select a file and a linearization file.")

        # Single file translation
        if len(self.entry_file_path.get().split(" ")) == 1:
            output_file_path = filedialog.asksaveasfilename(filetypes=[("Spectrum files", "*.MCA")],
                                                            initialfile=f"{self.entry_file_path.get().split('/')[-1].split('.Spe')[0]}")

            output_name = output_file_path.split("/")[-1]
            output_file_path = output_file_path.replace(output_name, "")
            output_names = [output_name]
        # Multiple file translation
        else:
            output_file_path = filedialog.askdirectory()

            # List of file names
            output_names = []
            for file_path in self.entry_file_path.get().split(" "):
                output_names.append(file_path.split("/")[-1].split(".Spe")[0])

            # Set default in case no directory was selected
            if output_file_path == "":
                output_file_path = os.path.dirname(os.path.abspath(__file__))

        for i, file_path in enumerate(self.entry_file_path.get().split(" ")):

            try:
                FileTranslator.translate_MAESTRO_file(file_path, self.entry_linearization_path.get(),
                                                        output_file_path, output_names[i], info_dict)
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred: {e}")
                return
        # Show success message (if no error occurred)
        messagebox.showinfo("Success", f"File translated successfully.\n\nOutput file: {output_file_path}")

    def reset(self):
        self.entry_file_path.delete(0, tk.END)
        self.entry_linearization_path.delete(0, tk.END)
        self.detector_var.set("silicon")
        self.particle_var.set("proton")
    
    def browse_spe_file(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("MAESTRO files", "*.Spe")], multiple=True)
        self.entry_file_path.delete(0, tk.END)
        self.entry_file_path.insert(0, file_path)
    
    def browse_linearization_file(self):
        script_directory = os.path.dirname(os.path.abspath(__file__))
        file_path = filedialog.askopenfilename(initialdir=script_directory, filetypes=[("Linearization files", "*.csv")])
        self.entry_linearization_path.delete(0, tk.END)
        self.entry_linearization_path.insert(0, file_path)

root = tk.Tk()
app = GUI_FileTranslator(root)

# Fixed window size
root.resizable(False, False)
root.mainloop()