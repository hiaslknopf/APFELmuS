import tkinter as tk
import subprocess

import ctypes

def launch_gui(script):
    subprocess.Popen(["python", script])

# Create a custom style for the buttons
button_style = {
    "font": ("Arial", 12),
    "bg": "lightblue",
    "fg": "black",
    "relief": "raised",
    "borderwidth": 2,
    "width": 30,
    "height": 2
}

def help_me():
    print("Help yourself!")

myappid = 'APFELMUS_Launcher' # arbitrary string
ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid)

root = tk.Tk()
root.title("APFELmuS Launcher")
root.geometry("300x300")

# Load the APFELmuS logo
root.iconbitmap("ressources/logo.ico")

# Buttons to launch the different GUIs
button_gui1 = tk.Button(root, text="Analyze Spectra\n(From Pulse Height -> uDos Spectrum)", command=lambda: launch_gui('guis/GUI_Analysis.py'), **button_style)
button_gui1.pack(pady=5)

button_gui2 = tk.Button(root, text="Edge Calibration\n(Get scaling factor or Edge position)", command=lambda: launch_gui('guis/GUI_Calibration.py'), **button_style)
button_gui2.pack(pady=5)

button_gui3 = tk.Button(root, text="Create Calibration File\n(From a pulse .Spe file)", command=lambda: launch_gui('guis/GUI_Linearization.py'), **button_style)
button_gui3.pack(pady=5)

button_gui4 = tk.Button(root, text="Attach Linearization\n(To a MAESTRO .Spe file)", command=lambda: launch_gui('guis/GUI_FileTranslator.py'), **button_style)
button_gui4.pack(pady=5)

button_gui5 = tk.Button(root, text="Help", command=help_me, **button_style)
button_gui5.pack(pady=5)

if __name__ == "__main__":
    root.mainloop()