import sys

# Whereever you stored the repository project folder
projekt_dir = 'C:/Users/knopf/Desktop/APFELmuS'
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import Linearization
import numpy as np

""" This is a template demonstrating the generation of a calibraiton file from pulser measurements (Channel vs mV). """

#----------------------------------------------------

pulser_calibration = ''

pulser_spectrum = ''
pulse_mV_list = np.linspace(20, 250, 24)

output_file_name = 'silicon_linearization_med_gain'
output_file_path = 'output'

#----------------------------------------------------

# This is it! The cutoff may be necessary to avoid fitting the initial noise peak (See testplot)
# At the moment: Point to point interpolation 'intrpol' or single linear fit 'linear' are available
Linearization.get_linearization(output_file_name, output_file_path, pulser_calibration, pulse_mV_list, pulser_spectrum, method='interpol', testplot=True, cutoff=75)