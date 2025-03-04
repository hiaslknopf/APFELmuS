import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from lib import Linearization

import numpy as np

""" This is a template demonstrating the generation of a linarization file from pulser measurements (Channel vs mV). """

#----------------------------------------------------

pulser_calibration = 'ressources/TGF4242_output_20231214.csv'

pulser_spectrum = 'tests/testdata/original_data/si_lin_spec_med_gain.Spe'
pulse_mV_list = np.linspace(20, 250, 24)

output_file_name = 'silicon_linearization_med_gain'
output_file_path = 'tests/testdata'

#----------------------------------------------------

# This is it! The cutoff may be necessary to avoid fitting the initial noise peak (See testplot)
# At the moment: Point to point interpolation 'intrpol' or single linear fit 'linear' are available
Linearization.get_linearization(output_file_name, output_file_path,
                                pulse_mV_list, pulser_spectrum,
                                pulser_calibration=pulser_calibration, method='interpol',
                                testplot=True, cutoff_front=75)