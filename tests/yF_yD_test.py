import os
import sys

import numpy as np

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" Calculate f(y) to yd(y) distribution and plot spectrum + mean values """

plot = True

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_file('tests/testdata/spectra/exp_carbon_data_1.MCA')

if plot:
    Output.plot_single(campaign1.measurements['exp_carbon_data_1'], 'Raw Spectrum')

#----------------------------------------------------
# I happen to know these values
edge_pos_mV = 416.4087
ymax = 148.33
Calibrate.get_chord_length(campaign1.measurements['exp_carbon_data_1'], 'slab', 10)
Calibrate.scale_energy_axis(campaign1.measurements['exp_carbon_data_1'], edge_pos_mV, ymax)

# Noise cutoff for proper means calculation
Spectrum.cutoff(campaign1.measurements['exp_carbon_data_1'], lineal_energy=4.0)

if plot:
    Output.plot_single(campaign1.measurements['exp_carbon_data_1'], 'Calibrated Spectrum')

#----------------------------------------------------

Spectrum.probability_function(campaign1.measurements['exp_carbon_data_1'], 'F') # F(y)
Spectrum.probability_density(campaign1.measurements['exp_carbon_data_1']) # f(y)
Spectrum.normalize_linear_spectrum(campaign1.measurements['exp_carbon_data_1']) # f(y) normalized to area 1

# Calculate means
y_F, y_D = Output.means_from_fy(campaign1.measurements['exp_carbon_data_1'])

print(f'Frequency mean y_F:  {y_F:.2f} keV/um')
print(f'Dose mean y_D: {y_D:.2f} keV/um')

if plot:
    Output.plot_single(campaign1.measurements['exp_carbon_data_1'], y_F=y_F, y_D=y_D)