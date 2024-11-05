import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" Read in a folder, calibrate, evaluate and normalize one spectrum """

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()

#Read in a whole folder --> campaign dict
campaign1.read_folder('tests/testdata/spectra')

exp_file = 'exp_carbon_data_2'
cal_file = 'exp_carbon_edge'

#----------------------------------------------------

# Calibrate
Calibrate.get_chord_length(campaign1.measurements[exp_file], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements[cal_file], 'slab', 10)

edge_pos_mV, _ = Calibrate.get_edge_pos(campaign1.measurements[cal_file], marker_point='hDD', fit_bounds=[1800,2400], check_plot=True)
ymax, Lmax, _ = Calibrate.get_stopping_power(campaign1.measurements[cal_file], 'mean', 'SRIM', 0.01, plot=True)

Calibrate.scale_energy_axis(campaign1.measurements[exp_file], edge_pos_mV, ymax, 'lineal', 'mean')

# Calculate uDos representation
Spectrum.cutoff(campaign1.measurements[exp_file], channels=14)

Spectrum.probability_function(campaign1.measurements[exp_file], 'F')
Spectrum.probability_density(campaign1.measurements[exp_file])
Spectrum.extrapolate(campaign1.measurements[exp_file], 'linear')
Spectrum.dose_density(campaign1.measurements[exp_file])
Spectrum.weighted_probability_density(campaign1.measurements[exp_file])
Spectrum.logarithmic_binning(campaign1.measurements[exp_file], 60)

# Plotting
Output.plot_single(campaign1.measurements[exp_file], 'Binned')

# Normalization
Spectrum.normalize_spectrum(campaign1.measurements[exp_file])
Output.plot_single(campaign1.measurements[exp_file], 'Normalized')