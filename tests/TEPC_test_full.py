import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" Demonstration of a full analysis workflow for a TEPC measurement campaign """

plot = True

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()

#Have a look at the stuff in the folder
campaign1.get_files_with_format('tests/testdata')

#Read in a whole folder --> campaign dict
campaign1.read_folder('tests/testdata')

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'Raw MCA spectrum: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'Raw MCA spectrum: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'Raw MCA spectrum: High gain')

#----------------------------------------------------

Calibrate.get_chord_length(campaign1.measurements['sandra_low'], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements['sandra_mid'], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements['sandra_high'], 'slab', 10)

edge_pos_mV, _ = Calibrate.get_edge_pos(campaign1.measurements['sandra_low'], marker_point='hTC', fit_bounds=[20,150], check_plot=True)
#ymax, Lmax, _ = Calibrate.get_stopping_power(campaign1.measurements['sandra_low'], 'mean', 'SRIM', 0.01, plot=True) #TEPC stopping powers not yet added to ressources
ymax = 200

Calibrate.scale_energy_axis(campaign1.measurements['sandra_low'], edge_pos_mV, ymax, 'lineal', 'mean')
Calibrate.scale_energy_axis(campaign1.measurements['sandra_mid'], edge_pos_mV, ymax, 'lineal', 'mean')
Calibrate.scale_energy_axis(campaign1.measurements['sandra_high'], edge_pos_mV, ymax, 'lineal', 'mean')

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'Lineal Energy axis: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'Lineal Energy axis: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'Lineal Energy axis: High gain')

Spectrum.cutoff(campaign1.measurements['sandra_low'], channels=14)
Spectrum.cutoff(campaign1.measurements['sandra_mid'], channels=18)
Spectrum.cutoff(campaign1.measurements['sandra_high'], channels=140)

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'Spectrum with cutoff: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'Spectrum with cutoff: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'Spectrum with cutoff: High gain')

# The highest gain may be extrapolated towards zero (only possible for f(y))
Spectrum.probability_function(campaign1.measurements['sandra_high'], 'F')
Spectrum.probability_density(campaign1.measurements['sandra_high'])
Spectrum.extrapolate(campaign1.measurements['sandra_high'], 'linear')
Spectrum.dose_density(campaign1.measurements['sandra_high']) #d(y) calculated after extrapolation
Spectrum.weighted_probability_density(campaign1.measurements['sandra_high'])

Spectrum.probability_function(campaign1.measurements['sandra_low'], 'D')
Spectrum.probability_function(campaign1.measurements['sandra_mid'], 'D')
Spectrum.probability_density(campaign1.measurements['sandra_low'])
Spectrum.probability_density(campaign1.measurements['sandra_mid'])
Spectrum.weighted_probability_density(campaign1.measurements['sandra_low'])
Spectrum.weighted_probability_density(campaign1.measurements['sandra_mid'])

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'yd(y): Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'yd(y): Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'yd(y): High gain')

Spectrum.logarithmic_binning(campaign1.measurements['sandra_low'], 60)
Spectrum.logarithmic_binning(campaign1.measurements['sandra_mid'], 60)
Spectrum.logarithmic_binning(campaign1.measurements['sandra_high'], 60)

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'After logarithmic binning: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'After logarithmic binning: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'After logarithmic binning: High gain')

#An exception is the merge function for which you have to create a new entry in the measurement dict for the returned merged spectrum
sp_to_merge = [campaign1.measurements['sandra_low'], campaign1.measurements['sandra_mid'], campaign1.measurements['sandra_high']]
merged_spectrum = Spectrum.merge_spectra(sp_to_merge)

#This is done with the add_measurement method
#You also have to assign an identifier for this new spectrum
campaign1.add_measurement('sandra_merged', merged_spectrum)
Output.plot_single(campaign1.measurements['sandra_merged'], 'Merged Spectrum')

Spectrum.normalize_spectrum(campaign1.measurements['sandra_merged'])
Output.plot_single(campaign1.measurements['sandra_merged'], 'Normalized Spectrum')

#-------------------------------------------------

#In the end you can either get some parameters from the spectra
y_D = Output.get_mean_value(campaign1.measurements['sandra_merged'])

#Or plot the current state of the data
Output.plot_single(campaign1.measurements['sandra_merged'], output_path='plots/single.png', mean=y_D)