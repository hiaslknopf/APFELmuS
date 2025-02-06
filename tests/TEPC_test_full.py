import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" Demonstration of a full analysis workflow for a TEPC measurement campaign """

plot = False

#----------------------------------------------------

#Set up a new measurement campaign
campaign1 = MicroDosimetry()

#Have a look at the stuff in the folder
campaign1.get_files_with_format('tests/testdata/spectra')

#Read in a whole folder --> campaign dict
campaign1.read_folder('tests/testdata/spectra')

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'Raw MCA spectrum: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'Raw MCA spectrum: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'Raw MCA spectrum: High gain')

#----------------------------------------------------

Calibrate.get_chord_length(campaign1.measurements['sandra_low'], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements['sandra_mid'], 'slab', 10)
Calibrate.get_chord_length(campaign1.measurements['sandra_high'], 'slab', 10)

# I happen to know these values
edge_pos_mV = 416.4087
ymax = 148.33

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
Spectrum.logarithmic_binning(campaign1.measurements['sandra_high'], 60)

Spectrum.probability_density(campaign1.measurements['sandra_high'])
#Spectrum.extrapolate(campaign1.measurements['sandra_high'], 'linear') # TODO: Extrapolation has to be fixed; Doesnt really work for log binned data
Spectrum.dose_density(campaign1.measurements['sandra_high']) #d(y) calculated after extrapolation

Spectrum.probability_function(campaign1.measurements['sandra_low'], 'D')
Spectrum.probability_function(campaign1.measurements['sandra_mid'], 'D')


Spectrum.logarithmic_binning(campaign1.measurements['sandra_low'], 60)
Spectrum.logarithmic_binning(campaign1.measurements['sandra_mid'], 60)

Spectrum.probability_density(campaign1.measurements['sandra_low'])
Spectrum.probability_density(campaign1.measurements['sandra_mid'])

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'Spectrum with cutoff: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'Spectrum with cutoff: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'Spectrum with cutoff: High gain')


Spectrum.weighted_probability_density(campaign1.measurements['sandra_low'])
Spectrum.weighted_probability_density(campaign1.measurements['sandra_mid'])
Spectrum.weighted_probability_density(campaign1.measurements['sandra_high'])

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'yd(y): Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'yd(y): Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'yd(y): High gain')

if plot:
    Output.plot_single(campaign1.measurements['sandra_low'], 'After logarithmic binning: Low gain')
    Output.plot_single(campaign1.measurements['sandra_mid'], 'After logarithmic binning: Mid gain')
    Output.plot_single(campaign1.measurements['sandra_high'], 'After logarithmic binning: High gain')

#Output.csv_output(campaign1.measurements['sandra_low'], 'tests/testdata')
#Output.csv_output(campaign1.measurements['sandra_mid'], 'tests/testdata')
#Output.csv_output(campaign1.measurements['sandra_high'], 'tests/testdata')

#An exception is the merge function for which you have to create a new entry in the measurement dict for the returned merged spectrum
# This can be done using all 3 gains or only 2 of them
sp_to_merge = [campaign1.measurements['sandra_mid'], campaign1.measurements['sandra_high'], campaign1.measurements['sandra_low']]
overlap_regions = [[2.5,7], [8.5,20]]
scaling_factors = [0.6, 0.35]
stitching_points = [5, 18]
merged_spectrum = Spectrum.merge_spectra(sp_to_merge, overlap_regions, scaling_factors, stitching_points, testplots=True)

sp_to_merge = [campaign1.measurements['sandra_high'], campaign1.measurements['sandra_low']]
overlap_regions = [[4.5, 8]]
scaling_factors = [0.6]
stitching_points = [6.5]
merged_spectrum, _ = Spectrum.merge_spectra(sp_to_merge, overlap_regions, scaling_factors, stitching_points, testplots=True)

#This is done with the add_measurement method
#You also have to assign an identifier for this new spectrum
campaign1.add_measurement('sandra_merged', merged_spectrum)
Output.plot_single(campaign1.measurements['sandra_merged'], 'Merged Spectrum')

#-------------------------------------------------

#In the end you can calculate mean values
Spectrum.normalize_log_spectrum(campaign1.measurements['sandra_merged'])

#And plot the final data
Output.plot_single(campaign1.measurements['sandra_merged'], output_path='plots/single.png')