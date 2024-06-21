import sys

# Whereever you stored the repository project folder
projekt_dir = 'C:/Users/knopf/Desktop/APFELmuS'
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Calibrate, Spectrum, Output

""" This is a template demonstrating how to obtain a yd(y) spectrum from a given MCA spectrum. """

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_folder('analysis')

file_name = ''
scaling_factor = 3.71

plots = True

#----------------------------------------------------
"""
#Have a look at the raw data
if plots:
    Output.plot_single(campaign1.measurements[file_name])

#Cut the noise
Spectrum.cutoff(campaign1.measurements[file_name], 25)
if plots:
    Output.plot_single(campaign1.measurements[file_name])

# The detector chord length needs to be known
Calibrate.get_chord_length(campaign1.measurements[file_name], 'slab', 10, False)

# Assuming you already know the edge position and the maximum energy loss for this detector
edge_pos_mV = 300 #mV
ymax = 1112 #keV

#In the end, every spectrum can be rescaled by knowing the edge position and the maximum energy loss
Calibrate.scale_energy_axis(campaign1.measurements[file_name], edge_pos_mV, ymax, 'lineal', 'mean')
if plots:
    Output.plot_single(campaign1.measurements[file_name])

#Calculate the microdosimetry spectrum
Spectrum.probability_function(campaign1.measurements[file_name], 'D')
if plots:
    Output.plot_single(campaign1.measurements[file_name])
Spectrum.probability_density(campaign1.measurements[file_name])
if plots:
    Output.plot_single(campaign1.measurements[file_name])
Spectrum.weighted_probability_density(campaign1.measurements[file_name])
Spectrum.logarithmic_binning(campaign1.measurements[file_name], 60)
Spectrum.normalize_spectrum(campaign1.measurements[file_name])

# Finally, the yd(y) spectrum is obtained
Output.plot_single(campaign1.measurements[file_name])"""

#----------------------------------------------------

# Even faster when a conversion factor is already known

Spectrum.cutoff(campaign1.measurements[file_name], 25)
Calibrate.scale_energy_axis_with_factor(campaign1.measurements[file_name], scaling_factor)
Spectrum.probability_function(campaign1.measurements[file_name], 'D')
Spectrum.probability_density(campaign1.measurements[file_name])
Spectrum.weighted_probability_density(campaign1.measurements[file_name])
Spectrum.logarithmic_binning(campaign1.measurements[file_name], 60)	
Spectrum.normalize_spectrum(campaign1.measurements[file_name])
Output.plot_single(campaign1.measurements[file_name])

Spectrum.retrieve_original_spectrum(campaign1.measurements[file_name])

#----------------------------------------------------

""" # Plot the whole campaign
for file in campaign1.measurements:
    Calibrate.scale_energy_axis_with_factor(campaign1.measurements[file], scaling_factor)
    Spectrum.probability_function(campaign1.measurements[file], 'D')
    Spectrum.probability_density(campaign1.measurements[file])
    Spectrum.weighted_probability_density(campaign1.measurements[file])
    Spectrum.logarithmic_binning(campaign1.measurements[file], 60)	
    Spectrum.normalize_spectrum(campaign1.measurements[file])

Output.plot_campaign(campaign1, step='True') """