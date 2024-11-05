import os
import sys

projekt_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output

""" Demonstration of the multiplot feature. Read in folder and show all ROOT files in one plot. """

#----------------------------------------------------

campaign1 = MicroDosimetry()
campaign1.read_folder('tests/testdata/spectra')
files = campaign1.get_files_with_format('tests/testdata/spectra')

#----------------------------------------------------

#Remove all MCA files and just show root (for example)
for file in files:
    if file[-3:] == 'MCA'or file[-3:] == 'csv':
        print(file[:-4])
        MicroDosimetry.delete_measurement(campaign1, file[:-4])

#----------------------------------------------------

Output.plot_campaign(campaign1, 'Multiplot Test', 'output')