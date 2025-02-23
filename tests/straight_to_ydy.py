import sys

projekt_dir = r"C:\Users\knopf\Desktop\APFELmuS"
if projekt_dir not in sys.path:
    sys.path.append(projekt_dir)

from MicroDosimetry import MicroDosimetry
from lib import Output
from lib.Spectrum import ROOT_straight_to_ydy

""" Calculate yd(y) representation from ROOT files """

path1 = r"C:\Users\knopf\Dropbox\GATE_Sims_12um_SiC_Beam_Model\148.2MeV\no_pe_beam_model_proton_12um_SiC_148.2MeV_14.43cmH2O.root"
file1 = 'no_pe_beam_model_proton_12um_SiC_148.2MeV_14.43cmH2O'

output_path = r"C:\Users\knopf\Desktop"

#Set up campaign and read files
campaign = MicroDosimetry()
campaign.read_file(path1, num_root_bins=4096)

# Do some magic
ROOT_straight_to_ydy(campaign.measurements[file1], num_log_bins=60, mean_chord_length=12)

# Output
Output.plot_single(campaign.measurements[file1])
Output.csv_output(campaign.measurements[file1], output_path=output_path, name=f'{file1}_ydy')