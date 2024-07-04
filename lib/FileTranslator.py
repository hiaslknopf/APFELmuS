__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import os
import pandas as pd
from dateutil import parser
from datetime import datetime

""" This is a translator tool that converts spectrum files (.Spe) and linearization files (.csv) into the current version of MCA files (.MCA) for MicroDosimetry analysis.
It is currently designed to work with ORTEC MAESTRO .Spe output files.
"""

def _read_MAESTRO_file(filepath):
    """ Read MAESTRO file into dataframes
        One for the metadata and one for the spectrum data

    Args:
        filepath
    Returns:
        header_dict: Dict containing header information (metadata)
        data_df: Pandas dataframe containing spectrum information   
    """

    #Strip info from MAESTRO file header
    header_df = pd.DataFrame(pd.read_csv(filepath, header=None, nrows=12).values.reshape(-1,2))
    
    date = parser.parse(header_df.iloc[3][1])
    num_channels = int((header_df.iloc[5][1])[2:]) + 1
    live_time = int(header_df.iloc[4][1].split(' ')[0])
    real_time = int(header_df.iloc[4][1].split(' ')[1])

    header_dict = {'DATE': date, 'NUM_CHANNELS': num_channels,
                    'LIVE_TIME': live_time, 'REAL_TIME': real_time}

    #Get spectrum information from MAESTRO file
    data_df = pd.read_csv(filepath, header=None, skiprows=12, nrows = num_channels, names=['COUNTS'])

    return data_df, header_dict, num_channels

def __read_Linearization_file(filepath, num_channels):
    """ Read linearization file into dataframe. Expects a tab separated csv with header + empty line with channel number as index (1 or 3 columns).

    Args:
        filepath
    Returns:
        linearization_df: Pandas dataframe containing linearization information   
    """

    #For comma separated csv with header + empty line
    linearization_df = pd.read_csv(filepath, sep=',')

    if len(linearization_df.columns) == 4:
        linearization_df.columns = ['CHANNEL', 'LOW [mV]', 'MID [mV]', 'HIGH [mV]']
    else:
        linearization_df.columns = ['CHANNEL', 'INPUT [mV]']

    return linearization_df

def translate_MAESTRO_file(input_MAESTRO:str, input_linearization:str, output_path:str, name:str, info_dict:dict):

    """ Translates a MAESTRO .Spe file together with a linearization file of an arbitrary format into an .MCA file for further use with the MicroDosimetry package.

    Args:
        input_MAESTRO: The .Spe file to translate
        input_linearization: The linearization file to be attached to this measurement
        output_path: The location, the .MCA file is to be stored
        name: Name of the file to be stored
        info_dict: Dictionary of information for the MCA file header {'DETECTOR', 'GAIN', 'PARTICLE'}; GAIN can be set 'None' if only one amplifier was used
    """

    #TODO: Rethink the attachment of linearization files to the measurements. Maybe theres a better way to do this.

    input_MAESTRO = os.path.normpath(input_MAESTRO)
    input_linearization = os.path.normpath(input_linearization)
    output_path = os.path.normpath(output_path)

    #Read in MAESTRO .Spe and linearization .csv
    MAESTRO_df, header_dict, num_channels = _read_MAESTRO_file(input_MAESTRO)
    linearization_df = __read_Linearization_file(input_linearization, num_channels)

    #Check if output path exists
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Remove folder from name
    name = name.split('\\')[-1]
    if name.split('.')[-1] == 'Spe':
        name = name.split('.Spe')[0]
    print(name)

    #Write an MCA file (csv, sep=';')
    with open(f'{output_path}/{name}.MCA', 'w') as f:

        #MCA file header
        l1 = f"DATE;{datetime.strftime(header_dict['DATE'], '%m/%d/%y %H:%M:%S')}\n"
        l2 = f"NUM_CHANNELS;{header_dict['NUM_CHANNELS']}\n"
        l3 = f"DETECTOR;{info_dict['DETECTOR']}\n"
        l4 = f"GAIN;{info_dict['GAIN']}\n"
        l5 = f"PARTICLE;{info_dict['PARTICLE']}\n"
        l6 = f"LIVE_TIME;{header_dict['LIVE_TIME']}\n"
        l7 = f"REAL_TIME;{header_dict['REAL_TIME']}\n"
        l8 = "\n"

        l9 = "CHANNEL;COUNTS;mV\n"

        f.writelines([l1, l2, l3, l4, l5,l6, l7, l8, l9])

        #Write the data
        for i in range(num_channels):
            if info_dict['GAIN'] == 'None':
                f.write(f"{i+1};{MAESTRO_df.iloc[i].values[0]};{linearization_df.iloc[i]['INPUT [mV]']}\n")
            else:
                f.write(f"{i+1};{MAESTRO_df.iloc[i].values[0]};{linearization_df.iloc[i][info_dict['GAIN'] + ' [mV]']}\n")