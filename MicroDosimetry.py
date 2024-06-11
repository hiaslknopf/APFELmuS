import pandas as pd
from glob import glob
import os
from datetime import datetime
import uproot
import math
import numpy as np

from lib import Output

__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

class MicroDosimetry():
    """ Class for analysing microdosmetry measurement campaigns (consisting of one or more spectra). """

    def __init__(self):
        """ Create a dictionary for storing an identifier (filename) with a measurement dataframe"""

        self.measurements = {} #Dict for filenames + data in analysis folder

    def get_files_with_format(self, folderpath):
        """ Console printout of folder content """

        #TODO: Hübscherer Output

        folderpath = os.path.normpath(folderpath)
        files = glob('{}/*.*'.format(folderpath))
        files = [file.split(os.path.sep)[-1] for file in files]

        for i, file in enumerate(files):
            print(i+1, ':\t', file)   
        
        return files

    def read_folder(self, folderpath):
        """ Read in a whole folder of measurement files

        Args:
            folderpath: Relative path of the folder to be analysed
        """
        folderpath = os.path.normpath(folderpath)
        files = glob('{}/*.*'.format(folderpath))

        # Remove all but .MCA, .root and .csv files
        files = [file for file in files if file.endswith('.MCA') or file.endswith('.root') or file.endswith('.csv')]

        for file in files:
            self.read_file(file)

    def read_file(self, file):
        """ Read in a measurement file with a read function from the Measurement class.
        Currently this is implemented for my version of MCA files (header + 3 columns).
        If another format is to be read: Set up a new read function and enable some kind of selection.

        Args:
            file: Relative path of the file to be read in
        """

        #Create a Measurement instance
        measurement = Measurement()
    
        #Read in the information
        file = os.path.normpath(file)
        filename, extension = os.path.splitext(file)
        
        if extension == '.MCA':
            measurement.read_MCA_file(file)
        elif extension == '.root':
            measurement.read_ROOT_file(file)
        elif extension == '.csv':

            #Check if it is a uDos file
            with open(file, 'r') as f:
                lines = f.readlines()
                if 'MicroDosimetry' in lines[0]:
                    measurement.read_uDos_analysis(file)
                else:
                    pass

        #Get rid of filepath and extension
        filename = file.split(os.path.sep)[-1].split(extension)[0]

        #Fill the measurement dict
        self.measurements[filename] = measurement
        self.measurements[filename]._name = filename #Attribute of the measurement object

    def add_measurement(self, identifier, measurement):
        """ Add a new measurement to the measurement dict

        Args:
            identifier: Unique identifier to act as the dict index
            measurement: Measurement object containing the new (merged) data
        """

        self.measurements[identifier] = measurement
        self.measurements[identifier]._name = identifier #Attribute of the measurement object
        print(f'Added new measurement: {identifier}')

    def delete_measurement(self, identifier):
        """ Remove a measurement from measurement dict
        
        Args:
            identifier: Measurement identifier (filename without ending if not set manually)
        """

        del self.measurements[identifier]

    def attach_info(self, measurement, info_dict):
        """ Attach information to a measurement object.
        This can be used to change several attributes at once or to add initial information to a ROOT data.
        
        Attributes:
            particle: Particle type (proton, carbon, helium)
            detector: Detector type (silicon, diamond, sic, TEPC)
            mean_chord_length: Mean chord length (in um)
            max_chord_length: Max chord length (in um)
            num_channels: Number of channels (only for ROOT files)
            gain: Gain (only for MCA files)

        Args:
            measurement: Measurement object to be altered
            info_dict: Dictionary containing the information to be attached (_particle, _detector, _num_channels)
        """

        if measurement.gain == 'SIMULATION':
            Measurement.attach_root_info(self, measurement, info_dict)
        else:
            Measurement.attach_mca_info(self, measurement, info_dict)

class Measurement():
    """ Actual analysis class. One instance per file in analysis folder aka one entry in measurement dict """

    def __init__(self, num_channels=0, date = datetime.now(), detector=None, gain=None) -> None: 
        """ Set up a measurement object containing all info for a given spectrum """

        #General information
        self._name: str = ''
        self._num_channels:int = num_channels
        self._date: datetime = date
        self._particle: str = ''
        self._detector: str = detector
        self._gain: str = gain

        #Data information
        self._x_axis: str = ''
        self._y_axis: str = ''
        
        self._data: pd.DataFrame = pd.DataFrame()
        self._original_data: pd.DataFrame = pd.DataFrame()
        self._abs_filepath: str = ''

        #Calibration information - Set by calibration functions, Required for analysis
        self._mean_chord_length: float = 0.0
        self._max_chord_length: float = 0.0

        #Live time information
        self._live_time: float = 0.0
        self._real_time: float = 0.0

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        raise ValueError('Name cannot be changed here')

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, dataframe):
        self._data = dataframe
        print(f'Overwritten entire data!')

    @property
    def original_data(self):
        return self._original_data

    @original_data.setter
    def original_data(self, dataframe):
        raise ValueError('Original data cannot be changed')
    
    @property
    def abs_filepath(self):
        return self._abs_filepath
    
    @abs_filepath.setter
    def abs_filepath(self):
        raise ValueError('The filepath cannot be changed')

    @property
    def x_axis(self):
        return self._x_axis

    @x_axis.setter
    def x_axis(self, new_x_axis):
        self._data.rename(columns={self._x_axis: new_x_axis}, inplace=True)
        self._x_axis = new_x_axis
        print(f'Overwritten x axis to {new_x_axis}')

    @property
    def y_axis(self):
        return self._y_axis

    @y_axis.setter
    def y_axis(self, new_y_axis):
        self._data.rename(columns={self._y_axis: new_y_axis}, inplace=True)
        self._y_axis = new_y_axis
        print(f'Overwritten y axis to {new_y_axis}')

    @property
    def num_channels(self):
        return self._num_channels

    @num_channels.setter
    def num_channels(self, num_channels):
        self._num_channels = num_channels

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, date):
        raise ValueError('Date cannot be set manually')

    @property
    def detector(self):
        return self._detector

    @detector.setter
    def detector(self, detector):
        self._detector = detector
    
    @property
    def particle(self):
        return self._particle

    @particle.setter
    def particle(self, particle):
        self._particle = particle

    @property
    def gain(self):
        return self._gain

    @gain.setter
    def gain(self, gain):
        self._gain = gain
    
    @property
    def mean_chord_length(self):
        return self._mean_chord_length
    
    @mean_chord_length.setter
    def mean_chord_length(self, mean_chord_length):
        self._mean_chord_length = float(mean_chord_length)
        #print(f'Overwritten mean chord length to {mean_chord_length} um')
    
    @property
    def max_chord_length(self):
        return self._max_chord_length
    
    @max_chord_length.setter
    def max_chord_length(self, max_chord_length):
        self._max_chord_length = float(max_chord_length)
        #print(f'Overwritten max chord length to {max_chord_length} um')
    
    @property
    def live_time(self):
        return self._live_time

    @live_time.setter
    def live_time(self, live_time):
        raise ValueError('Live time cannot be set manually')

    @property
    def real_time(self):
        return self._real_time
    
    @real_time.setter
    def real_time(self, real_time):
        raise ValueError('Real time cannot be set manually')

    def __get_info_from_MCA_file(self, file):
        """ Get information from MCA file header.

        Parameters:
            file: MCA file to be read
        """

        file = os.path.normpath(file)
        header = pd.read_csv(file, header=None, nrows=7, sep=';', index_col=0)
        self.__set_attributes(df=header, file=file)

    def __get_info_from_uDos_file(self, file):
        """ Get information from uDos file footer.
        
        Parameters:
            file: uDos file to be read (csv format)
        """

        label_dict = Output.labels_dict
        
        file = os.path.normpath(file)

        # Get number of channels
        num_lines = 0
        with open(file, 'r') as f:
            lines = f.readlines()[6:-10]
            for line in lines:
                num_lines += 1

        # Read footer
        footer = pd.read_csv(file, header=None, skiprows=7+num_lines, index_col=0, sep='\t', engine='python')
        self.__set_attributes(df=footer, file=file, num_channels=num_lines)

        # Get x_axis and y_axis
        with open(file, 'r') as f:
            lines = f.readlines()
            x_axis = lines[4].split(';')[0]
            y_axis = lines[4].split(';')[1].split('\n')[0]
        f.close()

        if '$-1$' in x_axis or 'um' in x_axis:
            x_axis = x_axis.replace('-1', '$^{-1}$')
            x_axis = x_axis.replace('um', 'µm')
        if '$-1$' in y_axis or 'um' in y_axis:
            y_axis = y_axis.replace('-1', '$^{-1}$')
            y_axis = y_axis.replace('um', 'µm')

        x_axis = next((key for key, value in label_dict['x_labels'].items() if value == x_axis), None)
        y_axis = next((key for key, value in label_dict['y_labels'].items() if value == y_axis), None)

        if x_axis == None:
            raise NotImplementedError(f'In {file}: x_axis {x_axis} is not a viable option. Check labels_dict in Output.py')
        if y_axis == None:
            raise NotImplementedError(f'In {file}: y_axis {y_axis} is not a viable option. Check labels_dict in Output.py')
        
        return x_axis, y_axis
  
    def __set_attributes(self, df, file, num_channels=0):
        """ Set attributes of the measurement object from a given header or footer dataframe
            To be used by the __get_info_from ... functions """
                
        # Check if all necessary information is here
        if not 'DATE' in df.index:
            raise IOError('Date not found in file header')
    
        if num_channels == 0:
            if not 'NUM_CHANNELS' in df.index:
                raise IOError('Number of channels not found in file header')

        if not 'GAIN' in df.index or not 'DETECTOR' in df.index:
            raise IOError('Gain and detector not found in file header')
        
        if not 'PARTICLE' in df.index:
            raise IOError('Particle type not found in file header')

        if not 'LIVE_TIME' in df.index or not 'REAL_TIME' in df.index:
            raise IOError('Live time and real time not found in file header')
        
        # Get information from header
        try:
            date = datetime.strptime(df.loc['DATE'].values[0], '%m/%d/%y %H:%M:%S')
        except:
            date = datetime.strptime(df.loc['DATE'].values[0], '%Y-%m-%d %H:%M:%S')

        if num_channels == 0:
            num_channels = int(df.loc['NUM_CHANNELS'].values[0])

        detector = df.loc['DETECTOR'].values[0]
        gain = df.loc['GAIN'].values[0]
        particle = df.loc['PARTICLE'].values[0]
        live_time = float(df.loc['LIVE_TIME'].values[0])
        real_time = float(df.loc['REAL_TIME'].values[0])

        # Set attributes
        if isinstance(date, datetime) == True:
            self._date = date
        else:
            raise NotImplementedError(f'In {file}: Check date format: {date}')
        
        if math.log2(num_channels).is_integer() == True:
            self._num_channels = num_channels
        else:
            raise NotImplementedError(f'In {file}: Number of channels {num_channels} is not a power of two')
        
        if detector in ['silicon', 'diamond', 'sic', 'TEPC']:
            self._detector = detector
        else:
            raise NotImplementedError(f'In {file}: Detector type {detector} is not a viable option (silicon, diamond, sic, TEPC)')
        
        if gain in ['LOW', 'MID', 'HIGH', '', 'None', math.nan, np.nan, 'SIMULATION']:
            self._gain = gain
        else:
            raise NotImplementedError(f'In {file}: Gain {gain} is not a viable option.')
        
        if particle in ['proton', 'carbon', 'helium']:
            self._particle = particle
        else:
            raise NotImplementedError(f'In {file}: Particle type {particle} is not a viable option.')

    def read_MCA_file(self, file):
        """ Read in information from a given MCA file. """

        #Get header information and set it in measurement object
        file = os.path.normpath(file)
        self.__get_info_from_MCA_file(file)

        #For some reason you need to create this twice so the df are not linked to each other
        self._data = pd.read_csv(file, skiprows=7, sep=';')
        self._original_data = pd.read_csv(file, skiprows=7, sep=';')
        self._abs_filepath = file

        # Add x and y axis
        self._x_axis = 'CHANNEL'
        self._y_axis = 'COUNTS'

    def read_ROOT_file(self, file, NUM_ROOT_CHANNELS=4096):
        """ Read in information from a given ROOT file (GATE, Geant4 Energy deposition simulation). """

        #TODO: Include reading and attaching Track length information

        if math.log2(NUM_ROOT_CHANNELS).is_integer() == True:
            self._num_channels = NUM_ROOT_CHANNELS
        else:
            raise NotImplementedError(f'In {file}: Number of channels {NUM_ROOT_CHANNELS} is not a power of two')

        self._date = datetime.fromtimestamp(os.path.getmtime(file)).strftime('%m/%d/%y %H:%M:%S')
        self._detector = 'None (Not yet attached)'
        self._gain = 'SIMULATION'
        self._particle = 'None (Not yet attached)'
        self._x_axis = 'ENERGY'
        self._y_axis = 'COUNTS'

        keys = uproot.open(file)['Hits'].keys()

        if any('edep' in s for s in keys) == False and any('Edep' in s for s in keys) == False and any('energy' in s for s in keys) == False:
            print(keys)
            raise IOError('ROOT file does not contain an energy deposition tree (Hits/edep)')
        else:
            if 'edep' in keys:
                root_data = np.asarray(uproot.open(file)['Hits']['edep'])
            elif 'Edep' in keys:
                root_data = np.asarray(uproot.open(file)['Hits']['Edep'])
            elif 'energy' in keys:
                root_data = np.asarray(uproot.open(file)['Hits']['energy'])

            hist_data, bins = np.histogram(root_data, bins=NUM_ROOT_CHANNELS, density=False)
            bins = bins[:-1]

            #For some reason you need to create this twice so the df are not linked to each other
            self._data = pd.DataFrame({self.x_axis: np.multiply(bins, 1000), self.y_axis: hist_data})
            self._original_data = pd.DataFrame({self.x_axis: np.multiply(bins, 1000), self.y_axis: hist_data})
            self._abs_filepath = file

        track_data = None
        if 'track' in keys:
            track_data = np.asarray(uproot.open(file)['Hits']['track'])
        elif 'length' in keys:
            track_data = np.asarray(uproot.open(file)['Hits']['length'])
        elif 'TrackLength' in keys:
            track_data = np.asarray(uproot.open(file)['Hits']['TrackLength'])

        if track_data != None:
            #Automatically attach chord length info if given in ROOT files
            print('ROOT file contains chord length information -> Attached to the measurement')
            self._mean_chord_length = np.mean(track_data)
            self._max_chord_length = np.max(track_data)

    def read_uDos_analysis(self, file):
        """ Read in information from an already analysed uDos file (*csv). """

        file = os.path.normpath(file)
        x_axis, y_axis = self.__get_info_from_uDos_file(file)
        
        self._data = pd.read_csv(file, sep=';', skiprows=4, skipfooter=11, engine='python')
        self._original_data = pd.DataFrame(data='This was created from an already analysed uDos file', index=[0], columns=['This was created from an already analysed uDos file'])
        self._abs_filepath = file

        # Set x and y axis to standard format
        self._x_axis = x_axis
        self._y_axis = y_axis

        # Change the column names to the standard format
        self._data.columns = [x_axis, y_axis]

    def attach_mca_info(self, measurement, info_dict):
        """ Attach information to a MCA file measurement object.
        
        Args:
            measurement: Measurement object to be altered
            info_dict: Dictionary containing the information to be attached (_particle, _detector, _num_channels)
        """

        if measurement.gain == 'SIMULATION':
            raise ValueError('This function is only for MCA files')
        
        print('-----------------------------------')
        print(f'ATTACHING INFO TO {measurement.name}')
        print('This might require a recalibration of the energy axis !!!')
        print('-----------------------------------')

        for key, value in info_dict.items():
            setattr(measurement, key, value)

    def attach_root_info(self, measurement, info_dict):
        """ Attach information to a ROOT file measurement object.
        This can also be used to resample the spectrum (num channels)
        
        Args:
            measurement: Measurement object to be altered
            info_dict: Dictionary containing the information to be attached (_particle, _detector, _num_channels)
        """

        if measurement.gain != 'SIMULATION':
            raise ValueError('This function is only for ROOT files (GATE, Geant4)')

        for key, value in info_dict.items():
            setattr(measurement, key, value)
        
        # Resample the ROOT data (Histogram bins = num channels)
        if 'num_channels' in info_dict.keys():

            if measurement.x_axis != 'ENERGY' and measurement.y_axis != 'COUNTS ':
                raise ValueError('Resampling is only possible if no data manipulation has occured yet (energy vs counts)')

            file = measurement.abs_filepath
            
            if not 'edep' in uproot.open(file)['Hits'].keys():
                raise IOError('ROOT file does not contain an energy deposition tree (Hits/edep)')
            else:
                root_data = np.asarray(uproot.open(file)['Hits']['edep'])
                hist_data, bins = np.histogram(root_data, bins=info_dict['num_channels'], density=False)
                bins = bins[:-1]

                #For some reason you need to create this twice so the df are not linked to each other
                self._data = pd.DataFrame({measurement.x_axis: np.multiply(bins, 1000), measurement.y_axis: hist_data})
                self._original_data = pd.DataFrame({measurement.x_axis: np.multiply(bins, 1000), measurement.y_axis: hist_data})