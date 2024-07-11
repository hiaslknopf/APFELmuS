__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy import interpolate
import os

from lib import FileTranslator

#TODO: Test this also for LO; MID; HIGH

""" This script contains a helper tool to generate a calibration file (Channel -> mV) needed for the linearity/calibration procedure.
Necessary inpus are a Pulser calibration (mV goal vs mV measured) and 1 or more MCA linearization spectra (Gaussian peaks for different mV values).
"""

# Include whereever this script is located
ressources_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) + '/ressources'

def __gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / stddev) ** 2 / 2)

def __piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def _plot_pulser_cal(pulser_df, label, ax):

    if ax is None:
        fig, ax = plt.subplots()

    ax.errorbar(pulser_df['mV_goal'], pulser_df['mV_measured'], yerr=pulser_df['sigma'], fmt='x', label=label)
    ax.plot(np.linspace(0, 1000, 1000), np.linspace(0, 1000, 1000), 'k--', label='Ideal response')

    #ax.set_title('INPUT - Pulser calibration: ' + pulser + '-' + date)
    ax.set_xlabel('Goal [mV]')
    ax.set_ylabel('Measured [mV]')
    ax.grid(True)
    ax.legend()
    plt.show()

def _plot_peaks(spectrum, positions, popt_list, ax=None):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        
    ax.step(spectrum.index, spectrum['COUNTS'], label='MCA Spectrum', where='mid')
    ax.fill_between(spectrum.index, spectrum['COUNTS'], 0, step='mid', color='cornflowerblue', alpha=0.5)
    ax.plot(0,0,color='r', label='Peak Channel')
    ax.plot(0,0,color='b', label='Gaussian fit')

    for i, pos in enumerate(positions):
        popt = popt_list[i]

        pos = int(pos)

        ax.axvline(x=spectrum.index[pos], color='r', linestyle='--', alpha=0.75)
        ax.plot(spectrum.index, __gaussian(spectrum.index, *popt), 'b--', alpha=0.5)

        ax.text(spectrum.index[pos], spectrum['COUNTS'][pos]+25, f'{i+1}', ha='center', va='bottom', color='black', fontsize=12, fontweight='bold')

    ax.set_title(f'INPUT - Linearization Spectrum:')
    ax.set_xlabel('Channel')
    ax.set_ylabel('Counts')
    ax.set_xlim(0, max(positions)*1.1)
    ax.set_ylim(0, 1.1 * max(spectrum['COUNTS']))
    ax.grid(True)
    ax.legend()

def _plot_linearization(channel, mV, lin_curve_channel, lin_curve_mV, method, num_channels, ax=None):

    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(channel, mV, 'x', label='Measured Peak Heights')
    ax.plot(lin_curve_channel, lin_curve_mV, 'r--', label=f'{method} fit')

    ax.set_xlabel('Channel')
    ax.set_ylabel('Test input voltage [mV]')
    ax.grid(which='both')
    ax.set_xlim(0, num_channels)
    ax.set_ylim(0, 1.1 * max(mV))
    ax.legend()

    plt.show()

def _read_spectrum(filename):
    """ Reads a spectrum file and returns a pandas dataframe with the data.

    Args:
        filename: The spectrum file to be read
    Returns:
        data: Dataframe containing the spectrum data
    """
    data = pd.read_csv(filename, sep='\t', header=None, names=['mV', 'counts'])

    return data

def _read_pulser(filename, pulse_mV_list, testplot=False, ax=None):
    """ Reads a pulser file and returns an array containing the data.
    
    Args:
        filename: The pulser file to be read
        testplot: If True, plots are generated and shown for each step
    Returns:
        pulser_df: Dataframe containing the pulser data
    """

    name = os.path.basename(filename)
    #pulser = name.split('_')[0]
    #date = name.split('_')[2].split('.')[0]

    pulser_df = pd.read_csv(filename, sep=',', skiprows=1, names=['mV_goal', 'mV_measured', 'sigma'], dtype=float)
    # Check if sigma is in uV or mV
    with open(filename, 'r') as f:
        line = f.readline()
        if 'uV' in line:
            pulser_df['sigma'] = pulser_df['sigma'] / 1000
        
    label='Measured values'

    # If any of the pulse_mV_list values are not in the pulser_df, interpolate them
    diff = list(set(pulse_mV_list).difference(pulser_df['mV_goal'].values))
    diff = np.setdiff1d(pulse_mV_list,pulser_df['mV_goal'].values)
    if diff.any():
        print('Interpolating Pulser calibration for missing values: ', diff)
        label = 'Interpolated values\nfrom measurement'

        # Interpolate for actual values in pulse_mV_list - Uncertainties are set to 0
        pulser_interp = interpolate.interp1d(pulser_df['mV_goal'], pulser_df['mV_measured'], kind='linear', bounds_error=False, fill_value='extrapolate')
        pulser_df = pd.DataFrame({'mV_goal': pulse_mV_list, 'mV_measured': pulser_interp(pulse_mV_list), 'sigma': 0})
    
    # Remove any values that are not in pulse_mV_list
    pulser_df = pulser_df[pulser_df['mV_goal'].isin(pulse_mV_list)]

    if testplot:
        _plot_pulser_cal(pulser_df, label, ax)

    return pulser_df

def _fit_peaks(spectrum, testplot=False, cutoff_front=False, cutoff_back=False):
    """ Fit the peaks in the MCA spectrum and return the positions.

    Args:
        spectrum: The spectrum to be fitted
        name: The name of the spectrum
        output_path: The path to save the plots
        testplot: If True, plots are generated and shown for each step
        cutoff_front: Low energy noise cutoff, so the noise peak is not fitted (may be necessary)
        cutoff_back: High energy noise cutoff, so the noise peak is not fitted (may be necessary)
    Returns:
        positions: The positions of the peaks in the spectrum (channels)    
    """

    #Noise cutoff, so initial slope is not recognized as peak
    if cutoff_front:
        spectrum['COUNTS'][:cutoff_front] = 0
    
    if cutoff_back:
        spectrum['COUNTS'][cutoff_back:] = 0

    positions_gauss = []
    #It may be necessary to tune the prominence to find all peaks
    positions, prop = find_peaks(spectrum['COUNTS'], height=50, prominence=120, distance=10) 

    # If peak found on the left edge, remove it (noise artifact)
    if positions[0] == 0 or positions[0] == cutoff_front:
        print('Peak found on the left edge, removing it')
        positions = positions[1:]

    print(len(positions), 'peaks found in the spectrum')

    popt_list = []
    for i, pos in enumerate(positions):
        initial_guesses = [spectrum['COUNTS'][pos], spectrum.index[pos], 1.0]
        popt, pcov = curve_fit(__gaussian, spectrum.index, spectrum['COUNTS'], p0=initial_guesses)

        popt_list.append(popt)
        positions_gauss.append(popt[1])        

    if testplot:
        _plot_peaks(spectrum, positions_gauss, popt_list, ax=None)

    positions = np.array(positions_gauss) #Take the gaussian fit positions (more accurate)

    return positions, popt_list

def _fit_linearization(channel, mV, num_channels, method):
    """ Fit the linearization curve and return the parameters

    Args:
        channel: The channel numbers of the peaks
        mV: The measured mV signal of the peaks
        num_channels: The number of channels of the MCA
        method: The method for the linearization curve fit (linear, piecewise, linear_force_zero)
    Returns:
        x, y of the fit function
    """

    if method == 'linear':
        # Single line linear fit

        popt = np.polyfit(channel, mV, 1)
        y_fit = np.polyval(popt, np.linspace(1, num_channels, num_channels))
    
    if method == 'interpol':
        # Point to point linear interpolation

        #Linear extrapolation
        new_mV = mV[-2] + (num_channels - channel[-2]) / (channel[-1] - channel[-2]) * (mV[-1] - mV[-2])
        channel = np.append(channel, num_channels)
        mV = np.append(mV, new_mV)

        #Extrapolate towards channel 0
        new_mV = mV[1] - (channel[1] - 1) / (channel[2] - channel[1]) * (mV[2] - mV[1])
        channel = np.insert(channel, 0, 0)
        mV = np.insert(mV, 0, new_mV)

        y_fit = interpolate.interp1d(channel, mV)(np.linspace(1, num_channels, num_channels))
    
    elif method == 'linear_force_zero':
        # Single linear fit with forced zero intercept
        popt, _ = curve_fit(lambda x, a: a*x, channel, mV)
        y_fit = np.array([popt[0]*x for x in np.linspace(1, num_channels, num_channels)])
    
    else:
        raise ValueError(f'Linearizaiton method {method} not implemented. (linear, interpol, linear_force_zero)')
        
    return np.linspace(1, num_channels, num_channels), y_fit

def get_linearization(name:str, output_path:str, pulse_mV_list:list, spectrum_file1:str, spectrum_file2:str=None, spectrum_file3:str=None,
                      pulser_calibration:str = f'{ressources_path}/1to1_response.csv', method:str='linear', testplot:bool=False,
                      cutoff_front:bool=False, cutoff_back:bool=False):
    """ Get a calibration file (csv) in 4 steps

    - Read the pulser and spectrum file
    - Fit the peaks (Gaussian)
    - Fit the linearization curve
    - Write the calibration file (csv)

    Args:
        name: The name of the calibration file
        output_path: The path to save the plots
        pulse_mV_list: The list of mV values to be used for the pulser calibration
        pulser_calibration: Calibration file of the pulser used (default: 1to1_response.csv)
        spectrum_file<n>: The spectrum file(s) to be read (max 3)
        method: The method for the linearization curve fit (linear, piecewise)
        testplot: If True, plots are generated and shown for each step
        cutoff_front: Low energy noise cutoff, so the noise peak is not fitted (may be necessary)
        cutoff_back: High energy noise cutoff, so the noise peak is not fitted (may be necessary)
    """

    #Read the pulser file
    pulser_df = _read_pulser(pulser_calibration, pulse_mV_list, testplot)
    print(f'Read Pulser calibration: {pulser_calibration}\n')

    #Read the spectrum file(s)
    spectra_list = []
    peaks_list = []

    if spectrum_file2:
        spectrum_df1, num_channels1 = _read_spectrum(spectrum_file1)
        spectrum_df2, num_channels2 = _read_spectrum(spectrum_file2)
        spectrum_df3, num_channels3 = _read_spectrum(spectrum_file3)

        spectra_list.append(spectrum_df1, spectrum_df2, spectrum_df3)

        if num_channels1 != num_channels2 or num_channels1 != num_channels3:
            raise ValueError('The number of channels in the spectrum files is not equal!')
        else:
            num_channels = num_channels1
    else:
        spectrum_df, _, num_channels = FileTranslator._read_MAESTRO_file(spectrum_file1)
        spectra_list.append(spectrum_df)
    
    print(f'Read MCA spectrum: {spectrum_file1}')
    print(f'Number of channels: {num_channels}\n')
    
    print(f'Gaussian peak fitting for {len(spectra_list)} spectra...')
    #Fit the peaks
    for i in range(len(spectra_list)):
        peaks, popt_list = _fit_peaks(spectra_list[i], testplot, cutoff_front, cutoff_back)
        peaks_list.append(peaks)

    # Get noise figure from average sigma of the gaussian fit
    sigma = np.mean([popt[2] for popt in popt_list])

    #Fit the linearization curve
    if spectrum_file2:
        columns = ['CHANNEL', 'LOW [mV]', 'MID [mV]', 'HIGH [mV]']
        linearization_df = pd.DataFrame(index=np.arange(num_channels), columns=columns)
    else:
        columns = ['CHANNEL', 'INPUT [mV]']
        linearization_df = pd.DataFrame(index=np.arange(num_channels), columns=columns)

    for i in range(len(peaks_list)):
        print(f'Fitting linearization curve for spectrum {i+1}...')

        if len(peaks_list[i]) != len(pulser_df['mV_goal']):
            raise ValueError(f'The number of peaks in the spectrum {len(peaks_list[i])} does not match the number of pulser calibration points in pulser_mV_list {len(pulser_df)}!')

        channel = peaks_list[i]
        mV = pulser_df['mV_measured'].values

        lin_curve_channel, lin_curve_mV = _fit_linearization(channel, mV, num_channels, method)

        if spectrum_file2:
            linearization_df['CHANNEL'] = lin_curve_channel.astype(int)
            linearization_df[columns[i]] = np.round(lin_curve_mV, 3)
        else:
            linearization_df['CHANNEL'] = lin_curve_channel.astype(int)
            linearization_df['INPUT [mV]'] = np.round(lin_curve_mV, 3)

        if testplot:
            _plot_linearization(channel, mV, lin_curve_channel, lin_curve_mV, method, num_channels, ax=None)

    #Write the calibration file
    pd.DataFrame.to_csv(linearization_df, f'{output_path}/{name}.csv', sep=',', index=False)
    print(f'Wrote linearization file: {output_path}/{name}')

    print('\n')
    print(f'Average noise figure from fits (FWHM): {sigma*2.355:.2f} channels')

def linearization_from_MAESTRO_rpt(maestro_rpt_file: str, pulse_mV_list:list, pulser_calibration:str = f'{ressources_path}/1to1_response.csv',
                                   output_path:str='', name:str='linearization', 
                                   method:str='interpol', testplot:bool=False, num_channels:int=4096):
    
    """ Get a calibration file (csv) from a MAESTRO report file in 4 steps

    - Read the pulser and spectrum file
    - Fit the peaks (Gaussian)
    - Fit the linearization curve
    - Write the calibration file (csv)

    Args:
        maestro_rpt_file: The MAESTRO report file to be read
        pulse_mV_list: The list of mV values to be used for the pulser calibration
        output_path: The path to save the plots
        method: The method for the linearization curve fit (linear, piecewise)
        testplot: If True, plots are generated and shown for each step
    """

    #Read the MAESTRO report file
    channel = FileTranslator._read_MAESTRO_rpt_file(maestro_rpt_file)
    print(f'Read MAESTRO report file: {maestro_rpt_file}\n')

    #Read the pulser file
    pulser_df = _read_pulser(pulser_calibration, pulse_mV_list, testplot)
    print(f'Read Pulser calibration: {pulser_calibration}\n')


    #Fit the linearization curve
    print(f'Fitting linearization curve...')
    mV = pulser_df['mV_measured'].values

    lin_curve_channel, lin_curve_mV = _fit_linearization(channel, mV, num_channels, method)

    #Write the calibration file
    columns = ['CHANNEL', 'INPUT [mV]']
    linearization_df = pd.DataFrame(index=np.arange(num_channels), columns=columns)

    linearization_df['CHANNEL'] = lin_curve_channel.astype(int)
    linearization_df['INPUT [mV]'] = np.round(lin_curve_mV, 3)

    pd.DataFrame.to_csv(linearization_df, f'{output_path}/{name}.csv', sep=',', index=False)
    print(f'Wrote linearization file: {output_path}/{name}')

    if testplot:
        _plot_linearization(channel, mV, lin_curve_channel, lin_curve_mV, method, num_channels, ax=None)
