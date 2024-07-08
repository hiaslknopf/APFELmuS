__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

from MicroDosimetry import Measurement
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d
from scipy.optimize import minimize

""" This module contains all the functions for the manipulation of spectra """

def _find_nearest_idx_and_value(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def _gain_match_chi2(k, y_ref, y_match):
    """ Adjust the gain of y_match with a factor k to match the reference spectrum y_ref"""
    a = y_ref - k*y_match
    
    # Remove all math.nan values
    a = a[~np.isnan(a)]
    
    chi2 = np.sum(np.square(a))
    return chi2

def _lin_2_log(measurement, n_bins_per_decade: int):
    """ Calculate n evenly spaced logarithmic bins per decade for a given x_axis

    Args:
        Measurement: the spectrum to be analyzed
        n_bins_per_decade: Number off log bins per decade
    Returns:
        lin_to_log_x: List of logarithmic bins
    """
    
    x_val = measurement.data[measurement.x_axis]

    #Get metrics of x-axis
    start = measurement.data[x_val >= 0.0][measurement.x_axis].values[0]
    stop = x_val.max()
    num_decades = np.log10(stop)-np.log10(start)

    if start == None:
        start = 0 #For extrapolated data
    
    #Create log bins
    log_bins_x = np.logspace(np.log10(start), np.log10(stop), num=int(np.ceil(num_decades*n_bins_per_decade)), endpoint=True)
    lin_to_log_x = [np.abs(x_val.to_numpy() - log_bins_x[i]).argmin() for i in range(len(log_bins_x))]

    return lin_to_log_x

def cutoff(measurement, channels:int=None, energy:float=None, lineal_energy:float=None):
    """ Low energy cutoff for noisy channels
    Give either channels or energy as cutoff value
    
    Args:
        measurement: The spectrum to be analyzed
        channels: Number of initial channels to be set zero
        energy: Energy value to be set as cutoff (keV)
        lineal_energy: Lineal energy value to be set as cutoff (keV/um)
    """

    if channels:        
        measurement.data.loc[:channels, measurement.y_axis] = 0
        print(f'Channels {channels} and below cut off: {measurement.name}')
    
    elif energy:
        if measurement.x_axis != 'ENERGY':
            raise ValueError('This manipulation only works for a calibrated energy axes')
        
        idx, value = _find_nearest_idx_and_value(measurement.data[measurement.x_axis].tolist(), energy)
        measurement.data.loc[:idx, measurement.y_axis] = 0
        print(f'Energy {value} (Channel {idx}) and below cut off: {measurement.name}')
    
    elif lineal_energy:
        if measurement.x_axis != 'LINEAL ENERGY':
            raise ValueError('This manipulation only works for a calibrated lineal energy axes')
        
        idx, value = _find_nearest_idx_and_value(measurement.data[measurement.x_axis].tolist(), lineal_energy)
        measurement.data.loc[:idx, measurement.y_axis] = 0
        print(f'Lineal energy {value} (Channel {idx}) and below cut off: {measurement.name}')

def logarithmic_binning(measurement, n_bins_per_decade:int):
    """ Rescale measurement dataframe to a (semi)logarithmic x-axis
    
    Args:
        measurement: the spectrum to be analyzed
        n_bins_per_decade: Number of log bins per decade
    """

    if measurement.x_axis != 'LINEAL ENERGY':
        raise ValueError('This manipulation only works for a calibrated lineal energy axes')

    x = measurement.data[measurement.x_axis].to_numpy()
    y = measurement.data[measurement.y_axis].to_numpy()

    lin_to_log_x = _lin_2_log(measurement, n_bins_per_decade)
    
    quasi_log_x = np.asarray([])
    quasi_log_y = np.asarray([])

    quasi_log_x = x[lin_to_log_x]
    quasi_log_y = y[lin_to_log_x]

    #get rid of old entries
    measurement.data.drop(columns=measurement.x_axis, inplace=True)
    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data.index = range(0, len(quasi_log_x))

    #Replace data (new length)
    measurement.data[measurement.x_axis] = quasi_log_x.tolist()
    measurement.data[measurement.y_axis] = quasi_log_y.tolist()

    print(f'Spectrum has been logarithmically binned with {n_bins_per_decade} bins per decade: {measurement.name}')

def probability_function(measurement, type_of_dist:str):
    """ To get from MCA spectrum counts(channel) to F(y) or D(y)
    
    The value of the distribution function F(y) is the probability that the lineal energy is equal to or less than y
    The same goes for the dose distribution function D(y)

    Args:
        measurement: Spectrum to be analyzed
        type_of_dist: 'F' for F(y), 'D' for D(y)
    """
    
    #TODO: Does this also work for logarithmically binned data?

    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis != 'COUNTS':
        raise ValueError('This manipulation only works for energy calibrated MCA spectra [lineal energy, counts]')

    num_p_above_cutoff = int(np.sum(measurement.data[measurement.y_axis].tolist()))
    weighted_num_p_above_cutoff = float(np.sum(np.multiply(measurement.data[measurement.x_axis].tolist(), measurement.data[measurement.y_axis].tolist())))

    F_temp = [0.0]
    D_temp = [0.0]

    if type_of_dist == 'F':
        for i in range(1,len(measurement.data)):
            F_temp.append(np.sum(measurement.data[measurement.y_axis].tolist()[:i]))
        F_norm = np.divide(F_temp, num_p_above_cutoff, out=np.zeros_like(F_temp), where=num_p_above_cutoff!=0)
 
        measurement.data.drop(columns=measurement.y_axis, inplace=True)
        measurement.data[measurement.y_axis] = F_norm
        
        if measurement.x_axis == 'LINEAL ENERGY':
            #measurement.x_axis = 'LINEAL ENERGY'
            measurement.y_axis = 'F(y)'
        elif measurement.x_axis == 'ENERGY':
            #measurement.x_axis = 'ENERGY'
            measurement.y_axis = 'F(E)'

        print(f'Probability function {type_of_dist}(y) calculated: {measurement.name}')

    elif type_of_dist == 'D':
        for i in range(1,len(measurement.data)):
            mul = np.multiply(measurement.data[measurement.x_axis].tolist()[:i], measurement.data[measurement.y_axis].tolist()[:i])
            D_temp.append(np.sum(mul))
        D_norm = np.divide(D_temp, weighted_num_p_above_cutoff, out=np.zeros_like(D_temp), where=weighted_num_p_above_cutoff!=0)

        measurement.data.drop(columns=measurement.y_axis, inplace=True)
        measurement.data[measurement.y_axis] = D_norm

        if measurement.x_axis == 'LINEAL ENERGY':
            #measurement.x_axis = 'LINEAL ENERGY'
            measurement.y_axis = 'D(y)'
        elif measurement.x_axis == 'ENERGY':
            #measurement.x_axis = 'ENERGY'
            measurement.y_axis = 'D(E)'
        
        print(f'Probability function {type_of_dist}(y) calculated: {measurement.name}')

    else:
        raise KeyError(f'{type_of_dist} is not a valuable option for the type of distribution')

def probability_density(measurement):
    """ To get from probability F(y) or D(y) to probability distribution f(y) or d(y)
    
    Args:
        measurement: Spectrum to be analyzed
    """

    #TODO: With logscale data you have to deal with delta=0 somehow
    
    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis not in ['F(y)', 'D(y)', 'F(E)', 'D(E)']:
        raise ValueError('This manipulation only works for probability functions [y, F(y)] or [y, D(y)]')

    y = measurement.data[measurement.x_axis].to_numpy()
    F = measurement.data[measurement.y_axis].to_numpy()

    pdf = [0.0]

    #This is just numerical derivation
    for i in range(1,len(measurement.data)):
        delta_y = y[i] - y[i-1]
        #Just a temporary sanity check
        if y[i] < y[i-1]:
            raise ValueError('Something is wrong with the data: y[i] < y[i-1]')
        else:
            delta_F = F[i] - F[i-1]
            pdf.append(np.divide(delta_F,delta_y, out=np.zeros_like(delta_F), where=delta_y!=0))
        
    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data[measurement.y_axis] = pdf

    if measurement.x_axis == 'LINEAL ENERGY':
        #measurement.x_axis = 'LINEAL ENERGY'
        if measurement.y_axis == 'F(y)':
            measurement.y_axis = 'f(y)'
        elif measurement.y_axis == 'D(y)':
            measurement.y_axis = 'd(y)'
    
    elif measurement.x_axis == 'ENERGY':
        #measurement.x_axis = 'ENERGY'
        if measurement.y_axis == 'F(E)':
            measurement.y_axis = 'f(E)'
        elif measurement.y_axis == 'D(E)':
            measurement.y_axis = 'd(E)'

    print(f'Probability density function {measurement.y_axis} calculated: {measurement.name}')

def dose_density(measurement):
    """ Calculate dose pdf d(y) from from frequency pdf f(y)
    
    Args:
        measurement: Spectrum to be analyzed
    """
    
    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis not in ['f(y)', 'f(E)']:
        raise ValueError('This manipulation only works for a pdf [y, f(y)]')

    lineal_energy = measurement.data[measurement.x_axis].to_numpy()
    pdf = measurement.data[measurement.y_axis].to_numpy()

    y_F = np.trapz(lineal_energy*pdf, x=lineal_energy, dx=1.0)
    d = measurement.data.loc[:, measurement.x_axis]/ y_F*measurement.data.loc[:, measurement.y_axis]
    #d = np.divide(measurement.data.loc[:, measurement.x_axis], np.multiply(y_F, measurement.data.loc[:, measurement.y_axis]))
    
    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data[measurement.y_axis] = d

    if measurement.x_axis == 'LINEAL ENERGY':
        #measurement.x_axis = 'LINEAL ENERGY'
        measurement.y_axis = 'd(y)'
    elif measurement.x_axis == 'ENERGY':
        #measurement.x_axis = 'ENERGY'
        measurement.y_axis = 'd(E)'

    print(f'Dose density function {measurement.y_axis} calculated: {measurement.name}')

def weighted_probability_density(measurement):
    """ To get from probability distribution f(y) or d(y) to weighted pdf yf(y) or yd(y)
    
    Args:
        measurement: Spectrum to be analyzed
    """

    #TODO: Does this also work for logarithmically binned data?
    
    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis not in ['f(y)', 'd(y)', 'f(E)', 'd(E)']:
        raise ValueError('This manipulation only works for probability density functions [y, f(y)] or [y, d(y)]')

    y = measurement.data[measurement.x_axis].to_numpy()
    pdf = measurement.data[measurement.y_axis].to_numpy()

    weighted_pdf = [0.0]

    weighted_pdf = np.multiply(y, pdf)
        
    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data[measurement.y_axis] = weighted_pdf

    if measurement.x_axis == 'LINEAL ENERGY':
        #measurement.x_axis = 'LINEAL ENERGY'
        if measurement.y_axis == 'f(y)':
            measurement.y_axis = 'yf(y)'
        elif measurement.y_axis == 'd(y)':
            measurement.y_axis = 'yd(y)'
    
    elif measurement.x_axis == 'ENERGY':
        #measurement.x_axis = 'ENERGY'
        if measurement.y_axis == 'f(E)':
            measurement.y_axis = 'Ef(E)'
        elif measurement.y_axis == 'd(E)':
            measurement.y_axis = 'Ed(E)'
    else:
        raise ValueError('This manipulation only works for probability density functions [y, f(y)] or [y, d(y)]')

    print(f'Weighted probability density function {measurement.y_axis} calculated: {measurement.name}')

def normalize_spectrum(measurement):
    """ Normalize a spectrum to A=1 --> Numerical integration of arbitrary binned data.
    See ICRU report 36, annex B. Implemented for pdf/y spectra.
    
    Args:
        measurement: Spectrum to be analyzed
    """
    #TODO: Something about this does not quite work yet --> A!=1 for some spectra
    
    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis not in ['f(y)', 'd(y)', 'yf(y)', 'yd(y)', 'f(E)', 'd(E)', 'Ef(E)', 'Ed(E)']:
        raise ValueError('Normalization is only implemented for lineal energy pdfs: f(y), d(y), yf(y), yd(y)')

    x = measurement.data[measurement.x_axis].tolist()
    pdf = measurement.data[measurement.y_axis].tolist()

    #TODO: Evtl. wär das sowieso kein schlechter sanity check iwo im setter?
    if all(i >= 0.0 for i in x) == False:
        raise ValueError('There are negative values on your (lineal) energy axis!')
    elif all(i >= 0.0 for i in pdf) == False:
        np.abs(pdf)
        print('Negative values on your pdf axis have been removed!')
        #raise ValueError('There are negative values on your pdf axis!')

    x.insert(0,0.0)
    x.insert(-1,0.0)

    geometrical_mean_p = [] 
    geometrical_mean_m = [] 

    for i in range(1,len(x)-1):
        geometrical_mean_p.append(np.sqrt(x[i]*x[i+1])) #y_{i+1/2}
    for i in range(1,len(x)-1):
        geometrical_mean_m.append(np.sqrt(x[i]*x[i-1])) #y_{i-1/2}

    x = x[1:len(x)-2]
    pdf = pdf[:len(pdf)-1] 

    sum_for_norm = []

    for i in range(len(pdf)):
        if pdf[i] < 0.0:
            pdf[i] = 0.0
            pdf[i] = 0.0
        sum_for_norm.append(pdf[i]*(geometrical_mean_p[i]-geometrical_mean_m[i]))

    norm_factor = np.sum(sum_for_norm)
    pdf_norm = np.divide(pdf, norm_factor)
    
    #TODO: Why does Sandra do this? Shouldn't every Spectrum be normalized to A=1
    #if measurement.y_axis in ['yf(y)', 'yd(y)']:
    #    pdf_norm = np.multiply(pdf_norm, x)
    
    pdf_norm = np.append(pdf_norm, 0.0)

    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data[measurement.y_axis] = np.abs(pdf_norm)
    #measurement.data = measurement.data

    #Sanity check
    area = np.trapz(measurement.data[measurement.y_axis], x=measurement.data[measurement.x_axis], dx=1.0)
    print(f"{measurement.name} normalized to A={area: .3f}")

def extrapolate(measurement, method:str):
    """ Append an extrapolation towards lower energies to the measurement dataframe
    
    Args:
        measurement: The spectrum to be analyzed
        method: Method of extrapolation ('linear', ...)
    """

    if measurement.x_axis not in ['ENERGY', 'LINEAL ENERGY'] or measurement.y_axis not in ['f(y)', 'f(E)']:
        raise ValueError('This manipulation only works for probability density functions [y, f(y)]')

    if measurement.gain not in ['HIGH', 'None', np.nan, None]:
        raise ValueError('An extrapolation should only be performed for the highest gain level in a measurement. Also not for simulations.')

    x_val = measurement.data[measurement.x_axis]
    y_val = measurement.data[measurement.y_axis]

    #Linear x-axis
    if len(y_val) == measurement.num_channels:

        start_pos = measurement.data[y_val > 0.0].index[0]

        #Replace all entries below start_pos --> aka all the zeroes in the df
        if method == 'linear':
            x = x_val.to_numpy()[:start_pos]
            k = y_val.to_numpy()[start_pos]/x_val.to_numpy()[start_pos]
            d = 0.0

            extrapol = np.concatenate((np.array(k*x + d), np.array(y_val[start_pos:])))
        else:
            raise ValueError(f"Extrapolation method {method} currently not implemented")

    else:
        #TODO: Implement this for log binned data?
        raise NotImplementedError('Currently only implemented for linear x-axes')
    
    #Replace data y axis
    measurement.data.drop(columns=measurement.y_axis, inplace=True)
    measurement.data[measurement.y_axis] = extrapol.tolist()
    
    print(f'Spectrum exttrapolated: {measurement.name}')

def merge_spectra(measurements: list, overlap_regions: list, scaling_factors: list, stitching_points: list, testplots: bool = False):
    """ Merge several spectra (df objects) into one dataframe object and match the gains
        Be sure to provide the parameters in the order, the data appears in the spectrum: HIGH -> MID -> LOW
        The highest gain means the lowest energies !!!
    
    Args:
        measurements: The list of spectra to be merged (two or three gains possible: LOW, MID, HIGH)
        overlap_regions: The overlap regions for the different gains: MIGH-MID, MID-LOW in (lineal) energy units (list of lists)
        scaling_factors: The scaling factors or an educated guess from the lowest blabla (list of floats)
        stitching_points: The points where the different gains are stitched together (list of floats)
    Returns:
        merged_spectrum: The merged spectrum (Measurement object)
        stitching_info: The information about final parameters (dict)
    """

    #TODO: Brauch ich eine gemeinsame x-Achse für alle 3? Also iwie interpolieren, dass gleich gebinnt

    color_map = {'HIGH': 'cornflowerblue', 'MID': 'orange', 'LOW': 'forestgreen',
                 'overlap_low': 'r', 'overlap_high': 'b'}
    GRANULARITY = 100 # For the overlap regions

    #---------Check if data matches-----------
    check_x = []
    check_y = []

    if len(measurements) < 2:
        raise ValueError('Please provide at least two spectra for merging')

    for i in range(len(measurements)):
        check_x.append(measurements[i].x_axis)
        check_y.append(measurements[i].y_axis)

    if(len(set(check_x)) == True):
        print('all x axes equal')
    else:
        raise ValueError('The spectra have to have the same x axes')

    if(len(set(check_y)) == True):
        print('all y axes equal')
    else:
        raise ValueError('The spectra have to have the same y axes')

    if measurements[0].x_axis not in ['ENERGY', 'LINEAL ENERGY']:
        raise ValueError('The spectra have to have a calibrated energy axis')
    
    if measurements[0].gain == 'SIMULATION':
        raise ValueError('The spectra have to be measured spectra, not simulation data')
    
    # Check is upper and lower bounds are provided and in the right order
    if len(overlap_regions) != len(measurements)-1:
        raise ValueError('Please provide sufficient overlap intervals for the different gains in the right order')
    if len(scaling_factors) != len(measurements)-1:
        raise ValueError('Please provide sufficient scaling factors for the different gains in the right order')
    if len(stitching_points) != len(measurements)-1:
        raise ValueError('Please provide sufficient stitching points for the different gains in the right order')
    
    for i in range(len(overlap_regions)):
        if overlap_regions[i][0] > overlap_regions[i][1]:
            raise ValueError('The overlap regions have to be provided in the right order')

    # Sort the measurements by gain (in the order they appear in the spectrum !!!)
    sorting_order = {'HIGH': 0, 'MID': 1, 'LOW': 2}
    measurements.sort(key=lambda x: sorting_order[x.gain])

    data_x = []
    data_y = []
    for meas in range(len(measurements)):
        data_x.append(measurements[meas].data[measurements[meas].x_axis].tolist())
        data_y.append(measurements[meas].data[measurements[meas].y_axis].tolist())

    # Overlap regions
    bounds_dict = {'lower_overlap' : {'lower_energy_guess': overlap_regions[0][0], 'upper_energy_guess': overlap_regions[0][1],
                                      'lower_idx_low': 0, 'upper_idx_low': 0, 'lower_idx_high': 0, 'upper_idx_high': 0,
                                      'lower_energy_low': 0, 'upper_energy_low': 0, 'lower_energy_high': 0, 'upper_energy_high': 0,
                                      'x_axis_low': [], 'x_axis_high': [], 'y_axis_low': [], 'y_axis_high': []}}
    if len(measurements) == 3:
        bounds_dict['upper_overlap'] = {'lower_energy_guess': overlap_regions[1][0], 'upper_energy_guess': overlap_regions[1][1],
                                        'lower_idx_low': 0, 'upper_idx_low': 0, 'lower_idx_high': 0, 'upper_idx_high': 0,
                                        'lower_energy_low': 0, 'upper_energy_low': 0, 'lower_energy_high': 0, 'upper_energy_high': 0,
                                        'x_axis_low': [], 'x_axis_high': [], 'y_axis_low': [], 'y_axis_high': []}
                                        
    # Fill the dictionary
    for i, meas in enumerate(measurements):
        if i == 0:
            # HIGHEST GAIN
            lower_idx_low, lower_energy_low = _find_nearest_idx_and_value(data_x[i], bounds_dict['lower_overlap']['lower_energy_guess'])
            upper_idx_low, upper_energy_low = _find_nearest_idx_and_value(data_x[i], bounds_dict['lower_overlap']['upper_energy_guess'])
            bounds_dict['lower_overlap']['lower_idx_low'] = lower_idx_low
            bounds_dict['lower_overlap']['lower_energy_low'] = lower_energy_low
            bounds_dict['lower_overlap']['upper_idx_low'] = upper_idx_low
            bounds_dict['lower_overlap']['upper_energy_low'] = upper_energy_low
        if i == 1:
            # SECOND HIGHEST GAIN
            lower_idx_high, lower_energy_high = _find_nearest_idx_and_value(data_x[i], bounds_dict['lower_overlap']['lower_energy_guess'])
            upper_idx_high, upper_energy_high = _find_nearest_idx_and_value(data_x[i], bounds_dict['lower_overlap']['upper_energy_guess'])
            bounds_dict['lower_overlap']['lower_idx_high'] = lower_idx_high
            bounds_dict['lower_overlap']['lower_energy_high'] = lower_energy_high
            bounds_dict['lower_overlap']['upper_idx_high'] = upper_idx_high
            bounds_dict['lower_overlap']['upper_energy_high'] = upper_energy_high

            if len(measurements) == 3:
                # OPTIONAL LOWEST GAIN (3 gains)
                lower_idx_low, lower_energy_low = _find_nearest_idx_and_value(data_x[i], bounds_dict['upper_overlap']['lower_energy_guess'])
                upper_idx_low, upper_energy_low = _find_nearest_idx_and_value(data_x[i], bounds_dict['upper_overlap']['upper_energy_guess'])
                bounds_dict['upper_overlap']['lower_idx_low'] = lower_idx_low
                bounds_dict['upper_overlap']['lower_energy_low'] = lower_energy_low
                bounds_dict['upper_overlap']['upper_idx_low'] = upper_idx_low
                bounds_dict['upper_overlap']['upper_energy_low'] = upper_energy_low

        if i == 2:
            # OPTIONAL LOWEST GAIN (3 gains)
            lower_idx_high, lower_energy_high = _find_nearest_idx_and_value(data_x[i], bounds_dict['upper_overlap']['lower_energy_guess'])
            upper_idx_high, upper_energy_high = _find_nearest_idx_and_value(data_x[i], bounds_dict['upper_overlap']['upper_energy_guess'])
            bounds_dict['upper_overlap']['lower_idx_high'] = lower_idx_high
            bounds_dict['upper_overlap']['lower_energy_high'] = lower_energy_high
            bounds_dict['upper_overlap']['upper_idx_high'] = upper_idx_high
            bounds_dict['upper_overlap']['upper_energy_high'] = upper_energy_high
   
    # Create x-axes and y-axes for the overlap regions
    bounds_dict['lower_overlap']['x_axis_low'] = np.linspace(bounds_dict['lower_overlap']['lower_energy_low'], bounds_dict['lower_overlap']['upper_energy_low'], GRANULARITY)
    bounds_dict['lower_overlap']['x_axis_high'] = np.linspace(bounds_dict['lower_overlap']['lower_energy_high'], bounds_dict['lower_overlap']['upper_energy_high'], GRANULARITY)
    bounds_dict['lower_overlap']['y_axis_low'] = interp1d(data_x[0][bounds_dict['lower_overlap']['lower_idx_low']:bounds_dict['lower_overlap']['upper_idx_low']],
                                                          data_y[0][bounds_dict['lower_overlap']['lower_idx_low']:bounds_dict['lower_overlap']['upper_idx_low']],
                                                          bounds_error=False)(bounds_dict['lower_overlap']['x_axis_low'])
    bounds_dict['lower_overlap']['y_axis_high'] = interp1d(data_x[1][bounds_dict['lower_overlap']['lower_idx_high']:bounds_dict['lower_overlap']['upper_idx_high']],
                                                           data_y[1][bounds_dict['lower_overlap']['lower_idx_high']:bounds_dict['lower_overlap']['upper_idx_high']],
                                                           bounds_error=False)(bounds_dict['lower_overlap']['x_axis_high'])
    
    if len(measurements) == 3:
        bounds_dict['upper_overlap']['x_axis_low'] = np.linspace(bounds_dict['upper_overlap']['lower_energy_low'],bounds_dict['upper_overlap']['upper_energy_low'], GRANULARITY)
        bounds_dict['upper_overlap']['x_axis_high'] = np.linspace(bounds_dict['upper_overlap']['lower_energy_high'], bounds_dict['upper_overlap']['upper_energy_high'], GRANULARITY)
        bounds_dict['upper_overlap']['y_axis_low'] = interp1d(data_x[1][bounds_dict['upper_overlap']['lower_idx_low']:bounds_dict['upper_overlap']['upper_idx_low']],
                                                              data_y[1][bounds_dict['upper_overlap']['lower_idx_low']:bounds_dict['upper_overlap']['upper_idx_low']],
                                                              bounds_error=False)(bounds_dict['upper_overlap']['x_axis_low'])
        bounds_dict['upper_overlap']['y_axis_high'] = interp1d(data_x[2][bounds_dict['upper_overlap']['lower_idx_high']:bounds_dict['upper_overlap']['upper_idx_high']],
                                                               data_y[2][bounds_dict['upper_overlap']['lower_idx_high']:bounds_dict['upper_overlap']['upper_idx_high']]
                                                               , bounds_error=False)(bounds_dict['upper_overlap']['x_axis_high'])

    #testplot
    if testplots:
        fig, ax = plt.subplots()

        for i in range(len(measurements)):
            ax.step(data_x[i], data_y[i], label = f"{measurements[i].gain}", color = color_map[measurements[i].gain])
        
        # Plot the overlap regions
        ax.axvline(x = bounds_dict['lower_overlap']['lower_energy_low'], color=color_map['overlap_low'], linestyle='--')
        ax.axvline(x = bounds_dict['lower_overlap']['upper_energy_low'], color=color_map['overlap_low'], linestyle='--')
        ax.axvline(x = bounds_dict['lower_overlap']['lower_energy_high'], color='darkorange', linestyle='--')
        ax.axvline(x = bounds_dict['lower_overlap']['upper_energy_high'], color='darkorange', linestyle='--')
        ax.fill_betweenx([0, 5*max(data_y[0])], bounds_dict['lower_overlap']['lower_energy_low'], bounds_dict['lower_overlap']['upper_energy_low'],
                         color='r', alpha=0.25, label='Overlap region 1') 

        if len(measurements) == 3:
            ax.axvline(x = bounds_dict['upper_overlap']['lower_energy_low'], color=color_map['overlap_high'], linestyle='--')
            ax.axvline(x = bounds_dict['upper_overlap']['upper_energy_low'], color=color_map['overlap_high'], linestyle='--')
            ax.axvline(x = bounds_dict['upper_overlap']['lower_energy_high'], color='darkblue', linestyle='--')
            ax.axvline(x = bounds_dict['upper_overlap']['upper_energy_high'], color='darkblue', linestyle='--')
            ax.fill_betweenx([0, 5*max(data_y[0])], bounds_dict['upper_overlap']['lower_energy_low'], bounds_dict['upper_overlap']['upper_energy_low'],
                            color='b', alpha=0.25, label='Overlap region 2')

        ax.legend()
        ax.set_xscale('log')
        ax.set_xlim(0.01,1000)
        ax.set_ylim(0, 1.1*max(np.max(data_y[i]) for i in range(len(data_y))))

        ax.set_title('Unscaled spectra - Overlap regions')

        plt.show()

    # -----------------Gain matching-----------------

    # First region: Match from left to right
    # Scale HIGH to MID or MID to LOW
    scaling_factor_low = float(minimize(_gain_match_chi2, scaling_factors[0],
                                  args=(bounds_dict['lower_overlap']['y_axis_high'], bounds_dict['lower_overlap']['y_axis_low'])).x)
    # Second region: Match from right to left
    # Scale LOW to MID or MID to HIGH
    if len(measurements) == 3:
        scaling_factor_high = float(minimize(_gain_match_chi2, scaling_factors[1],
                                       args=(bounds_dict['upper_overlap']['y_axis_low'], bounds_dict['upper_overlap']['y_axis_high'])).x)

    if testplots:
        fig, ax = plt.subplots()

        ax.step(data_x[0], np.multiply(data_y[0],scaling_factor_low), label = f"{measurements[0].gain}*{scaling_factor_low:.3f}", color = color_map[measurements[0].gain])
        ax.step(data_x[1], data_y[1], label = f"{measurements[1].gain}", color = color_map[measurements[1].gain])
        if len(measurements) == 3:
            ax.step(data_x[2], np.multiply(data_y[2], scaling_factor_high), label = f"{measurements[2].gain}*{scaling_factor_high:.3f}", color = color_map[measurements[2].gain])

        ax.axvline(x = bounds_dict['lower_overlap']['lower_energy_low'], color=color_map['overlap_low'], linestyle='--', label = 'Overlap region 1 - low')
        ax.axvline(x = bounds_dict['lower_overlap']['upper_energy_low'], color=color_map['overlap_low'], linestyle='--')
        ax.fill_betweenx([0, 5*max(data_y[0])], bounds_dict['lower_overlap']['lower_energy_low'], bounds_dict['lower_overlap']['upper_energy_low'], color='r', alpha=0.25)

        if len(measurements) == 3:
            ax.axvline(x = bounds_dict['upper_overlap']['lower_energy_low'], color=color_map['overlap_high'], linestyle='--', label = 'Overlap region 2 - low')
            ax.axvline(x = bounds_dict['upper_overlap']['upper_energy_low'], color=color_map['overlap_high'], linestyle='--')
            ax.fill_betweenx([0, 5*max(data_y[0])], bounds_dict['upper_overlap']['lower_energy_low'], bounds_dict['upper_overlap']['upper_energy_low'], color='b', alpha=0.25)
        
        ax.legend()
        ax.set_xscale('log')
        ax.set_xlim(0.01,1000)
        ax.set_ylim(0, 1.5)
        plt.show()

    #-----------------Stitch the spectra together-----------------

    stitching_idx_low, _ = _find_nearest_idx_and_value(data_x[0], stitching_points[0])
    data_x[0] = data_x[0][:stitching_idx_low]
    data_y[0] = np.multiply(data_y[0][:stitching_idx_low], scaling_factor_low)

    if len(measurements) == 2:
        stitching_idx_high, _ = _find_nearest_idx_and_value(data_x[1], stitching_points[0])
        data_x[1] = data_x[1][stitching_idx_high:]
        data_y[1] = data_y[1][stitching_idx_high:]
    
    if len(measurements) == 3:
        stitching_idx_med_1, _ = _find_nearest_idx_and_value(data_x[1], stitching_points[0])
        stitching_idx_med_2, _ = _find_nearest_idx_and_value(data_x[1], stitching_points[1])   
        data_x[1] = data_x[1][stitching_idx_med_1:stitching_idx_med_2]
        data_y[1] = data_y[1][stitching_idx_med_1:stitching_idx_med_2]

        stitching_idx_high, _ = _find_nearest_idx_and_value(data_x[2], stitching_points[1])
        data_x[2] = data_x[2][stitching_idx_high:]
        data_y[2] = np.multiply(data_y[2][stitching_idx_high:], scaling_factor_high)

    if testplots:
        fig, ax = plt.subplots()

        ax.step(data_x[0],data_y[0], label = f"{measurements[0].gain}*{scaling_factor_low:.3f}", color = color_map[measurements[0].gain])
        ax.step(data_x[1], data_y[1], label = f"{measurements[1].gain}", color = color_map[measurements[1].gain])
        if len(measurements) == 3:
            ax.step(data_x[2], data_y[2], label = f"{measurements[2].gain}*{scaling_factor_high:.3f}", color = color_map[measurements[2].gain])

        ax.axvline(x = stitching_points[0], color=color_map['overlap_low'], linestyle='--', label = 'Stitching point 1', lw=2)
        if len(measurements) == 3:
            ax.axvline(x = stitching_points[1], color=color_map['overlap_high'], linestyle='--', label = 'Stitching point 2', lw=2)
        
        ax.legend()
        ax.set_xscale('log')
        ax.set_xlim(0.01,1000)
        ax.set_ylim(0, 1.1*max(np.max(data_y[i]) for i in range(len(data_y))))
        plt.show()

    #-----------------Merge the spectra-----------------    
    merged_spectrum_x = np.concatenate((data_x[0], data_x[1]))
    if len(measurements) == 3:
        merged_spectrum_x = np.concatenate((merged_spectrum_x, data_x[2]))

    merged_spectrum_y = np.concatenate((data_y[0], data_y[1]))
    if len(measurements) == 3:
        merged_spectrum_y = np.concatenate((merged_spectrum_y, data_y[2]))

    #Create new Dataframe
    merged_df = pd.DataFrame()
    merged_df[measurements[0].x_axis] = merged_spectrum_x
    merged_df[measurements[0].y_axis] = merged_spectrum_y

    #Create new measurement object
    merged_spectrum = Measurement()
    merged_spectrum._data = merged_df

    merged_spectrum._x_axis = measurements[0].x_axis
    merged_spectrum._y_axis = measurements[0].y_axis
    merged_spectrum._num_channels = len(merged_spectrum_x)
    merged_spectrum._detector = measurements[0].detector
    merged_spectrum._particle = measurements[0].particle
    merged_spectrum._gain = 'merged'

    # What about date, live time, dead time, etc.?

    print(f'Spectra have been merged')

    return merged_spectrum 

def retrieve_original_spectrum(measurement):
    """ Retrieve the original spectrum again (from the original_data attribute)
    
    Args:
        measurement: The spectrum to be analyzed
    """

    # New dataframe
    measurement.data = pd.DataFrame()

    if measurement.gain == 'SIMULATION': #for ROOT files
        measurement.x_axis = 'ENERGY'
        measurement.y_axis = 'COUNTS'

        # Create new columns first (Else it will not work with the ROOT files)
        measurement.data[measurement.x_axis] = np.zeros(len(measurement.original_data[measurement.x_axis]))
        measurement.data[measurement.y_axis] = np.zeros(len(measurement.original_data[measurement.y_axis]))

        measurement.data[measurement.x_axis].update(measurement.original_data[measurement.x_axis])
        measurement.data[measurement.y_axis].update(measurement.original_data[measurement.y_axis])
    else:
        measurement.x_axis = 'CHANNEL'
        measurement.y_axis = 'COUNTS'

        measurement.data[measurement.x_axis] = measurement.original_data[measurement.x_axis].tolist()
        measurement.data[measurement.y_axis] = measurement.original_data[measurement.y_axis].tolist()
        measurement.data['mV'] = measurement.original_data['mV']
   
def add_spectra(measurements_list:list):
    """ Add the counts of an arbitrary number of spectra
        This can be used to obtain a more distinct particle edge for the calibration procedure
        
        Args:
            measurements_list: List of spectra to be added
    """

    #Check if all spectra have the same length
    if all(len(measurements_list[i].data) == len(measurements_list[0].data) for i in range(1,len(measurements_list))):
        print('All spectra have the same length')
    else:
        raise ValueError('The spectra have to have the same length')

    #Check if all spectra have the same x-axis
    if all(measurements_list[i].x_axis == measurements_list[0].x_axis for i in range(1,len(measurements_list))):
        print('All spectra have the same x-axis')
    else:
        raise ValueError('The spectra have to have the same x-axis')

    #Check if all spectra have the same y-axis
    if all(measurements_list[i].y_axis == measurements_list[0].y_axis for i in range(1,len(measurements_list))):
        print('All spectra have the same y-axis')
    else:
        raise ValueError('The spectra have to have the same y-axis')
    
    #Add the spectra
    x = measurements_list[0].data[measurements_list[0].x_axis].to_numpy()
    tot_sum = np.zeros(len(x), dtype=float)

    # Keep mV axis if it exists
    if 'mV' in measurements_list[0].data.columns:
        mV = measurements_list[0].data['mV'].to_numpy()

    for i in range(len(measurements_list)):
        print(f'Adding spectrum {i+1} of {len(measurements_list)}')
        data = pd.to_numeric(measurements_list[i].data[measurements_list[i].y_axis], errors='coerce')
        tot_sum += data

    #Create new Dataframe
    sum_df = pd.DataFrame()
    sum_df[measurements_list[0].x_axis] = x
    sum_df[measurements_list[0].y_axis] = tot_sum

    if 'mV' in measurements_list[0].data.columns:
        sum_df['mV'] = mV

    #Create new measurement object
    sum_spectrum = Measurement()
    sum_spectrum._data = sum_df

    sum_spectrum._x_axis = measurements_list[0].x_axis
    sum_spectrum._y_axis = measurements_list[0].y_axis
    sum_spectrum._num_channels = measurements_list[0].num_channels
    sum_spectrum._date = measurements_list[0].date
    sum_spectrum._detector = measurements_list[0].detector
    sum_spectrum._gain = measurements_list[0].gain
    sum_spectrum._particle = measurements_list[0].particle

    print(f'Spectra have been added')

    return sum_spectrum