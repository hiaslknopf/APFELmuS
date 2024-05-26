__author__ = "Matthias Knopf"
__email__ = "matthias.knopf@tuwien.ac.at"
__date__ = "2024"

from MicroDosimetry import Measurement
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft
from sklearn.preprocessing import normalize

""" This module contains all the functions for the manipulation of spectra """

def _lin_2_log(measurement, n_bins_per_decade):
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

def cutoff(measurement, channels):
    """ Low energy cutoff for noisy channels
    
    Args:
        measurement: The spectrum to be analyzed
        channels: Number of channels/values to be set zero
    """

    measurement.data.loc[:channels, measurement.y_axis] = 0
    print(f'Channels {channels} and below cut off: {measurement.name}')

def logarithmic_binning(measurement, n_bins_per_decade):
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

def probability_function(measurement, type_of_dist):
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

    #TODO: Why does Sandra do this? Shouldn't every Spectrum be normalized to A=1
    #if measurement.y_axis in ['yf(y)', 'yd(y)']:
    #    #Extra y needs to be stripped away
    #    pdf = np.divide(pdf,x)

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

def extrapolate(measurement, method):
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

def merge_spectra(measurements):
    """ Merge several spectra (df objects) into one dataframe object and match the gains
    
    Args:
        measurements: The list of spectra to be merged
    """

    #---------Check if data matches-----------
    check_x = []
    check_y = []

    if len(measurements) != 3:
        raise ValueError('You can only merge 3 spectra: LOW, MID, HIGH')

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
        raise ValueError('The spectra have to have a calibrated lineal energy axis')
    
    if measurements[0].gain == 'SIMULATION':
        raise ValueError('The spectra have to be measured spectra, not simulations')

    #---------- Find best fit for a common x-axis ----------
    #TODO: Linear fit siehe Sandra (dass gleiches sampling für alle 3)



    #---------- Algorithm for overlap and gain match -----------
    #TODO: Implement algorithms for finding these or choose manually

    #Overlap regions
    LM_overlap = [8,5]
    MH_overlap = [20,10]

    #Gain of spectra = scaled y_axis
    scale_LM = 0.4642104662836921 #Low gain gets scaled
    scale_MH = 1.0282117319899817 #Mid gain gets scaled

    #-------------Have a look at the data-----------------
    fig, ax = plt.subplots()

    for i in range(len(measurements)):

        x = measurements[i].data[measurements[i].x_axis]
        y = measurements[i].data[measurements[i].y_axis]

        if measurements[i].gain == 'LOW':
            y *= scale_LM
            ax.plot(x, y, label = f"{measurements[i].gain} * {scale_LM: .2f}")
        elif measurements[i].gain == 'HIGH':
            y *= scale_MH
            ax.plot(x, y, label = f"{measurements[i].gain} * {scale_MH: .2f}")
        else:
            ax.plot(x, y, label = f"{measurements[i].gain}")

    ax.set_xlim(0.01,1000)
    ax.set_ylim(0.0,1.1*np.max(y))

    if measurements[0].x_axis == 'LINEAL ENERGY':
        ax.set_xlabel('y [keV/um]')
    elif measurements[0].x_axis == 'ENERGY':
        ax.set_xlabel('E [keV]')

    ax.set_ylabel(measurements[0].y_axis)
    ax.set_title('Merged spectra - overlap')

    #Current overlap regions
    ax.axvline(x = LM_overlap[0], color='r', linestyle='--', label = 'LM overlap region')
    ax.axvline(x = LM_overlap[1], color='r', linestyle='--')
    ax.axvline(x = MH_overlap[0], color='g', linestyle='--', label='MH overlap region')
    ax.axvline(x = MH_overlap[1], color='g', linestyle='--')

    ax.legend()
    ax.set_xscale('log')

    plt.show()

    #-----------------Merge the data into a new measurement-------------------

    #TODO: Take boundaries from algorithm
    upper_HG = 175
    lower_MG = 20
    upper_MG = 80
    lower_LG = 100

    x_list = []
    y_list = []


    fig, ax = plt.subplots()

    for i in range(len(measurements)):

        x = measurements[i].data[measurements[i].x_axis].tolist()
        y = measurements[i].data[measurements[i].y_axis].tolist()

        if measurements[i].gain == 'LOW':
            ax.plot(x[lower_LG:], y[lower_LG:], label = f"{measurements[i].gain}")
            ax.axvline(x = x[lower_LG], color='b', linestyle='--')
        elif measurements[i].gain == 'HIGH':
            ax.plot(x[:upper_HG], y[:upper_HG], label = f"{measurements[i].gain}")
            ax.axvline(x = x[upper_HG-1], color='y', linestyle='--')
        else:
            ax.plot(x[lower_MG:upper_MG], y[lower_MG:upper_MG], label = f"{measurements[i].gain}")
            ax.axvline(x = x[lower_MG], color='r', linestyle='--')
            ax.axvline(x = x[upper_MG-1], color='r', linestyle='--')

    ax.set_xlim(0.01,1000)
    ax.set_ylim(0.0,1.1*np.max(y))

    if measurements[0].x_axis == 'LINEAL ENERGY':
        ax.set_xlabel('y [keV/um]')
    elif measurements[0].x_axis == 'ENERGY':
        ax.set_xlabel('E [keV]')

    ax.set_ylabel(measurements[0].y_axis)
    ax.set_title('Merged spectra - Cutoff')

    ax.legend()
    ax.set_xscale('log')

    plt.show()


    for i in range(len(measurements)):

        x = measurements[i].data[measurements[i].x_axis].tolist()
        y = measurements[i].data[measurements[i].y_axis].tolist()

        if measurements[i].gain == 'LOW':
            x_list.append(x[lower_LG:])
            y_list.append(y[lower_LG:])
        elif measurements[i].gain == 'HIGH':
            x_list.append(x[:upper_HG])
            y_list.append(y[:upper_HG])
        else:
            x_list.append(x[lower_MG:upper_MG])
            y_list.append(y[lower_MG:upper_MG])

    x_list = [item for sublist in x_list for item in sublist]
    y_list = [item for sublist in y_list for item in sublist]

    #Create new Dataframe
    merged_df = pd.DataFrame()
    merged_df[measurements[0].x_axis] = x_list
    merged_df[measurements[0].y_axis] = y_list
    
    merged_df = merged_df.sort_values(measurements[0].x_axis)
    merged_df = merged_df.reset_index(drop=True)

    #Create new measurement object
    merged_spectrum = Measurement()
    merged_spectrum._data = merged_df

    merged_spectrum._x_axis = measurements[0].x_axis
    merged_spectrum._y_axis = measurements[0].y_axis
    merged_spectrum._num_channels = measurements[0].num_channels
    merged_spectrum._date = measurements[0].date
    merged_spectrum._detector = measurements[0].detector
    merged_spectrum._gain = 'merged'

    print(f'Spectra have been merged')

    return merged_spectrum 

def retrieve_original_spectrum(measurement):
    """ Retrieve the original spectrum again (from the original_data attribute)
    
    Args:
        measurement: The spectrum to be analyzed
    """

    measurement.data = pd.DataFrame()

    if measurement.gain == 'SIMULATION': #for ROOT files
        measurement.x_axis = 'ENERGY'
        measurement.y_axis = 'COUNTS'

        measurement.data[measurement.x_axis].update(measurement.original_data[measurement.x_axis])
        measurement.data[measurement.y_axis].update(measurement.original_data[measurement.y_axis])
    else:
        measurement.x_axis = 'CHANNEL'
        measurement.y_axis = 'COUNTS'

        measurement.data[measurement.x_axis] = measurement.original_data[measurement.x_axis].tolist()
        measurement.data[measurement.y_axis] = measurement.original_data[measurement.y_axis].tolist()
        measurement.data['mV'] = measurement.original_data['mV']
   
def get_LET_distribution(measurement, chord_length_dist, LET_spec='pdf'):
    """ Calculate an LET distribution from a calibrated f(E) pdf
    Following an algorithm according to Kellerer 1972, Phys.Med.Biol., 17, NO.2, 232-240
    
    Args:
        measurement: The spectrum to be analyzed
        chord_length_dict: Chord length distribution information
        LET_spec: 'pdf' for t(LET) or 'weighted' for LET*t(LET)
    """

    testplot = False

    #TODO: Check if chord_length is already attached

    if measurement.x_axis != 'ENERGY' or measurement.y_axis not in ['COUNTS', 'F(E)', 'f(E)', 'Ef(E)']:
        raise ValueError('This calculation only works for energy calibrated spectra')
    if measurement.mean_chord_length == 0.0:
        raise ValueError('The chord length function has to be called first')

    #Transform pdf data into right format Ef(E)
    if measurement.y_axis == 'COUNTS':
        probability_function(measurement, 'F')
        probability_density(measurement)
        weighted_probability_density(measurement)
    elif measurement.y_axis == 'F(E)':
        probability_density(measurement)
        weighted_probability_density(measurement)
    elif measurement.y_axis == 'f(E)':
        weighted_probability_density(measurement)
    
    #Get pdf data
    E = measurement.data[measurement.x_axis].to_numpy()
    phi = measurement.data[measurement.y_axis].to_numpy()

    #Normalize distributions for Fourier Transform
    gamma = np.multiply(chord_length_dist['pdf'], chord_length_dist['chord_len'])
    norm_gamma = normalize(gamma[:,np.newaxis], axis=0).ravel()
    norm_phi = normalize(phi[:,np.newaxis], axis=0).ravel()

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(chord_length_dist['chord_len'], norm_gamma, label='chord length dist g(l)')
        ax.set_xlabel('Chord length [µm]')
        ax.legend()
        plt.show()

        fig, ax = plt.subplots()
        ax.plot(E, norm_phi, label='normalize phi = Ef(E)')
        ax.set_xlabel('Energy [keV]')
        ax.legend()
        plt.show()

    #Fourier transform
    phi_ft = rfft(norm_phi)
    gamma_ft = rfft(norm_gamma)
    rfftx = rfftfreq(measurement.num_channels, 1 / 2*measurement.num_channels)

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(rfftx, np.abs(phi_ft), label='phi*')
        ax.plot(rfftx, np.abs(gamma_ft), label='gamma*')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    #Calculate tau
    tau_ft = np.divide(phi_ft, gamma_ft, out=np.zeros_like(phi_ft), where=gamma_ft!=0)

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(rfftx, np.abs(tau_ft), label='tau*')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    #Backtrafo
    tau = irfft(tau_ft)
    LET_pdf = np.divide(tau, chord_length_dist['chord_len'], out=np.zeros_like(tau), where=chord_length_dist['chord_len']!=0)

    #Get rid of the negative values
    LET_pdf[LET_pdf < 0] = 0
    tau[tau < 0] = 0

    if testplot:
        fig, ax = plt.subplots()
        ax.plot(LET_pdf, label='LET_pdf = t(LET)')
        ax.plot(tau, label='tau = LET*t(LET)')
        ax.set_xlabel('Discrete step [a.u.]')
        ax.legend()
        plt.show()

    measurement.data.drop(columns=measurement.x_axis, inplace=True)
    measurement.data.drop(columns=measurement.y_axis, inplace=True)

    #TODO: Something is not right about the x-axis
    #TODO: Kellerer mentions this: Intervall of periodiciallity has to be chosen accordingly, that all spectra fit
    measurement.data[measurement.x_axis] = np.divide(E, chord_length_dist['chord_len'], out=np.zeros_like(E), where=chord_length_dist['chord_len']!=0)
    measurement.x_axis = 'LET'

    if LET_spec == 'pdf':
        measurement.data[measurement.y_axis] = LET_pdf
        measurement.y_axis = 't(LET)'
    elif LET_spec == 'weighted':
        measurement.data[measurement.y_axis] = tau
        measurement.y_axis = 'LETt(LET)'
    
    print(f'LET distribution calculated from energy calibrated spectrum: {measurement.name}')