import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy import signal
from . import exp_decay, exp_decay_sin, get_calibration_dict, get_title, \
    save_fig, smooth_data_butter, smooth_data_SG, plot_cf_data, \
    get_sample_name, g_from_qubit

# TODO: write fit functions: qubit_from_ssb_measure,
#                            qubit_from_ssb_power_sweep,
#                            qubit_from_ssb_volt_sweep
# TODO: docstrings
# TODO: remove hard coding of linear_magnitude


###########################
# VNA
###########################


def find_peaks(dataset, fs, cutoff=0.2e-6, order=5,
               subplot=None, widths=np.linspace(1, 150)):
    """
    Function which given a 1d array smoothes the data and finds resonances
    and plots results

    Args:
        dataset (qcodes DataSet)
        fs (float): frequency of sampling, passed to smoothed_data_butter
        cutoff (float): used for smoothing passed to smooth_data_butter
                    default 30e9
        order (int): used for smoothing passed to smooth_data_butter, default 5
        subplot (matplotlib AxesSubplot): subplot which this data should be
                    plotted on, default None will create new one
        widths (array): array peak widths to search for, passed to
                    signal.find_peaks_cwt, default np.linspace(1, 150)

    Returns:
        peakind (array): indices of resonances found
        frequencies (array): frequencies of resonances found
        subplot (matplotlib AxesSubplot): plot of results
    """
    setpoints = next(getattr(dataset, key)
                     for key in dataset.arrays.keys() if "set" in key)

    # smooth data
    unsmoothed_data = dataset.vna_linear_magnitude
    smoothed_data = smooth_data_butter(
        unsmoothed_data, fs, cutoff=cutoff, order=order)

    # find peak indices
    peakind = signal.find_peaks_cwt(np.multiply(smoothed_data, -1), widths)

    try:
        num = dataset.data_num
    except AttributeError:
        num = dataset.location_provider.counter
        print('warning: check title, could be wrong datanum')

    # plot: unsmoothed data, smoothed data and add peak estimate values
    fig, subplot = plot_cf_data([unsmoothed_data, smoothed_data],
                                xdata=setpoints,
                                data_num=num,
                                axes_labels=['frequency(Hz)', 'S21'])
    subplot.plot(setpoints[peakind], smoothed_data[peakind], 'gs')
    txt = '{} resonances found at {}'.format(len(peakind), setpoints[peakind])

    fig.suptitle('{}_{}_find_peaks'.format(num, get_sample_name()),
                 fontsize=12)

    # TODO save this info!!
    print(txt)
    return peakind, setpoints[peakind], subplot


def get_resonator_push(dataset):
    """
    Function which gets the change in resonance frequency from a power
    sweep dataset.

    Args:
        dataset (qcodes DataSet)

    Returns:
        low_res (float): calculated resonance freq at low power
        high_res (float): calculated resonance freq at low power
        dif (float): value of push in Hz
        axarr (numpy.ndarray): subplot array
    """
    # get data for high and low power from dataset
    mag_arrays = dataset.vna_linear_magnitude
    freq_array = dataset.frequency_set[0]
    pow_array = dataset.vna_power_set
    mag_high = mag_arrays[0]
    mag_low = mag_arrays[-1]
    smoothed_mag_low = smooth_data_SG(mag_low, 15, 6)
    smoothed_mag_high = smooth_data_SG(mag_high, 15, 6)

    # get freqeuncy of resonances for low and high power and power values
    low_res = freq_array[smoothed_mag_low.argmin()]
    high_res = freq_array[smoothed_mag_high.argmin()]
    low_pow = pow_array[-1]
    high_pow = pow_array[0]
    dif = low_res - high_res

    # for all pow sweeps smooth and get resonance
    res_data = np.zeros(len(mag_arrays))
    for i, sweep in enumerate(mag_arrays):
        smoothed_data = smooth_data_SG(sweep, 15, 6)
        res_data[i] = freq_array[smoothed_data.argmin()]

    # plot
    fig, axarr = plt.subplots(2)

    # subplot 1: high and low power cuts, smoothed and unsmoothed
    plot_cf_data(smoothed_mag_high, mag_high, smoothed_mag_low, mag_low,
                 subplot=axarr[0], xdata=freq_array,
                 legend_labels=['pow={}'.format(high_pow),
                                'pow={},smoothed'.format(high_pow),
                                'pow={}'.format(low_pow),
                                'pow={}, smoothed'.format(low_pow)],
                 axes_labels=['frequency (Hz)', 'S21'])

    # subplot 2: resonance for all power sweep with max and min lines plotted
    axarr[1].plot(pow_array, res_data, 'k-')
    axarr[1].plot([high_pow, low_pow], [high_res, high_res], 'r', lw=2,
                  label='high power res = {}'.format(high_res))
    axarr[1].plot([high_pow, low_pow], [low_res, low_res], 'b', lw=2,
                  label='low power res = {}'.format(low_res))
    axarr[1].set_xlabel('power (dBm)')
    axarr[1].set_ylabel('frequency (Hz)')
    axarr[1].legend(loc='upper right', fontsize=10)

    plt.tight_layout()

    try:
        # TODO: replace with data.location_provider.counter?
        fig.data_num = dataset.data_num
        fig.suptitle('dataset {}'.format(fig.data_num), fontsize=12)
        fig.text(0, 0, 'bare res: {}, pushed res: {}, push: {}'.format(
            high_res, low_res, dif))
    except AttributeError as e:
        fig.data_num = dataset.location_provider.counter
        print('dataset has no data_num set: {}'.format(e))

    return [low_res, high_res, dif], fig


###########################
# Alazar
###########################

def find_extreme(data, x_key="freq", y_key="mag", extr="min"):
    try:
        x_key_array_name = [v for v in data.arrays.keys() if x_key in v][0]
    except IndexError:
        raise KeyError('keys: {} not in data array '
                       'names: {}'.format(x_key,
                                          list(data.arrays.keys())))
    try:
        y_key_array_name = [v for v in data.arrays.keys() if y_key in v][0]
    except IndexError:
        raise KeyError('keys: {} not in data array '
                       'names: {}'.format(y_key,
                                          list(data.arrays.keys())))

    x_data = getattr(data, x_key_array_name)
    y_data = getattr(data, y_key_array_name)
    if extr is "min":
        index = np.argmin(y_data)
        val = np.amin(y_data)
    elif extr is "max":
        index = np.argmax(y_data)
        val = np.amax(y_data)
    else:
        raise ValueError('extr must be set to "min" or "max", given'
                         ' {}'.format(extr))
    extr_freq = x_data[index]
    return extr_freq, val


def recalculate_g(dec_chans=None):
    """
    Function which uses the values in the calibration dictionary for expected
    qubit position, actual position, resonator push and g value to recalculate
    the g value for the current qubit and compare it to the old value.

    Args:
        dec_chans (default None): if dec_chans given, gate value of current
            qubit is compared to value at which resonator data was taken to
            check validity of recalculated g
    """
    c_dict = get_calibration_dict()
    expected = c_dict['expected_qubit_positions'][c_dict['current_qubit']]
    actual = c_dict['actual_qubit_positions'][c_dict['current_qubit']]
    res_data = c_dict['resonator_pushes'][c_dict['current_qubit']]
    old_g = c_dict['g_values'][c_dict['current_qubit']]
    new_g = g_from_qubit(actual, res_data[0], res_data[2])
    if dec_chans is not None:
        current_voltage = dec_chans[c_dict['current_qubit']].get_latest()
        if (c_dict['gatability'][c_dict['current_qubit']] and
                (current_voltage !=
                    c_dict['gate_volts'][c_dict['current_qubit']])):
            print('New g factor calculated will not be a good estimate '
                  'as current gate value is not the same as the value when '
                  'the push on the resonator was measured.')
    print('expected qubit freq: {}\n (from g of {}, push on resonator {})\n'
          'actual qubit freq: {}\n (for same push gives g of {}'.format(
              expected, old_g, res_data[2], actual, new_g))
    return new_g


def qubit_from_ssb_measure(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError


def qubit_from_ssb_power_sweep(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError


def qubit_from_ssb_volt_sweep(dataset, gradient_sign=1, min_res_width=4e6,
                              high_voltage=True):
    raise NotImplementedError


def get_t2(data, x_name='delay', y_name='magnitude',
           counter=None, plot=True, subplot=None):
    """
    Function which fits results of a data set to a sine wave modulated
    by an exponential decay and returns the fit parameters and the standard
    deviation errors on them.

    Args:
        data (qcodes dataset): 1d sweep to be fit to
        x_name (str) (default 'delay'): x axis key used to search data.arrays
            for corresponding data
        y_name (str) (default 'magnitude'): y axis key
        counter (int) (default None): used to set title and name for saving if
            dataset does not have a data_num (which it should!)
        plot (default True)
        subplot (default None): subplot to plot in otherwise makes new figure
    """
    x_data = getattr(getattr(data, x_name), 'ndarray')
    y_data = getattr(getattr(data, y_name), 'ndarray')
    x_units = getattr(getattr(data, x_name), 'unit')
    y_units = getattr(getattr(data, y_name), 'unit')
    popt, pcov = curve_fit(exp_decay_sin, x_data, y_data,
                           p0=[0.003, 1e-7, 10e7, 0, 0.01])
    errors = np.sqrt(np.diag(pcov))
    print('fit to equation of form y = a * exp(-x / b) * sin(c * x + d) + e'
          'gives:\na {}, b {}, c {}, d {}, e{}\n'
          'with one standard deviation errors:\n'
          'a {}, b {}, c {}, d {}, e{}'.format(popt[0], popt[1], popt[2],
                                               popt[3], popt[4], errors[0],
                                               errors[1], errors[2],
                                               errors[3], errors[4]))
    if plot:
        if subplot is None:
            fig, ax = plt.subplots()
        else:
            ax = subplot
            fig = ax.figure
        try:
            num = data.data_num
        except AttributeError:
            num = counter
        try:
            qubit = get_calibration_dict()['current_qubit']
            title = '{}_{}_T2'.format(get_title(num), qubit)
            name = '{}_{}_T2'.format(num, qubit)
        except Exception:
            title = '{}_T2'.format(get_title(num))
            name = '{}_T2'.format(num)

        if (not hasattr(fig, "data_num")) and (counter is not None):
            fig.data_num = num
        ax.plot(x_data,
                exp_decay_sin(x_data, *popt),
                label='fit: T2 {}{}'.format(popt[1],
                                            x_units))
        ax.plot(x_data, y_data, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        save_fig(ax, name=name)
        return ax, popt, errors
    else:
        return popt, errors


def get_t1(data, x_name='delay', y_name='magnitude',
           counter=None, plot=True, subplot=None):
    """
    Function which fits results of a data set to an exponential decay and
    returns the fit parameters and the standard deviation errors on them.

    Args:
        data (qcodes dataset): 1d sweep to be fit to
        x_name (str) (default 'delay'): x axis key used to search data.arrays
            for corresponding data
        y_name (str) (default 'magnitude'): y axis key
        counter (int) (default None): used to set title and name for saving if
            dataset does not have a data_num (which it should!)
        plot (default True)
        subplot (default None): subplot to plot in otherwise makes new figure
    """
    x_data = getattr(getattr(data, x_name), 'ndarray')
    y_data = getattr(getattr(data, y_name), 'ndarray')
    x_units = getattr(getattr(data, x_name), 'unit')
    y_units = getattr(getattr(data, y_name), 'unit')
    popt, pcov = curve_fit(exp_decay, x_data, y_data, p0=[0.05, 1e-6, 0.01])
    errors = np.sqrt(np.diag(pcov))
    print('fit to equation of form y = a * exp(-x / b) + c gives:\n'
          'a {}, b {}, c {}\n'
          'with one standard deviation errors:\n'
          'a {}, b {}, c {}'.format(popt[0], popt[1], popt[2],
                                    errors[0], errors[1], errors[2]))
    if plot:
        if subplot is None:
            fig, ax = plt.subplots()
        else:
            ax = subplot
            fig = ax.figure
        try:
            num = data.data_num
        except AttributeError:
            num = counter
        try:
            qubit = get_calibration_dict()['current_qubit']
            title = '{}_{}_T1'.format(get_title(num), qubit)
            name = '{}_{}_T1'.format(num, qubit)
        except Exception:
            title = '{}_T1'.format(get_title(num))
            name = '{}_T1'.format(num)

        if (not hasattr(fig, "data_num")) and (counter is not None):
            fig.data_num = num
        ax.plot(x_data,
                exp_decay(x_data, *popt),
                label='fit: T1 {}{}'.format(popt[1],
                                            x_units))
        ax.plot(x_data, y_data, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        save_fig(ax, name=name)
        return ax, popt, errors
    else:
        return popt, errors
