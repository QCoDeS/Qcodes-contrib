import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt
from math import factorial, sqrt
from scipy import signal
from . import plot_cf_data, get_latest_counter

# TODO: get vna_helper funtions to 'inherit' from general


def resonator_sweep_setup(v1, power=-30, pm_range=200e6, avg=5,
                          bw=1000, npts=2001, centre=7.25e9):
    """
    Function which sets up the VNA with default settings
    for full resonator sweep.

    Args:
        v1 (instrument): VNA instrument to set params on
        power (int): default -30
        pm_range (float): default -30
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 2001
        centre (float): default 7.25e9
    """
    v1.rf_on()
    v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.power(power)
    v1.start(centre - pm_range)
    v1.stop(centre + pm_range)
    v1.npts(npts)


def power_sweep_setup(v1, avg=5, bw=1000, npts=201):
    """
    Function which sets up the VNA with default settings for
    power sweep on a resonator.

    Args:
        v1 (instrument): VNA instrument to set params on
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 201
    """
    v1.rf_on()
    v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.npts(npts)


def do_power_sweep(v1, centre, pm_range=10e6,
                   pow_start=-10, pow_stop=-50, pow_step=1):
    """
    Function which sweeps power in a loop and takes vna trace data

    Args:
        v1 (instrument): instrument to sweep power and frequency on
        centre (float): centre of frequency scan
        pm_range (float): plus minus range
        ie scan will be from centre-pm_range to centre+pm_range
        default 10e6
        pow_start (int): dBm starting power, default -10
        pow_stop (int): dBm finishing power default -60
        pow_step (int): dBm powe step, default 1

    Returns:
        dataset (qcodes DataSet)
        plot (QtPlot)
    """
    v1.start(centre - pm_range)
    v1.stop(centre + pm_range)
    loop = qc.Loop(v1.power.sweep(
        pow_start, pow_stop, pow_step)).each(v1.trace)
    dataset = loop.get_data_set()
    plot = qc.QtPlot(figsize=(700, 500))
    plot.add(dataset.VNA_magnitude)
    _ = loop.with_bg_task(plot.update, plot.save).run()
    dataset.data_num = get_latest_counter()
    plot.data_num = get_latest_counter()
    return dataset, plot


def gate_sweep_setup(v1, avg=5, bw=1000, npts=201, power=-40):
    """
    Function which sets up the VNA with default settings for power sweep
    on a resonator.

    Args:
        v1 (instrument): VNA instrument to set params on
        avg (int): default 5
        bw (int): default 1000
        npts (int): default 201
        power (int): default -40
    """
    v1.rf_on()
    v1.spec_mode('off')
    v1.bandwidth(bw)
    v1.avg(avg)
    v1.npts(npts)
    v1.power(power)


def gates_to_zero(dec_chans):
    """
    Sets all gate channels to 0
    """
    for chan in dec_chans:
        chan(0)


def do_gate_sweep(v1, centre, chan, reset_after=True, pm_range=10e6,
                  gate_start=0, gate_stop=-0.1, gate_step=0.01):
    """
    Function which sweeps gate voltage in a loop and takes vna trace data

    Args:
        v1 (instrument): VNA instrument to 'do frequency sweep' on
                         and get data from
        centre (float): centre of frequency scan
        chan (int): decadac channel controlling gate to sweep
        reset_after (bool): reset gate to gate_start after?
        pm_range (float): plus minus range
        ie scan will be from centre-pm_range to centre+pm_range
        default 10e6
        gate_start (float): V starting gate voltage, default 0
        gate_stop (float): V finishing gate voltaget -0.1 volts
        gate_step (float): V powe step, default 10 milivolt

    Returns:
        dataset (qcodes DataSet)
        plot (QtPlot)
    """
    v1.start(centre - pm_range)
    v1.stop(centre + pm_range)
    loop = qc.Loop(chan.sweep(gate_start, gate_stop, gate_step)).each(v1.trace)
    dataset = loop.get_data_set()
    plot = qc.QtPlot(figsize=(700, 500))
    plot.add(dataset.VNA_magnitude)
    _ = loop.with_bg_task(plot.update, plot.save).run()
    if reset_after:
        chan(gate_start)
    dataset.data_num = get_latest_counter()
    plot.data_num = get_latest_counter()
    return dataset, plot


def smooth_data_SG(y, window_size, order, deriv=0, rate=1):
    """
    Function which does savitzky_golay data smoothing method

    Args:
        y (numpy array): array for filtering
        window_size (int): size of window in points for smoothing
        order (int): polynomial order to use for smoothing

    Returns:
        filtered array
    """
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half_window = (window_size - 1) // 2

    # precompute coefficients
    b = np.mat([[k**i for i in order_range]
                for k in range(-half_window, half_window + 1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)

    # pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))

    return np.convolve(m[::-1], y, mode='valid')


def butter_lowpass(cutoff, fs, order):
    """
    Function which generates scipy butterworth filter to data

    Args:
        cutoff (float): frequency above which components should be removed
        order (int): filter order

    Returns:
        filter coefs: (numerator and denominator polynomials of IIR filter)
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def smooth_data_butter(data, fs, cutoff, order):
    """
    Function which applies butterwoth filter to data

    Args:
        data (numpy array)
        cutoff (float): passed to butter_lowpass

    Returns:
        filtered data
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


def find_peaks(dataset, fs, cutoff=30e9, order=5,
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
    setpoints = dataset.frequency

    # smooth data
    unsmoothed_data = dataset.VNA_magnitude
    smoothed_data = smooth_data_butter(
        unsmoothed_data, fs, cutoff=cutoff, order=order)

    # find peak indices
    peakind = signal.find_peaks_cwt(np.multiply(smoothed_data, -1), widths)

    # plot: unsmoothed data, smoothed data and add peak estimate values
    subplot = plot_cf_data(unsmoothed_data,
                           smoothed_data,
                           xdata=setpoints,
                           datanum=dataset.data_num,
                           axes_labels=['frequency(Hz)', 'S21'])
    subplot.plot(setpoints[peakind], smoothed_data[peakind], 'gs')
    txt = '{} resonances found at {}'.format(len(peakind), setpoints[peakind])

    fig = subplot.figure
    try:
        fig.suptitle('dataset {}'.format(dataset.data_num), fontsize=12)
    except AttributeError as e:
        print('cannot set title, dataset has no data_num: {}'.format(e))

    # TODO save this info!!
    print(txt)
    return peakind, setpoints[peakind], subplot


def plot_resonances(dataset, indices, subplot=None):
    """
    Function which does simple plot of data and resonance points

    Args:
        dataset (qcodes DataSet)
        indices (array): array of data indices of resonances

    Returns:
        subplot (matplotlib AxesSubplot): plot of results
    """
    if subplot is None:
        fig = plt.figure()
        subplot = plt.subplot(111)
        try:
            fig.data_num = dataset.data_num
        except AttributeError as e:
            print('dataset has no data_num set: {}'.format(e))

    setpoints = dataset.frequency
    magnitude = dataset.VNA_magnitude
    subplot.plot(setpoints, magnitude, 'b')
    subplot.plot(setpoints[indices], magnitude[indices], 'gs')
    subplot.set_xlabel('frequency(Hz)')
    subplot.set_ylabel('S21')
    subplot.figure.suptitle('dataset {}'.format(fig.data_num), fontsize=12)
    # TODO save this info!!
    return subplot


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
    mag_arrays = dataset.VNA_magnitude
    freq_array = dataset.frequency[0]
    pow_array = dataset.VNA_power_set
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
        fig.data_num = dataset.data_num
        fig.suptitle('dataset {}'.format(fig.data_num), fontsize=12)
        fig.text(0, 0, 'bare res: {}, pushed res: {}, push: {}'.format(
            high_res, low_res, dif))
    except AttributeError as e:
        print('dataset has no data_num set: {}'.format(e))

    return [low_res, high_res, dif], fig


def qubit_from_push(g, push, bare_res):
    """
    Get estimated qubit location given coupling, push on resonator
    and bare resonator position.

    Args:
        g: coupling in Hz
        push: push in Hz
        bare_res: resonator position in Hz

    Returns:
        estimated qubit frequency
    """
    delta = g**2 / push
    return bare_res - delta


def g_from_qubit(qubit, bare_res, push):
    """
    Get estimated coupling strength given qubit frequency, push
    on resonator and bare resonator position.

    Args:
        g: coupling in Hz
        push: push in Hz
        bare_res: resonator position in Hz

    Returns:
        estimated qubit frequency
    """
    delta = bare_res - qubit
    return sqrt(delta * push)
