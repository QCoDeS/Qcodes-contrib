import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from . import plot_data_single_window, plot_data, sweep1d, exp_decay, exp_decay_sin, \
    get_title, measure, save_fig, get_calibration_dict, get_calibration_val, set_calibration_val

# TODO: init decision for differentiating between vna functions and alazar
#   functions
# TODO: exception types
# TODO: TWPA settings
# TODO: remove do_ssb_pow_sweep, qubit_sweep_setup, do_ssb_time_sweep
#   do_ssb_gate_sweep do_ssb_gate_sweep
# TODO: write fit functions: qubit_from_ssb_measure, qubit_from_ssb_power_sweep,
#    qubit_from_ssb_volt_sweep, find_qubit, find_cavity_val
# TODO: _ ?
# TODO: Remove do_rabi_freq_sweep, do_rabi_pow_sweep
# TODO: find qubit to compare to average


def config_alazar(alazar, seq_mode=0, clock_source='EXTERNAL_CLOCK_10MHz_REF'):
    """
    Function which puts alazar in sequency mode, configures and sets clock
    source

    Args:
        alazar instrument
        seq mode (0 or 1) (default 0)
        clock_source (default 'EXTERNAL_CLOCK_10MHz_REF')
    """
    if seq_mode not in [0, 1]:
        raise ValueError('must set seq mode to 0 or 1')
    if seq_mode:
        io_mode = 'AUX_IN_TRIGGER_ENABLE'
        io_param = 'TRIG_SLOPE_POSITIVE'
    else:
        io_mode = 'AUX_IN_AUXILIARY'
        io_param = 'NONE'
    alazar.config(clock_source=clock_source,
                  sample_rate=500000000,
                  clock_edge='CLOCK_EDGE_RISING',
                  decimation=1,
                  coupling=['DC', 'DC'],
                  channel_range=[.4, .4],
                  impedance=[50, 50],
                  trigger_operation='TRIG_ENGINE_OP_J',
                  trigger_engine1='TRIG_ENGINE_J',
                  trigger_source1='EXTERNAL',
                  trigger_slope1='TRIG_SLOPE_POSITIVE',
                  trigger_level1=140,
                  trigger_engine2='TRIG_ENGINE_K',
                  trigger_source2='DISABLE',
                  trigger_slope2='TRIG_SLOPE_POSITIVE',
                  trigger_level2=128,
                  external_trigger_coupling='DC',
                  external_trigger_range='ETR_2V5',
                  trigger_delay=0,
                  timeout_ticks=0,
                  aux_io_mode=io_mode,
                  aux_io_param=io_param
                  )
    print('Configured alazar, seq mode {}'.format(seq_mode))
    print('-------------------------')


def get_alazar_seq_mode(alazar):
    """
    Gets the sequencing mode of alazar based on AUX_IN_TRIGGER_ENABLE and
    TRIG_SLOPE_POSITIVE

    Args:
        alazar instrument
    """
    if (alazar.aux_io_mode() is 'AUX_IN_TRIGGER_ENABLE' and
            alazar.aux_io_param() is 'TRIG_SLOPE_POSITIVE'):
        return 1
    elif (alazar.aux_io_mode() is 'AUX_IN_AUXILIARY' and
          alazar.aux_io_param() is 'NONE'):
        return 0
    else:
        raise ValueError('aux_io_mode: {}, aux_io_param: {} '
                         'do not correspond to seq_mode on or off')


def set_alazar_seq_mode(alazar, mode):
    """
    Sets the sequence mode of the alazar

    Args:
        alazar instrument
        mode (0 or 1)
    """
    if mode == 1:
        alazar.config(sample_rate=alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_TRIGGER_ENABLE',
                      aux_io_param='TRIG_SLOPE_POSITIVE')
    elif mode == 0:
        alazar.config(sample_rate=alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_AUXILIARY',
                      aux_io_param='NONE')
    else:
        raise ValueError('must set seq mode to 0 or 1')


def get_demod_freq(cavity, localos, acq_ctrl):
    """
    Gets demodulation frequency based on cavity and
    local oscilator and checks that aqc controller is demodulating
    at this frequency.

    Args:
        cavity instrument (r&s SGS100)
        localos instrument (r&s SGS100)
        acq_ctrl instrument (alazar acq controller)

    Returns:
        demod freq
    """
    lo, cav = localos.frequency(), cavity.frequency()
    demod = lo - cav
    acq_freqs = acq_ctrl.demod_freqs()
    if demod not in acq_freqs:
        raise Exception('demod freq {} (from cavity freq {} and localos '
                        'freq {}) not in acq controller demodulation '
                        'frequencies: {}.'.format(demod, cav, lo, acq_freqs))
    else:
        return demod


def set_single_demod_freq(cavity, localos, acq_ctrls, demod_freq,
                          cav_freq=None):
    """
    Optionally sets cavity frequency. Then sets detuning of localos and checks
    that acq controllers demodulate at this detuning

    Args:
        cavity instrument (r&s SGS100)
        localos instrument (r&s SGS100)
        acq_ctrls list of instruments (alazar acq controllers)
        demod freq (float) to be set
        cav freq (optional float) freq to set cavity to
    """
    if cav_freq is not None:
        cavity.frequency(cav_freq)
    else:
        cav_freq = cavity.frequency()
    localos.frequency(cav_freq + demod_freq)
    for ctrl in acq_ctrls:
        remove_demod_freqs(ctrl)
        ctrl.demod_freqs.add_demodulator(demod_freq)


def remove_demod_freqs(acq_ctrl):
    """
    Strips demod freqs from acquisition controller

    Args:
        acq_ctrl (alazar acq controller)
    """
    freqs = acq_ctrl.demod_freqs.get()
    for freq in freqs:
        acq_ctrl.demod_freqs.remove_demodulator(freq)


def do_cavity_freq_sweep(cavity, localos, cavity_freq, acq_ctrl,
                         cavity_pm=10e6, freq_step=1e6, demod_freq=None, live_plot=True,
                         key=None, save=True):
    """
    Function which sweeps the cavity frequency around central value by pm_range
    and measures using given acq_ctr, also steps local os to keep same demod
    freq.

    Args:
        cavity instrument (r&s SGS100)
        localos instrument (r&s SGS100)
        cavity_freq: central cavity drive freq
        acq_ctrl instrument (alazar acq controller)
        cavity_pm (float) (default 10e6): sweep range will be cavity_freq +-
            this value
        freq_step (float) (default 1e6)
        demod_freq (float) (default None): default uses the current value
        live_plot (default true)
        key (default None): string specifying specific parameter array to be
            plotted, default is to plot all
        save (default True): whether to save png on completion, nb if you
            choose not to live plot and the measured parameter returns
            multiple values then unless you specify a specific one to plot
            via 'key' then none will be saved.

    Args
    """
    if demod_freq is None:
        demod_freq = get_demod_freq(cavity, localos, acq_ctrl)
    loop = qc.Loop(cavity.frequency.sweep(cavity_freq - cavity_pm,
                                          cavity_freq + cavity_pm,
                                          freq_step)).each(
        qc.Task(localos.frequency.set,
                (cavity.frequency + demod_freq)),
        acq_ctrl.acquisition)
    if live_plot:
        dataset = loop.get_data_set()
        dataset.data_num = dataset.location_provider.counter
        plot = plot_data_single_window(dataset, acq_ctrl.acquisition, key=key)
        try:
            if save:
                _ = loop.with_bg_task(plot.update, plot.save).run()
            else:
                _ = loop.with_bg_task(plot.update).run()
                print('warning: plots not saved, if you want to save this'
                      'plot run plot.save()')
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        data.data_num = data.location_provider.counter
        plots = plot_data(data, key=key)
        if (key is not None) and save:
            plots.save()
        else:
            print('warning: plots not saved. To save one choose '
                  'and run plots[i].save()')
        return data, plots


def calibrate_cavity(cavity, localos, acq_ctrl, alazar, centre_freq=None, demod_freq=None,
                     calib_update=True, cavity_power=-35, lo_power=15, live_plot=True):
    if centre_freq is None:
        centre_freq = cavity.frequency()
    if demod_freq is None:
        demod_freq = get_demod_freq(cavity, localos, acq_ctrl)
    alazar.seq_mode(0)
    cavity.status('on')
    cavity.power(cavity_power)
    localos.power(lo_power)
    data, plot = do_cavity_freq_sweep(cavity, localos, centre_freq, acq_ctrl,
                                       cavity_pm=3e6, freq_step=0.1e6, demod_freq=demod_freq,
                                    live_plot=True, key="mag", save=True)
    good_cavity_freq = find_extreme(data, x_key="frequency_set", y_key="mag", extr="min")[0] + 3e5
    if calib_update:
        set_calibration_val('cavity_freqs', good_cavity_freq)
        set_calibration_val('cavity_pows', cavity_power)
        set_calibration_val('demod_freqs', demod_freq)
    set_single_demod_freq(cavity, localos, [acq_ctrl], demod_freq,
                          cav_freq=good_cavity_freq)
    print('cavity_freq set to {}'.format(good_cavity_freq))


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


def set_cavity_from_calib_dict(cavity, localos, acq_ctrls, num_avg=1000):
    for acq_ctrl in acq_ctrls:
        acq_ctrl.int_time(get_calibration_val('int_times'))
        acq_ctrl.int_delay(get_calibration_val('int_delays'))
        acq_ctrl.num_avg(num_avg)
    cavity.power(get_calibration_val('cavity_pows'))
    set_single_demod_freq(cavity, localos, acq_ctrls,
                          get_calibration_val('demod_freqs'),
                          cav_freq=get_calibration_val('cavity_freqs'))


def sweep_2d_ssb(qubit, acq_ctrl, centre_freq, sweep_param,
                 start, stop, step, delay=0.01, live_plot=True,
                 key=None, save=True):
    qubit.frequency(centre_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=centre_freq + 100e6,
                                            setpoints_stop=centre_freq - 100e6)
    return sweep1d(acq_ctrl.acquisition, sweep_param, start,
                   stop, step, delay=delay, live_plot=live_plot, key=key, save=save)


def qubit_from_ssb_measure(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError


def qubit_from_ssb_power_sweep(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError


def qubit_from_ssb_volt_sweep(dataset, gradient_sign=1, min_res_width=4e6,
                              high_voltage=True):
    raise NotImplementedError
    # voltage_array_name = [n for n in dataset.arrays.keys() if "voltage" in n][
    #     0]
    # magnitude_array_name = [n for n in dataset.arrays.keys() if "mag" in n][0]
    # frequency_array_name = [n for n in dataset.arrays.keys() if "ssb" in n][0]

    # if high_voltage:
    #     voltage_index = np.argmax(getattr(dataset, voltage_array_name))
    # else:
    #     voltage_index = np.argmin(getattr(dataset, voltage_array_name))

    # mag_array = getattr(dataset, magnitude_array_name)[voltage_index]
    # freq_array = getattr(dataset, frequency_array_name)[0]
    # return find_qubit(freq_array, mag_array,
    #                   gradient_sign=gradient_sign,
    #                   min_res_width=min_res_width)


def find_qubit(awg, alazar, acq_ctrl, qubit, start_freq=4e9, stop_freq=6e9, qubit_power=None, calib_update=True):
    if 'ssb_qubit' not in acq_ctrl.acquisition.setpoint_names[0][0]:
        ssb_seq = make_ssb_qubit_seq()
        set_up_sequence(awg1, alazar, [acq_ctrl], ssb_seq, seq_mode=1)
    else:
        alazar.seq_mode(1)
    if qubit_power is None:
        qubit_power = get_calibration_val('spec_powers')
    qubit.status('on')
    qubit.power(qubit_power)
    qubit_freq = None
    qubit_mag = 0
    for centre in np.linspace(start_freq+100e6, stop_freq-100e6, num=(stop_freq-start_freq)/200e6):
        ssb_centre = centre
        data, pl = measure_ssb(qubit, acq_ctrl, ssb_centre, live_plot=True, key="mag")
        freq, maximum = find_extreme(data, x_key="set", extr="max")
        if maximum > qubit_mag:
            qubit_freq = freq
            qubit_mag = maximum
    if calib_update:
        set_calibration_val('actual_qubit_positions', qubit_freq)
        set_calibration_val('spec_powers', qubit_power)
    print('qubit found at {}, mag {}'.format(qubit_freq, qubit_mag))
    return qubit_freq, qubit_mag

# def find_qubit(freq_array, mag_array, gradient_sign=1, min_res_width=4e6):
    # raise NotImplementedError
    # max_freq = np.amax(freq_array)
    # min_freq = np.amin(freq_array)
    # sampling_rate = len(freq_array) / (max_freq - min_freq)
    # cutoff = 2 / min_res_width
    # smoothed_data = smooth_data_butter(mag_array, sampling_rate, cutoff, 5)
    # if gradient_sign > 1:
    #     qubit_freq_index = np.argmax(smoothed_data)
    # else:
    #     qubit_freq_index = np.argmin(smoothed_data)
    # qubit_freq = freq_data_array[qubit_freq_index]
    # return qubit_freq


def find_cavity_val(acq_controller, cavity, localos, old_position=None,
                    demod_freq=None):
    raise NotImplementedError
    # if old_position is None:
    #     old_position = cavity.frequency()
    # mag_vals = np.zeros(5)
    # freq_vals = old_position + 0.25e6 * (np.arange(5) - 2)
    # demod_freq = demod_freq or acq_controller.demod_freqs()[0]
    # magnitude_array_name = [
    #     n for n in acq_controller.acquisition.names if "mag" in n][0]
    # mag_array_index = acq_controller.acquisition.names.index(
    #     magnitude_array_name)
    # for i, f in enumerate(freq_vals):
    #     set_single_demod_freq(cavity, localos, [acq_controller], demod_freq,
    #                           cav_freq=f)
    #     mag_vals[i] = acq_controller.acquisition()[mag_array_index]
    # # coeffs, stats = np.polynomial.polynomial.polyfit(freq_vals, mag_vals, 1)
    # # new_slope = coeffs[1]
    # gradients_array = np.gradient(mag_vals)
    # new_index = np.argmax(np.absolute(gradients_array))
    # new_position = freq_vals[new_index]
    # set_single_demod_freq(cavity, localos, [acq_controller], demod_freq,
    #                       cav_freq=new_position)
    # cavity.frequency(new_position)
    # # return new_position


def do_tracking_ssb_gate_sweep(qubit, cavity, localos, rec_acq_ctrl,
                               ave_acq_ctrl, gate,
                               initial_qubit_freq, initial_cavity_freq,
                               demod_freq, gate_start, gate_stop, gate_step=0.01,
                               live_plot=True):
    raise NotImplementedError
#     set_single_demod_freq(cavity, localos, [ave_acq_ctrl], demod_freq,
#                           cav_freq=initial_cavity_freq)
#     qubit.frequency(initial_qubit_freq + 100e6)
#     rec_acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
#                                                 base_label='Qubit Drive Frequency',
#                                                 base_unit='Hz',
#                                                 setpoints_start=initial_qubit_freq + 100e6,
#                                                 setpoints_stop=initial_qubit_freq - 100e6)
#     loop = qc.Loop(gate.sweep(gate_start, gate_stop, gate_step)).each(
#         qc.Task(find_cavity_val, ave_acq_ctrl, cavity, localos),
#         rec_acq_ctrl.acquisition,
#         cavity.frequency)
#     if live_plot:
#         dataset = loop.get_data_set()
#         plot = plot_data_live(dataset, rec_acq_ctrl.acquisition)
#         try:
#             _ = loop.with_bg_task(plot.update).run()  # TODO
#         except KeyboardInterrupt:
#             print("Measurement Interrupted")
#         return dataset, plot
#     else:
#         data = loop.run()
#         plots = plot_data(data)
#         return data, plots


def measure_ssb(qubit, acq_ctrl, centre_freq, live_plot=True,
  key=None, save=True):
    qubit.frequency(centre_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=centre_freq + 100e6,
                                            setpoints_stop=centre_freq - 100e6)
    return measure(acq_ctrl.acquisition, key=key, save=save)


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
