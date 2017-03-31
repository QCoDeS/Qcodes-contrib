import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from . import plot_data_live, plot_data, sweep1d, exp_decay, exp_decay_sin, get_title, measure, save_fig

# TODO: init decision for differentiating between vna functions and alazar
# functions
# TODO: exception types
# TWPA settings
# rabi_setup, ssb setup duplications
# TODO: docstrings


def config_alazar(alazar, seq_mode=0, clock_source='EXTERNAL_CLOCK_10MHz_REF'):
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
    if mode == 1:
        alazar.config(sample_rate = alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_TRIGGER_ENABLE',
                      aux_io_param='TRIG_SLOPE_POSITIVE')
    elif mode == 0:
        alazar.config(sample_rate = alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_AUXILIARY',
                      aux_io_param='NONE')
    else:
        raise ValueError('must set seq mode to 0 or 1')


def get_demod_freq(cavity, localos, acq_ctrl):
    lo, cav = localos.frequency(), cavity.frequency()
    demod = lo - cav
    acq_freqs = acq_ctrl.demod_freqs()
    if demod not in acq_freqs:
        raise Exception('demod freq {} (from cavity freq {} and localos '
                        'freq {}) not in acq controller demodulation '
                        'frequencies: {}.'.format(demod, cav, lo, acq_freqs))
    else:
        return demod


def set_single_demod_freq(cavity, localos, acq_controllers, demod_freq,
                          cav_freq=None):
    if cav_freq is not None:
        cavity.frequency(cav_freq)
    else:
        cav_freq = cavity.frequency()
    localos.frequency(cav_freq + demod_freq)
    for ctrl in acq_controllers:
        remove_demod_freqs(ctrl)
        ctrl.demod_freqs.add_demodulator(demod_freq)


def remove_demod_freqs(acq_controller):
    freqs = acq_controller.demod_freqs.get()
    for freq in freqs:
        acq_controller.demod_freqs.remove_demodulator(freq)


def cavity_sweep_setup(cavity, localos, qubit=None, twpa=None,
                       cav_pow=-35, locos_pow=15):
    cavity.power(cav_pow)
    localos.power(locos_pow)
    cavity.status('on')
    localos.status('on')
    if qubit is not None:
        qubit.status('off')
    if twpa is not None:
        twpa.status('on')


def do_cavity_freq_sweep(cavity, localos, cavity_freq, acq_ctrl,
                         cavity_pm=10e6, freq_step=1e6, live_plot=True):
    demod_freq = get_demod_freq(cavity, localos, acq_ctrl)
    loop = qc.Loop(cavity.frequency.sweep(cavity_freq - cavity_pm,
                                          cavity_freq + cavity_pm,
                                          freq_step)).each(
        qc.Task(localos.frequency.set,
                (cavity.frequency + demod_freq)),
        acq_ctrl.acquisition)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, acq_ctrl.acquisition)
        try:
            _ = loop.with_bg_task(plot.update).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        return data, plots


def ssb_setup(qubit, cavity=None, localos=None, twpa=None, qubit_pow=-25):
    qubit.power(qubit_pow)
    qubit.status('on')
    if cavity is not None:
        cavity.status('on')
    if localos is not None:
        localos.status('on')
    if twpa is not None:
        twpa.status('on')


def rabi_setup(qubit, cavity=None, localos=None, twpa=None,
               qubit_pow=-5, qubit_freq=None):
    qubit.power(qubit_pow)
    qubit.status('on')
    if cavity is not None:
        cavity.status('on')
    if localos is not None:
        localos.status('on')
    if twpa is not None:
        twpa.status('on')
    if qubit_freq is not None:
        qubit.frequency(qubit_freq)


def do_demod_freq_sweep(cavity, localos, acq_ctrl, manual_param,
                        demod_centre=15e6, demod_pm=5e6, demod_step=1e6,
                        live_plot=True):
    cav_freq = cavity.frequency()
    loop = qc.Loop(localos.frequency.sweep(cav_freq + demod_centre - demod_pm,
                                           cav_freq + demod_centre + demod_pm,
                                           demod_step)).each(
        qc.Task(remove_demod_freqs, acq_ctrl),
        qc.Task(acq_ctrl.demod_freqs.add_demodulator,
                (localos.frequency - cav_freq)),
        qc.Task(manual_param.set, (localos.frequency - cav_freq)),
        acq_ctrl.acquisition,
        manual_param)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, acq_ctrl.acquisition)
        try:
            _ = loop.with_bg_task(plot.update).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        return data, plots

def do_ssb_pow_sweep(qubit, acq_ctrl, qubit_freq, qubit_pow_start=-5,
                     qubit_pow_stop=-25, qubit_pow_step=2, live_plot=True):
    qubit.frequency(qubit_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=qubit_freq + 100e6,
                                            setpoints_stop=qubit_freq - 100e6)
    return sweep1d(acq_ctrl.acquisition, qubit.power, qubit_pow_start,
                   qubit_pow_stop, qubit_pow_step, live_plot=live_plot)

def do_ssb_time_sweep(qubit, acq_ctrl, manual_param, qubit_freq, time=500,
                     live_plot=True):
    qubit.frequency(qubit_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=qubit_freq + 100e6,
                                            setpoints_stop=qubit_freq - 100e6)
    return sweep1d(acq_ctrl.acquisition, manual_param, 0,
                   time, 1, live_plot=live_plot)

def do_ssb_gate_sweep(qubit, acq_ctrl, gate, qubit_freq, gate_start, gate_stop,
                      gate_step, live_plot=True):
    qubit.frequency(qubit_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=qubit_freq + 100e6,
                                            setpoints_stop=qubit_freq - 100e6)
    return sweep1d(acq_ctrl.acquisition, gate, gate_start,
                   gate_stop, gate_step, live_plot=live_plot)

def qubit_from_ssb_measure(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError

def qubit_from_ssb_power_sweep(dataset, gradient_sign=1, min_res_width=4e6):
    raise NotImplementedError

def qubit_from_ssb_volt_sweep(dataset, gradient_sign=1, min_res_width=4e6, high_voltage=True):
    voltage_array_name = [n for n in dataset.arrays.keys() if "voltage" in n][0]
    magnitude_array_name = [n for n in dataset.arrays.keys() if "mag" in n][0]
    frequency_array_name = [n for n in dataset.arrays.keys() if "ssb" in n][0]

    if high_voltage:
            voltage_index = np.argmax(getattr(dataset, voltage_array_name))
    else:
            voltage_index = np.argmin(getattr(dataset, voltage_array_name))

    mag_array = getattr(dataset, magnitude_array_name)[voltage_index]
    freq_array = getattr(dataset, frequency_array_name)[0]
    return find_qubit(freq_array, mag_array,
        gradient_sign=gradient_sign,
        min_res_width=min_res_width)


def find_qubit(freq_array, mag_array, gradient_sign=1, min_res_width=4e6):
    max_freq = np.amax(freq_array)
    min_freq = np.amin(freq_array)
    sampling_rate = len(freq_array) / (max_freq - min_freq) 
    cutoff = 2 / min_res_width
    smoothed_data = smooth_data_butter(mag_array, sampling_rate, cutoff, 5)
    if gradient_sign  > 1:
        qubit_freq_index = np.argmax(smoothed_data)
    else:
        qubit_freq_index = np.argmin(smoothed_data)
    qubit_freq = freq_data_array[qubit_freq_index]
    return qubit_freq


def find_cavity_val(acq_controller, cavity, localos, old_position=None, demod_freq=None):
    if old_position is None:
        old_position = cavity.frequency()
    mag_vals = np.zeros(5)
    freq_vals = old_position + 0.25e6 * (np.arange(5) - 2)
    demod_freq = demod_freq or acq_controller.demod_freqs()[0]
    magnitude_array_name = [n for n in acq_controller.acquisition.names if "mag" in n][0]
    mag_array_index = acq_controller.acquisition.names.index(magnitude_array_name)
    for i, f in enumerate(freq_vals):
        set_single_demod_freq(cavity, localos, [acq_controller], demod_freq,
                              cav_freq=f)
        mag_vals[i] = acq_controller.acquisition()[mag_array_index]
    # coeffs, stats = np.polynomial.polynomial.polyfit(freq_vals, mag_vals, 1)
    # new_slope = coeffs[1]
    gradients_array = np.gradient(mag_vals)
    new_index = np.argmax(np.absolute(gradients_array))
    new_position = freq_vals[new_index]
    set_single_demod_freq(cavity, localos, [acq_controller], demod_freq,
                              cav_freq=new_position)
    cavity.frequency(new_position)
    # return new_position


def do_tracking_ssb_gate_sweep(qubit, cavity, localos, rec_acq_ctrl, ave_acq_ctrl, gate,
                                    initial_qubit_freq, initial_cavity_freq, demod_freq,
                                    gate_start, gate_stop, gate_step=0.01,
                                    live_plot=True):
    set_single_demod_freq(cavity, localos, [ave_acq_ctrl], demod_freq,
                              cav_freq=initial_cavity_freq)
    qubit.frequency(initial_qubit_freq + 100e6)
    rec_acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=initial_qubit_freq + 100e6,
                                            setpoints_stop=initial_qubit_freq - 100e6)
    loop = qc.Loop(gate.sweep(gate_start, gate_stop, gate_step)).each(
        qc.Task(find_cavity_val, ave_acq_ctrl, cavity, localos),
        rec_acq_ctrl.acquisition,
        cavity.frequency)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, rec_acq_ctrl.acquisition)
        try:
            _ = loop.with_bg_task(plot.update).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        return data, plots



def measure_ssb(qubit, acq_ctrl, centre_freq, live_plot=True):
    qubit.frequency(centre_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_qubit_drive_freq',
                                            base_label='Qubit Drive Frequency',
                                            base_unit='Hz',
                                            setpoints_start=centre_freq + 100e6,
                                            setpoints_stop=centre_freq - 100e6)
    return measure(acq_ctrl.acquisition)


def do_rabi_freq_sweep(qubit, acq_ctrl, qubit_freq, qubit_freq_pm=4e6,
                       qubit_freq_step=0.5e6, live_plot=True):
    qubit_freq_start = qubit_freq - qubit_freq_pm
    qubit_freq_stop = qubit_freq + qubit_freq_pm
    return sweep1d(acq_ctrl.acquisition, qubit.frequency, qubit_freq_start,
                   qubit_freq_stop, qubit_freq_step, live_plot=live_plot)


def do_rabi_pow_sweep(qubit, acq_ctrl, qubit_pow_start=-3,
                      qubit_pow_stop=-10, qubit_pow_step=1, live_plot=True):
    return sweep1d(acq_ctrl.acquisition, qubit.power, qubit_pow_start,
                   qubit_pow_stop, qubit_pow_step, live_plot=live_plot)


def do_ramsey_freq_sweep(qubit, acq_ctrl, qubit_freq, qubit_freq_pm=4e6,
                         qubit_freq_step=0.5e6, live_plot=True):
    qubit_freq_start = qubit_freq - qubit_freq_pm
    qubit_freq_stop = qubit_freq + qubit_freq_pm
    return sweep1d(acq_ctrl.acquisition, qubit.frequency, qubit_freq_start,
                   qubit_freq_stop, qubit_freq_step, live_plot=live_plot)


def get_t1(data, x_name='delay', y_name='magnitude',
           counter=None, plot=True, subplot=None):
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
        if counter is None:
            try:
                counter = data.location_provider.counter
                title = get_title(counter) + '_T1'
                if not hasattr(fig, 'counter'):
                    fig.counter = counter
            except AttributeError:
                title = 'T1'
        else:
            title = get_title(counter) + '_T1'
        fig.counter = counter
        ax.plot(x_data, exp_decay(x_data, *popt), label='fit: T1 {}{}'.format(popt[1],
                                                                     x_units))
        ax.plot(x_data, y_data, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        try:
          qubit = get_calibration_dict()['current_qubit']
          name = 'qubit{}_t1'.format(qubit)
        except Exception:
          name='t1'
        save_fig(ax, name=name)
        return ax, popt, errors
    else:
        return popt, errors


def get_t2(data, x_name='delay', y_name='magnitude',
           counter=None, plot=True, subplot=None):
    x_data = getattr(getattr(data, x_name), 'ndarray')
    y_data = getattr(getattr(data, y_name), 'ndarray')
    x_units = getattr(getattr(data, x_name), 'unit')
    y_units = getattr(getattr(data, y_name), 'unit')
    popt, pcov = curve_fit(exp_decay_sin, x_data, y_data, p0=[0.003, 1e-7, 10e7, 0, 0.01])
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
        if counter is None:
            try:
                counter = data.location_provider.counter
                title = get_title(counter) + '_T1'
                if not hasattr(fig, 'counter'):
                    fig.counter = counter
            except AttributeError:
                title = 'T2*'
        else:
            title = get_title(counter) + '_T2*'
        fig.counter = counter
        ax.plot(x_data, exp_decay_sin(x_data, *popt), label='fit: T2 {}{}'.format(popt[1],
                                                                     x_units))
        ax.plot(x_data, y_data, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        try:
          qubit = get_calibration_dict()['current_qubit']
          name = 'qubit{}_t2'.format(qubit)
        except Exception:
          name='t1'
        save_fig(ax, name=name)
        return ax, popt, errors
    else:
        return popt, errors