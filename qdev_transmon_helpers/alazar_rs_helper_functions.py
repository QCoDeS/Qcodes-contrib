import qcodes as qc
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from . import get_latest_counter, save_fig, get_title, plot_data_live, \
    plot_data, sweep1d, exp_decay_sin, exp_decay

# TODO: init decision for differentiating between vna functions and alazar
# functions
# TODO: exception types
# TWPA settings
# rabi_setup, ssb setup duplications
# TODO: docstrings


def config_alazar(alazar, seq_mode=0):
    if seq_mode not in [0, 1]:
        raise ValueError('must set seq mode to 0 or 1')
    if seq_mode:
        io_mode = 'AUX_IN_TRIGGER_ENABLE'
        io_param = 'TRIG_SLOPE_POSITIVE'
    else:
        io_mode = 'AUX_IN_AUXILIARY'
        io_param = 'NONE'
    alazar.config(clock_source='EXTERNAL_CLOCK_10MHz_REF',
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
        alazar.config(aux_io_mode='AUX_IN_TRIGGER_ENABLE',
                      aux_io_param='TRIG_SLOPE_POSITIVE')
    elif mode == 0:
        alazar.config(aux_io_mode='AUX_IN_AUXILIARY',
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
        acq_controller.remove_demodulator(freq)


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
        plot = plot_data_live(dataset, acq_ctrl)
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
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
        acq_ctrl.acquisition)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, acq_ctrl)
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        return data, plots


def do_ssb_pow_sweep(qubit, acq_ctrl, qubit_freq, qubit_pow_start=-5,
                     qubit_pow_stop=-25, qubit_pow_step=2, live_plot=True):
    qubit.frequency(qubit_freq + 10e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_drive',
                                            base_label='qubit drive freq',
                                            base_unit='Hz',
                                            setpoints_start=qubit_freq + 10e6,
                                            setpoints_stop=qubit_freq - 10e6)
    return sweep1d(acq_ctrl.acquisition, qubit.power, qubit_pow_start,
                   qubit_pow_stop, qubit_pow_step, live_plot=live_plot)


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
           data_num=None, plot=True, subplot=None):
    xdata = getattr(data, x_name)
    ydata = getattr(data, y_name)
    popt, pcov = curve_fit(exp_decay, xdata, ydata)
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
            x_units = xdata.units
            y_units = ydata.units
        except AttributeError:
            x_units = None
            y_units = None
        if data_num is None:
            try:
                data_num = data.data_num
                title = get_title(data_num) + '_T1'
                if not hasattr(fig, 'data_num'):
                    fig.data_num = data_num
            except AttributeError:
                title = 'T1'
        else:
            title = get_title(data_num) + '_T1'
        fig.data_num = data_num
        ax.plot(exp_decay(xdata, *popt), label='fit: T1 {}{}'.format(popt[1],
                                                                     x_units))
        ax.plot(ydata, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        save_fig(ax, name='t1_fit')
        return ax, popt, errors
    else:
        return popt, errors


def get_t2(data, x_name='delay', y_name='magnitude',
           data_num=None, plot=True, subplot=None):
    xdata = getattr(data, x_name)
    ydata = getattr(data, y_name)
    popt, pcov = curve_fit(exp_decay_sin, xdata, ydata)
    errors = np.sqrt(np.diag(pcov))
    print('fit to equation of form y = a * exp(-x / b) * sin(c * x + d) + e'
          'gives:\na {}, b {}, c {}, d {}, e{}\n'
          'with one standard deviation errors:\n'
          'a {}, b {}, c {}, d {}, e{}'.format(popt[0], popt[1], popt[2],
                                               popt[4], popt[5], errors[0],
                                               errors[1], errors[2],
                                               errors[4], errors[4]))
    if plot:
        if subplot is None:
            fig, ax = plt.subplots()
        else:
            ax = subplot
            fig = ax.figure
        try:
            x_units = xdata.units
            y_units = ydata.units
        except AttributeError:
            x_units = None
            y_units = None
        if data_num is None:
            try:
                data_num = data.data_num
                title = get_title(data_num) + '_T1'
                if not hasattr(fig, 'data_num'):
                    fig.data_num = data_num
            except AttributeError:
                title = 'T2*'
        else:
            title = get_title(data_num) + '_T2*'
        fig.data_num = data_num
        ax.plot(exp_decay(xdata, *popt), label='fit: T1 {}{}'.format(popt[1],
                                                                     x_units))
        ax.plot(ydata, label='data')
        ax.set_xlabel('{} ({})'.format(x_name, x_units))
        ax.set_ylabel('{} ({})'.format(y_name, y_units))
        ax.set_title(title)
        ax.legend(loc='upper right', fontsize=10)
        save_fig(ax, name='t2_fit')
        return ax, popt, errors
    else:
        return popt, errors
