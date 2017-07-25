import qcodes as qc
from . import plot_data_single_window, plot_data, sweep1d, measure, \
    get_calibration_val

# TODO: exception types
# TODO: TWPA settings
# TODO: _ ?
# TODO: write set_demod_freqs


def config_alazar(alazar, seq_mode='off',
                  clock_source='EXTERNAL_CLOCK_10MHz_REF'):
    """
    Function which puts alazar in sequency mode, configures and sets clock
    source

    Args:
        alazar instrument
        seq mode ('on' or 'off') (default 'off')
        clock_source (default 'EXTERNAL_CLOCK_10MHz_REF')
    """
    if seq_mode not in ['on', 'off']:
        raise ValueError('must set seq mode to "on" or "off"')
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
        return 'on'
    elif (alazar.aux_io_mode() is 'AUX_IN_AUXILIARY' and
          alazar.aux_io_param() is 'NONE'):
        return 'off'
    else:
        raise ValueError('aux_io_mode: {}, aux_io_param: {} '
                         'do not correspond to seq_mode on or off')


def set_alazar_seq_mode(alazar, mode):
    """
    Sets the sequence mode of the alazar

    Args:
        alazar instrument
        mode ('off' or 'on')
    """
    if mode is 'on':
        alazar.config(sample_rate=alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_TRIGGER_ENABLE',
                      aux_io_param='TRIG_SLOPE_POSITIVE')
    elif mode is 'off':
        alazar.config(sample_rate=alazar.sample_rate(),
                      clock_edge=alazar.clock_edge(),
                      clock_source=alazar.clock_source(),
                      aux_io_mode='AUX_IN_AUXILIARY',
                      aux_io_param='NONE')
    else:
        raise ValueError('must set seq mode to "on" or "off"')


def get_demod_freq(cavity, localos, acq_ctrl, SSBfreq=0):
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
    demod = lo - (cav - SSBfreq)
    acq_freqs = acq_ctrl.demod_freqs()
    if demod not in acq_freqs:
        raise Exception('demod freq {} (from cavity freq {} and localos '
                        'freq {}) not in acq controller demodulation '
                        'frequencies: {}.'.format(demod, cav, lo, acq_freqs))
    else:
        return demod


def set_single_demod_freq(cavity, localos, acq_ctrls, demod_freq=None,
                          cavity_freq=None, localos_freq=None, SSBfreq=0):
    if all([i is not None for i in [demod_freq, cavity_freq, localos_freq]]):
        raise Exception('Set up demodulation by setting max 2 out of '
                        '[demod_freq, cavity_freq, localos_freq], not all '
                        'three')
    cavity.status('on')
    localos.status('on')
    if localos_freq is None:
        if cavity_freq is None:
            cavity_freq = cavity.frequency()
        else:
            cavity.frequency(cavity_freq)
        if demod_freq is None:
            demod_freq = localos.frequency() - (cavity_freq - SSBfreq)
        else:
            localos.frequency(cavity_freq - SSBfreq + demod_freq)
    elif cavity_freq is None:
        localos.frequency(localos_freq)
        if demod_freq is None:
            demod_freq = localos_freq - (cavity.frequency() - SSBfreq)
        else:
            cavity.frequency(localos_freq - demod_freq + SSBfreq)
    else:
        localos.frequency(localos_freq)
        cavity.frequency(cavity_freq)
        demod_freq = localos_freq - (cavity_freq - SSBfreq)
    for ctrl in acq_ctrls:
        remove_demod_freqs(ctrl)
        ctrl.demod_freqs.add_demodulator(demod_freq)


def set_demod_freqs(cavity_list, localos, acq_ctrls,
                    localos_freq=None, cavity_freqs=None,
                    demod_freqs=None):
    if all([i is not None for
            i in [demod_freqs, cavity_freqs, localos_freq]]):
        raise Exception('Set up demodulation by settin max 2 out of '
                        '[demod_freqs, cavity_freqs, localos_freq], not all '
                        'three')
    if any(i is not None for i in [demod_freqs, cavity_freqs]):
        length = len(cavity_list)
        if not all(len(x) == length for x in [cavity_freqs, demod_freqs]):
            raise Exception(
                'cavity list length not equal to cavity freqs list length '
                'and demod_freqs lengths: {}'.format(
                    [length, len(cavity_freqs), len(demod_freqs)]))

    localos.status('on')
    if all([i is not None for i in [cavity_freqs, demod_freqs]]):
        localos_freq = cavity_freqs[0] + demod_freqs[0]
        for i, cav_f in enumerate(cavity_freqs):
            if cav_f[i] + demod_freqs[i] != localos_freq:
                raise Exception(
                    'Cavity frequencies and demod frequencies set'
                    ' but these do not all correspond to the same value for '
                    'the local oscillator. This is unphsyical')
            cavity_list[i].frequency(cav_f)
        localos.frequency(localos_freq)
    elif localos_freq is None:
        localos_freq = localos.frequency()
        if cavity_freqs is not None:
            for i, cavity in enumerate(cavity_list):
                cavity.frequency(cavity_freqs[i])
                demod_freqs[i] = localos_freq - cavity_freqs[i]
        elif demod_freqs is not None:
            for i, cavity in enumerate(cavity_list):
                cavity.frequency(localos_freq + demod_freqs[i])
        else:
            demod_freqs = [localos_freq - c.frequency() for c in cavity_list]

    for ctrl in acq_ctrls:
        remove_demod_freqs(ctrl)

    for i, cav in enumerate(cavity_list):
        cav.status('on')
        for ctrl in acq_ctrls:
            ctrl.demod_freqs.add_demodulator(demod_freqs[i])


# def set_demod_freqs_SSB(cavity, localos, acq_ctrls,
#                         localos_freq=None, cavity_freq=None,
#                         demod_freqs=None, SSBfreqs=None):
#     if all([i is not None for
#             i in [demod_freqs, cavity_freq, localos_freq, SSBfreqs]]):
#         raise Exception('Set up demodulation by settin max 2 out of '
#                         '[demod_freqs, SSBfreqs, cavity_freq, localos_freq], '
#                         'not all three')
#     elif not any([i is not None for i in demod_freqs, SSBfreqs]):
#         raise Exception('Must set one from demod_freqs and SSBfreqs')

#     localos.status('on')
#     cavity.status('on')

#     if localos_freq is None:
#         if cavity_freq is None:
#             cavity_freq = cavity.frequency()
#         if demod_freqs is None:
#             localos_freq = localos.frequency()
#             demod_freqs = [(localos_freq - cavity_freq + SSBfreqs[i])
#                            for i in range(freq_count)]
#         elif SSBfreqs is None:
#             SSBfreqs = [0] * freq_count
#         else:
#             localos.frequency()

#     if localos_freq is not None:
#         localos.frequency(localos_freq)
#     else:
#         localos_freq = localos.frequency()

#     if cavity_freqs is not None:
#         if len(cavity_freqs) != len(cavity_list):
#             raise Exception(
#                 'cavity freq list length not equal to cavity list length')
#         demod_freqs = [localos_freq - c for c in cavity_freqs]
#     else:
#         demod_freqs = [localos_freq - c.frequency() for c in cavity_list]
#         cavity_freqs = [None] * len(cavity_list)
#     for ctrl in acq_ctrls:
#         remove_demod_freqs(ctrl)

#     for i, cav in enumerate(cavity_list):
#         cav.status('on')
#         if cavity_freqs[i] is not None:
#             cav.frequency(cavity_freqs[i])
#         for ctrl in acq_ctrls:
#             ctrl.demod_freqs.add_demodulator(demod_freqs[i])


def remove_demod_freqs(acq_ctrl):
    """
    Function whish removes demod freqs from acquisition controller

    Args:
        acq_ctrl (alazar acq controller)
    """
    freqs = acq_ctrl.demod_freqs.get()
    for freq in freqs:
        acq_ctrl.demod_freqs.remove_demodulator(freq)


def do_cavity_freq_sweep(cavity, localos, cavity_freq, acq_ctrl,
                         cavity_pm=10e6, freq_step=1e6, demod_freq=None,
                         live_plot=True, key=None, save=True):
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

    Returns:
        data, plot(s)
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


def set_cavity_from_calib_dict(cavity, localos, acq_ctrls, num_avg=1000):
    """
    Funtion which sets the cavity, local oscillator and the acq controllers to
    have correct demodulation settings for single qubit readout as well as
    averaging settings and int_time, int_delay and cavity power, localos_power

    Args:
        cavity (R&S instrument)
        localos (R&S instrument)
        acq_ctrls (list of alazar acq controller instruments)
        num_avg (int): num of averages for acq controller
    """
    for acq_ctrl in acq_ctrls:
        acq_ctrl.int_time(get_calibration_val('int_time'))
        acq_ctrl.int_delay(get_calibration_val('int_delay'))
        acq_ctrl.num_avg(num_avg)
    cavity.power(get_calibration_val('cavity_pow'))
    localos.power(get_calibration_val('localos_pow'))
    set_single_demod_freq(cavity, localos, acq_ctrls,
                          get_calibration_val('demod_freq'),
                          cavity_freq=get_calibration_val('cavity_freq'))


def sweep2d_ssb(qubit, acq_ctrl, centre_freq, sweep_param,
                start, stop, step, delay=0.01, live_plot=True,
                key=None, save=True):
    """
    Function which sets up a ssb spectroscopy 'hardware controlled sweep'
    +-100MHz around cenre_freq on one axis and sweeps another parameter on
    the other axis. Assumes correct awg upload. Produces a 2d plot.

    Args:
        qubit (R&S instrument)
        acq_ctrls(alazar acq controller instrument)
        centre_freq (float): freq to centre ssb spectoscopy around
        sweep_param (qcodes parameter): param top sweep on y axis
        start: start value for sweep_param
        stop: stop value for sweep param
        step: step value for sweep param
        delay (default 0.01): delay value beween step of sweep param
        live_plot (bool) (default True)
        key (string) (default None): string key used to search for data
            array to plot and save (otherwise plots all)
        save (bool) (default True): whether to save png

    Returns:
        sweep1d result
    """
    qubit.frequency(centre_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(
        base_name='ssb_qubit_drive_freq',
        base_label='Qubit Drive Frequency',
        base_unit='Hz',
        setpoints_start=centre_freq + 100e6,
        setpoints_stop=centre_freq - 100e6)
    return sweep1d(acq_ctrl.acquisition, sweep_param, start,
                   stop, step, delay=delay, live_plot=live_plot,
                   key=key, save=save)


def measure_ssb(qubit, acq_ctrl, centre_freq,
                key=None, save=True):
    """
    Function which does a 'hardware controlled sweep' for single sideband
    spectroscopy uploaded to the awg +-100MHz around the centre_freq.
    Produces a 3d plot.

    Args:
        qubit (R&S instrument)
        acq_ctrls(alazar acq controller instrument)
        centre_freq (float): freq to centre ssb spectoscopy around
        live_plot (bool) (default True)
        key (string) (default None): string key used to search for data
            array to plot and save (otherwise plots all)
        save (bool) (default True): whether to save png

    Returns:
        measure result
    """
    qubit.frequency(centre_freq + 100e6)
    acq_ctrl.acquisition.set_base_setpoints(
        base_name='ssb_qubit_drive_freq',
        base_label='Qubit Drive Frequency',
        base_unit='Hz',
        setpoints_start=centre_freq + 100e6,
        setpoints_stop=centre_freq - 100e6)
    return measure(acq_ctrl.acquisition, key=key, save=save)
