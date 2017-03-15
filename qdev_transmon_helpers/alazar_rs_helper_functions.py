import qcodes as qc
from . import get_latest_counter, get_sample_name, get_data_file_format

# TODO: init decision for differentiating between vna functions and alazar
# functions
# TODO: exception types
# getlatest vs counter thing
# TWPA settings
# rabi_setup, ssb setup duplications


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


def remove_demod_freqs(acq_controller):
    freqs = acq_controller.demod_freqs.get()
    for freq in freqs:
        acq_controller.remove_demodulator(freq)


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


def do_cavity_freq_sweep(cavity, localos, cavity_freq, acq_ctrl,
                         cavity_pm=10e6, freq_step=1e6, live_plot=True):
    data_num = get_latest_counter() + 1  # TODO: replcae with dataset counter
    str_data_num = '{0:03d}'.format(data_num)
    title = get_data_file_format().format(
        sample_name=get_sample_name(),
        counter=str_data_num)
    demod_freq = get_demod_freq(cavity, localos, acq_ctrl)
    loop = qc.Loop(cavity.frequency.sweep(cavity_freq - cavity_pm,
                                          cavity_freq + cavity_pm,
                                          freq_step)).each(
        qc.Task(localos.frequency.set,
                (cavity.frequency + demod_freq)),
        acq_ctrl.acquisition)
    dataset = loop.get_data_set()
    dataset.data_num = data_num
    if live_plot:
        plot = qc.QtPlot(figsize=(700, 500))
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
            if i == 0:
                # TODO: test if you can get away with only doing this once
                plot.subplots[0].setTitle(title)
            else:
                plot.subplots[i].setTitle("")
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = [plot]
        plot.data_num = data_num
    else:
        try:
            _ = loop.run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = []
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            pl = qc.QtPlot(getattr(dataset, inst_meas_name))
            pl.subplots[0].showGrid(True, True)
            pl.setTitle(title)
            pl.data_num = data_num
            plots.append(pl)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
        print('plots not saved, choose a plot to save and run plots[i].save()')
    return dataset, plots


def do_demod_freq_sweep(cavity, localos, acq_ctrl, manual_param,
                        demod_centre=15e6, demod_pm=5e6, demod_step=1e6,
                        live_plot=True):
    data_num = get_latest_counter() + 1  # TODO: replcae with dataset counter
    str_data_num = '{0:03d}'.format(data_num)
    title = get_data_file_format().format(
        sample_name=get_sample_name(),
        counter=str_data_num)
    cav_freq = cavity.frequency()
    loop = qc.Loop(localos.frequency.sweep(cav_freq + demod_centre - demod_pm,
                                           cav_freq + demod_centre + demod_pm,
                                           demod_step)).each(
        qc.Task(remove_demod_freqs, acq_ctrl),
        qc.Task(acq_ctrl.demod_freqs.add_demodulator,
                (localos.frequency - cav_freq)),
        qc.Task(manual_param.set, (localos.frequency - cav_freq)),
        acq_ctrl.acquisition)
    dataset = loop.get_data_set()
    dataset.data_num = data_num
    if live_plot:
        plot = qc.QtPlot(figsize=(700, 500))
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
            if i == 0:
                # TODO: test if you can get away with only doing this once
                plot.subplots[0].setTitle(title)
            else:
                plot.subplots[i].setTitle("")
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = [plot]
        plot.data_num = data_num
    else:
        try:
            _ = loop.run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = []
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            pl = qc.QtPlot(getattr(dataset, inst_meas_name))
            pl.subplots[0].showGrid(True, True)
            pl.setTitle(title)
            pl.data_num = data_num
            plots.append(pl)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
        print('plots not saved, choose a plot to save and run plots[i].save()')
    return dataset, plots


def ssb_setup(qubit, cavity=None, localos=None, twpa=None, qubit_pow=-25):
    qubit.power(qubit_pow)
    qubit.status('on')
    if cavity is not None:
        cavity.status('on')
    if localos is not None:
        localos.status('on')
    if twpa is not None:
        twpa.status('on')


def do_ssb_pow_sweep(qubit, acq_ctrl, qubit_freq, qubit_pow_start=-5,
                     qubit_pow_stop=-25, qubit_pow_step=2, live_plot=True):
    data_num = get_latest_counter() + 1  # TODO: replcae with dataset counter
    str_data_num = '{0:03d}'.format(data_num)
    title = get_data_file_format().format(
        sample_name=get_sample_name(),
        counter=str_data_num)
    qubit.frequency(qubit_freq + 10e6)
    acq_ctrl.acquisition.set_base_setpoints(base_name='ssb_drive',
                                            base_label='qubit drive freq',
                                            base_unit='Hz',
                                            setpoints_start=qubit_freq + 10e6,
                                            setpoints_stop=qubit_freq - 10e6)
    loop = qc.Loop(qubit.power.sweep(qubit_pow_start, qubit_pow_stop,
                                     qubit_pow_step)).each(
        acq_ctrl.acquisition)
    dataset = loop.get_data_set()
    dataset.data_num = data_num
    if live_plot:
        plot = qc.QtPlot(figsize=(700, 500))
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
            if i == 0:
                # TODO: test if you can get away with only doing this once
                plot.subplots[0].setTitle(title)
            else:
                plot.subplots[i].setTitle("")
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = [plot]
        plot.data_num = data_num
    else:
        try:
            _ = loop.run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        plots = []
        for i, name in enumerate(acq_ctrl.names):
            inst_meas_name = "{}_{}".format(acq_ctrl._instrument.name, name)
            pl = qc.QtPlot(getattr(dataset, inst_meas_name))
            pl.subplots[0].showGrid(True, True)
            pl.setTitle(title)
            pl.data_num = data_num
            plots.append(pl)
            plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
            plot.subplots[i].showGrid(True, True)
        print('plots not saved, choose a plot to save and run plots[i].save()')
    return dataset, plots


def rabi_setup(qubit, cavity=None, localos=None, twpa=None, qubit_pow=-5):
    qubit.power(qubit_pow)
    qubit.status('on')
    if cavity is not None:
        cavity.status('on')
    if localos is not None:
        localos.status('on')
    if twpa is not None:
        twpa.status('on')

