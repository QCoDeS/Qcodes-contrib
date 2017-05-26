import qcodes as qc
from . import plot_data_single_window, plot_data


def sweep_awg_amp(meas_param, awg_ch_1, awg_ch_2, chan_amp, start, stop, step,
                  delay=0.01, live_plot=True, key=None, save=True):
    loop = qc.Loop(chan_amp.sweep(start, stop, step)).each(
        qc.Task(awg_ch_1.set, chan_amp),
        qc.Task(awg_ch_2.set, chan_amp),
        meas_param)
    if live_plot:
        dataset = loop.get_data_set()
        dataset.data_num = dataset.location_provider.counter
        plot = plot_data_single_window(dataset, meas_param, key=key)
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
