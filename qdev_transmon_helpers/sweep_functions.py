import qcodes as qc
from . import plot_data_single_window, plot_data

# TODO: _?
# TODO: check working of counter/data_num


def sweep2d(meas_param, sweep_param1, start1, stop1, step1,
            sweep_param2, start2, stop2, step2, delay=0.01,
            live_plot=True, key=None, save=True):
    """
    Function which does a 2 dimensional sweep and optionally plots the results.

    Args:
        meas_param: parameter which we want the value of at each point
        sweep_param1: parameter to be swept in outer loop (default on y axis)
        start1: starting value for sweep_param1
        stop1: final value for sweep_param1
        step1: value to step sweep_param1 by
        sweep_param2: parameter to be swept in inner loop (default on x axis)
        start2: starting value for sweep_param2
        stop2: final value for sweep_param2
        step2: value to step sweep_param2 by
        delay (default 0.01): mimimum time to spend on each point
        live_plot (default True): live plot or not, currently live plot means
            plotting all meas_param arrays (eg magnitude and phase) in one
            window
        key (default None): string specifying specific parameter array to be
            plotted, default is to plot all
        save (default True): whether to save png on completion, nb if you
            choose not to live plot and the measured parameter returns
            multiple values then unless you specify a specific one to plot
            via 'key' then none will be saved.

    Returns:
        data (qcodes dataset)
        plot(s): depending on live_plot value returns single window with
            meas_param values or one per meas_param array
    """
    loop = qc.Loop(sweep_param2.sweep(
        start2, stop2, step2), delay).loop(sweep_param1.sweep(
            start1, stop1, step1), delay).each(meas_param)
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


def sweep1d(meas_param, sweep_param, start, stop, step,
            delay=0.01, live_plot=True, key=None, save=True):
    """
    Function which does a 1 dimensional sweep and optionally plots the results.

    Args:
        meas_param: parameter which we want the value of at each point
        sweep_param: parameter to be swept in outer loop (default on y axis)
        start: starting value for sweep_param1
        stop: final value for sweep_param1
        step: value to step sweep_param1 by
        delay (default 0.01): mimimum time to spend on each point
        live_plot (default True): live plot or not, currently live plot means
            plotting all meas_param attributes (eg magnitude and phase) in one
            window
        key (default None): string specifying specific parameter array to be
            plotted, default is to plot all
        save (default True): whether to save png on completion, nb if you
            choose not to live plot and the measured parameter returns
            multiple values then unless you specify a specific one to plot
            via 'key' then none will be saved.

    Returns:
        data (qcodes dataset)
        plot(s): depending on live_plot value returns single window with
            meas_param values or one per meas_param attribute
    """
    loop = qc.Loop(sweep_param.sweep(
        start, stop, step), delay).each(meas_param)
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


def measure(meas_param, plot=True, key=None, save=True):
    """
    Function which does a 1 dimensional sweep and optionally plots the results.

    Args:
        meas_param: parameter which we want the value of at each point
        plot (default True): live plot not possible but will plot the
            meas_param attributes in seperate windows and return them
            if this is True
        key (default None): string specifying specific parameter array to be
            plotted, default is to plot all
        save (default True): whether to save png on completion

    Returns:
        data (qcodes dataset)
        plots: optionally depending on plot value
    """
    data = qc.Measure(meas_param).run()
    data.data_num = data.location_provider.counter
    if plot:
        plots = plot_data(data, key=key)
        if (key is not None) and save:
            plots.save()
        else:
            print('warning: plots not saved. To save one choose '
                  'and run plots[i].save()')
        return data, plots
    else:
        return data
