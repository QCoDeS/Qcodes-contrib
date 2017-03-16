import qcodes as qc
from . import get_latest_counter, get_title, plot_data_live, plot_data

# TODO: if counter thing works change data_num for couter everywhere
# TODO: delete first measure if do1d hard works
# TODO: docstrings


def measure(param, plot=True, plot_variable=None):
    """
    Function which does qc.Measure with live plot and plot.save in fewer lines

    Args:
        Parameter to measure
        plot (default=True): do you want a plot?
        plot_variable (default=None): what do you want to plot?
            None -> deafultParameter_array

    Returns:
        data (qcodes DataSet)
        plot (QtPlot)
    """
    data = qc.Measure(param).run()
    data_num = get_latest_counter()
    data.data_num = data_num
    if plot:
        title = get_title(data_num)
        plots = []
        # TODO does this actually check for being a setpoint?
        for value in data.arrays.keys():
            if "set" not in value:
                pl = qc.QtPlot(getattr(data, value))
                pl.subplots[0].setTitle(title)
                pl.subplots[0].showGrid(True, True)
                pl.data_num = data_num
                plots.append(pl)
        print('warning: none of these plots are saved, '
              'choose one and run plot.save()')
        return data, plots
    else:
        return data


def sweep2d(meas_param, sweep_param1, start1, stop1, step1,
            sweep_param2, start2, stop2, step2, delay=0.01, live_plot=True):
    loop = qc.Loop(sweep_param2.sweep(
        start2, stop2, step2), delay).loop(sweep_param1.sweep(
            start1, stop1, step1), delay).each(meas_param)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, meas_param)
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        print('warning: plots not saved, choose a plot to save '
              'and run plots[i].save()')
        return data, plots


def sweep1d(meas_param, sweep_param, start, stop, step,
            delay=0.01, live_plot=True):
    loop = qc.Loop(sweep_param.sweep(
        start, stop, step), delay).each(meas_param)
    if live_plot:
        dataset = loop.get_data_set()
        plot = plot_data_live(dataset, meas_param)
        try:
            _ = loop.with_bg_task(plot.update, plot.save).run()  # TODO
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        plots = plot_data(data)
        print('warning: plots not saved, choose a plot to save '
              'and run plots[i].save()')
        return data, plots


def measure2(meas_param):
    data = qc.Measure(meas_param).run()
    plots = plot_data(data)
    print('warning: plots not saved, choose a plot to save '
          'and run plots[i].save()')
    return data, plots
