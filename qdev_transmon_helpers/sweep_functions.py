import qcodes as qc
from . import plot_data_live, plot_data

# TODO: if counter thing works change data_num for couter everywhere
# TODO: delete first measure if do1d hard works
# TODO: docstrings


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
    else: # TODO: switch this back when plot_data fixed
        # data = loop.run()
        # plots = plot_data(data)
        # print('warning: plots not saved, choose a plot to save '
        #       'and run plots[i].save()')
        # return data, plots
        return data


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
    else: # TODO: switch this back when plot_data fixed
        # data = loop.run()
        # plots = plot_data(data)
        # print('warning: plots not saved, choose a plot to save '
        #       'and run plots[i].save()')
        # return data, plots
        return data


def measure(meas_param, plot=True):
    data = qc.Measure(meas_param).run()
    if plot:
        plots = plot_data_live(data, meas_param) # TODO: switch this back when plot_data fixed
        print('warning: plots not saved, choose a plot to save '
              'and run plots[i].save()')
        return data, plots
    else:
        return data
