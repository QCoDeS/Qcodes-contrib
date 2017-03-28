import qcodes as qc
import matplotlib.pyplot as plt
import numpy as np

# TODO: check x and y behaviour of plotting if not specified
# TODO: remove commented out code!
# TODO: docstrings

def get_title(counter, current_qubit=None):
    str_counter = '{0:03d}'.format(counter)
    if current_qubit is None:
        title = "{counter}_{sample_name}".format(
            sample_name=get_sample_name(),
            counter=str_counter)
    else:
        try:
            title = "{counter}_{sample_name}_qubit{current_qubit}".format(
                sample_name=get_sample_name(),
                counter=str_counter,
                current_qubit=current_qubit)
        except Exception:
            title = "{counter}_{sample_name}".format(
                sample_name=get_sample_name(),
                counter=str_counter)
    return title

##########################
# QtPlot - data plotting #
##########################


def plot_data(data, current_qubit=None, with_title=True):
    counter = data.location_provider.counter
    if with_title is True:
        title = get_title(counter, current_qubit=current_qubit)
    else:
        title = ''
    plots = []
    for value in data.arrays.keys():
        if "set" not in value:
            pl = qc.QtPlot(getattr(data, value))
            pl.subplots[0].setTitle(title)
            pl.subplots[0].showGrid(True, True)
            # pl.counter = counter
            plots.append(pl)
    return plots


def plot_data_live(dataset, meas_param, current_qubit=None, with_title=True):
    counter = dataset.location_provider.counter
    if with_title is True:
        title = get_title(counter, current_qubit=current_qubit)
    else:
        title = ''
    plot = qc.QtPlot()
    for i, name in enumerate(meas_param.names):
        inst_meas_name = "{}_{}".format(meas_param._instrument.name, name)
        plot.add(getattr(dataset, inst_meas_name), subplot=i + 1)
        plot.subplots[i].showGrid(True, True)
        if i == 0:
                # TODO: test if you can get away with only doing this once
            plot.subplots[0].setTitle(title)
        else:
            plot.subplots[i].setTitle("")
    plot.counter = counter
    return plot


###############################
# MatPlot - analysis plotting #
###############################


def plot_cf_data(data1, data2, data3=None, data4=None, datanum=None,
                 subplot=None, xdata=None,
                 legend_labels=[[]] * 4, axes_labels=[[]] * 2):
    """
    Function to plot multiple arrays (of same length) on one axis

    Args:
        data1 (array)
        data2 (array)
        data3 (array): optional
        data4 (array): optional
        datanum (int): number to ascribe to the data optional, should
            match the name under which the
            dataset to reference is saved
        subplot (matplotlib AxesSubplot): optional subplot which this data
            should be plotted on default None will create new one
        xdata (array): optional x axis data, default None results in indices
            of data1 being used
        legend_labels ['d1 label', ..]: optional labels for data
        axes_labels ['xlabel', 'ylabel']: optional labels for axes

    Returns:
        fig, sub (matplotlib.figure.Figure, matplotlib AxesSubplot) if
            subplot kwarg not None
    """
    # set up plotting values
    if subplot is None:
        fig, sub = plt.subplots()
    if datanum is not None:
        sub.figure.counter = datanum
        sub.set_title(get_title(datanum))
    if len(legend_labels) < 4:
        legend_labels.extend([[]] * (4 - len(legend_labels)))
    if xdata is None:
        xdata = np.arange(len(data1))

    # plot data
    sub.plot(xdata, data1, 'r', linewidth=1.0, label=legend_labels[0])
    sub.plot(xdata, data2, 'b', linewidth=1.0, label=legend_labels[1])
    if data3 is not None:
        sub.plot(xdata, data3, 'g', linewidth=1.0, label=legend_labels[2])
    if data4 is not None:
        sub.plot(xdata, data4, 'y', linewidth=1.0, label=legend_labels[3])

    # apply plotting values
    if any(legend_labels):
        sub.legend(loc='upper right', fontsize=10)
    if any(axes_labels):
        sub.set_xlabel(axes_labels[0])
        sub.set_ylabel(axes_labels[1])

    if subplot is None:
        return fig, sub
