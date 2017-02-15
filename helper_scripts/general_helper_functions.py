import qcodes as qc
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from time import localtime, strftime
import logging
from IPython import get_ipython
EXPERIMENT_VARS = {'analysis_location_status': False,
                   'data_location_status': False,
                   'python_log_location_status': False,
                   'jupyter_log_location_status': False}

# TODO: redo load to be less hacky, apparenlty counter is accessible now
#       so data_num can be replaces
# TODO do1d, do2d, plot wrapper


def set_file_locations():
    """
    Wrapper function which calls set_log_location,
    set_data_location and set_analysis_location.
    """
    set_log_locations()
    set_data_location()
    set_analysis_location()


def in_ipynb():
    """
    Tests whether code is being run in an ipython
    (or jupyter) notebook

    Returns:
        bool
    """
    try:
        get_ipython().config
        return True
    except (NameError, AttributeError):
        return False


def set_qubit_count(count):
    """
    Sets qubit_count in EXPERIMENT_VARS dictionary

    Args: count
    """
    EXPERIMENT_VARS['qubit_count'] = count


def get_qubit_count():
    """
    Returns: value of qubit_count in local dict
    """
    try:
        return EXPERIMENT_VARS['qubit_count']
    except KeyError:
        raise KeyError('qubit_count not set, please call qubit_count')


def set_sample_name(name):
    """
    Sets sample_name in EXPERIMENT_VARS dictionary

    Args: sample_name
    """
    EXPERIMENT_VARS['sample_name'] = name


def get_sample_name():
    """
    Returns: value of sample_name in local dict
    """
    try:
        return EXPERIMENT_VARS['sample_name']
    except KeyError:
        raise KeyError('sample_name not set, please call set_sample_name')


def set_data_location():
    """
    Sets location for qcodes to save data based on
    the qcodes.config.user.data_location value and sets global
    variable data_location.
    """
    sample_name = get_sample_name()

    try:
        data_location_format = qc.config.user.data_location.format(sample_name=sample_name)
    except KeyError:
        raise KeyError('data_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    loc_provider = qc.data.location.FormatLocation(
        fmt=data_location_format)
    qc.data.data_set.DataSet.location_provider = loc_provider
    data_location = data_location_format.replace('{counter}', '')
    EXPERIMENT_VARS['data_location_status'] = True
    logging.info('Set data location: {}'.format(data_location))
    print('Set data location: {}'.format(data_location))
    print('-------------------------')


def get_data_location():
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['data_location_status']:
        data_location_format = qc.config.user.data_location.format(sample_name=sample_name)
        data_location = data_location_format.replace('{counter}', '')
        return data_location
    else:
        raise Exception('data_location not set, please call set_data_location')


def set_analysis_location():
    """
    Sets location for qcodes to save analysis data based on
    the qcodes.config.user.analysis_location value and sets global
    variable analysis_location.
    """
    sample_name = get_sample_name()

    try:
        analysis_location = qc.config.user.analysis_location.format(
            sample_name=sample_name)
    except KeyError:
        raise KeyError('analysis_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if not os.path.exists(analysis_location):
        os.makedirs(analysis_location)
    EXPERIMENT_VARS['analysis_location_status'] = True
    logging.info('Set up analysis location: {}'.format(analysis_location))
    print('Set up analysis location: {}'.format(analysis_location))
    print('-------------------------')


def get_analysis_location():
    """
    Returns: analysis_location based on config and local sample_name
    """

    sample_name = get_sample_name()
    if EXPERIMENT_VARS['analysis_location_status']:
        return qc.config.user.analysis_location.format(
            sample_name=sample_name)
    else:
        raise Exception(
            'analysis_location not set, please call set_analysis_location')


def set_log_locations():
    """
    Sets location for qcodes to save log files based on
    the qcodes.config.user.log_location value. Within this folder
    creates python_logs file and (if in notebook) creates jupyter_logs
    file, starts jupyter log. Sets global variables python_log_location,
    jupyter_log_location.
    """
    sample_name = get_sample_name()
    try:
        log_location = qc.config.user.log_location.format(
            sample_name=sample_name)
    except KeyError:
        raise KeyError('log_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if not os.path.exists(log_location):
        os.makedirs(log_location)

    python_log_location = log_location + 'python_logs/'
    if not os.path.exists(python_log_location):
        os.makedirs(python_log_location)
    python_logfile_name = str(python_log_location +
                              strftime('%Y-%m-%d_%H-%M', localtime()) +
                              '_pythonlogfile.log')
    logging.basicConfig(filename=python_logfile_name,
                        level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d_%H-%M-%S')
    EXPERIMENT_VARS['python_log_location_status'] = True
    print('Set up python log location: {}'.format(python_log_location))
    print('-------------------------')

    if in_ipynb():
        jupyter_log_location = log_location + 'jupyter_logs/'
        if not os.path.exists(jupyter_log_location):
            os.makedirs(jupyter_log_location)
        jupyter_logfile_name = str(jupyter_log_location +
                                   strftime('%Y-%m-%d_%H-%M', localtime()) +
                                   '_ipythonlogfile.txt')
        get_ipython().magic("logstart -t '%s'" % jupyter_logfile_name)
        EXPERIMENT_VARS['jupyter_log_location_status'] = True
        print('Set up jupyter log location: {}'.format(jupyter_log_location))
        print('-------------------------')


def get_log_locations():
    if (EXPERIMENT_VARS['python_log_location_status'] or
            EXPERIMENT_VARS['jupyter_log_location_status']):
        logs = {}
        sample_name = EXPERIMENT_VARS['sample_name']
        log_location = qc.config.user.log_location.format(
            sample_name=sample_name)
        if EXPERIMENT_VARS['python_log_location_status']:
            logs['python_log'] = log_location + 'python_logs/'
        if EXPERIMENT_VARS['jupyter_log_location_status']:
            logs['jupyter_log'] = log_location + 'jupyter_logs/'
        return logs
    else:
        raise Exception('no logs set, please run set_log_locations to start'
                        ' logs')


def get_latest_counter():
    """
    Really hacky function which looks in the folder the data is
    being stored in (given that is has the format currently used by the
    transmon team) and looks for the highest number file name and
    returns that number, this should be replaced by having a unique id
    for a dataset which can be linked to the plot, analysis,
    other data sets etc.

    Returns:
        latest counter int
    """
    path = qc.data.data_set.DataSet.location_provider.fmt[:-9]
    files = [re.sub("[^0-9]", "", f) for f in os.listdir(path)]
    int_files = [int(file_name) for file_name in files]
    return max(int_files)


def load(num, plot=True, plot_all=False):
    """
    Function like shownum which loads dataset and QtPlot from default location

    Args:
    num (int)
    plot(default=True): do you want to return a plot as well as dataset?
    plot_all (bool): should all params be plotted or just default?

    Returns:
    dataset (qcodes DataSet)
    plot (QtPlot)
    """
    # format num to have 3 digits if necessary
    if len(str(num)) > 3:
        formatted_num = str(num)
    else:
        formatted_num = "{0:0=3d}".format(num)

    # load data from default location
    loaded_data = qc.load_data(
        qc.data.data_set.DataSet.location_provider.formatter.format(
            qc.data.data_set.DataSet.location_provider.fmt,
            counter=formatted_num))

    # set num to be the same as that in the folder
    loaded_data.data_num = num

    if plot:
        # plot default parameter or all paraters which aren't setpoint
        # arrays if plot_all is True
        plot = qc.QtPlot()
        if not plot_all:
            plot.add(loaded_data.default_parameter_array())
        else:
            for i, array_key in enumerate(loaded_data.arrays.keys()):
                data_array = loaded_data.arrays[array_key]
            if not data_array.is_setpoint:
                plot.add(data_array, subplot=i + 1)
            plot.data_num = num
        return loaded_data, plot
    else:
        return loaded_data


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
        pl (QtPlot)
    """
    dataset = qc.Measure(param).run()
    dataset.data_num = get_latest_counter()
    if plot:
        pl_var = plot_variable or dataset.default_parameter_array()
        pl = qc.QtPlot(pl_var, figsize=(700, 500))
        pl.save()
        pl.data_num = get_latest_counter()
        return dataset, pl
    else:
        return dataset


def save_plot(plot_to_save, name=None):
    """
    Function which saves a plot in qc config analysis_location

    Args:
        plot_to_save (matplotlib AxesSubplot or Figure)
        name (default None): name to save  if none tries
                to save it using plot_to_save.data_num
    """

    fig = getattr(plot_to_save, 'figure', None)

    if fig is None:
        fig = plot_to_save

    if name is None:
        try:
            name = 'analysis_{}.png'.format(fig.data_num)
        except AttributeError:
            raise AttributeError('No name specified and fig has no data_num'
                                 ': please specify a name for the plot')
    else:
        name = name + '_{}.png'.format(fig.data_num)

    analysis_location = get_analysis_location()
    fig.savefig(analysis_location + name)


def plot_cf_data(data1, data2, data3=None, data4=None, datanum=None,
                 subplot=None, xdata=None,
                 legend_labels=None, axes_labels=None):
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
        subplot (matplotlib AxesSubplot): subplot which this data should
                    be plotted on default None will create new one
        xdata (array): x axis data, default None results in indices
                    of data1 being used
        legend_labels ['d1 label', ..]: labels for data, default None
        axes_labels ['xlabel', 'ylabel']: labels for axes

    Returns:
        subplot (matplotlib AxesSubplot)
    """
    # set up plotting values
    if subplot is None:
        fig = plt.figure()
        subplot = plt.subplot(111)
    if datanum is not None:
        subplot.figure.data_num = datanum
    if legend_labels is None:
        legend_labels = [None] * 4
    if xdata is None:
        xdata = np.arange(len(data1))

    # plot data
    subplot.plot(xdata, data1, 'r', linewidth=1.0, label=legend_labels[0])
    subplot.plot(xdata, data2, 'b', linewidth=1.0, label=legend_labels[1])
    if data3 is not None:
        subplot.plot(xdata, data3, 'g', linewidth=1.0, label=legend_labels[2])
    if data4 is not None:
        subplot.plot(xdata, data4, 'y', linewidth=1.0, label=legend_labels[3])

    # apply plotting values
    if legend_labels[0] is not None:
        subplot.legend(loc='upper right', fontsize=10)
    if axes_labels is not None:
        subplot.set_xlabel(axes_labels[0])
        subplot.set_ylabel(axes_labels[1])

    return subplot
