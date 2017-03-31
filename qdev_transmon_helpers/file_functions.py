import qcodes as qc
import re
import os
from qcodes.plots.pyqtgraph import QtPlot
import matplotlib.pyplot as plt
from time import localtime, strftime
import logging
from IPython import get_ipython
EXPERIMENT_VARS = {'analysis_loc': False,
                   'data_loc_fmt': False,
                   'python_log_loc': False,
                   'jupyter_log_loc': False,
                   'pulse_loc': False}


# TODO: test in_ipynb vs qc.in_notebook (and potentially replace)
# TODO: test data.location_provider.counter
# TODO: remove commented out code!
# TODO: docstrings


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

    Args:
        count
    """
    EXPERIMENT_VARS['qubit_count'] = count


def get_qubit_count():
    """
    Returns:
        value of qubit_count in EXPERIMENT_VARS dictionary
    """
    try:
        return EXPERIMENT_VARS['qubit_count']
    except KeyError:
        raise KeyError('qubit_count not set, please call qubit_count')


def set_sample_name(name):
    """
    Sets sample_name in EXPERIMENT_VARS dictionary

    Args:
        sample_name
    """
    EXPERIMENT_VARS['sample_name'] = name


def get_sample_name():
    """
    Returns:
        value of sample_name in EXPERIMENT_VARS dictionary.
    """
    try:
        return EXPERIMENT_VARS['sample_name']
    except KeyError:
        raise KeyError('sample_name not set, please call set_sample_name')


def set_data_location():
    """
    Sets location for qcodes to save data based on the
    qcodes.config.user.data_location value, the
    qc.config.user.data_format and the sample name and sets
    data_loc_fmt in EXPERIMENT_VARS dictionary to True.
    """
    sample_name = get_sample_name()
    try:
        dat_loc = qc.config.user.data_location.format(sample_name=sample_name)
        dat_fmt = qc.config.user.data_format
        dat_loc_fmt = dat_loc + dat_fmt
    except KeyError:
        raise KeyError('data_location or data_format not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if EXPERIMENT_VARS['data_loc_fmt']:
        print('Data location already set at {}.\n'
              'Data file format already set as {}.'.format(dat_loc, dat_fmt))
        print('-------------------------')
    else:
        if not os.path.exists(dat_loc):
            os.makedirs(dat_loc)
        loc_provider = qc.FormatLocation(fmt=dat_loc_fmt)
        qc.data.data_set.DataSet.location_provider = loc_provider
        EXPERIMENT_VARS['data_loc_fmt'] = True
        logging.info('Set data location: {}'.format(dat_loc))
        print('Set data location: {}'.format(dat_loc))
        print('-------------------------')
        logging.info('Set data file format: {}'.format(dat_fmt))
        print('Set data file format: {}'.format(dat_fmt))
        print('-------------------------')


def get_data_location():
    """
    Returns:
        location of data and pngs from qc.config.user.data_location and
        sample_name if data_location set by calling set_data_location()
    """
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['data_loc_fmt']:
        dat_loc = qc.config.user.data_location.format(sample_name=sample_name)
        return dat_loc
    else:
        raise Exception('data_location not set, please call set_data_location')


def get_data_file_format():
    """
    Returns
        format of name given to data files from qc.config.user.data_format
        if data_location set by calling set_data_location()
    """
    if EXPERIMENT_VARS['data_loc_fmt']:
        dat_fmt = qc.config.user.data_format
        return dat_fmt
    else:
        raise Exception('data_location not set, please call set_data_location')


def set_analysis_location():
    """
    Sets location for qcodes to save analysis data based on
    the qcodes.config.user.analysis_location value and sets analysis_loc
    in EXPERIMENT_VARS dictionary to True.
    """
    sample_name = get_sample_name()

    try:
        analysis_location = qc.config.user.analysis_location.format(
            sample_name=sample_name)
    except KeyError:
        raise KeyError('analysis_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if EXPERIMENT_VARS['analysis_loc']:
        print('Analysis location already set at {}.'.format(analysis_location))
        print('-------------------------')
    else:
        if not os.path.exists(analysis_location):
            os.makedirs(analysis_location)
        EXPERIMENT_VARS['analysis_loc'] = True
        logging.info('Set up analysis location: {}'.format(analysis_location))
        print('Set up analysis location: {}'.format(analysis_location))
        print('-------------------------')


def get_analysis_location():
    """
    Returns:
        analysis_location based on qc.config.user.analysis_location
        and sample_name if analysis_location set by calling set_analysis_location()
    """
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['analysis_loc']:
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
    file, starts jupyter log. Sets python_log_loc and jupyter_log_loc
    in EXPERIMENT_VARS dictionary to True.
    """
    sample_name = get_sample_name()
    try:
        log_location = qc.config.user.log_location.format(
            sample_name=sample_name)
        python_log_location = log_location + 'python_logs/'
        jupyter_log_location = log_location + 'jupyter_logs/'
    except KeyError:
        raise KeyError('log_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if EXPERIMENT_VARS['python_log_loc']:
        print('Python log already started at {}.'.format(python_log_location))
        print('-------------------------')
    else:
        if not os.path.exists(python_log_location):
            os.makedirs(python_log_location)
        python_logfile_name = str(python_log_location +
                                  strftime('%Y-%m-%d_%H-%M', localtime()) +
                                  '_pythonlogfile.log')
        logging.basicConfig(filename=python_logfile_name,
                            level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s %(message)s',
                            datefmt='%Y-%m-%d_%H-%M-%S')
        EXPERIMENT_VARS['python_log_loc'] = True
        print('Set up python log location: {}'.format(python_log_location))
        print('-------------------------')
    if EXPERIMENT_VARS['jupyter_log_loc']:
        print('Jupyter log already started at {}.'.format(
            jupyter_log_location))
        print('-------------------------')
    else:
        if in_ipynb():
            if not os.path.exists(jupyter_log_location):
                os.makedirs(jupyter_log_location)
            jupyter_logfile_name = str(jupyter_log_location +
                                       strftime('%Y-%m-%d_%H-%M',
                                                localtime()) +
                                       '_ipythonlogfile.txt')
            get_ipython().magic("logstart -t %s" % jupyter_logfile_name)
            EXPERIMENT_VARS['jupyter_log_loc'] = True
            print('Set up jupyter log location: {}'.format(
                jupyter_log_location))
            print('-------------------------')


def get_log_locations():
    """
    Returns:
        dictionary of 'python_log' and 'jupyter_log' locations where set by
        calling get_log_locations() based on qc.config.user.log_location and
        sample_name.
    """
    if (EXPERIMENT_VARS['python_log_loc'] or
            EXPERIMENT_VARS['jupyter_log_loc']):
        logs = {}
        sample_name = EXPERIMENT_VARS['sample_name']
        log_location = qc.config.user.log_location.format(
            sample_name=sample_name)
        if EXPERIMENT_VARS['python_log_loc']:
            logs['python_log'] = log_location + 'python_logs/'
        if EXPERIMENT_VARS['jupyter_log_loc']:
            logs['jupyter_log'] = log_location + 'jupyter_logs/'
        return logs
    else:
        raise Exception('no logs set, please run set_log_locations to start'
                        ' logs')


def set_pulse_location():
    """
    Sets location for qcodes to save awg pulses on
    the qcodes.config.user.pulse_location value and sets pulse_loc
    in EXPERIMENT_VARS dictionary to True.
    """
    sample_name = get_sample_name()

    try:
        pulse_location = qc.config.user.pulse_location.format(
            sample_name=sample_name)
    except KeyError:
        raise KeyError('pulse_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if EXPERIMENT_VARS['pulse_loc']:
        print('Pulse location already set at {}.'.format(pulse_location))
        print('-------------------------')
    else:
        if not os.path.exists(pulse_location):
            os.makedirs(pulse_location)
        EXPERIMENT_VARS['pulse_loc'] = True
        logging.info('Set up pulse location: {}'.format(pulse_location))
        print('Set up pulse location: {}'.format(pulse_location))
        print('-------------------------')


def get_pulse_location():
    """
    Returns:
        pulse_location based on qc.config.user.pulse_location
        and sample_name if pulse_location set by calling set_pulse_location()
    """
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['pulse_loc']:
        return qc.config.user.pulse_location.format(
            sample_name=sample_name)
    else:
        raise Exception(
            'pulse_location not set, please call set_pulse_location')


def set_file_locations():
    """
    Wrapper function which calls set_log_location,
    set_data_location and set_analysis_location.
    """
    set_log_locations()
    set_data_location()
    set_analysis_location()
    set_pulse_location()


def get_latest_counter(path=None):
    """
    Really hacky function which looks in the folder the data is
    being stored in and looks for the highest number file name and
    returns that number, this should be replaced by having a unique id
    for a dataset which can be linked to the plot, analysis,
    other data sets etc.

    Returns:
        latest counter int
    """
    if path is None:
        path = get_data_location()
    try:
        file_names = [re.sub("[^0-9]", "", f) for f in os.listdir(path)]
    except FileNotFoundError:
        raise FileNotFoundError('No numbered files in ' + path)
    file_ints = [int(f) for f in file_names if f]
    if not file_ints:
        raise FileNotFoundError('No numbered files in ' + path)
    return max(file_ints)


def get_title(counter):
    str_counter = '{0:03d}'.format(counter)
    title = "{counter}_{sample_name}".format(
        sample_name=get_sample_name(),
        counter=str_counter)
    return title



def load(counter, plot=True):
    """
    Function like shownum which loads dataset and (optionally)
    QtPlot from default location

    Args:
        num (int)
        plot(default=True): do you want to return plots as well as dataset?

    Returns:
        dataset (qcodes DataSet)
        plot (QtPlot)
    """
    # load data from default location
    str_counter = '{0:03d}'.format(counter)
    data = qc.load_data(
        qc.DataSet.location_provider.fmt.format(sample_name=get_sample_name(),
                                                counter=str_counter))

    if plot:
        # TODO: replace when counter works: plots = plot_data(data)
        title = get_title(counter) 
        plots = []
        for value in data.arrays.keys():
            if "set" not in value:
                pl = qc.QtPlot(getattr(data, value), figsize=(700, 500))
                pl.subplots[0].setTitle(title)
                pl.subplots[0].showGrid(True, True)
                plots.append(pl)
        return data, plots
    else:
        return data


    # if plot:
    #     # plot default parameter or all paraters which aren't setpoint
    #     # arrays if plot_all is True
    #     plot = qc.QtPlot()
    #     if not plot_all:
    #         plot.add(loaded_data.default_parameter_array())
    #     else:
    #         for i, array_key in enumerate(loaded_data.arrays.keys()):
    #             data_array = loaded_data.arrays[array_key]
    #         if not data_array.is_setpoint:
    #             plot.add(data_array, subplot=i + 1)
    #         plot.data_num = num
    #     return loaded_data, plot
    # else:
    #     return loaded_data

    # if plot:
    #     pl_var = plot_variable or dataset.default_parameter_array()
    #     pl = qc.QtPlot(pl_var, figsize=(700, 500))
    #     pl.save()
    #     pl.data_num = get_latest_counter()
    #     return dataset, pl
    # else:
    #     return dataset


def save_fig(plot_to_save, name='analysis', counter=None, pulse=False):
    """
    Function which saves a matplot figure in analysis_location from
    get_analysis_location()

    Args:
        plot_to_save (matplotlib AxesSubplot or Figure)
        name  (str): plot will be saved with '{data_num}_{name}.png'
            so data_num and/or name must be unique, default 'analysis'
        data_num (int): data_num for fig naming as above, if not specified
            will try to use one from the plot.
    """

    fig = getattr(plot_to_save, 'figure', None)

    if fig is None:
        fig = plot_to_save

    if counter is None:
        try:
            str_counter = '{0:03d}'.format(fig.counter)
        except AttributeError:
            str_counter = ''
            if name is 'analysis':
                raise AttributeError('No name specified and fig has '
                                     'no counter: please specify a '
                                     'name for the plot')
    else:
        str_counter = '{0:03d}'.format(counter)

    full_name = str_counter + '_' + name + '.png'

    if pulse:
        location = get_pulse_location()
    else:
        location = get_analysis_location()
    fig.savefig(location + full_name)


# TODO: check x and y behaviour of plotting if not specified
# TODO: remove commented out code!
# TODO: docstrings

##########################
# QtPlot - data plotting #
##########################


def plot_data(data, with_title=True, key=None):
    counter = data.location_provider.counter
    if with_title is True:
        title = get_title(counter)
    else:
        title = ''
    if key is None:
        plots = []
        for value in data.arrays.keys():
            if "set" not in value:
                pl = qc.QtPlot(getattr(data, value), figsize=(700, 500))
                pl.subplots[0].setTitle(title)
                pl.subplots[0].showGrid(True, True)
                plots.append(pl)
        return plots
    else:
        for value in data.arrays.keys():
            if key in value:
                pl = qc.QtPlot(getattr(data, value), figsize=(700, 500))
                pl.subplots[0].setTitle(title)
                pl.subplots[0].showGrid(True, True)
        return pl

def plot_data_live(dataset, meas_param, with_title=True):
    counter = dataset.location_provider.counter
    if with_title is True:
        title = get_title(counter)
    else:
        title = ''
    plot = qc.QtPlot(figsize=(700 * len(meas_param.names), 500))
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
    else:
        fig, sub = subplot.figure, subplot
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
