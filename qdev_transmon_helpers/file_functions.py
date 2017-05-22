import qcodes as qc
import re
import os
import numpy as np
import matplotlib.pyplot as plt
from time import localtime, strftime
import logging
from IPython import get_ipython
from os.path import abspath
from os.path import sep
from collections import defaultdict
from functools import reduce
import operator
from dateutil import parser
import warnings

EXPERIMENT_VARS = {'analysis_loc': False,
                   'data_loc_fmt': False,
                   'python_log_loc': False,
                   'jupyter_log_loc': False,
                   'pulse_loc': False,
                   'metadata_list': []}


# TODO: test in_ipynb vs qc.in_notebook (and potentially replace)
# TODO: replace plot_cf_data with a less rubbish version
# TODO: either get 'file format' to behave as wanted or remove freedom
# TODO: docstrings
# TODO: plot_subset to work with both dimensions soft?

def getFromDict(dataDict, mapList):
    """
    Basic helper function which given a mapList (list) as a path gets the value
    from dataDict (nested dictionary structure)
    """
    return reduce(operator.getitem, mapList, dataDict)


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
        count (int)
    """
    EXPERIMENT_VARS['qubit_count'] = count


def get_qubit_count():
    """
    Returns:
        value of qubit_count in EXPERIMENT_VARS dictionary (int)
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
        value of sample_name in EXPERIMENT_VARS dictionary (str).
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
        location of data and pngs (str): from qc.config.user.data_location and
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
        format of name given to data files (str): from
        qc.config.user.data_format if data_location set by calling
        set_data_location()
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
        analysis_location (str) based on qc.config.user.analysis_location
        and sample_name if analysis_location set by calling
        set_analysis_location()
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
    warnings.simplefilter('error', UserWarning)

    sample_name = get_sample_name()
    try:
        log_location = abspath(qc.config.user.log_location.format(
            sample_name=sample_name))
        python_log_location = sep.join([log_location, 'python_logs', ""])
        jupyter_log_location = sep.join([log_location, 'jupyter_logs', ""])
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
        python_logfile_name = "{}{}{}".format(python_log_location,
                                  strftime('%Y-%m-%d_%H-%M-%S', localtime()),
                                  '_pythonlogfile.log')
        logging.basicConfig(filename=python_logfile_name,
                            level=logging.INFO,
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
            jupyter_logfile_name = "{}{}{}".format(jupyter_log_location,
                                       strftime('%Y-%m-%d_%H-%M-%S',
                                                localtime()),
                                       '_ipythonlogfile.txt')
            try:
                get_ipython().magic("logstart -t {} append".format(jupyter_logfile_name))
                EXPERIMENT_VARS['jupyter_log_loc'] = True
                print('Set up jupyter log location: {}'.format(
                jupyter_log_location))
                print('-------------------------')
            except Warning as w:
                print('Could not set up jupyter log: ', w)
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
        log_location = abspath(qc.config.user.log_location.format(
            sample_name=sample_name))
        if EXPERIMENT_VARS['python_log_loc']:
            logs['python_log'] = sep.join([log_location, 'python_logs', ""])
        if EXPERIMENT_VARS['jupyter_log_loc']:
            logs['jupyter_log'] = sep.join([log_location, 'jupyter_logs', ""])
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
        pulse_location (str): based on qc.config.user.pulse_location
        and sample_name if pulse_location set by calling set_pulse_location()
    """
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['pulse_loc']:
        return qc.config.user.pulse_location.format(
            sample_name=sample_name)
    else:
        raise Exception(
            'pulse_location not set, please call set_pulse_location')


def add_to_metadata_list(*args):
    """
    Args:
        qcodes parameters to be added to EXPERIMENT_VARS['metadata_list']
        which flags them as parameters of interest. The values of these are
        printed when a dataset is loaded.
    """
    for param in args:
        inst_param_list = [param._instrument.name, param.name]
        if inst_param_list not in EXPERIMENT_VARS['metadata_list']:
            EXPERIMENT_VARS['metadata_list'].append(inst_param_list)


def remove_from_metadata_list(*args):
    """
    Args:
        qcodes parameters to be removed from EXPERIMENT_VARS['metadata_list']
        so that their values are no longer printed when a dataset is loaded.
    """
    for param in args:
        inst_param_list = [param._instrument.name, param.name]
    if inst_param_list in EXPERIMENT_VARS['metadata_list']:
        EXPERIMENT_VARS['metadata_list'].remove(inst_param_list)


def set_file_locations():
    """
    Wrapper function which calls set_log_location,
    set_data_location, set_analysis_location and set_pulse_location.
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
        latest counter (int)
    """
    if path is None:
        path = get_data_location()
    try:
        file_names = [re.sub("[^0-9]", "", f) for f in os.listdir(path)]
    except OSError as e:
        raise OSError('Error looking for numbered files in {}:'
                      ''.format(path, e))
    file_ints = [int(f) for f in file_names if f]
    if not file_ints:
        raise OSError('No numbered files in ' + path)
    return max(file_ints)


def get_title(counter):
    """
    Returns:
        title (str): creates a title based on the counter and sample name to
        be used for figures
    """
    if counter is None:
        return get_sample_name()
    else:
        str_counter = '{0:03d}'.format(counter)
        return "{counter}_{sample_name}".format(
            sample_name=get_sample_name(),
            counter=str_counter)


def get_metadata(dataset, display=True, specific_list=None):
    """
    Function which gets the metadata dictionary for a dataset.

    Args:
        dataset(qcodes dataset ot int): dataset or int counter of dataset for
            which we want metadata
        display (default True): should metadata be printed as well as returned?
        specific_list (default None):

    Returns:
        meta_dict dictionary of the form:
            {'instr': {'param':
                            {'value': val,
                            'unit': u}}}
            for a parameter with name 'param' belonging to instrument 'instr'
            with value 'val' and unit 'u'
    """
    missing_keys = []
    if isinstance(dataset, int):
        dataset = load(dataset, plot=False)
    snapshot = dataset.snapshot()
    meta_dict = defaultdict(dict)
    for instr, param in specific_list or EXPERIMENT_VARS['metadata_list']:
        try:
            unit = getFromDict(snapshot, ["station",
                                          "instruments",
                                          instr,
                                          "parameters",
                                          param,
                                          "unit"])
            value = getFromDict(snapshot, ["station",
                                           "instruments",
                                           instr,
                                           "parameters",
                                           param,
                                           "value"])
            meta_dict[instr][param] = {}
            meta_dict[instr][param]['value'] = value
            meta_dict[instr][param]['unit'] = unit
        except KeyError:
            missing_keys.append([instr, param])
    if display:
        print_metadata(meta_dict)
    if len(missing_keys) > 0:
        print('\nSome of the specified parameters were not found in the '
              'snapshot metadata for this dataset:')
        for missing in missing_keys:
            print(missing[0] + ' ' + missing[1])
    return meta_dict


def print_metadata(meta_dict):
    """
    Function which given a metadata dictionary generated by get_metadata
    prints it.
    """
    for instr in meta_dict:
        print(instr)
        for param in meta_dict[instr]:
            print('\t{} : {} {}'.format(param,
                                        meta_dict[instr][param]['value'],
                                        meta_dict[instr][param]['unit']))


def load(counter, plot=True, metadata=True, matplot=False):
    """
    Function like shownum which loads dataset, (optionally)
    QtPlot from default location and (optionally) prints metadata for that
    measurement given current EXPERIMENT_VARS['metadata_list'] entries

    Args:
        counter (int)
        plot(default True): do you want to return plots as well as dataset?
        metadata (default True): do you want to print the metadata?
        matplot (bool) (default False): default is to QtPlot the data
    Returns:
        dataset (qcodes DataSet)
        plot (QtPlot): optional
    """
    str_counter = '{0:03d}'.format(counter)
    data = qc.load_data(
        qc.DataSet.location_provider.fmt.format(sample_name=get_sample_name(),
                                                counter=str_counter))
    data.data_num = counter

    if metadata:
        get_metadata(data, display=True)
        get_data_duration(data)
    if plot:
        plots = plot_data(data, matplot=matplot)
        return data, plots
    else:
        return data


def get_data_duration(dataset):
    if 'loop' in dataset.metadata.keys():
        start = parser.parse(dataset.metadata['loop']['ts_start'])
        end = parser.parse(dataset.metadata['loop']['ts_end'])
    elif 'measurement' in dataset.metadata.keys():
        start = parser.parse(dataset.metadata['measurement']['ts_start'])
        end = parser.parse(dataset.metadata['measurement']['ts_end'])
    else:
        raise KeyError('Could not find "loop" or "measurement" in dataset.metadata.keys()')
    dur = end - start
    m, s = divmod(dur.total_seconds(), 60)
    h, m = divmod(m, 60)
    print("data taking duration %d:%02d:%02d" % (h, m, s))
    return dur.total_seconds()


def save_plot(dataset, key):
    """
    Function for saving one of the subplots from a dataset based on a
    given key.

    Args:
        dataset (int or qcodes dataset): dataset to plot or counter
            from which dataset is loaded.
        key (str): string specifying parameter array/values to plot
        (eg key="mag" will search arrays of the dataset for any with "mag"
         in the name).
    """
    if isinstance(dataset, int):
        dataset = load(dataset, plot=False)
    plot = plot_data(dataset, key=key)
    plot.save()
    return plot


def save_fig(plot_to_save, name='analysis', counter=None, pulse=False):
    """
    Function which saves a matplot figure in analysis_location from
    get_analysis_location()

    Args:
        plot_to_save (matplotlib AxesSubplot or Figure)
        name  (str): plot will be saved with '{data_num}_{name}.png'
            so data_num and/or name must be unique, default 'analysis'
        counter (int): counter for fig naming as above, if not specified
            will try to use one from the plot.
        pulse (bool): if true saves fig in pulse_lib folder from config,
            otherwise save in analysis folder from config, default False.
    """

    fig = getattr(plot_to_save, 'figure', plot_to_save) or plot_to_save

    if counter is None:
        try:
            str_counter = '{0:03d}'.format(fig.data_num)
        except AttributeError:
            str_counter = ''
            if name is 'analysis':
                raise AttributeError('No name specified and fig has '
                                     'no data_num: please specify a '
                                     'name for the plot')
    else:
        str_counter = '{0:03d}'.format(counter)

    full_name = str_counter + '_' + name + '.png'

    if pulse:
        location = get_pulse_location()
    else:
        location = get_analysis_location()
    fig.savefig(location + full_name)


##########################
# QtPlot - data plotting #
##########################


def plot_data(data, key=None, matplot=False):
    """
    Plotting function for plotting arrays of a dataset in seperate
    QtPlots, cannot be used with live_plot if there is more than one
    subplot. If key is specified returns only the plot for the array
    with the key in the name.
    Args:
        data (qcodes dataset): dataset to be plotted
        key (str): key which if specified is used to select the first array
            from the dataset for plotting with a name which contains this key.
        matplot (bool) (default False): default is to QtPlot the data
    """
    if hasattr(data, "data_num"):
        title = title = get_title(data.data_num)
    else:
        title = ""
    if key is None:
        plots = []
        for value in data.arrays.keys():
            if "set" not in value:
                if matplot:
                    pl = qc.MatPlot(getattr(data, value))
                else:
                    pl = qc.QtPlot(getattr(data, value), figsize=(700, 500))
                    pl.subplots[0].setTitle(title)
                    pl.subplots[0].showGrid(True, True)
                plots.append(pl)
        return plots
    else:
        try:
            key_array_name = [v for v in data.arrays.keys() if key in v][0]
        except IndexError:
            raise KeyError('key: {} not in data array '
                           'names: {}'.format(key,
                                              list(data.arrays.keys())))
        if matplot:
            pl = qc.MatPlot(getattr(data, value))
        else:
            pl = qc.QtPlot(getattr(data, key_array_name), figsize=(700, 500))
            pl.subplots[0].setTitle(title)
            pl.subplots[0].showGrid(True, True)
        return pl


def plot_data_single_window(dataset, meas_param, key=None):
    """
    Plotting function for plotting arrays of a dataset in a single window
    (works with live plot but not great if you want to save a png)

    Args:
        dataset (qcodes dataset): dataset to be plotted
        meas_param: parameter being measured
        key (str) (default None): string to search the array names of the
            measured param for, all arrays with key in name will be added
            as subplots. Default is to plot all.
    """
    if hasattr(dataset, "data_num"):
        title = title = get_title(dataset.data_num)
    else:
        title = ""
    plot_array_names = []
    for array_name in meas_param.names:
        if (key is None) or (key in array_name):
            plot_array_names.append("{}_{}".format(meas_param._instrument.name,
                                                   array_name))
    if len(plot_array_names) == 0:
        raise KeyError('key: {} not in parameter array '
                       'names: {}'.format(key,
                                          list(meas_param.names)))
    plot = qc.QtPlot(figsize=(700 * len(plot_array_names), 500))
    for i, plot_array_name in enumerate(plot_array_names):
        plot.add(getattr(dataset, plot_array_name), subplot=i)
        plot.subplots[i].showGrid(True, True)
    plot.subplots[0].setTitle(title)
    return plot


###############################
# MatPlot - analysis plotting #
###############################


def plot_cf_data(data_list, data_num=None,
                 subplot=None, xdata=None,
                 legend_labels=[], axes_labels=[]):
    """
    Function to plot multiple arrays (of same length) on one axis

    Args:
        data_list: list of arrays to be compared
        data_num (int): number to ascribe to the data optional, should
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
    if data_num is not None:
        fig.data_num = data_num
        sub.set_title(get_title(data_num))
    if (len(legend_labels) == 0) or (len(legend_labels) != len(data_list)):
        legend_labels = [[]]*len(data_list)
    if (len(axes_labels) == 0) or (len(axes_labels) != 2):
        axes_labels = [[]]*2
    if xdata is None:
        xdata = np.arange(len(data_list[0]))

    # plot data
    for i, data in enumerate(data_list):
        sub.plot(xdata, data, linewidth=1.0, label=legend_labels[i])

    # apply plotting values
    if any(legend_labels):
        box = sub.get_position()
        sub.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        sub.legend(fontsize=10, bbox_to_anchor=(1, 1))
    if any(axes_labels):
        sub.set_xlabel(axes_labels[0])
        sub.set_ylabel(axes_labels[1])

    if subplot is None:
        return fig, sub

def line_cut(array, vals, axis='y', data_num=None):
    x_data = np.array(getattr(array, "set_arrays")[1][0])
    y_data = np.array(getattr(array, "set_arrays")[0])
    x_label = '{} ({})'.format(getattr(array, "set_arrays")[1].label, getattr(array, "set_arrays")[1].unit)
    y_label = '{} ({})'.format(getattr(array, "set_arrays")[0].label,getattr(array, "set_arrays")[0].unit)
    z_label = array.name
    values = np.zeros(len(vals))
    if axis is 'x':
        z_data = np.zeros((len(vals), len(y_data)))
        for i, v in enumerate(vals):
            x_index = np.where(x_data == v)
            z_data[i] = array[:, x_index]
        fig, sub = plot_cf_data(z_data,
                                xdata=y_data,
                                data_num=data_num,
                                legend_labels=["{} {}".format(v,x_label) for v in vals],
                                axes_labels=[y_label, z_label])
    
    elif axis is 'y':
        z_data = np.zeros((len(vals), len(x_data)))
        for i, v in enumerate(vals):
            y_index = np.where(y_data == v)
            z_data[i] = array[y_index, :]
        fig, sub = plot_cf_data(z_data,
                                xdata=x_data,
                                data_num=data_num,
                                legend_labels=[str(v)+y_label for v in vals],
                                axes_labels=[x_label, z_label])
    return fig, sub

def plot_subset(array, x_start=None, x_stop=None, y_start=None, y_stop=None):
    x_data = np.array(getattr(array, "set_arrays")[1][0])
    y_data = np.array(getattr(array, "set_arrays")[0])
    x_label = '{} ({})'.format(getattr(array, "set_arrays")[1].label, getattr(array, "set_arrays")[1].unit)
    y_label = '{} ({})'.format(getattr(array, "set_arrays")[0].label,getattr(array, "set_arrays")[0].unit)
    x_indices = np.where((x_data >= (x_start or -1*np.inf)) &  (x_data <= (x_stop or np.inf)))[0]
    y_indices = np.where((y_data >= (y_start or -1*np.inf)) &  (y_data <= (y_stop or np.inf)))[0]
    pl = qc.MatPlot(x_data[x_indices[0]:x_indices[-1]],
                    y_data[y_indices[0]:y_indices[-1]],
                    array[y_indices[0]:y_indices[-1], x_indices[0]:x_indices[-1]])
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(array.name)
    return pl