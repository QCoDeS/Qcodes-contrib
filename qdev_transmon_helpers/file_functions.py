import qcodes as qc
import re
import os
from time import localtime, strftime
import logging
from IPython import get_ipython
from os.path import abspath
from os.path import sep
import warnings

EXPERIMENT_VARS = {'analysis_loc': False,
                   'data_loc_fmt': False,
                   'python_log_loc': False,
                   'ipython_log_loc': False,
                   'temp_dict_loc': False,
                   'pulse_loc': False}


# TODO: either get 'file format' to behave as wanted or remove freedom


def in_ipynb():
    """
    Tests whether code is being run in an ipython session

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
    Sets qubit_count in EXPERIMENT_VARS dictionary and the global variable
    (to avoid confusion when using spyder)

    Args:
        count (int)
    """
    global qubit_count
    qubit_count = count
    EXPERIMENT_VARS['qubit_count'] = count


def get_qubit_count():
    """
    Returns:
        value of qubit_count in EXPERIMENT_VARS dictionary (int)
    """
    global qubit_count
    try:
        count = EXPERIMENT_VARS['qubit_count']
        qubit_count = count
        return count
    except KeyError:
        raise KeyError('qubit_count not set, please call set_qubit_count')


def set_sample_name(name):
    """
    Sets sample_name in EXPERIMENT_VARS dictionary and the global variable
    (to avoid confusion when using spyder)

    Args:
        sample_name
    """
    global sample_name
    sample_name = name
    EXPERIMENT_VARS['sample_name'] = name


def get_sample_name():
    """
    Returns:
        value of sample_name in EXPERIMENT_VARS dictionary (str).
    """
    global sample_name
    try:
        name = EXPERIMENT_VARS['sample_name']
        sample_name = name
        return name
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


def set_temp_dict_location():
    """
    Sets location for qcodes to save analysis data based on
    the qcodes.config.user.analysis_location value and sets analysis_loc
    in EXPERIMENT_VARS dictionary to True.
    """
    sample_name = get_sample_name()

    try:
        temp_dict_location = qc.config.user.temp_dict_location.format(
            sample_name=sample_name)
    except KeyError:
        raise KeyError('temp_dict_location not set in config, see '
                       '"https://github.com/QCoDeS/Qcodes/blob/master'
                       '/docs/examples/Configuring_QCoDeS.ipynb"')
    if EXPERIMENT_VARS['temp_dict_loc']:
        print('Temp dictionary location already set at {}.'.format(
            temp_dict_location))
        print('-------------------------')
    else:
        if not os.path.exists(temp_dict_location):
            os.makedirs(temp_dict_location)
        EXPERIMENT_VARS['temp_dict_loc'] = True
        logging.info('Set up temp dict location: {}'.format(
            temp_dict_location))
        print('Set up temp dict location: {}'.format(temp_dict_location))
        print('-------------------------')


def get_temp_dict_location():
    """
    Returns:
        analysis_location (str) based on qc.config.user.analysis_location
        and sample_name if analysis_location set by calling
        set_analysis_location()
    """
    sample_name = get_sample_name()
    if EXPERIMENT_VARS['temp_dict_loc']:
        return qc.config.user.temp_dict_location.format(
            sample_name=sample_name)
    else:
        raise Exception(
            'temp_dict_location not set, please call set_temp_dict_location')


def set_log_locations():
    """
    Sets location for qcodes to save log files based on
    the qcodes.config.user.log_location value. Within this folder
    creates python_logs file and (if in notebook) creates ipython_logs
    file, starts ipython log. Sets python_log_loc and ipython_log_loc
    in EXPERIMENT_VARS dictionary to True.
    """
    warnings.simplefilter('error', UserWarning)

    sample_name = get_sample_name()
    try:
        log_location = abspath(qc.config.user.log_location.format(
            sample_name=sample_name))
        python_log_location = sep.join([log_location, 'python_logs', ""])
        ipython_log_location = sep.join([log_location, 'ipython_logs', ""])
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
                                              strftime('%Y-%m-%d_%H-%M-%S',
                                                       localtime()),
                                              '_pythonlogfile.log')
        logging.basicConfig(filename=python_logfile_name,
                            level=logging.INFO,
                            format='%(asctime)s %(levelname)s %(message)s',
                            datefmt='%Y-%m-%d_%H-%M-%S')
        EXPERIMENT_VARS['python_log_loc'] = True
        print('Set up python log location: {}'.format(python_log_location))
        print('-------------------------')
    if EXPERIMENT_VARS['ipython_log_loc']:
        print('ipython log already started at {}.'.format(
            ipython_log_location))
        print('-------------------------')
    else:
        if in_ipynb():
            if not os.path.exists(ipython_log_location):
                os.makedirs(ipython_log_location)
            ipython_logfile_name = "{}{}{}".format(
                ipython_log_location,
                strftime('%Y-%m-%d_%H-%M-%S',
                         localtime()),
                '_ipythonlogfile.txt')
            try:
                get_ipython().magic("logstart -t {} append".format(
                    ipython_logfile_name))
                EXPERIMENT_VARS['ipython_log_loc'] = True
                print('Set up ipython log location: {}'.format(
                    ipython_log_location))
                print('-------------------------')
            except Warning as w:
                print('Could not set up ipython log: ', w)
                print('-------------------------')


def get_log_locations():
    """
    Returns:
        dictionary of 'python_log' and 'ipython_log' locations where set by
        calling get_log_locations() based on qc.config.user.log_location and
        sample_name.
    """
    if (EXPERIMENT_VARS['python_log_loc'] or
            EXPERIMENT_VARS['ipython_log_loc']):
        logs = {}
        sample_name = EXPERIMENT_VARS['sample_name']
        log_location = abspath(qc.config.user.log_location.format(
            sample_name=sample_name))
        if EXPERIMENT_VARS['python_log_loc']:
            logs['python_log'] = sep.join([log_location, 'python_logs', ""])
        if EXPERIMENT_VARS['ipython_log_loc']:
            logs['ipython_log'] = sep.join([log_location, 'ipython_logs', ""])
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


def set_file_locations():
    """
    Wrapper function which calls set_log_location,
    set_data_location, set_analysis_location and set_pulse_location.
    """
    set_log_locations()
    set_data_location()
    set_analysis_location()
    set_pulse_location()
    set_temp_dict_location()


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
