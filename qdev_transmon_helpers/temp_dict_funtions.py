from collections import defaultdict
from functools import reduce
import operator
import pprint

from . import get_temp_dict_location

import pickle
import numpy as np
import FileNotFoundError
from . import get_qubit_count

# TODO: make dict_keys better
# docstrings


def validate_temp_dicts():
    validate_calibration_dictionary()
    validate_metadata_dictionary()


################################
# Calibration Dictionary
################################
vna_dict_keys = ['expected_qubit_freqs', 'g_values',
                 'resonances', 'resonator_pushes', 'gatability', 'gate_volts']
alazar_dict_keys = ['current_qubit', 'int_times', 'int_delays', 'cavity_freqs',
                    'cavity_pows', 'demod_freqs', 'pi_pulse_amplitudes',
                    'qubit_freqs', 'pi_pulse_durations', 'pi_pulse_powers',
                    'spec_powers', 'pi_pulse_sigma_cutoff', 'spec_awg_amps',
                    'pi_pulse_awg_amps']


def get_calibration_dict():
    """
    Function which gets the calibration_dict.p saved in get_analysis_location()
    Returns:
        calibration_dict
    """
    path = get_temp_dict_location()
    file_name = 'calibration_dict.p'
    try:
        calibration_dict = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        print('calib dict not found, making one')
        pickle.dump({}, open(path + file_name, 'wb'))
        calibration_dict = {}
    return calibration_dict


def update_calibration_dict(update_dict):
    """
    Function which updates (or creates) a pickled python dictionary saved
    in the same folder as from get_analysis_location() called
    calibration_dict.p

    Args:
        update_dict
    """
    path = get_temp_dict_location()
    file_name = 'calibration_dict.p'
    try:
        dict_to_update = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        dict_to_update = {}

    dict_to_update.update(update_dict)
    pickle.dump(dict_to_update, open(path + file_name, 'wb'))


def set_current_qubit(index, vna=True, alazar=True):
    """
    Sets the value of the 'current_qubit' in the calibration dictionary. This
    can then be used to access index of this qubit when getting or setting
    other values in the calibration dictionary.

    Args:
        index (int) in range of qubit count (starts at 0)
        vna (default True), alazar (default True): keys to check are in the
            calibration dictionary based on default lists
    """
    qubit_count = get_qubit_count()
    if index >= qubit_count:
        raise ValueError('Expects qubit index less than qubit count: {}. '
                         'Received {}'.format(qubit_count, index))
    update_calibration_dict({'current_qubit': index})
    print('current_qubit set to {}'.format(index))


def get_current_qubit():
    """
    Gets the value of the current qubit as set in the calibration dictionary

    Returns:
        qubit index
    """
    c_dict = get_calibration_dict()
    return c_dict['current_qubit']


def set_calibration_val(key, qubit_value, qubit_index=None):
    """
    Sets the relevant qubit_value found at the index of a
    value in the calibration dictionary based on a given key

    Args:
        key (str): key name in calibration dictionary
        qubit_value (float): value for particular qubit
        qubit_index (int) (default None): index of qubit to which val
            corresponds. Default is to use 'current_qubit' value.
    """
    c_dict = get_calibration_dict()
    vals = c_dict[key].copy()
    if qubit_index is None:
        vals[c_dict['current_qubit']] = qubit_value
    else:
        vals[qubit_index] = qubit_value
    update_calibration_dict({key: vals})


def get_calibration_val(key, qubit_index=None):
    """
    Gets the relevant index of a value in the calibration
    dictionary based on a given key

    Args:
        key (str): key name in calibration dictionary
        qubit_index (int) (default None): index of qubit you want the value of
            the key for. Default is to use 'current_qubit' value.

    Returns:
        qubit_value
    """
    c_dict = get_calibration_dict()
    if qubit_index is None:
        return c_dict[key][c_dict['current_qubit']]
    else:
        return c_dict[key][c_dict[qubit_index]]


def validate_calibration_dictionary(vna=True, alazar=True):
    """
    Function which checks that the calibration dictionary contains the
    keys specified in the default lists for vna and alazar measurements
    and that values correspond to lists of the same length as th number of
    qubits. Populates them with 0s if not

    Args:
        vna (default True): check for vna measurement keys?
        alazar (default True): check for alazar measurement keys?
    """
    c_dict = get_calibration_dict()
    qubit_num = get_qubit_count()
    qubit_length_list = np.zeros(get_qubit_count())
    missing_keys = []
    wrong_length_keys = []
    required_keys = alazar_dict_keys * alazar + vna_dict_keys * vna
    for k in required_keys:
        if k not in c_dict:
            missing_keys.append(k)
        elif not (k is 'current_qubit' or len(c_dict[k]) == qubit_num):
            wrong_length_keys.append(k)
    for k in missing_keys:
        c_dict[k] = qubit_length_list
    for k in wrong_length_keys:
        c_dict[k] = qubit_length_list
    update_calibration_dict(c_dict)
    if missing_keys:
        print('{} added to calibration_dictionary'.format(missing_keys))
    if wrong_length_keys:
        print('{} were reset to correct length'.format(wrong_length_keys))


#########################################
# Metadata dictionary
#########################################


def getFromDict(dataDict, mapList):
    """
    Basic helper function which given a mapList (list) as a path gets the value
    from dataDict (nested dictionary structure)
    """
    return reduce(operator.getitem, mapList, dataDict)


def get_metadata_list():
    path = get_temp_dict_location()
    file_name = 'metadata_list.p'
    try:
        metadata_list = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        print('metadata list not found, making one')
        pickle.dump({}, open(path + file_name, 'wb'))
        metadata_list = {}
    return metadata_list


def add_to_metadata_list(*args):
    """
    Args:
        qcodes parameters to be added to EXPERIMENT_VARS['metadata_list']
        which flags them as parameters of interest. The values of these are
        printed when a dataset is loaded.
    """
    metadata_list = get_metadata_list()
    for param in args:
        inst_param_list = [param._instrument.name, param.name]
        if inst_param_list not in metadata_list:
            metadata_list.append(inst_param_list)
    set_metadata_list(metadata_list)


def add_to_metadata_list_manual(instr_name, param_name):
    metadata_list = get_metadata_list()
    inst_param_list = [instr_name, param_name]
    if inst_param_list not in metadata_list:
        metadata_list.append(inst_param_list)
    set_metadata_list(metadata_list)


def set_metadata_list(updated_list):
    path = get_temp_dict_location()
    file_name = 'metadata_list.p'
    pickle.dump(updated_list, open(path + file_name, 'wb'))


def remove_from_metadata_list(*args):
    metadata_list = get_metadata_list()
    for param in args:
        inst_param_list = [param._instrument.name, param.name]
    if inst_param_list in metadata_list:
        metadata_list.remove(inst_param_list)
    set_metadata_list(metadata_list)


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


def validate_metadata_dictionary():
    c_dict = get_metadata_list()
    if c_dict is []:
        print('Metadata list is empty')


##############################################
# Pulse dictionary
##############################################


default_pulse_dict = {'cycle_duration': 20e-6, 'sample_rate': 1e9,
                      'pulse_end': 10e-6, 'pulse_readout_delay': 30e-9,
                      'readout_time': 4e-6, 'marker_time': 500e-9,
                      'marker_readout_delay': 0, 'qubit_time': 1e-6}


def get_pulse_dict():
    """
    Gets the pulse dictionary stored as 'pulse_dict.p' in
    pulse location file or makes an empty one if none found.

    Returns:
        pulse dictionary
    """
    path = get_temp_dict_location()
    file_name = 'pulse_dict.p'
    try:
        pulse_dict = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        print('pulse_dict not found, making one from default values')
        pulse_dict = default_pulse_dict
        pickle.dump(pulse_dict, open(path + file_name, 'wb'))
        pprint.pprint(pulse_dict, width=1)
    return pulse_dict


def update_pulse_dict(update_dict):
    """
    Updates the pulse dictionary 'pulse_dict.p' in
    pulse location file with given dict.

    Args:
        update_dict (dictionary)
    """
    path = get_temp_dict_location()
    file_name = 'pulse_dict.p'
    dict_to_update = get_pulse_dict()
    dict_to_update.update(update_dict)
    pickle.dump(dict_to_update, open(path + file_name, 'wb'))


def update_pulse_val(key, val):
    """
    Sets a value in the pulse dictionary 'pulse_dict.p' in
    pulse location file.

    Args:
        key, val
    """
    p_dict = get_pulse_dict()
    p_dict[key] = val
    update_pulse_dict(p_dict)


def get_pulse_val(key):
    """
    Gets a value from the pulse dictionary 'pulse_dict.p' in
    pulse location file.

    Args:
        key

    Returns:
        val
    """
    p_dict = get_pulse_dict()
    return p_dict[key]


def validate_pulse_dictionary():
    """
    Function which checks that the pulse dictionary contains the
    keys specified in the default list.
    """
    p_dict = get_pulse_dict()
    missing_keys = []
    for k in default_pulse_dict:
        if k not in p_dict:
            missing_keys.append(k)
    for k in missing_keys:
        p_dict[k] = default_pulse_dict[k]
    update_pulse_dict(p_dict)
    if missing_keys:
        print('{} added to pulse_dictionary'.format(missing_keys))
