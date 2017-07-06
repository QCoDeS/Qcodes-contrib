import pprint
import collections
from . import get_temp_dict_location
import os
import pickle
import copy
import numpy as np
from . import get_qubit_count


################################
# Calibration Dictionary
################################


vna_dict_keys = ['expected_qubit_freq', 'g_value', 'bare_res_freq',
                 'pushed_res_freq', 'gate_volt']

alazar_dict_keys = ['int_time', 'int_delay', 'cavity_freq',
                    'cavity_pow', 'localos_pow', 'demod_freq',
                    'qubit_freq', 'pi_pulse_pow', 'spec_pow']

default_pulse_dict = {'cycle_time': 20e-6, 'sample_rate': 1e9,
                      'pulse_end': 10e-6, 'pulse_readout_delay': 30e-9,
                      'readout_time': 4e-6, 'readout_amp': 1,
                      'marker_time': 500e-9,
                      'marker_readout_delay': 0, 'qubit_spec_time': 1e-6,
                      'pulse_mod_time': 1.5e-6, 'pi_pulse_amp': 1,
                      'pi_half_pulse_amp': 0.5, 'pi_pulse_sigma': None,
                      'pi_pulse_dur': None, 'sigma_cutoff': 4,
                      'z_pulse_amp': None, 'z_pulse_dur': None,
                      'z_half_pulse_amp': None, 'drag_coef': 0.5}


def dd_f():
    return None


default_calib_dict = collections.defaultdict(dd_f)


def get_calibration_dict():
    """
    Function which gets the calibration_dict.p saved in
    get_temp_dict_location()

    Returns:
        calibration_dict
    """
    path = get_temp_dict_location()
    file_name = 'calibration_dict.p'
    if os.path.exists(path + file_name):
        calibration_dict = pickle.load(open(path + file_name, "rb"))
    else:
        print('calib dict not found, making one')
        calibration_dict = default_calib_dict
        for k in default_pulse_dict:
            v = np.array([default_pulse_dict[k]] * get_qubit_count())
            calibration_dict[k] = v
        pickle.dump(calibration_dict, open(path + file_name, 'wb'))
    return calibration_dict


def _update_calibration_dict(update_dict):
    """
    Function which updates (or creates) a pickled python dictionary saved
    in the same folder as from get_temp_dict_location() called
    calibration_dict.p

    Args:
        update_dict
    """
    path = get_temp_dict_location()
    file_name = 'calibration_dict.p'
    if os.path.exists(path + file_name):
        calibration_dict = pickle.load(open(path + file_name, "rb"))
    else:
        print('calib dict not found, making one')
        calibration_dict = default_calib_dict
        for k in default_pulse_dict:
            v = np.array([default_pulse_dict[k]] * get_qubit_count())
            calibration_dict[k] = v
    calibration_dict.update(update_dict)
    pickle.dump(calibration_dict, open(path + file_name, 'wb'))


def set_current_qubit(index):
    """
    Sets the value of the 'current_qubit' in the calibration dictionary. This
    can then be used to access index of this qubit when getting or setting
    other values in the calibration dictionary.

    Args:
        index (int) in range of qubit count (starts at 0)
    """
    qubit_count = get_qubit_count()
    if index >= qubit_count:
        raise ValueError('Expects qubit index less than qubit count: {}. '
                         'Received {}'.format(qubit_count, index))
    _update_calibration_dict({'current_qubit': index})


def get_current_qubit():
    """
    Gets the value of the current qubit as set in the calibration dictionary

    Returns:
        qubit index
    """
    c_dict = get_calibration_dict()
    return c_dict['current_qubit']


def get_allowed_keys(vna=True, alazar=True, pulse=True):
    """
    Gets the keys allowed in the calibration dictionary. These are hard coded.

    Returns:
        set of keys
    """
    keys_list = []
    if vna:
        keys_list.extend(vna_dict_keys)
    if alazar:
        keys_list.extend(alazar_dict_keys)
    if pulse:
        keys_list.extend(default_pulse_dict.keys())
    return keys_list


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
    allowed_keys = get_allowed_keys()
    if key not in allowed_keys:
        raise Exception('{} not in allowed calibration dict keys: {}'
                        ''.format(key, allowed_keys))
    c_dict = get_calibration_dict()
    if key in c_dict:
        vals = copy.copy(c_dict[key])
    else:
        vals = np.zeros(get_qubit_count())
    vals[qubit_index or c_dict['current_qubit']] = qubit_value
    _update_calibration_dict({key: vals})


def set_calibration_array(key, array):
    """
    Sets values for all qubits in the calibration dictionary

    Args:
        key (str) of the values to be set
        array (list or numpy array) of the values for each qubit
    """
    qubit_count = get_qubit_count()
    if len(array) != qubit_count:
        raise Exception('array given must be the same length as the number'
                        ' of qubits: {}'.fromat(qubit_count))
    allowed_keys = get_allowed_keys()
    if key not in allowed_keys:
        raise Exception('{} not in allowed calibration dict keys: {}'
                        ''.format(key, allowed_keys))
    _update_calibration_dict({key: np.array(array)})


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
    return c_dict[key][qubit_index or c_dict['current_qubit']]


def get_calibration_array(key):
    """
    Gets the values from the calibration dictionary for all qubits.

    Args:
        key (str): key name in calibration dictionary

    Returns:
        qubit_values (array)
    """
    c_dict = get_calibration_dict()
    return c_dict[key]


def print_pulse_settings():
    """
    Pretty prints pulse settings
    """
    calib_dict = get_calibration_dict()
    pulse_dict = {}
    for k in default_pulse_dict:
        pulse_dict[k] = calib_dict[k]
    pprint.pprint(pulse_dict, width=1)


#########################################
# Metadata List
#########################################


def get_metadata_list():
    path = get_temp_dict_location()
    file_name = 'metadata_list.p'
    try:
        metadata_list = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        print('metadata list not found, making one')
        metadata_list = []
        pickle.dump(metadata_list, open(path + file_name, 'wb'))
    return metadata_list


def add_to_metadata_list(*args):
    """
    Args:
        qcodes parameters to be added to metadata_list.p
        which flags them as parameters of interest. The values of these are
        printed when a dataset is loaded.
    """
    metadata_list = get_metadata_list()
    for param in args:
        inst_param_tuple = (param._instrument.name, param.name)
        if inst_param_tuple not in metadata_list:
            metadata_list.append(inst_param_tuple)
    _set_metadata_list(metadata_list)


def add_to_metadata_list_manual(instr_name, param_name):
    """
    Args:
        intr_name (str): name of instrument which param belongs to
        param_name (str): name of param
    """
    metadata_list = get_metadata_list()
    inst_param_tuple = (instr_name, param_name)
    if inst_param_tuple not in metadata_list:
        metadata_list.append(inst_param_tuple)
    _set_metadata_list(metadata_list)


def _set_metadata_list(updated_list):
    """
    Finction which creates 'metadata_list.p' in get_temp_dict_location
    and dumps picled list there (overwrites and existing list)

    Args:
        list to write
    """
    path = get_temp_dict_location()
    file_name = 'metadata_list.p'
    pickle.dump(updated_list, open(path + file_name, 'wb'))


def remove_from_metadata_list(*args):
    """
    Function which removes the instr_param_tuples from the matadata_list.p
    if present

    Args:
        qcodes parameters
    """
    metadata_list = get_metadata_list()
    for param in args:
        inst_param_tuple = (param._instrument.name, param.name)
    if inst_param_tuple in metadata_list:
        metadata_list.remove(inst_param_tuple)
    _set_metadata_list(metadata_list)


def remove_from_metadata_list_manual(instr_name, param_name):
    """
    Function which removes the instr_param_tuple from the matadata_list.p
    if present

    Args:
        intr_name (str): name of instrument which param belongs to
        param_name (str): name of param
    """
    metadata_list = get_metadata_list()
    inst_param_tuple = (instr_name, param_name)
    if inst_param_tuple in metadata_list:
        metadata_list.remove(inst_param_tuple)
    _set_metadata_list(metadata_list)
