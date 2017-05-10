import pickle
import numpy as np
from . import get_analysis_location, g_from_qubit, get_qubit_count

# TODO: make dict_keys better

vna_dict_keys = ['expected_qubit_positions', 'g_values',
                 'resonances', 'resonator_pushes', 'gatability', 'gate_volts']
alazar_dict_keys = ['current_qubit', 'int_times', 'int_delays', 'cavity_freqs',
                    'cavity_pows', 'demod_freqs', 'pi_pulse_amplitudes', 't1s',
                    't1_errors', 't2s', 't2_errrors', 'actual_qubit_positions',
                    'pi_pulse_durations', 'pi_pulse_powers', 'spec_powers']


def get_calibration_dict():
    """
    Function which gets the calibration_dict.p saved in get_analysis_location()
    Returns:
        calibration_dict
    """
    path = get_analysis_location()
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
    path = get_analysis_location()
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
    validate_calibration_dictionary(vna, alazar)
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


def recalculate_g(dec_chans=None):
    """
    Function which uses the values in the calibration dictionary for expected
    qubit position, actual position, resonator push and g value to recalculate
    the g value for the current qubit and compare it to the old value.

    Args:
        dec_chans (default None): if dec_chans given, gate value of current
            qubit is compared to value at which resonator data was taken to
            check validity of recalculated g
    """
    c_dict = get_calibration_dict()
    expected = c_dict['expected_qubit_positions'][c_dict['current_qubit']]
    actual = c_dict['actual_qubit_positions'][c_dict['current_qubit']]
    res_data = c_dict['resonator_pushes'][c_dict['current_qubit']]
    old_g = c_dict['g_values'][c_dict['current_qubit']]
    new_g = g_from_qubit(actual, res_data[0], res_data[2])
    if dec_chans is not None:
        current_voltage = dec_chans[c_dict['current_qubit']].get_latest()
        if (c_dict['gatability'][c_dict['current_qubit']] and
                (current_voltage !=
                    c_dict['gate_volts'][c_dict['current_qubit']])):
            print('New g factor calculated will not be a good estimate '
                  'as current gate value is not the same as the value when '
                  'the push on the resonator was measured.')
    print('expected qubit freq: {}\n (from g of {}, push on resonator {})\n'
          'actual qubit freq: {}\n (for same push gives g of {}'.format(
              expected, old_g, res_data[2], actual, new_g))
    return new_g


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
