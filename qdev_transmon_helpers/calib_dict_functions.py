import pickle
import numpy as np
from . import get_analysis_location, g_from_qubit, get_qubit_count
vna_dict_keys = ['expected_qubit_positions', 'g_values',
                 'resonances', 'resonator_pushes', 'gatability', 'gate_volts']
alazar_dict_keys = ['current_qubit', 'int_times', 'int_delays', 'cavity_freqs',
                    'cavity_pows', 'demod_freqs', 'pi_pulse_amplitudes', 't1s',
                    't1_errors', 'actual_qubit_positions',
                    'pi_pulse_durations', 'pi_pulse_powers']


def get_calibration_dict():
    """
    Function which gets the calibration_dict.p saved in get_analysis_location()
    Returns:
        calibration_dict
    """
    path = get_analysis_location()
    file_name = 'analysis_dict.p'
    calibration_dict = pickle.load(open(path + file_name, "rb"))
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
    update_calibration_dict({'current_qubit': index})
    print('current_qubit set to {}'.format(index))
    validate_calibration_dictionary(vna, alazar)


def update_calibration_val(key, val, qubit_index=None):
    c_dict = get_calibration_dict()
    if qubit_index is None:
        c_dict[key][c_dict['current_qubit']] = val
    else:
        c_dict[key][qubit_index] = val
    update_calibration_dict(c_dict)


def get_calibration_val(key):
    c_dict = get_calibration_dict()
    return c_dict[key][c_dict['current_qubit']]


def recalculate_g(dec_chans=None):
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
    print('expeced qubit freq: {}\n (from g of {}, push on resonator {})'
          'actual qubit freq: {}\n (for same push gives g of {}'.format(
              expected, old_g, res_data[2], actual, new_g))
    return new_g


def validate_calibration_dictionary(vna=True, alazar=True):
    c_dict = get_calibration_dict()
    qubit_num = get_qubit_count()
    qubit_length_list = np.zeros(get_qubit_count())
    missing_keys = []
    wrong_length_keys = []
    if vna:
        missing_keys.extend([k for k in vna_dict_keys if k not in c_dict])
        wrong_length_keys.extend(
            [k for k in vna_dict_keys if len(c_dict[k]) != qubit_num])
    if alazar:
        missing_keys.extend([k for k in alazar_dict_keys if k not in c_dict])
        wrong_length_keys.extend(
            [k for k in alazar_dict_keys if len(c_dict[k]) != qubit_num])
    for k in missing_keys:
        c_dict[k] = qubit_length_list
    for k in wrong_length_keys:
        c_dict[k] = qubit_length_list
    print('{} added to calibration_dictionary\n'
          '{} were reset to correct length'.format(missing_keys,
                                                   wrong_length_keys))
