# import numpy as np
import pickle
from . import get_latest_counter, \
    get_calibration_val, get_pulse_location

# TODO: rount -> int/ciel
# TODO: test segments
# TODO: make_save_send_load_awg_file -> not doing same thing twice!
# TODO: docstrings
# TODO: checks
# TODO: floquet

#################################################################
# General AWG and Sequence functions
#################################################################


def get_current_seq(awg):
    seq_name = awg.current_seq()
    path = get_pulse_location()
    seq = pickle.load(open(path + seq_name, "rb"))
    return seq


def make_save_send_load_awg_file(awg, sequence, file_name):
    """
    WYSIYWYG

    Args:
        awg instrument for upload
        sequence to be uploaded
    """
    awg.make_and_save_awg_file(*sequence.unwrap(), filename=file_name)
    awg.make_send_and_load_awg_file(*sequence.unwrap())


def check_sample_rate(awg):
    """
    Checks sample rate in pulse dict against that on awg

    Args:
        awg instrument for checking
    """
    sr = get_calibration_val('sample_rate')
    if sr != awg.clock_freq():
        awg.clock_freq(sr)
    print('awg clock freq set to {}'.format(sr))


def check_seq_uploaded(awg, seq_type, dict_to_check,
                       start=None, stop=None, step=None):
    uploaded_seq = get_current_seq(awg)
    if uploaded_seq.labels['seq_type'] is not seq_type:
        return False
    for k in dict_to_check:
        if k not in uploaded_seq.labels:
            return False
        elif uploaded_seq.labels[k] != dict_to_check[k]:
            return False
    try:
        if ([start, stop, step] !=
                [uploaded_seq.start, uploaded_seq.stop, uploaded_seq.step]):
            return False
    except AttributeError:
        return False
    return True
