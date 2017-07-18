import numpy as np
import pickle
import os
from . import get_calibration_dict, get_allowed_keys, gaussian_array, \
    gaussian_derivative_array, flat_array, cos_gaussian_array, \
    sin_gaussian_array, cos_array, sin_array, cos_multi_array, sin_multi_array

from . import Segment, Waveform, Element, Sequence


def make_readout_wf(first_in_seq=False, channel=4):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    measurement_segment = Segment(
        name='cavity_measurement',
        gen_func=flat_array,
        func_args={'amp': 1, 'dur': calib_dict['readout_time'][qubit]},
        time_markers={
            1: {'delay_time': [calib_dict['marker_readout_delay'][qubit]],
                'duration_time': [calib_dict['marker_time'][qubit]]}})

    time_before_readout = (calib_dict['pulse_end'][qubit] +
                           calib_dict['pulse_readout_delay'][qubit])
    time_after_readout = (calib_dict['cycle_time'][qubit] -
                          calib_dict['pulse_end'][qubit] -
                          calib_dict['pulse_readout_delay'][qubit] -
                          calib_dict['readout_time'][qubit])
    wait1_segment = Segment(
        name='wait_before_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_before_readout})
    wait2_segment = Segment(
        name='wait_after_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_readout})
    readout_wf = Waveform(
        channel=channel,
        segment_list=[wait1_segment, measurement_segment, wait2_segment],
        sample_rate=calib_dict['sample_rate'][qubit])
    if first_in_seq:
        readout_wf.add_marker(
            2, 0, int(calib_dict['marker_time'][qubit] *
                      calib_dict['sample_rate'][qubit]))
    return readout_wf


def make_readout_ssb_wf_I(freq_list, first_in_seq=False, channel=3):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    measurement_segment = Segment(
        name='cavity_measurement_i',
        gen_func=cos_multi_array,
        func_args={'amp': 1, 'dur': calib_dict['readout_time'][qubit],
                   'freq_list':freq_list},
        time_markers={
            1: {'delay_time': [calib_dict['marker_readout_delay'][qubit]],
                'duration_time': [calib_dict['marker_time'][qubit]]}})
    time_before_readout = (calib_dict['pulse_end'][qubit] +
                           calib_dict['pulse_readout_delay'][qubit])
    time_after_readout = (calib_dict['cycle_time'][qubit] -
                          calib_dict['pulse_end'][qubit] -
                          calib_dict['pulse_readout_delay'][qubit] -
                          calib_dict['readout_time'][qubit])
    wait1_segment = Segment(
        name='wait_before_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_before_readout})
    wait2_segment = Segment(
        name='wait_after_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_readout})
    readout_wf = Waveform(
        channel=channel,
        segment_list=[wait1_segment, measurement_segment, wait2_segment],
        sample_rate=calib_dict['sample_rate'][qubit])
    if first_in_seq:
        readout_wf.add_marker(
            2, 0, int(calib_dict['marker_time'][qubit] *
                      calib_dict['sample_rate'][qubit]))
    return readout_wf


def make_readout_ssb_wf_Q(freq_list, channel=4):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    measurement_segment = Segment(
        name='cavity_measurement_q',
        gen_func=sin_multi_array,
        func_args={'amp': 1, 'dur': calib_dict['readout_time'][qubit],
                   'freq_list':freq_list, 'positive': False})
    time_before_readout = (calib_dict['pulse_end'][qubit] +
                           calib_dict['pulse_readout_delay'][qubit])
    time_after_readout = (calib_dict['cycle_time'][qubit] -
                          calib_dict['pulse_end'][qubit] -
                          calib_dict['pulse_readout_delay'][qubit] -
                          calib_dict['readout_time'][qubit])
    wait1_segment = Segment(
        name='wait_before_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_before_readout})
    wait2_segment = Segment(
        name='wait_after_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_readout})
    readout_wf = Waveform(
        channel=channel,
        segment_list=[wait1_segment, measurement_segment, wait2_segment],
        sample_rate=calib_dict['sample_rate'][qubit])
    return readout_wf