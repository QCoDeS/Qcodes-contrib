import numpy as np
import pickle
import os
from . import get_calibration_dict, get_allowed_keys, gaussian_array, \
    gaussian_derivative_array, flat_array, cos_gaussian_array, \
    sin_gaussian_array, cos_array, sin_array

from . import Segment, Waveform, Element, Sequence

# TODO: marker on gate list sequence
# TODO: refactor sequence making functions to have less than a million args
# TODO: square pulses
# TODO: docstrings


def save_sequence(sequence, file_name):
    if os.path.exists(file_name):
        raise Exception('File already exists at this location with this '
                        'name: {}'.format(file_name))
    else:
        pickle.dump(sequence, open(file_name, 'wb'))


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


####################################################################
# Sequence building functions (vary param over sequence)
####################################################################

def make_varying_sequence(element_template, vary_ch, vary_seg,
                          arg_to_vary, start, stop, step, name=None,
                          variable_name=None,
                          variable_label=None, variable_unit=None,
                          readout_ch=4, marker_points=100):
    var_name = variable_name or arg_to_vary
    seq_name = name or arg_to_vary + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=var_name, start=start, stop=stop, step=step,
        variable_label=variable_label, variable_unit=variable_unit)
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch].segment_list[vary_seg].func_args[arg_to_vary] = val
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


def make_time_varying_sequence(element_template, vary_ch, vary_seg,
                               arg_to_vary, start, stop, step,
                               compensate_seg, total_time, name=None,
                               variable_name=None,
                               variable_label=None, variable_unit=None,
                               readout_ch=4,
                               marker_points=100):
    var_name = variable_name or arg_to_vary
    seq_name = name or arg_to_vary + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=var_name, start=start, stop=stop, step=step,
        variable_label=variable_label, variable_unit=variable_unit)
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch].segment_list[vary_seg].func_args[arg_to_vary] = val
        c_dur = total_time - elem[vary_ch].duration
        elem[vary_ch].segment_list[compensate_seg].func_args["dur"] = c_dur
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


def make_multi_varying_sequence(element_template,
                                vary_ch_1, vary_seg_1, arg_to_vary_1,
                                start1, stop1, step1,
                                vary_ch_2, vary_seg_2, arg_to_vary_2,
                                start2, stop2, step2,
                                name=None, variable_name=None,
                                variable_label=None,
                                variable_unit=None,
                                readout_ch=4, marker_points=100):
    var_name = variable_name or arg_to_vary_1 + '_' + arg_to_vary_2
    seq_name = name or variable_name + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=var_name, start=start1, stop=stop1,
        step=step1, variable_label=variable_label, variable_unit=variable_unit)
    variable_array_2 = np.linspace(
        start2, stop2, num=round(abs(stop2 - start2) / step2 + 1))
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch_1].segment_list[vary_seg_1].func_args[
            arg_to_vary_1] = val
        elem[vary_ch_2].segment_list[vary_seg_2].func_args[
            arg_to_vary_2] = variable_array_2[i]
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


def make_time_multi_varying_sequence(element_template,
                                     vary_ch_1, vary_seg_1, arg_to_vary_1,
                                     start1, stop1, step1,
                                     vary_ch_2, vary_seg_2, arg_to_vary_2,
                                     start2, stop2, step2,
                                     compensate_seg_1, compensate_seg_2,
                                     total_time, name=None,
                                     variable_name=None,
                                     variable_label=None, variable_unit=None,
                                     readout_ch=4, marker_points=100):
    var_name = variable_name or arg_to_vary_1 + '_' + arg_to_vary_2
    seq_name = name or variable_name + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=var_name, start=start1, stop=stop1,
        step=step1, variable_label=variable_label, variable_unit=variable_unit)
    variable_array_2 = np.linspace(
        start2, stop2, num=round(abs(stop2 - start2) / step2 + 1))
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch_1].segment_list[vary_seg_1].func_args[
            arg_to_vary_1] = val
        c1_dur = total_time - elem[vary_ch_1].duration
        elem[vary_ch_1].segment_list[
            compensate_seg_1].func_args["dur"] = c1_dur
        elem[vary_ch_2].segment_list[vary_seg_2].func_args[
            arg_to_vary_2] = variable_array_2[i]
        c2_dur = total_time - elem[vary_ch_2].duration
        elem[vary_ch_2].segment_list[
            compensate_seg_2].func_args["dur"] = c2_dur
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


#################################################################
# Element building functions (for use executing gates)
#################################################################

def make_x_y_carrier_gaussian_pulses(params, drag=False):
    pulse_dict = {}
    x_y_wait = Segment(
        name='XY_wait', gen_func=flat_array,
        func_args={
            'amp': 0,
            'dur': 2 * params['pi_pulse_sigma'] * params['sigma_cutoff']})

    x_y_pi = Segment(
        name='pi', gen_func=gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'amp': params['pi_pulse_amp']})

    x_y_pi_half = Segment(
        name='pi/2', gen_func=gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'amp': params['pi_half_pulse_amp']})

    x_y_pi_half_neg = Segment(
        name='-pi/2', gen_func=gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})

    pulse_dict['XY_wait'] = x_y_wait
    pulse_dict['X_I'] = x_y_pi
    pulse_dict['X/2_I'] = x_y_pi_half
    pulse_dict['-X/2_I'] = x_y_pi_half_neg
    pulse_dict['Y_Q'] = x_y_pi
    pulse_dict['Y/2_Q'] = x_y_pi_half
    pulse_dict['-Y/2_Q'] = x_y_pi_half_neg
    if drag:
        x_y_pi_drag = Segment(
            name='pi_drag',
            gen_func=gaussian_derivative_array,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'amp': (params['pi_pulse_amp'] *
                               params['drag_coef'])})

        x_y_pi_half_drag = Segment(
            name='pi/2_drag',
            gen_func=gaussian_derivative_array,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'amp': (params['pi_half_pulse_amp'] *
                               params['drag_coef'])})

        x_y_pi_half_neg_drag = Segment(
            name='-pi/2_drag',
            gen_func=gaussian_derivative_array,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'amp': (params['pi_half_pulse_amp'] *
                               params['drag_coef']),
                       'positive': False})

        pulse_dict['X_Q'] = x_y_pi_drag
        pulse_dict['X/2_Q'] = x_y_pi_half_drag
        pulse_dict['-X/2_Q'] = x_y_pi_half_neg_drag
        pulse_dict['Y_I'] = x_y_pi_drag
        pulse_dict['Y/2_I'] = x_y_pi_half_drag
        pulse_dict['-Y/2_I'] = x_y_pi_half_neg_drag
    else:
        pulse_dict['X_Q'] = x_y_wait
        pulse_dict['X/2_Q'] = x_y_wait
        pulse_dict['-X/2_Q'] = x_y_wait
        pulse_dict['Y_I'] = x_y_wait
        pulse_dict['Y/2_I'] = x_y_wait
        pulse_dict['-Y/2_I'] = x_y_wait
    return pulse_dict


def make_x_y_carrier_flat_pulses(params):
    pulse_dict = {}
    x_y_wait = Segment(
        name='XY_wait', gen_func=flat_array,
        func_args={
            'amp': 0,
            'dur': 2 * params['pi_pulse_sigma'] * params['sigma_cutoff']})

    x_y_pi = Segment(
        name='pi', gen_func=flat_array,
        func_args={'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp']})

    x_y_pi_half = Segment(
        name='pi/2', gen_func=flat_array,
        func_args={'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp']})

    x_y_pi_half_neg = Segment(
        name='-pi/2', gen_func=flat_array,
        func_args={'dur': params['pi_pulse_dur'],
                   'amp': -1 * params['pi_half_pulse_amp']})

    pulse_dict['XY_wait'] = x_y_wait
    pulse_dict['X_I'] = x_y_pi
    pulse_dict['X_Q'] = x_y_wait
    pulse_dict['X/2_I'] = x_y_pi_half
    pulse_dict['X/2_Q'] = x_y_wait
    pulse_dict['-X/2_I'] = x_y_pi_half_neg
    pulse_dict['-X/2_Q'] = x_y_wait
    pulse_dict['Y_I'] = x_y_wait
    pulse_dict['Y_Q'] = x_y_pi
    pulse_dict['Y/2_I'] = x_y_wait
    pulse_dict['Y/2_Q'] = x_y_pi_half
    pulse_dict['-Y/2_I'] = x_y_wait
    pulse_dict['-Y/2_Q'] = x_y_pi_half_neg
    return pulse_dict


def make_x_y_ssb_gaussian_pulses(params, SSBfreq):
    pulse_dict = {}
    x_y_wait = Segment(
        name='XY_wait', gen_func=flat_array,
        func_args={
            'amp': 0,
            'dur': 2 * params['pi_pulse_sigma'] * params['sigma_cutoff']})
    ssb_pi_x_I = Segment(
        name='X_I', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_pulse_amp']})
    ssb_pi_x_Q = Segment(
        name='X_Q', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_pulse_amp'],
                   'positive': False})
    ssb_pi_half_x_I = Segment(
        name='X/2_I', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_half_x_Q = Segment(
        name='X/2_Q', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_x_I = Segment(
        name='-X/2_I', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_x_Q = Segment(
        name='-X/2_Q', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_y_I = Segment(
        name='Y_I', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_pulse_amp']})
    ssb_pi_y_Q = Segment(
        name='Y_Q', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_pulse_amp']})
    ssb_pi_half_y_I = Segment(
        name='Y/2_I', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_half_y_Q = Segment(
        name='Y/2_Q', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_half_neg_y_I = Segment(
        name='-Y/2_I', gen_func=sin_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_y_Q = Segment(
        name='-Y/2_Q', gen_func=cos_gaussian_array,
        func_args={'sigma': params['pi_pulse_sigma'],
                   'sigma_cutoff': params['sigma_cutoff'],
                   'SSBfreq': SSBfreq,
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    pulse_dict['XY_wait'] = x_y_wait
    pulse_dict['X_I'] = ssb_pi_x_I
    pulse_dict['X_Q'] = ssb_pi_x_Q
    pulse_dict['X/2_I'] = ssb_pi_half_x_I
    pulse_dict['X/2_Q'] = ssb_pi_half_x_Q
    pulse_dict['-X/2_I'] = ssb_pi_half_neg_x_I
    pulse_dict['-X/2_Q'] = ssb_pi_half_neg_x_Q
    pulse_dict['Y_I'] = ssb_pi_y_I
    pulse_dict['Y_Q'] = ssb_pi_y_Q
    pulse_dict['Y/2_I'] = ssb_pi_half_y_I
    pulse_dict['Y/2_Q'] = ssb_pi_half_y_Q
    pulse_dict['-Y/2_I'] = ssb_pi_half_neg_y_I
    pulse_dict['-Y/2_Q'] = ssb_pi_half_neg_y_Q

    return pulse_dict


def make_x_y_ssb_flat_pulses(params, SSBfreq):
    pulse_dict = {}
    x_y_wait = Segment(
        name='XY_wait', gen_func=flat_array,
        func_args={
            'amp': 0,
            'dur': 2 * params['pi_pulse_sigma'] * params['sigma_cutoff']})
    ssb_pi_x_I = Segment(
        name='X_I', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp']})
    ssb_pi_x_Q = Segment(
        name='X_Q', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp'],
                   'positive': False})
    ssb_pi_half_x_I = Segment(
        name='X/2_I', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_half_x_Q = Segment(
        name='X/2_Q', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_x_I = Segment(
        name='-X/2_I', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_x_Q = Segment(
        name='-X/2_Q', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_y_I = Segment(
        name='Y_I', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp']})
    ssb_pi_y_Q = Segment(
        name='Y_Q', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp']})
    ssb_pi_half_y_I = Segment(
        name='Y/2_I', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp']})
    ssb_pi_half_y_Q = Segment(
        name='Y/2_Q', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp']})
    ssb_pi_half_neg_y_I = Segment(
        name='-Y/2_I', gen_func=sin_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_pulse_amp'],
                   'positive': False})
    ssb_pi_half_neg_y_Q = Segment(
        name='-Y/2_Q', gen_func=cos_array,
        func_args={'freq': SSBfreq,
                   'dur': params['pi_pulse_dur'],
                   'amp': params['pi_half_pulse_amp'],
                   'positive': False})
    pulse_dict['XY_wait'] = x_y_wait
    pulse_dict['X_I'] = ssb_pi_x_I
    pulse_dict['X_Q'] = ssb_pi_x_Q
    pulse_dict['X/2_I'] = ssb_pi_half_x_I
    pulse_dict['X/2_Q'] = ssb_pi_half_x_Q
    pulse_dict['-X/2_I'] = ssb_pi_half_neg_x_I
    pulse_dict['-X/2_Q'] = ssb_pi_half_neg_x_Q
    pulse_dict['Y_I'] = ssb_pi_y_I
    pulse_dict['Y_Q'] = ssb_pi_y_Q
    pulse_dict['Y/2_I'] = ssb_pi_half_y_I
    pulse_dict['Y/2_Q'] = ssb_pi_half_y_Q
    pulse_dict['-Y/2_I'] = ssb_pi_half_neg_y_I
    pulse_dict['-Y/2_Q'] = ssb_pi_half_neg_y_Q

    return pulse_dict


def make_z_pulses(params):
    pulse_dict = {}
    z_pi = Segment(
        name='Z', gen_func=flat_array,
        func_args={'dur': params['z_pulse_dur'],
                   'amp': params['z_pulse_amp']})

    z_pi_half = Segment(
        name='Z/2', gen_func=flat_array,
        func_args={'dur': params['z_pulse_dur'],
                   'amp': params['z_half_pulse_amp']})
    z_pi_half_neg = Segment(
        name='-Z/2', gen_func=flat_array,
        func_args={'dur': params['z_pulse_dur'],
                   'amp': -1 * params['z_half_pulse_amp']})
    z_wait = Segment(
        name='Z_wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': params['z_pulse_dur']})
    pulse_dict['Z'] = z_pi
    pulse_dict['Z/2'] = z_pi_half
    pulse_dict['-Z/2'] = z_pi_half_neg
    pulse_dict['Z_wait'] = z_wait

    return pulse_dict


def make_pulse_dict(pulse_params={}, qubit_num=None,
                    SSBfreq=None, gaussian=True,
                    drag=False, z_gates=False):
    """
    Function which returns a dictionary of pulses for gates based on the
    current calibration dictionary values (or optionally those specified in a
    dictionary given).

    Args:
        pulse_params (dict) (default {}): values for the pulse params to
            override those in the calibration dictionary
        SSBfreq (float) (default None): sideband frequency if the carrier
            is not being used to drive the qubit. Lower sideband by default

    Returns:
        pulse_dict (dict): dictionary containing segments of pulses
    """
    if drag and SSBfreq is not None:
        raise Exception('Cannot use drag pulse and ssb pulses.')
    if drag and not gaussian:
        raise Exception('Cannot use drag and square pulses.')

    calib_dict = get_calibration_dict()
    qubit = qubit_num or calib_dict['current_qubit']
    params = {}
    pulse_keys = get_allowed_keys(vna=False, alazar=False, pulse=True)
    for k in pulse_keys:
        if k in pulse_params:
            params[k] = pulse_params[k]
        else:
            params[k] = calib_dict[k][qubit]

    pulse_dict = {}

    measurement = Segment(
        name='cavity_measurement', gen_func=flat_array,
        func_args={'amp': params['readout_amp'],
                   'dur': params['readout_time']},
        time_markers={1: {'delay_time': [-1 * params['marker_readout_delay']],
                          'duration_time': [params['marker_time']]}})

    measurement_wait = Segment(
        name='wait_measurement', gen_func=flat_array,
        func_args={'amp': 0, 'dur': params['readout_time']})

    identity = Segment(
        name='identity', gen_func=flat_array,
        func_args={
            'amp': 0,
            'dur': 2 * params['pi_pulse_sigma'] * params['sigma_cutoff']})

    wait = Segment(name='wait', gen_func=flat_array, func_args={'amp': 0})

    pulse_dict['measurement'] = measurement
    pulse_dict['measurement_wait'] = measurement_wait
    pulse_dict['wait_template'] = wait
    pulse_dict['identity'] = identity

    if SSBfreq is not None and gaussian:
        pulses = make_x_y_ssb_gaussian_pulses(
            params, SSBfreq)
    elif SSBfreq is not None:
        pulses = make_x_y_ssb_flat_pulses(
            params, SSBfreq)
    elif gaussian:
        pulses = make_x_y_carrier_gaussian_pulses(
            params, drag=drag)
    else:
        pulses = make_x_y_carrier_flat_pulses(
            params)
    pulse_dict.update(pulses)

    if z_gates:
        z_pulses = make_z_pulses(params)
        pulse_dict.update(z_pulses)

    return pulse_dict


def do_x_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['X_I']
    q_pulse = pulse_dict['X_Q']
    identity = pulse_dict['XY_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_x_pi_half(element, pulse_dict, channels=[1, 2, 3, 4], positive=True):
    identity = pulse_dict['XY_wait']
    if positive:
        i_pulse = pulse_dict['X/2_I']
        q_pulse = pulse_dict['X/2_Q']
    else:
        i_pulse = pulse_dict['-X/2_I']
        q_pulse = pulse_dict['-X/2_Q']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_y_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['Y_I']
    q_pulse = pulse_dict['Y_Q']
    identity = pulse_dict['XY_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_y_pi_half(element, pulse_dict, channels=[1, 2, 3, 4], positive=True):
    identity = pulse_dict['XY_wait']
    if positive:
        i_pulse = pulse_dict['Y/2_I']
        q_pulse = pulse_dict['Y/2_Q']
    else:
        i_pulse = pulse_dict['-Y/2_I']
        q_pulse = pulse_dict['-Y/2_Q']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_z_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['Z_wait']
    z_pulse = pulse_dict['Z']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(z_pulse)
    element[channels[3]].add_segment(identity)


def do_z_pi_half(element, pulse_dict, channels=[1, 2, 3, 4], positive=True):
    identity = pulse_dict['Z_wait']
    if positive:
        z_pulse = pulse_dict['Z/2']
    else:
        z_pulse = pulse_dict['-Z/2']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(z_pulse)
    element[channels[3]].add_segment(identity)


def do_identity(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['XY_wait']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def measure(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['measurement_wait']
    measurement = pulse_dict['measurement']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(measurement)


def wait(element, pulse_dict, dur, channels=[1, 2, 3, 4]):
    identity = pulse_dict['wait_template'].copy()
    identity.func_args['dur'] = dur
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def prepend_compensating_wait_to_element(element, pulse_dict, total_time):
    identity = pulse_dict['wait_template'].copy()
    identity.func_args['dur'] = total_time - element.duration
    for w in element:
        element[w].add_segment(identity, position=0)


def execute_gates(element, pulse_dict, gate_list,
                  channels=[1, 2, 3, 4], spacing=None):
    for i in gate_list:
        if i == 'I':
            do_identity(element, pulse_dict, channels=channels)
        if i == 'X':
            do_x_pi(element, pulse_dict, channels=channels)
        elif i == 'X/2':
            do_x_pi_half(element, pulse_dict, channels=channels)
        elif i == '-X/2':
            do_x_pi_half(element, pulse_dict, channels=channels,
                         positive=False)
        elif i == 'Y':
            do_y_pi(element, pulse_dict, channels=channels)
        elif i == 'Y/2':
            do_y_pi_half(element, pulse_dict, channels=channels)
        elif i == '-Y/2':
            do_y_pi_half(element, pulse_dict, channels=channels,
                         positive=False)
        elif i == 'Z':
            do_z_pi(element, pulse_dict, channels=channels)
        elif i == 'Z/2':
            do_z_pi_half(element, pulse_dict, channels=channels)
        if spacing is not None:
            wait(element, pulse_dict, spacing, channels=channels)
        elif i == '-Z/2':
            do_z_pi_half(element, pulse_dict, channels=channels,
                         positive=False)


###########################################################
# Sequence building function (with gates)
###########################################################

def make_element_from_gate_list(
        gate_list, SSBfreq=None, drag=False, gaussian=True,
        channels=[1, 2, 3, 4], spacing=None, calib_dict=None,
        pulse_dict=None):
    calib_dict = calib_dict or get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pulse_dict = pulse_dict or make_pulse_dict(
        SSBfreq=SSBfreq, drag=drag, gaussian=gaussian)
    element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    i_wf = Waveform(channel=channels[0])
    q_wf = Waveform(channel=channels[1])
    z_wf = Waveform(channel=channels[2])
    measure_wf = Waveform(channel=channels[3])
    element.add_waveform(i_wf)
    element.add_waveform(q_wf)
    element.add_waveform(z_wf)
    element.add_waveform(measure_wf)
    execute_gates(element, pulse_dict, gate_list,
                  spacing=spacing)
    wait(element, pulse_dict, calib_dict['pulse_readout_delay'][qubit],
         channels=channels)
    measure(element, pulse_dict, channels=channels)
    time_after_readout = (calib_dict['cycle_time'][qubit] -
                          calib_dict['pulse_end'][qubit] -
                          calib_dict['pulse_readout_delay'][qubit] -
                          calib_dict['readout_time'][qubit])
    wait(element, pulse_dict, time_after_readout, channels=channels)
    prepend_compensating_wait_to_element(element, pulse_dict,
                                         calib_dict['cycle_time'][qubit])
    return element


def make_sequence_from_gate_lists(
        gate_list, SSBfreq=None, drag=False, gaussian=True,
        channels=[1, 2, 3, 4], name=None, variable_label=None, spacing=None):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pulse_dict = make_pulse_dict(SSBfreq=SSBfreq, drag=drag, gaussian=gaussian)
    seq = Sequence(name=name or 'seq_from_gates',
                   variable_label=variable_label)

    for i, gate_list in enumerate(gate_list):
        element = make_element_from_gate_list(
            gate_list, SSBfreq=SSBfreq, drag=drag, gaussian=gaussian,
            channels=channels, spacing=spacing, pulse_dict=pulse_dict,
            calib_dict=calib_dict)
        if i == 0:
            marker_points = int(calib_dict['marker_time'][qubit] *
                                calib_dict['sample_rate'][qubit])
            element[channels[3]].add_marker(2, 0, marker_points)
        seq.add_element(element)
    return seq
