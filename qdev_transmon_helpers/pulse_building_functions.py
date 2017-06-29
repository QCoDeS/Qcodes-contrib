import numpy as np
from . import get_calibration_dict, get_allowed_keys, gaussian, \
    gaussian_derivative, flat, cos_gaussian, sin_gaussian

from . import Segment, Waveform, Element, Sequence

# TODO: marker on gate list sequence
# TODO: refactor sequence making functions to have less than a million args
# TODO: square pulses
# TODO: docstrings


def make_readout_wf(first_in_seq=False, channel=4):
    p_dict = get_calibration_dict()
    measurement_segment = Segment(
        name='cavity_measurement',
        gen_fuvnc=flat,
        func_args={'amp': 1, 'dur': p_dict['readout_time']},
        time_markers={1: {'delay_time': [p_dict['marker_readout_delay']],
                          'duration_time': [p_dict['marker_time']]}})

    time_before_readout = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    time_after_readout = (p_dict['cycle_time'] -
                          p_dict['pulse_end'] -
                          p_dict['pulse_readout_delay'] -
                          p_dict['readout_time'])
    wait1_segment = Segment(name='wait_before_readout', gen_func=flat,
                            func_args={'amp': 0, 'dur': time_before_readout})
    wait2_segment = Segment(name='wait_after_readout', gen_func=flat,
                            func_args={'amp': 0, 'dur': time_after_readout})
    readout_wf = Waveform(
        channel=channel,
        segment_list=[wait1_segment, measurement_segment, wait2_segment],
        sample_rate=p_dict['sample_rate'])
    if first_in_seq:
        readout_wf.add_marker(
            2, 0, int(p_dict['marker_time'] * p_dict['sample_rate']))
    return readout_wf


####################################################################
# Sequence building functions (vary param over sequence)
####################################################################

def make_varying_sequence(element_template, vary_ch, vary_seg,
                          arg_to_vary, start, stop, step, name=None,
                          variable_label=None, variable_unit=None,
                          readout_ch=4, marker_points=100):
    seq_name = name or arg_to_vary + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=arg_to_vary, start=start, stop=stop, step=step,
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
                               variable_label=None, variable_unit=None,
                               readout_ch=4,
                               marker_points=100):
    seq_name = name or arg_to_vary + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=arg_to_vary, start=start, stop=stop, step=step,
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
                                name=None, variable_label=None,
                                variable_unit=None,
                                readout_ch=4, marker_points=100):
    variable_name = arg_to_vary_1 + '_' + arg_to_vary_2
    seq_name = name or variable_name + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=variable_name, start=start1, stop=stop1,
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
                                     variable_label=None, variable_unit=None,
                                     readout_ch=4, marker_points=100):
    variable_name = arg_to_vary_1 + '_' + arg_to_vary_2
    seq_name = name or variable_name + '_varying_seq'
    sequence = Sequence(
        name=seq_name, variable=variable_name, start=start1, stop=stop1,
        step=step1, variable_label=variable_label, variable_unit=variable_unit)
    variable_array_2 = np.linspace(
        start2, stop2, num=round(abs(stop2 - start2) / step2 + 1))
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch_1].segment_list[vary_seg_1].func_args[
            arg_to_vary_1] = val
        elem[vary_ch_1].segment_list[compensate_seg_1].func_args[
            "dur"] = total_time - elem[vary_ch_1].duration
        elem[vary_ch_2].segment_list[vary_seg_2].func_args[
            arg_to_vary_2] = variable_array_2[i]
        elem[vary_ch_2].segment_list[compensate_seg_2].func_args[
            "dur"] = total_time - elem[vary_ch_2].duration
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


#################################################################
# Element building functions (for use executing gates)
#################################################################


def make_pulse_dict(pulse_params={}, SSBfreq=None, square=False,
                    drag=False):
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

    calib_dict = get_calibration_dict()
    params = {}
    pulse_keys = get_allowed_keys(vna=False, alazar=False, pulse=True)
    for k in pulse_keys:
        params[k] = pulse_params[k] if k in pulse_params else calib_dict[k]

    pulse_dict = {}

    measurement = Segment(
        name='cavity_measurement', gen_func=flat,
        func_args={'amp': params['readout_amp'],
                   'dur': params['readout_time']},
        time_markers={1: {'delay_time': [-1 * params['marker_readout_delay']],
                          'duration_time': [params['marker_time']]}})

    measurement_wait = measurement.copy()
    measurement_wait.func_args['amp'] = 0
    measurement_wait.name = 'wait_for_measurement'

    identity_template = Segment(name='wait', gen_func=flat,
                                func_args={'amp': 0})

    identity = identity_template.copy()
    identity = 'identity'
    identity.func_args['dur'] = 2 * params['pi_pulse_sigma'] * \
        params['sigma_cutoff']

    pulse_dict['measurement'] = measurement
    pulse_dict['measurement_wait'] = measurement_wait
    pulse_dict['wait_template'] = identity_template
    pulse_dict['identity'] = identity

    x_y_wait_pulse = identity_template.copy()
    x_y_wait_pulse.name = 'wait_for_x_y_pulse'
    x_y_wait_pulse.func_args['dur'] = params['pi_pulse_sigma'] * \
        params['sigma_cutoff']

    pulse_dict['x_y_wait'] = x_y_wait_pulse

    if SSBfreq is None:
        x_y_pi_pulse = Segment(
            name='pi_pulse', gen_func=gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'amp': params['pi_pulse_amp']})

        x_y_pi_half_pulse = Segment(
            name='pi_half_pulse', gen_func=gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'amp': params['pi_half_pulse_amp']})
        pulse_dict['pi_x_I'] = x_y_pi_pulse
        pulse_dict['pi_x_Q'] = x_y_wait_pulse
        pulse_dict['pi_half_I'] = x_y_pi_half_pulse
        pulse_dict['pi_x_half_Q'] = x_y_wait_pulse
        pulse_dict['pi_y_I'] = x_y_wait_pulse
        pulse_dict['pi_y_Q'] = x_y_pi_pulse
        pulse_dict['pi_y_half_I'] = x_y_wait_pulse
        pulse_dict['pi_half_y_Q'] = x_y_pi_half_pulse
        if drag:
            x_y_pi_pulse_drag = Segment(
                name='pi_pulse_drag',
                gen_func=gaussian_derivative,
                func_args={'sigma': params['pi_pulse_sigma'],
                           'sigma_cutoff': params['sigma_cutoff'],
                           'amp': (params['pi_pulse_amp'] *
                                   params['drag_coef'])})

            x_y_pi_half_pulse_drag = Segment(
                name='pi_half_pulse_drag',
                gen_func=gaussian_derivative,
                func_args={'sigma': params['pi_pulse_sigma'],
                           'sigma_cutoff': params['sigma_cutoff'],
                           'amp': (params['pi_half_pulse_amp'] *
                                   params['drag_coef'])})
            pulse_dict['pi_drag'] = x_y_pi_pulse_drag
            pulse_dict['pi_half_drag'] = x_y_pi_half_pulse_drag
    else:
        ssb_pi_pulse_x_I = Segment(
            name='pi_pulse_x_I', gen_func=cos_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_pulse_amp']})
        ssb_pi_pulse_x_Q = Segment(
            name='pi_pulse_x_Q', gen_func=sin_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_pulse_amp'],
                       'positive': False})
        ssb_pi_half_pulse_x_I = Segment(
            name='pi_half_pulse_x_I', gen_func=cos_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_half_pulse_amp']})
        ssb_pi_half_pulse_x_Q = Segment(
            name='pi_pulse_x_Q', gen_func=sin_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_half_pulse_amp'],
                       'positive': False})
        ssb_pi_pulse_y_I = Segment(
            name='pi_pulse_y_I', gen_func=sin_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_pulse_amp']})
        ssb_pi_pulse_y_Q = Segment(
            name='pi_pulse_y_Q', gen_func=cos_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_pulse_amp']})
        ssb_pi_half_pulse_y_I = Segment(
            name='pi_half_pulse_y_I', gen_func=sin_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_half_pulse_amp']})
        ssb_pi_half_pulse_y_Q = Segment(
            name='pi_half_pulse_y_Q', gen_func=cos_gaussian,
            func_args={'sigma': params['pi_pulse_sigma'],
                       'sigma_cutoff': params['sigma_cutoff'],
                       'SSBfreq': SSBfreq,
                       'amp': params['pi_half_pulse_amp']})
        pulse_dict['pi_x_I'] = ssb_pi_pulse_x_I
        pulse_dict['pi_x_Q'] = ssb_pi_pulse_x_Q
        pulse_dict['pi_half_x_I'] = ssb_pi_half_pulse_x_I
        pulse_dict['pi_half_x_Q'] = ssb_pi_half_pulse_x_Q
        pulse_dict['pi_y_I'] = ssb_pi_pulse_y_I
        pulse_dict['pi_y_Q'] = ssb_pi_pulse_y_Q
        pulse_dict['pi_half_x_I'] = ssb_pi_half_pulse_y_I
        pulse_dict['pi_half_x_Q'] = ssb_pi_half_pulse_y_Q

    z_pi_pulse = Segment(name='pi_pulse_z', gen_func=flat,
                         func_args={'dur': params['z_pulse_dur'],
                                    'amp': params['z_pulse_amp']})

    z_pi_half_pulse = Segment(name='pi_half_pulse_z', gen_func=flat,
                              func_args={'dur': params['z_pulse_dur'],
                                         'amp': params['z_half_pulse_amp']})
    z_wait_pulse = identity_template.copy()
    z_wait_pulse.name = 'wait_for_z_pulse'
    z_wait_pulse.func_args['dur'] = params['z_pulse_dur']

    pulse_dict['pi_z'] = z_pi_pulse
    pulse_dict['pi_half_z'] = z_pi_half_pulse
    pulse_dict['z_wait'] = z_wait_pulse

    return pulse_dict


def do_x_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['pi_x_I']
    q_pulse = pulse_dict['pi_x_Q']
    identity = pulse_dict['x_y_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_x_pi_half(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['pi_half_x_I']
    q_pulse = pulse_dict['pi_half_x_Q']
    identity = pulse_dict['x_y_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_y_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['pi_y_I']
    q_pulse = pulse_dict['pi_y_Q']
    identity = pulse_dict['x_y_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_y_pi_half(element, pulse_dict, channels=[1, 2, 3, 4]):
    i_pulse = pulse_dict['pi_half_y_I']
    q_pulse = pulse_dict['pi_half_y_Q']
    identity = pulse_dict['x_y_wait']
    element[channels[0]].add_segment(i_pulse)
    element[channels[1]].add_segment(q_pulse)
    element[channels[2]].add_segment(identity)
    element[channels[3]].add_segment(identity)


def do_z_pi(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['z_wait']
    z_pulse = pulse_dict['pi_z']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(z_pulse)
    element[channels[3]].add_segment(identity)


def do_z_pi_half(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['z_wait']
    z_pulse = pulse_dict['pi_half_z']
    element[channels[0]].add_segment(identity)
    element[channels[1]].add_segment(identity)
    element[channels[2]].add_segment(z_pulse)
    element[channels[3]].add_segment(identity)


def do_identity(element, pulse_dict, channels=[1, 2, 3, 4]):
    identity = pulse_dict['x_y_wait']
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
        w.add_segment(identity, position=0)


def execute_gates(element, pulse_dict, gate_list,
                  channels=[1, 2, 3, 4], spacing=None):
    for i in gate_list:
        if i is 'id':
            do_identity(element, pulse_dict, channels=channels)
        if i is 'x_pi':
            do_x_pi(element, pulse_dict, channels=channels)
        elif i is 'x_pi_half':
            do_x_pi_half(element, pulse_dict, channels=channels)
        elif i is 'y_pi':
            do_y_pi(element, pulse_dict, channels=channels)
        elif i is 'y_pi_half':
            do_y_pi_half(element, pulse_dict, channels=channels)
        elif i is 'z_pi':
            do_z_pi(element, pulse_dict, channels=channels)
        elif i is 'z_pi_half':
            do_z_pi_half(element, pulse_dict, channels=channels)
        if spacing is not None:
            wait(element, pulse_dict, spacing, channels=channels)


###########################################################
# Sequence building function (with gates)
###########################################################


def make_sequence_from_gate_lists(gate_list, SSBfreq=None, drag=False,
                                  channels=[1, 2, 3, 4], name=None,
                                  variable_label=None, spacing=None):
    calib_dict = get_calibration_dict()
    pulse_dict = make_pulse_dict(SSBfreq=SSBfreq, drag=drag)
    seq = Sequence(name=name, variable_label=variable_label)

    for i, gate_list in enumerate(gate_list):
        i_wf = Waveform(channel=channels[0])
        q_wf = Waveform(channel=channels[1])
        z_wf = Waveform(channel=channels[2])
        measure_wf = Waveform(channel=channels[3])
        element = Element()
        element.add_waveform(i_wf)
        element.add_waveform(q_wf)
        element.add_waveform(z_wf)
        element.add_waveform(measure_wf)
        element.sample_rate = calib_dict['sample_rate']
        execute_gates(element, pulse_dict, gate_list,
                      spacing=spacing)
        wait(element, pulse_dict, calib_dict['pulse_readout_delay'],
             channels=channels)
        measure(element, pulse_dict, channels=channels)
        prepend_compensating_wait_to_element(element, pulse_dict,
                                             calib_dict['cycle_time'])
        if i == 0:
            measure_wf.add_marker(2, 0, calib_dict['marker_time'])
        seq.add_element(element)
    return seq
