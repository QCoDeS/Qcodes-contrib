from . import make_readout_wf, get_calibration_dict, \
    make_time_varying_sequence, make_multi_varying_sequence, \
    make_time_multi_varying_sequence, \
    cos_array, sin_array, flat_array, gaussian_array, cos_gaussian_array, \
    sin_gaussian_array
from . import Segment, Waveform, Element, Sequence


def make_floquet_dur_sequence(start, stop, step, amp=1, floquet_freq=1e6,
    channels=[1, 4], form='cos'):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                    calib_dict['pulse_end'][qubit])


    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    if form == 'cos':
        floquet_drive = Segment(
            name='floquet_drive', gen_func=cos_array,
            func_args={'amp': amp, 'freq': floquet_freq})
    elif form == 'sin':
        floquet_drive = Segment(
            name='floquet_drive', gen_func=sin_array,
            func_args={'amp': amp, 'freq': floquet_freq})
    else:
        raise Exception('unrecognised form, should be cos or sin')
    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit})

    qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, floquet_drive, wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    floquet_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    floquet_element.add_waveform(qubit_wf)
    floquet_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    floquet_sequence = make_time_varying_sequence(
        floquet_element, channels[0], 1, 'dur', start, stop, step,
        0, calib_dict['cycle_time'][qubit],
        name='floquet_seq',
        variable_name='floquet_drive_dur', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    floquet_sequence.labels = {'seq_type': 'floquet'}
    return floquet_sequence


def make_floquet_dur_seq_gated(start, stop, step, amp=1, floquet_freq=1e6,
    channels=[1, 4], form='cos', pi_half_before=True, pi_half_after=False,
    gaussian=True, pi_half_after_neg=False):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                    calib_dict['pulse_end'][qubit])
    pi_half_amp = calib_dict['pi_half_pulse_amp'][qubit]


    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    if form == 'cos':
        floquet_drive = Segment(
            name='floquet_drive', gen_func=cos_array,
            func_args={'amp': amp, 'freq': floquet_freq})
    elif form == 'sin':
        floquet_drive = Segment(
            name='floquet_drive', gen_func=sin_array,
            func_args={'amp': amp, 'freq': floquet_freq})
    else:
        raise Exception('unrecognised form, should be cos or sin')
    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit})

    if gaussian:
        pi_half_sigma = calib_dict['pi_pulse_sigma'][qubit]
        pi_half_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_half_amp, 'sigma': pi_half_sigma})
        pi_half_segment_neg = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_half_amp, 'sigma': pi_half_sigma,
                       'positive': False})
    else:
        pi_half_dur = calib_dict['pi_half_pulse_dur'][qubit]
        pi_half_segment = Segment(
            name='square_pi_pulse', gen_func=flat_array,
            func_args={'amp': pi_half_amp, 'dur': pi_half_dur})
        pi_half_segment_neg = Segment(
            name='square_pi_pulse', gen_func=flat_array,
            func_args={'amp': -1 * pi_half_amp, 'dur': pi_half_dur})

    if pi_half_before and pi_half_after and not pi_half_after_neg:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, pi_half_segment, floquet_drive,
                                    pi_half_segment, wait_segment])
        floquet_seg_index = 2
    elif pi_half_before and pi_half_after:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, pi_half_segment, floquet_drive,
                                    pi_half_segment_neg, wait_segment])
        floquet_seg_index = 2
    elif pi_half_before:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, pi_half_segment, floquet_drive,
                                    wait_segment])
        floquet_seg_index = 2
    elif pi_half_after and not pi_half_after_neg:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, floquet_drive,
                                    pi_half_segment, wait_segment])
        floquet_seg_index = 1
    elif pi_half_after:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, floquet_drive,
                                    pi_half_segment_neg, wait_segment])
        floquet_seg_index = 1
    else:
        qubit_wf = Waveform(channel=channels[0],
                    segment_list=[compensating_wait_segment, floquet_drive,
                                    wait_segment])
        floquet_seg_index = 1

    readout_wf = make_readout_wf(channel=channels[-1])

    floquet_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    floquet_element.add_waveform(qubit_wf)
    floquet_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    floquet_sequence = make_time_varying_sequence(
        floquet_element, channels[0], floquet_seg_index, 'dur', start, stop, step,
        0, calib_dict['cycle_time'][qubit],
        name='floquet_seq',
        variable_name='floquet_drive_dur', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    floquet_sequence.labels = {'seq_type': 'floquet', 'pi_half_before': pi_half_before,
                                'pi_half_after': pi_half_after}
    return floquet_sequence