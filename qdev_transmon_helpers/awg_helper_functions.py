import numpy as np
from . import get_latest_counter, \
    make_SSB_I_gaussian, make_SSB_Q_gaussian, get_calibration_dict, \
    get_calibration_val, get_pulse_location,

from . import Segment, Waveform, Element, Sequence

# TODO: rount -> int/ciel
# TODO: test segments
# TODO: make_save_send_load_awg_file -> not doing same thing twice!
# TODO: allxy clean up
# TODO: docstrings
# TODO: checks
# TODO: T2 echo
# TODO: refactor so it works with new sequences etc

#################################################################
# General AWG and Sequence functions
#################################################################


def make_save_send_load_awg_file(awg, sequence):
    """
    WYSIYWYG

    Args:
        awg instrument for upload
        sequence to be uploaded
    """
    pulse_location = get_pulse_location()
    try:
        num = get_latest_counter(pulse_location) + 1
    except (FileNotFoundError, OSError):
        num = 1
    name = '{0:03d}_{name}.awg'.format(num, name=sequence.name)
    file_location = pulse_location + name
    awg.make_and_save_awg_file(*sequence.unwrap(), filename=file_location)
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


#######################################
# Some handy functions with segments
#######################################


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
    wait1_segment = Segment(name='wait', gen_func=flat,
                            func_args={'amp': 0, 'dur': time_before_readout})
    wait2_segment = Segment(name='wait', gen_func=flat,
                            func_args={'amp': 0, 'dur': time_after_readout})
    readout_wf = Waveform(
        channel=channel,
        segment_list=[wait1_segment, measurement_segment, wait2_segment],
        sample_rate=p_dict['sample_rate'])
    if first_in_seq:
        readout_wf.add_marker(
            2, 0, int(p_dict['marker_time'] * p_dict['sample_rate']))
    return readout_wf


def make_varying_sequence(element_template, vary_ch, vary_seg,
                          arg_to_vary, start, stop, step,
                          readout_ch=4, marker_points=100):
    sequence = Sequence(variable=arg_to_vary,
                        start=start, stop=stop, step=step)
    for i, val in enumerate(sequence.variable_array):
        elem = element_template.copy()
        elem[vary_ch].segment_list[vary_seg].func_args[arg_to_vary] = val
        if i == 0:
            elem[readout_ch].add_marker(2, 0, marker_points)
        sequence.add_element(elem)
    return sequence


def make_time_varying_sequence(element_template, vary_ch, vary_seg,
                               arg_to_vary, start, stop, step,
                               compensate_seg, total_time, readout_ch=4,
                               marker_points=100):
    sequence = Sequence(variable=arg_to_vary, start=start,
                        stop=stop, step=step)
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
                                readout_ch=4, marker_points=100):
    variable_name = arg_to_vary_1 + '_' + arg_to_vary_2
    sequence = Sequence(
        variable=variable_name, start=start1, stop=stop1, step=step1)
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
                                     total_time,
                                     readout_ch=4, marker_points=100):
    variable_name = arg_to_vary_1 + '_' + arg_to_vary_2
    sequence = Sequence(
        variable=variable_name, start=start1, stop=stop1, step=step1)
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


##################################################################
# AllXY
##################################################################


def make_allxySSB_seq(pi_duration, pi_amp, SSBfreq, sigma_cuttoff, channels=[1,2,4], pulse_mod=False):
    """
    Oh dear.
    """

    validate_pulse_dictionary()

    p_dict = get_pulse_dict()

    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    pi_pulseI_x = make_SSB_I_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp, SSBfreq)
    pi_pulseQ_x = make_SSB_Q_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp, SSBfreq)
    pi_half_pulseI_x = make_SSB_I_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp/2, SSBfreq)
    pi_half_pulseQ_x = make_SSB_Q_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp/2, SSBfreq)
    pi_pulseI_y = -1*make_SSB_Q_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp, SSBfreq)
    pi_pulseQ_y = make_SSB_I_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp, SSBfreq)
    pi_half_pulseI_y = -1*make_SSB_Q_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp/2, SSBfreq)
    pi_half_pulseQ_y = make_SSB_I_gaussian(p_dict['sample_rate'], pi_duration, sigma_cuttoff, pi_amp/2, SSBfreq)

    readout_waveform =  Waveform(length=total_points, channel=channels[2])
    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    pulse_points = len(pi_pulseI_x)
    pulse_mod_points = int(pi_duration * sigma_cuttoff * 4 / resolution)

    seq = Sequence(name='allxy', variable='operation_combination', variable_label='Operation Combination Id')

    elem_0 = Element()
    x_waveform_0 = Waveform(length=total_points, channel=1)
    y_waveform_0 = Waveform(length=total_points, channel=2)
    readout_first = readout_waveform.copy()
    readout_first.marker_2[10:10 + marker_points] = 1
    if pulse_mod:
        x_waveform_0.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_0.add_waveform(x_waveform_0)
    elem_0.add_waveform(y_waveform_0)
    elem_0.add_waveform(readout_first)
    seq.add_element(elem_0)

    elem_1 = Element()
    x_waveform_1 = Waveform(length=total_points, channel=channels[0])
    y_waveform_1 = Waveform(length=total_points, channel=channels[1])
    x_waveform_1.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    y_waveform_1.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    x_waveform_1.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    y_waveform_1.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_x
    if pulse_mod:
        x_waveform_1.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_1.add_waveform(x_waveform_1)
    elem_1.add_waveform(y_waveform_1)
    elem_1.add_waveform(readout_first)
    seq.add_element(elem_1)

    elem_2 = Element()
    x_waveform_2 = Waveform(length=total_points, channel=channels[0])
    y_waveform_2 = Waveform(length=total_points, channel=channels[1])
    x_waveform_2.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    y_waveform_2.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    x_waveform_2.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    y_waveform_2.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_y
    if pulse_mod:
        x_waveform_2.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_2.add_waveform(x_waveform_2)
    elem_2.add_waveform(y_waveform_2)
    elem_2.add_waveform(readout_waveform)
    seq.add_element(elem_2)

    elem_3 = Element()
    x_waveform_3 = Waveform(length=total_points, channel=channels[0])
    y_waveform_3 = Waveform(length=total_points, channel=channels[1])
    x_waveform_3.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    y_waveform_3.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    x_waveform_3.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    y_waveform_3.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_y
    if pulse_mod:
        x_waveform_3.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_3.add_waveform(x_waveform_3)
    elem_3.add_waveform(y_waveform_3)
    elem_3.add_waveform(readout_waveform)
    seq.add_element(elem_3)

    elem_4 = Element()
    x_waveform_4 = Waveform(length=total_points, channel=channels[0])
    y_waveform_4 = Waveform(length=total_points, channel=channels[1])
    x_waveform_4.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    y_waveform_4.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    x_waveform_4.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    y_waveform_4.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_x
    if pulse_mod:
        x_waveform_4.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_4.add_waveform(x_waveform_4)
    elem_4.add_waveform(y_waveform_4)
    elem_4.add_waveform(readout_waveform)
    seq.add_element(elem_4)

    elem_5 = Element()
    x_waveform_5 = Waveform(length=total_points, channel=channels[0])
    y_waveform_5 = Waveform(length=total_points, channel=channels[1])
    x_waveform_5.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    y_waveform_5.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    if pulse_mod:
        x_waveform_5.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_5.add_waveform(x_waveform_5)
    elem_5.add_waveform(y_waveform_5)
    elem_5.add_waveform(readout_waveform)
    seq.add_element(elem_5)

    elem_6 = Element()
    x_waveform_6 = Waveform(length=total_points, channel=channels[0])
    y_waveform_6 = Waveform(length=total_points, channel=channels[1])
    x_waveform_6.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    y_waveform_6.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    if pulse_mod:
        x_waveform_6.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_6.add_waveform(x_waveform_6)
    elem_6.add_waveform(y_waveform_6)
    elem_6.add_waveform(readout_waveform)
    seq.add_element(elem_6)

    elem_7 = Element()
    x_waveform_7 = Waveform(length=total_points, channel=channels[0])
    y_waveform_7 = Waveform(length=total_points, channel=channels[1])
    x_waveform_7.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    y_waveform_7.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    x_waveform_7.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    y_waveform_7.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    if pulse_mod:
        x_waveform_7.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_7.add_waveform(x_waveform_7)
    elem_7.add_waveform(y_waveform_7)
    elem_7.add_waveform(readout_waveform)
    seq.add_element(elem_7)

    elem_8 = Element()
    x_waveform_8 = Waveform(length=total_points, channel=channels[0])
    y_waveform_8 = Waveform(length=total_points, channel=channels[1])
    x_waveform_8.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    y_waveform_8.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    x_waveform_8.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    y_waveform_8.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    if pulse_mod:
        x_waveform_8.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_8.add_waveform(x_waveform_8)
    elem_8.add_waveform(y_waveform_8)
    elem_8.add_waveform(readout_waveform)
    seq.add_element(elem_8)

    elem_9 = Element()
    x_waveform_9 = Waveform(length=total_points, channel=channels[0])
    y_waveform_9 = Waveform(length=total_points, channel=channels[1])
    x_waveform_9.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    y_waveform_9.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    x_waveform_9.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    y_waveform_9.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    if pulse_mod:
        x_waveform_9.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_9.add_waveform(x_waveform_9)
    elem_9.add_waveform(y_waveform_9)
    elem_9.add_waveform(readout_waveform)
    seq.add_element(elem_9)

    elem_10 = Element()
    x_waveform_10 = Waveform(length=total_points, channel=channels[0])
    y_waveform_10 = Waveform(length=total_points, channel=channels[1])
    x_waveform_10.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    y_waveform_10.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    x_waveform_10.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    y_waveform_10.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    if pulse_mod:
        x_waveform_10.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_10.add_waveform(x_waveform_10)
    elem_10.add_waveform(y_waveform_10)
    elem_10.add_waveform(readout_waveform)
    seq.add_element(elem_10)

    elem_11 = Element()
    x_waveform_11 = Waveform(length=total_points, channel=channels[0])
    y_waveform_11 = Waveform(length=total_points, channel=channels[1])
    x_waveform_11.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    y_waveform_11.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    x_waveform_11.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    y_waveform_11.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    if pulse_mod:
        x_waveform_11.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_11.add_waveform(x_waveform_11)
    elem_11.add_waveform(y_waveform_11)
    elem_11.add_waveform(readout_waveform)
    seq.add_element(elem_11)

    elem_12 = Element()
    x_waveform_12 = Waveform(length=total_points, channel=channels[0])
    y_waveform_12 = Waveform(length=total_points, channel=channels[1])
    x_waveform_12.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    y_waveform_12.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    x_waveform_12.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    y_waveform_12.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    if pulse_mod:
        x_waveform_12.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_12.add_waveform(x_waveform_12)
    elem_12.add_waveform(y_waveform_12)
    elem_12.add_waveform(readout_waveform)
    seq.add_element(elem_12)

    elem_13 = Element()
    x_waveform_13 = Waveform(length=total_points, channel=channels[0])
    y_waveform_13 = Waveform(length=total_points, channel=channels[1])
    x_waveform_13.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    y_waveform_13.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    x_waveform_13.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    y_waveform_13.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    if pulse_mod:
        x_waveform_13.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_13.add_waveform(x_waveform_13)
    elem_13.add_waveform(y_waveform_13)
    elem_13.add_waveform(readout_waveform)
    seq.add_element(elem_13)

    elem_14 = Element()
    x_waveform_14 = Waveform(length=total_points, channel=channels[0])
    y_waveform_14 = Waveform(length=total_points, channel=channels[1])
    x_waveform_14.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    y_waveform_14.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    x_waveform_14.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    y_waveform_14.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    if pulse_mod:
        x_waveform_14.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_14.add_waveform(x_waveform_14)
    elem_14.add_waveform(y_waveform_14)
    elem_14.add_waveform(readout_waveform)
    seq.add_element(elem_14)

    elem_15 = Element()
    x_waveform_15 = Waveform(length=total_points, channel=channels[0])
    y_waveform_15 = Waveform(length=total_points, channel=channels[1])
    x_waveform_15.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    y_waveform_15.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    x_waveform_15.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    y_waveform_15.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    if pulse_mod:
        x_waveform_15.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_15.add_waveform(x_waveform_15)
    elem_15.add_waveform(y_waveform_15)
    elem_15.add_waveform(readout_waveform)
    seq.add_element(elem_15)

    elem_16 = Element()
    x_waveform_16 = Waveform(length=total_points, channel=channels[0])
    y_waveform_16 = Waveform(length=total_points, channel=channels[1])
    x_waveform_16.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    y_waveform_16.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    x_waveform_16.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    y_waveform_16.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    if pulse_mod:
        x_waveform_16.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_16.add_waveform(x_waveform_16)
    elem_16.add_waveform(y_waveform_16)
    elem_16.add_waveform(readout_waveform)
    seq.add_element(elem_16)

    elem_17 = Element()
    x_waveform_17 = Waveform(length=total_points, channel=channels[0])
    y_waveform_17 = Waveform(length=total_points, channel=channels[1])
    x_waveform_17.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    y_waveform_17.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    if pulse_mod:
        x_waveform_17.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_17.add_waveform(x_waveform_17)
    elem_17.add_waveform(y_waveform_17)
    elem_17.add_waveform(readout_waveform)
    seq.add_element(elem_17)

    elem_18 = Element()
    x_waveform_18 = Waveform(length=total_points, channel=channels[0])
    y_waveform_18 = Waveform(length=total_points, channel=channels[1])
    x_waveform_18.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    y_waveform_18.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    if pulse_mod:
        x_waveform_18.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_18.add_waveform(x_waveform_18)
    elem_18.add_waveform(y_waveform_18)
    elem_18.add_waveform(readout_waveform)
    seq.add_element(elem_18)

    elem_19 = Element()
    x_waveform_19 = Waveform(length=total_points, channel=channels[0])
    y_waveform_19 = Waveform(length=total_points, channel=channels[1])
    x_waveform_19.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    y_waveform_19.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    x_waveform_19.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    y_waveform_19.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    if pulse_mod:
        x_waveform_19.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_19.add_waveform(x_waveform_19)
    elem_19.add_waveform(y_waveform_19)
    elem_19.add_waveform(readout_waveform)
    seq.add_element(elem_19)

    elem_20 = Element()
    x_waveform_20 = Waveform(length=total_points, channel=channels[0])
    y_waveform_20 = Waveform(length=total_points, channel=channels[1])
    x_waveform_20.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    y_waveform_20.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    x_waveform_20.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    y_waveform_20.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    if pulse_mod:
        x_waveform_20.marker_1[pulse_end_points-pulse_mod_points:pulse_end_points] = 1
    elem_20.add_waveform(x_waveform_20)
    elem_20.add_waveform(y_waveform_20)
    elem_20.add_waveform(readout_waveform)
    seq.add_element(elem_20)

    return seq

def make_floquet_dur_sequence(floquet_freq, start, stop, step,
                              channels=[1, 4]):
    """
    Cosine modulated qubit drive at floquet frequency for varying durations at
    floquet amplitude. Square readout pulse with markers (1 for readout start,
    2 for seq start)
    """
    validate_pulse_dictionary()
    duration_floquet_sequence = Sequence(name='floquet',
                                         variable='floquet_dur',
                                         variable_label='time',
                                         variable_unit='S',
                                         step=step,
                                         start=start,
                                         stop=stop)
    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(p_dict['readout_time'] / resolution))
    pulse_end_points = int(np.ceil(p_dict['pulse_end'] / resolution))
    marker_points = int(np.ceil(p_dict['marker_time'] / resolution))
    total_points = int(np.ceil(p_dict['cycle_duration'] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    floquet_durations = duration_floquet_sequence.variable_array
    for i, floquet_duration in enumerate(floquet_durations):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1

        floquet_points = int(np.ceil(floquet_duration / resolution))
        floquet_time_array = np.arange(floquet_points) * resolution
        angle = floquet_time_array * floquet_freq * 2 * np.pi
        cos_array = np.cos(angle)
        floquet_start = pulse_end_points - floquet_points
        qubit_waveform.wave[floquet_start:pulse_end_points] = cos_array
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        duration_floquet_sequence.add_element(element)
    duration_floquet_sequence.check()
    return duration_floquet_sequence


def make_floquet_freq_sequence(floquet_dur, floquet_amp, start, stop, step,
                               channels=[1, 4]):
    """
    Cosine modulated qubit drive duation floquet_dur for varying frequencies
    at floquet amplitude. Square readout pulse with markers (1 for readout
    start, 2 for seq start)
    """
    validate_pulse_dictionary()
    frequency_floquet_sequence = Sequence(name='floquet',
                                          variable='floquet_freq',
                                          variable_label='floquet_freq',
                                          variable_unit='Hz',
                                          step=step,
                                          start=start,
                                          stop=stop)
    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(p_dict['readout_time'] / resolution))
    pulse_end_points = int(np.ceil(p_dict['pulse_end'] / resolution))
    marker_points = int(np.ceil(p_dict['marker_time'] / resolution))
    floquet_points = int(np.ceil(floquet_dur / resolution))
    total_points = int(np.ceil(p_dict['cycle_duration'] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    floquet_freqs = frequency_floquet_sequence.variable_array
    for i, floquet_freq in enumerate(floquet_freqs):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1

        floquet_time_array = np.arange(floquet_points) * resolution
        angle = floquet_time_array * floquet_freq * 2 * np.pi
        cos_array = floquet_amp * np.cos(angle)
        floquet_start = pulse_end_points - floquet_points
        qubit_waveform.wave[floquet_start:pulse_end_points] = cos_array
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        frequency_floquet_sequence.add_element(element)
    frequency_floquet_sequence.check()
    return frequency_floquet_sequence


def make_floquet_amp_sequence(floquet_dur, floquet_freq, start, stop, step,
                              channels=[1, 4]):
    """
    Cosine modulated qubit drive at floquet frequency for floquet_dur at
    varying amplitude. Square readout pulse with markers (1 for readout start,
    2 for seq start)
    """
    validate_pulse_dictionary()
    amplitude_floquet_sequence = Sequence(name='floquet',
                                          variable='floquet_amp',
                                          variable_label='floquet_amp',
                                          variable_unit='V',
                                          step=step,
                                          start=start,
                                          stop=stop)
    p_dict = get_pulse_dict()

    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(p_dict['readout_time'] / resolution))
    pulse_end_points = int(np.ceil(p_dict['pulse_end'] / resolution))
    marker_points = int(np.ceil(p_dict['marker_time'] / resolution))
    floquet_points = int(np.ceil(floquet_dur / resolution))
    total_points = int(np.ceil(p_dict['cycle_duration'] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    floquet_amps = amplitude_floquet_sequence.variable_array
    for i, floquet_amp in enumerate(floquet_amps):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1

        floquet_time_array = np.arange(floquet_points) * resolution
        angle = floquet_time_array * floquet_freq * 2 * np.pi
        cos_array = floquet_amp * np.cos(angle)
        floquet_start = pulse_end_points - floquet_points
        qubit_waveform.wave[floquet_start:pulse_end_points] = cos_array
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        amplitude_floquet_sequence.add_element(element)
    amplitude_floquet_sequence.check()
    return amplitude_floquet_sequence


##################################################################
# AllXY
##################################################################
