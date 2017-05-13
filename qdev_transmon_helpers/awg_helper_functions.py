import pickle
import numpy as np
import pprint
import matplotlib.pyplot as plt
from . import get_pulse_location, get_latest_counter, make_gaussian

from . import Sequence, Waveform, Element

# TODO: gaussian
# TODO: rount -> int/ciel
# TODO: test segments
# TODO: plot sequence -> plot variable channel count
# TODO: make_save_send_load_awg_file -> not doing same thing twice!
# TODO: allxy clean up


default_pulse_dict = {'cycle_duration': 20e-6, 'sample_rate': 1e9,
                      'pulse_end': 10e-6, 'pulse_readout_delay': 30e-9,
                      'readout_time': 4e-6, 'marker_time': 500e-9,
                      'marker_readout_delay': 0, 'qubit_time': 1e-6}


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
    print('pulse config:')
    pprint.pprint(p_dict, width=1)


def get_pulse_dict():
    """
    Gets the pulse dictionary stored as 'pulse_dict.p' in
    pulse location file or makes an empty one if none found.

    Returns:
        pulse dictionary
    """
    path = get_pulse_location()
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
    path = get_pulse_location()
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


def make_readout_seq(channels=[4]):
    """
    Square pulse duting readout time with one marker at start
    """
    validate_pulse_dictionary()
    readout_sequence = Sequence(name='plain_readout', variable='')
    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[0])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    element = Element()
    element.add_waveform(readout_waveform)
    readout_sequence.add_element(element)
    readout_sequence.check()
    return readout_sequence


def make_full_ssb_wave(freq=8e6, duration=20e-6, channels=[1, 2]):
    """
    Cosine and sine waves of given frequency for cycle duration, no markers
    """
    seq = Sequence(name='ssb_seq')
    element = Element()
    resolution = 1 / 1e9
    total_points = duration / resolution
    waveform_i = Waveform(channel=channels[0])
    waveform_q = Waveform(channel=channels[1])
    time_array = np.arange(total_points) * resolution
    angle = time_array * freq * 2 * np.pi
    cos_array = np.cos(angle)
    sin_array = np.sin(angle)
    waveform_i.wave = cos_array
    waveform_q.wave = -1 * sin_array
    element.add_waveform(waveform_i)
    element.add_waveform(waveform_q)
    seq.add_element(element)
    seq.check()
    return seq


def make_ssb_qubit_seq(start=0, stop=200e6, step=1e6, channels=[1, 2, 4]):
    """
    Cosine and sine waves for qubit time with range of frequencies, square
    readout wave for readout time. Markers on readout channel (1 for readout
    start, 2 for seq start)
    """
    validate_pulse_dictionary()
    ssb_sequence = Sequence(name='qubit_ssb',
                            variable='ssb_qubit_modulation_freq',
                            variable_unit='Hz',
                            step=step,
                            start=start,
                            stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = round(p_dict['qubit_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_template = Waveform(length=total_points, channel=channels[2])
    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    qubit_time_array = np.arange(qubit_points) * resolution
    freq_array = ssb_sequence.variable_array

    for i, freq in enumerate(freq_array):
        element = Element()
        qubit_i = Waveform(length=total_points, channel=channels[0])
        qubit_q = Waveform(length=total_points, channel=channels[1])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = pulse_end_points - qubit_points
        qubit_end = pulse_end_points
        angle = qubit_time_array * freq * 2 * np.pi
        cos_array = np.cos(angle)
        sin_array = np.sin(angle)
        qubit_i.wave[qubit_start:qubit_end] = cos_array
        qubit_q.wave[qubit_start:qubit_end] = -1 * sin_array
        element.add_waveform(qubit_i)
        element.add_waveform(qubit_q)
        element.add_waveform(readout_waveform)
        ssb_sequence.add_element(element)
    ssb_sequence.check()
    return ssb_sequence


def make_t1_seq(pi_duration, pi_amp, start=0, stop=5e-6, step=50e-9,
                channels=[1, 4]):
    """
    Square qubit drive for pi duration at pi amplitude on qubit with varying
    wait time before readout (square pulse for readout time). Markers on
    readout channel (1 for readout start, 2 for seq start)
    """
    validate_pulse_dictionary()
    t1_sequence = Sequence(name='t1',
                           variable='drive_readout_delay',
                           variable_label='Delay',
                           variable_unit='s',
                           step=step,
                           start=start,
                           stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = pi_duration / resolution
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    delay_array_points = np.round(
        t1_sequence.variable_array / resolution).astype(np.int)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = int(pulse_end_points - delay_points - qubit_points)
        qubit_end = int(qubit_start + qubit_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        t1_sequence.add_element(element)
    t1_sequence.check()
    return t1_sequence


def make_rabi_sequence(pi_amp, start=0, stop=200e-9, step=2e-9,
                       channels=[1, 4]):
    """
    Square qubit drive of pi amplitude of varying duration, square readout
    drive. Markers on readout channel (1 for readout start, 2 for seq start)
    """
    validate_pulse_dictionary()
    rabi_sequence = Sequence(name='rabi',
                             variable='drive_duration',
                             variable_label='Drive Duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

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

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    qubit_duration_array_points = np.round(
        rabi_sequence.variable_array / resolution).astype(int)

    for i, qubit_points in enumerate(qubit_duration_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = int(pulse_end_points - qubit_points)
        qubit_end = int(pulse_end_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        rabi_sequence.add_element(element)

    rabi_sequence.check()

    return rabi_sequence


def make_ramsey_sequence(pi_duration, pi_amp, start=0, stop=200e-9, step=2e-9,
                         channels=[1, 4]):
    """
    Two square pulses on qubit of pi duration of half pi amplitude separated
    by varying duration. Square readout with markers (1 for readout start,
    2 for seq start)
    """
    validate_pulse_dictionary()
    ramsey_sequence = Sequence(name='ramsey',
                               variable='drive_drive_delay',
                               variable_label='Delay',
                               variable_unit='s',
                               step=step,
                               start=start,
                               stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = round(pi_duration / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    delay_array_points = np.round(
        ramsey_sequence.variable_array / resolution).astype(np.int)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start_first = int(pulse_end_points - delay_points - 2 * qubit_points)
        qubit_end_first = int(qubit_start_first + qubit_points)
        qubit_start_second = int(pulse_end_points - qubit_points)
        qubit_end_second = int(qubit_start_second + qubit_points)
        qubit_waveform.wave[qubit_start_first:qubit_end_first] = pi_amp / 2
        qubit_waveform.wave[qubit_start_second:qubit_end_second] = pi_amp / 2
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        ramsey_sequence.add_element(element)
    ramsey_sequence.check()
    return ramsey_sequence


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

def make_rabi_gaussian_sequence(pi_amp, pulse_sigmas_number, start=0, stop=200e-9, step=2e-9,
                       channels=[1, 4]):
    """
    Square qubit drive of pi amplitude of varying duration, square readout
    drive. Markers on readout channel (1 for readout start, 2 for seq start)
    """
    validate_pulse_dictionary()
    rabi_sequence = Sequence(name='rabi',
                             variable='drive_duration',
                             variable_label='Drive Duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

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

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1

    for i, pi_duration in enumerate(rabi_sequence.variable_array):
        pi_pulse = make_gaussian(p_dict['sample_rate'], pi_duration, pulse_sigmas_number, pi_amp)
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = int(pulse_end_points - len(pi_pulse))
        qubit_end = int(pulse_end_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_pulse
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        rabi_sequence.add_element(element)

    rabi_sequence.check()

    return rabi_sequence


def make_allxy_seq(pi_duration, total_sigmas, channels=[1,2,4]):
    """
    Fuck this fucking function
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

    pi_pulse = make_gaussian(p_dict['sample_rate'], pi_duration, total_sigmas, 1)
    pi_half_pulse = make_gaussian(p_dict['sample_rate'], pi_duration, total_sigmas, 0.5)

    readout_waveform =  Waveform(length=total_points, channel=channels[2])
    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start_points:readout_marker_start_points + marker_points] = 1


    pulse_points = len(pi_pulse)

    seq = Sequence(name='allxy', variable='operation_combination', variable_label='Operation Combination Id')

    elem_0 = Element()
    x_waveform_0 = Waveform(length=total_points, channel=1)
    y_waveform_0 = Waveform(length=total_points, channel=2)
    readout_first = readout_waveform.copy()
    readout_first.marker_2[10:10 + marker_points] = 1
    elem_0.add_waveform(x_waveform_0)
    elem_0.add_waveform(y_waveform_0)
    elem_0.add_waveform(readout_first)
    seq.add_element(elem_0)

    elem_1 = Element()
    x_waveform_1 = Waveform(length=total_points, channel=channels[0])
    y_waveform_1 = Waveform(length=total_points, channel=channels[1])
    x_waveform_1.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    x_waveform_1.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_1.add_waveform(x_waveform_1)
    elem_1.add_waveform(y_waveform_1)
    elem_1.add_waveform(readout_first)
    seq.add_element(elem_1)

    elem_2 = Element()
    x_waveform_2 = Waveform(length=total_points, channel=channels[0])
    y_waveform_2 = Waveform(length=total_points, channel=channels[1])
    y_waveform_2.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    y_waveform_2.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_2.add_waveform(x_waveform_2)
    elem_2.add_waveform(y_waveform_2)
    elem_2.add_waveform(readout_waveform)
    seq.add_element(elem_2)

    elem_3 = Element()
    x_waveform_3 = Waveform(length=total_points, channel=channels[0])
    y_waveform_3 = Waveform(length=total_points, channel=channels[1])
    x_waveform_3.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    y_waveform_3.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_3.add_waveform(x_waveform_3)
    elem_3.add_waveform(y_waveform_3)
    elem_3.add_waveform(readout_waveform)
    seq.add_element(elem_3)

    elem_4 = Element()
    x_waveform_4 = Waveform(length=total_points, channel=channels[0])
    y_waveform_4 = Waveform(length=total_points, channel=channels[1])
    y_waveform_4.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    x_waveform_4.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_4.add_waveform(x_waveform_4)
    elem_4.add_waveform(y_waveform_4)
    elem_4.add_waveform(readout_waveform)
    seq.add_element(elem_4)

    elem_5 = Element()
    x_waveform_5 = Waveform(length=total_points, channel=channels[0])
    y_waveform_5 = Waveform(length=total_points, channel=channels[1])
    x_waveform_5.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    elem_5.add_waveform(x_waveform_5)
    elem_5.add_waveform(y_waveform_5)
    elem_5.add_waveform(readout_waveform)
    seq.add_element(elem_5)

    elem_6 = Element()
    x_waveform_6 = Waveform(length=total_points, channel=channels[0])
    y_waveform_6 = Waveform(length=total_points, channel=channels[1])
    y_waveform_6.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    elem_6.add_waveform(x_waveform_6)
    elem_6.add_waveform(y_waveform_6)
    elem_6.add_waveform(readout_waveform)
    seq.add_element(elem_6)

    elem_7 = Element()
    x_waveform_7 = Waveform(length=total_points, channel=channels[0])
    y_waveform_7 = Waveform(length=total_points, channel=channels[1])
    x_waveform_7.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    y_waveform_7.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_7.add_waveform(x_waveform_7)
    elem_7.add_waveform(y_waveform_7)
    elem_7.add_waveform(readout_waveform)
    seq.add_element(elem_7)

    elem_8 = Element()
    x_waveform_8 = Waveform(length=total_points, channel=channels[0])
    y_waveform_8 = Waveform(length=total_points, channel=channels[1])
    y_waveform_8.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    x_waveform_8.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_8.add_waveform(x_waveform_8)
    elem_8.add_waveform(y_waveform_8)
    elem_8.add_waveform(readout_waveform)
    seq.add_element(elem_8)

    elem_9 = Element()
    x_waveform_9 = Waveform(length=total_points, channel=channels[0])
    y_waveform_9 = Waveform(length=total_points, channel=channels[1])
    x_waveform_9.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    y_waveform_9.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_9.add_waveform(x_waveform_9)
    elem_9.add_waveform(y_waveform_9)
    elem_9.add_waveform(readout_waveform)
    seq.add_element(elem_9)

    elem_10 = Element()
    x_waveform_10 = Waveform(length=total_points, channel=channels[0])
    y_waveform_10 = Waveform(length=total_points, channel=channels[1])
    y_waveform_10.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    x_waveform_10.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_10.add_waveform(x_waveform_10)
    elem_10.add_waveform(y_waveform_10)
    elem_10.add_waveform(readout_waveform)
    seq.add_element(elem_10)

    elem_11 = Element()
    x_waveform_11 = Waveform(length=total_points, channel=channels[0])
    y_waveform_11 = Waveform(length=total_points, channel=channels[1])
    x_waveform_11.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    y_waveform_11.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_11.add_waveform(x_waveform_11)
    elem_11.add_waveform(y_waveform_11)
    elem_11.add_waveform(readout_waveform)
    seq.add_element(elem_11)

    elem_12 = Element()
    x_waveform_12 = Waveform(length=total_points, channel=channels[0])
    y_waveform_12 = Waveform(length=total_points, channel=channels[1])
    y_waveform_12.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    x_waveform_12.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_12.add_waveform(x_waveform_12)
    elem_12.add_waveform(y_waveform_12)
    elem_12.add_waveform(readout_waveform)
    seq.add_element(elem_12)

    elem_13 = Element()
    x_waveform_13 = Waveform(length=total_points, channel=channels[0])
    y_waveform_13 = Waveform(length=total_points, channel=channels[1])
    x_waveform_13.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    x_waveform_13.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_13.add_waveform(x_waveform_13)
    elem_13.add_waveform(y_waveform_13)
    elem_13.add_waveform(readout_waveform)
    seq.add_element(elem_13)

    elem_14 = Element()
    x_waveform_14 = Waveform(length=total_points, channel=channels[0])
    y_waveform_14 = Waveform(length=total_points, channel=channels[1])
    x_waveform_14.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    x_waveform_14.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_14.add_waveform(x_waveform_14)
    elem_14.add_waveform(y_waveform_14)
    elem_14.add_waveform(readout_waveform)
    seq.add_element(elem_14)

    elem_15 = Element()
    x_waveform_15 = Waveform(length=total_points, channel=channels[0])
    y_waveform_15 = Waveform(length=total_points, channel=channels[1])
    y_waveform_15.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    y_waveform_15.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulse
    elem_15.add_waveform(x_waveform_15)
    elem_15.add_waveform(y_waveform_15)
    elem_15.add_waveform(readout_waveform)
    seq.add_element(elem_15)

    elem_16 = Element()
    x_waveform_16 = Waveform(length=total_points, channel=channels[0])
    y_waveform_16 = Waveform(length=total_points, channel=channels[1])
    y_waveform_16.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    y_waveform_16.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_16.add_waveform(x_waveform_16)
    elem_16.add_waveform(y_waveform_16)
    elem_16.add_waveform(readout_waveform)
    seq.add_element(elem_16)

    elem_17 = Element()
    x_waveform_17 = Waveform(length=total_points, channel=channels[0])
    y_waveform_17 = Waveform(length=total_points, channel=channels[1])
    x_waveform_17.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    elem_17.add_waveform(x_waveform_17)
    elem_17.add_waveform(y_waveform_17)
    elem_17.add_waveform(readout_waveform)
    seq.add_element(elem_17)

    elem_18 = Element()
    x_waveform_18 = Waveform(length=total_points, channel=channels[0])
    y_waveform_18 = Waveform(length=total_points, channel=channels[1])
    y_waveform_18.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulse
    elem_18.add_waveform(x_waveform_18)
    elem_18.add_waveform(y_waveform_18)
    elem_18.add_waveform(readout_waveform)
    seq.add_element(elem_18)

    elem_19 = Element()
    x_waveform_19 = Waveform(length=total_points, channel=channels[0])
    y_waveform_19 = Waveform(length=total_points, channel=channels[1])
    x_waveform_19.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    x_waveform_19.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_19.add_waveform(x_waveform_19)
    elem_19.add_waveform(y_waveform_19)
    elem_19.add_waveform(readout_waveform)
    seq.add_element(elem_19)

    elem_20 = Element()
    x_waveform_20 = Waveform(length=total_points, channel=channels[0])
    y_waveform_20 = Waveform(length=total_points, channel=channels[1])
    y_waveform_20.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulse
    y_waveform_20.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulse
    elem_20.add_waveform(x_waveform_20)
    elem_20.add_waveform(y_waveform_20)
    elem_20.add_waveform(readout_waveform)
    seq.add_element(elem_20)

    return seq

def plot_sequence(sequence, elemnum=0, channels=[1, 2]):
    """
    Function which plots 2 channels and markers

    Args:
        sequence to plot
        elemnum to plot
        channels (list of two) to plot

    Returns:
        matplotlib fic
    """
    fig = plt.figure()
    plt_count = len(channels) * 2
    for i, chan in enumerate(channels):
        index_a = (plt_count * 100) + 10 + (2 * i) + 1
        index_b = index_a + 1
        ax_w = fig.add_subplot(index_a)
        ax_w.set_title('Channel {} waveform'.format(chan))
        ax_w.set_ylim([-1.1, 1.1])
        ax_m = fig.add_subplot(index_b)
        ax_m.set_title('Channel {} markers'.format(chan))
        ax_m.set_ylim([-0.1, 1.1])
        ax_w.plot(sequence[elemnum][chan].wave, lw=1, color='#009FFF')
        ax_m.plot(sequence[elemnum][chan].marker_1, lw=1,
                  color='#008B45', alpha=0.6, label='m1')
        ax_m.plot(sequence[elemnum][chan].marker_2, lw=1,
                  color='#FE6447', alpha=0.6, label='m2')
        ax_m.legend(loc='upper right', fontsize=10)
    plt.tight_layout()
    return fig


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
    sr = get_pulse_val('sample_rate')
    if sr != awg.clock_freq():
        awg.clock_freq(sr)
    print('awg clock freq set to {}'.format(sr))
