import pickle
import numpy as np
import pprint
import matplotlib.pyplot as plt
from . import get_pulse_location, get_latest_counter

from . import Sequence, Waveform, Element


default_pulse_dict = {'cycle_duration': 20e-6, 'sampling_rate': 1e9,
                      'pulse_end': 10e-6, 'pulse_readout_delay': 30e-9,
                      'readout_time': 4e-6, 'marker_time': 500e-9,
                      'marker_readout_delay': 0, 'qubit_time': 1e-6}


def validate_pulse_dictionary():
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
    path = get_pulse_location()
    file_name = 'pulse_dict.p'
    pulse_dict = pickle.load(open(path + file_name, "rb"))
    return pulse_dict


def update_pulse_dict(update_dict):
    path = get_pulse_location()
    file_name = 'pulse_dict.p'
    try:
        dict_to_update = pickle.load(open(path + file_name, "rb"))
    except FileNotFoundError:
        dict_to_update = {}

    dict_to_update.update(update_dict)
    pickle.dump(dict_to_update, open(path + file_name, 'wb'))


def update_pulse_val(key, val):
    p_dict = get_pulse_dict()
    p_dict[key] = val
    update_pulse_dict(p_dict)


def get_pulse_val(key):
    p_dict = get_pulse_dict()
    return p_dict[key]


def make_readout_seq(channels=[4]):
    validate_pulse_dictionary()
    readout_sequence = Sequence(name='plain_readout', variable='')
    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sampling_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = p_dict['readout_start'] - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[0])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start:readout_marker_start + marker_points] = 1
    readout_waveform.marker_2[
        readout_marker_start:readout_marker_start + marker_points] = 1

    element = Element()
    element.add_waveform(readout_waveform)

    readout_sequence.add_element(element)
    readout_sequence.check()

    return readout_sequence


def make_ssb_qubit_seq(start=0, stop=200e6, step=1e6, channels=[1, 2, 4]):
    validate_pulse_dictionary()
    ssb_sequence = Sequence(name='qubit_ssb',
                            variable='sideband_modulation_freq',
                            variable_unit='Hz',
                            step=step,
                            start=start,
                            stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sampling_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = p_dict['readout_start'] - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / p_dict['resolution'])
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = round(p_dict['qubit_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[2])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start:readout_marker_start + marker_points] = 1
    readout_waveform.marker_2[
        readout_marker_start:readout_marker_start + marker_points] = 1

    qubit_time_array = np.arange(qubit_points) * resolution
    freq_array = ssb_sequence.variable_array

    for i, freq in enumerate(freq_array):
        element = Element()
        element.add_waveform(readout_waveform)
        qubit_i = Waveform(length=total_points, channel=channels[0])
        qubit_q = Waveform(length=total_points, channel=channels[1])
        if i == 0:
            qubit_i.marker_1[10:10 + marker_points] = 1
            qubit_i.marker_2[10:10 + marker_points] = 1
        qubit_start = pulse_end_points - qubit_points
        qubit_end = pulse_end_points
        angle = qubit_time_array * freq * 2 * np.pi
        cos_array = np.cos(angle)
        sin_array = np.sin(angle)
        qubit_i.wave[qubit_start:qubit_end] = cos_array
        qubit_q.wave[qubit_start:qubit_end] = -1 * sin_array
        element.add_waveform(qubit_i)
        element.add_waveform(qubit_q)
        ssb_sequence.add_element(element)

    ssb_sequence.check()

    return ssb_sequence


def make_t1_seq(pi_duration, pi_amp, start=0, stop=5e-6, step=50e-9, channels=[1, 4]):
    validate_pulse_dictionary()
    t1_sequence = Sequence(name='t1',
                           variable='drive_readout_delay',
                           variable_unit='s',
                           step=step,
                           start=start,
                           stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sampling_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = p_dict['readout_start'] - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / p_dict['resolution'])
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = pi_duration / resolution
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[1])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start:readout_marker_start + marker_points] = 1
    readout_waveform.marker_2[
        readout_marker_start:readout_marker_start + marker_points] = 1

    delay_array_points = np.round(
        t1_sequence.variable_array / resolution).astype(np.int)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        element.add_waveform(readout_waveform)
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        if i == 0:
            qubit_waveform.marker_1[10:10 + marker_points] = 1
            qubit_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = pulse_end_points - delay_points - qubit_points
        qubit_end = qubit_start + qubit_points
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        element.add_waveform(qubit_waveform)
        t1_sequence.add_element(element)

    t1_sequence.check()

    return t1_sequence


def make_rabi_sequence(pi_amp, start=0, stop=200e-9, step=2e-9, channels=[1, 4]):
    validate_pulse_dictionary()
    rabi_sequence = Sequence(name='rabi',
                             variable='drive_duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sampling_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = p_dict['readout_start'] - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / p_dict['resolution'])
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[1])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start:readout_marker_start + marker_points] = 1
    readout_waveform.marker_2[
        readout_marker_start:readout_marker_start + marker_points] = 1

    qubit_duration_array_points = np.round(
        rabi_sequence.variable_array / resolution).astype(int)

    for i, qubit_points in enumerate(qubit_duration_array_points):
        element = Element()
        element.add_waveform(readout_waveform)
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        if i == 0:
            qubit_waveform.marker_1[10:10 + marker_points] = 1
            qubit_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start = pulse_end_points - qubit_points
        qubit_end = pulse_end_points
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        element.add_waveform(qubit_waveform)
        rabi_sequence.add_element(element) / 2

    rabi_sequence.check()

    return rabi_sequence


def make_ramsey_sequence(pi_duration, pi_amp, start=0, stop=200e-9, step=2e-9, channels=[1, 4]):
    validate_pulse_dictionary()
    ramsey_sequence = Sequence(name='ramsey',
                             variable='drive_drive_delay',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

    p_dict = get_pulse_dict()
    resolution = 1 / p_dict['sampling_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = p_dict['readout_start'] - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / p_dict['resolution'])
    marker_points = round(p_dict['marker_time'] / resolution)
    qubit_points = pi_duration / resoultion
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[1])

    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.marker_1[
        readout_marker_start:readout_marker_start + marker_points] = 1
    readout_waveform.marker_2[
        readout_marker_start:readout_marker_start + marker_points] = 1

    delay_array_points = np.round(
        ramsey_sequence.variable_array / resolution).astype(np.int)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        element.add_waveform(readout_waveform)
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        if i == 0:
            qubit_waveform.marker_1[10:10 + marker_points] = 1
            qubit_waveform.marker_2[10:10 + marker_points] = 1
        qubit_start_first = pulse_end_points - delay_points - 2 * qubit_points
        qubit_end_first = qubit_start_first + qubit_points
        qubit_start_second = pulse_end_points - delay_points - qubit_points
        qubit_end_second = qubit_start_second + qubit_points
        qubit_waveform.wave[qubit_start_first:qubit_end_first] = pi_amp / 2
        qubit_waveform.wave[qubit_start_second:qubit_end_second] = pi_amp / 2
        element.add_waveform(qubit_waveform)
        ramsey_sequence.add_element(element)

    ramsey_sequence.check()

    return ramsey_sequence

def plot_sequence(sequence, elemnum=0, chan_a=1, chan_b=2):
    fig = plt.figure()

    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412)
    ax3 = fig.add_subplot(413)
    ax4 = fig.add_subplot(414)
    ax1.set_title('Channel {} waveform'.format(chan_a))
    ax1.set_ylim([-1.1, 1.1])
    ax2.set_title('Channel {} markers'.format(chan_a))
    ax2.set_ylim([-0.1, 1.1])
    ax3.set_title('Channel {} waveform'.format(chan_b))
    ax3.set_ylim([-1.1, 1.1])
    ax4.set_title('Channel {} markers'.format(chan_b))
    ax4.set_ylim([-0.1, 1.1])


    ax1.plot(sequence[elemnum][chan_a].wave, lw=1, color='#e1cb66')
    ax2.plot(sequence[elemnum][chan_a].marker_1, lw=1, color='#FF4500', alpha=0.6)
    ax2.plot(sequence[elemnum][chan_a].marker_2, lw=1, color='#FF8C00', alpha=0.6)

    ax3.plot(sequence[elemnum][chan_b].wave, lw=1, color='#6689e1')
    ax4.plot(sequence[elemnum][chan_b].marker_1, lw=1, color='#6A5ACD', alpha=0.6)
    ax4.plot(sequence[elemnum][chan_b].marker_2, lw=1, color='#EE82EE', alpha=0.6)

    plt.tight_layout()
    return fig


def make_save_send_load_awg_file(awg, sequence):
    pulse_location = get_pulse_location()
    try:
        num = get_latest_counter(pulse_location) + 1
    except FileNotFoundError:
        num = 1
    name = '{0:03d}_{}.awg'.format(num, sequence.name)
    file_location = pulse_location + name
    awg.make_and_save_awg_file(*sequence.unwrap(), filename=file_location)
    awg.make_send_and_load_awg_file(*sequence.unwrap())


def check_sample_rate(awg):
    sr = get_pulse_val('sample_rate')
    if sr != awg.clock_freq():
        awg.clock_freq(sr)
    print('awg clock freq set to {}'.format(sr))
