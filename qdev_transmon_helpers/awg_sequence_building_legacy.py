from . import get_calibration_dict, gaussian, SSB_Q_gaussian, SSB_I_gaussian
from . import Segment, Waveform, Element, Sequence
import numpy as np


def make_readout_seq_points(channels=[4]):
    """
    Square pulse duting readout time with one marker at start
    """
    readout_sequence = Sequence(name='plain_readout', variable='')
    p_dict = get_calibration_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    readout_waveform = Waveform(length=total_points, channel=channels[0])
    readout_waveform.wave[readout_start_points:
                          readout_start_points + readout_points] = 1
    readout_waveform.add_marker(1, readout_marker_start_points, marker_points)

    element = Element()
    element.add_waveform(readout_waveform)
    readout_sequence.add_element(element)
    readout_sequence.check()
    return readout_sequence


################################################################
# Spectroscopy
################################################################


def make_full_ssb_wave_points(freq=8e6, duration=20e-6, channels=[1, 2]):
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


def make_ssb_qubit_seq(start=0, stop=200e6, step=1e6, channels=[1, 2, 4],
                       pulse_mod=False):
    """
    Cosine and sine waves for qubit time with range of frequencies, square
    readout wave for readout time. Markers on readout channel (1 for readout
    start, 2 for seq start)
    """
    ssb_sequence = Sequence(name='qubit_ssb',
                            variable='ssb_qubit_modulation_freq',
                            variable_unit='Hz',
                            step=step,
                            start=start,
                            stop=stop)

    p_dict = get_calibration_dict()
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
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    qubit_time_array = np.arange(qubit_points) * resolution
    freq_array = ssb_sequence.variable_array
    pulse_mod_points = qubit_points * 2

    for i, freq in enumerate(freq_array):
        element = Element()
        qubit_i = Waveform(length=total_points, channel=channels[0])
        qubit_q = Waveform(length=total_points, channel=channels[1])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = pulse_end_points - qubit_points
        qubit_end = pulse_end_points
        angle = qubit_time_array * freq * 2 * np.pi
        cos_array = np.cos(angle)
        sin_array = np.sin(angle)
        qubit_i.wave[qubit_start:qubit_end] = cos_array
        qubit_q.wave[qubit_start:qubit_end] = -1 * sin_array
        if pulse_mod:
            qubit_i.add_marker(1, pulse_end_points - pulse_mod_points,
                               pulse_mod_points)
            qubit_q.add_marker(1, pulse_end_points - pulse_mod_points,
                               pulse_mod_points)
        element.add_waveform(qubit_i)
        element.add_waveform(qubit_q)
        element.add_waveform(readout_waveform)
        ssb_sequence.add_element(element)
    ssb_sequence.check()
    return ssb_sequence


################################################################
# Rabi and T1
################################################################


def make_rabi_square_sequence_points(pi_amp=1, start=0, stop=200e-9, step=2e-9,
                                     channels=[1, 4], pulse_mod=False):
    """
    Square qubit drive of pi amplitude of varying duration, square readout
    drive. Markers on readout channel (1 for readout start, 2 for seq start)
    """
    rabi_sequence = Sequence(name='rabi',
                             variable='drive_duration',
                             variable_label='Drive Duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

    p_dict = get_calibration_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_time'] / resolution)

    pulse_mod_points = int(stop * 2 / resolution)

    readout_template = Waveform(length=total_points, channel=channels[1])

    wave = readout_template.wave
    wave[readout_start_points:readout_start_points + readout_points] = 1
    readout_template.wave = wave

    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    qubit_duration_array_points = np.round(
        rabi_sequence.variable_array / resolution).astype(int)

    for i, qubit_points in enumerate(qubit_duration_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = int(pulse_end_points - qubit_points)
        qubit_end = int(pulse_end_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        rabi_sequence.add_element(element)

    rabi_sequence.check()

    return rabi_sequence


def make_rabi_gaussian_sequence_points(sigma_cuttoff, pi_amp=1, start=0,
                                       stop=200e-9, step=2e-9,
                                       channels=[1, 4], pulse_mod=False):
    """
    Square qubit drive of pi amplitude of varying duration, square readout
    drive. Markers on readout channel (1 for readout start, 2 for seq start)
    """
    rabi_sequence = Sequence(name='rabi',
                             variable='gaussian_drive_duration',
                             variable_label='Drive Duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

    p_dict = get_calibration_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    pulse_mod_points = int(sigma_cuttoff * 4 * stop / resolution)

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    for i, pi_duration in enumerate(rabi_sequence.variable_array):
        pi_pulse = gaussian(pi_duration, sigma_cuttoff, pi_amp,
                            p_dict['sample_rate'])
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = int(pulse_end_points - len(pi_pulse))
        qubit_end = int(pulse_end_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_pulse
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        rabi_sequence.add_element(element)

    rabi_sequence.check()

    return rabi_sequence


def make_rabi_gaussianSSB_sequence_points(sigma_cuttoff, pi_amp=1, start=0,
                                          stop=200e-9, step=2e-9,
                                          SSBfreq=100e6, channels=[1, 2, 4],
                                          pulse_mod=True):
    """
    Square qubit drive of pi amplitude of varying duration, square readout
    drive. Markers on readout channel (1 for readout start, 2 for seq start)
    """
    rabi_sequence = Sequence(name='rabi',
                             variable='gaussian_drive_duration',
                             variable_label='Drive Duration',
                             variable_unit='s',
                             step=step,
                             start=start,
                             stop=stop)

    p_dict = get_calibration_dict()
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    pulse_mod_points = int(sigma_cuttoff * 4 * stop / resolution)

    readout_template = Waveform(length=total_points, channel=channels[2])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    for i, pi_duration in enumerate(rabi_sequence.variable_array):
        pi_pulseI = SSB_I_gaussian(pi_duration, sigma_cuttoff,
                                   SSBfreq, pi_amp, p_dict['sample_rate'])
        pi_pulseQ = SSB_Q_gaussian(pi_duration, sigma_cuttoff,
                                   SSBfreq, pi_amp, p_dict['sample_rate'])
        element = Element()
        qubit_waveformI = Waveform(length=total_points, channel=channels[0])
        qubit_waveformQ = Waveform(length=total_points, channel=channels[1])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = int(pulse_end_points - len(pi_pulseI))
        qubit_end = int(pulse_end_points)
        qubit_waveformI.wave[qubit_start:qubit_end] = pi_pulseI
        qubit_waveformQ.wave[qubit_start:qubit_end] = pi_pulseQ
        if pulse_mod:
            qubit_waveformI.add_marker(1, pulse_end_points - pulse_mod_points,
                                       pulse_mod_points)
        element.add_waveform(qubit_waveformI)
        element.add_waveform(qubit_waveformQ)
        element.add_waveform(readout_waveform)
        rabi_sequence.add_element(element)

    rabi_sequence.check()

    return rabi_sequence


def make_t1_square_seq_points(pi_duration, pi_amp=1, start=0, stop=5e-6,
                              step=50e-9, channels=[1, 4], pulse_mod=False):
    """
    Square qubit drive for pi duration at pi amplitude on qubit with varying
    wait time before readout (square pulse for readout time). Markers on
    readout channel (1 for readout start, 2 for seq start)
    """
    t1_sequence = Sequence(name='t1',
                           variable='drive_readout_delay',
                           variable_label='Delay',
                           variable_unit='s',
                           step=step,
                           start=start,
                           stop=stop)

    p_dict = get_calibration_dict()
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
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    delay_array_points = np.round(
        t1_sequence.variable_array / resolution).astype(np.int)

    pulse_mod_points = int((qubit_points * 4 + stop) / resolution)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = int(pulse_end_points - delay_points - qubit_points)
        qubit_end = int(qubit_start + qubit_points)
        qubit_waveform.wave[qubit_start:qubit_end] = pi_amp
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        t1_sequence.add_element(element)
    t1_sequence.check()
    return t1_sequence


def make_t1_gaussian_seq_points(pi_duration, sigma_cuttoff, pi_amp=1, start=0,
                                stop=5e-6, step=50e-9, channels=[1, 4],
                                pulse_mod=False):
    """
    Gaussian qubit drive for pi duration at pi amplitude on qubit with varying
    wait time before readout (square pulse for readout time). Markers on
    readout channel (1 for readout start, 2 for seq start)
    """
    t1_sequence = Sequence(name='t1',
                           variable='drive_readout_delay',
                           variable_label='Delay',
                           variable_unit='s',
                           step=step,
                           start=start,
                           stop=stop)

    p_dict = get_calibration_dict()
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
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    delay_array_points = np.round(
        t1_sequence.variable_array / resolution).astype(np.int)

    pulse_mod_points = int((sigma_cuttoff * 4 + stop) / resolution)

    for i, delay_points in enumerate(delay_array_points):
        pi_pulse = gaussian(pi_duration, sigma_cuttoff, pi_amp,
                            p_dict['sample_rate'])
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start = int(pulse_end_points - delay_points - len(pi_pulse))
        qubit_end = int(qubit_start + len(pi_pulse))
        qubit_waveform.wave[qubit_start:qubit_end] = pi_pulse
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        t1_sequence.add_element(element)
    t1_sequence.check()
    return t1_sequence


################################################################
# Ramsey and T2
################################################################

def make_ramsey_square_sequence_points(pi_duration, pi_half_amp=1 / 2, start=0,
                                       stop=200e-9, step=2e-9, channels=[1, 4],
                                       pulse_mod=False):
    """
    Two square pulses on qubit of pi duration of half pi amplitude separated
    by varying duration. Square readout with markers (1 for readout start,
    2 for seq start)
    """
    ramsey_sequence = Sequence(name='ramsey',
                               variable='drive_drive_delay',
                               variable_label='Delay',
                               variable_unit='s',
                               step=step,
                               start=start,
                               stop=stop)

    p_dict = get_calibration_dict()
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
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    delay_array_points = np.round(
        ramsey_sequence.variable_array / resolution).astype(np.int)

    pulse_mod_points = int((stop + 4 * pi_duration) / resolution)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start_first = int(
            pulse_end_points - delay_points - 2 * qubit_points)
        qubit_end_first = int(qubit_start_first + qubit_points)
        qubit_start_second = int(pulse_end_points - qubit_points)
        qubit_end_second = int(qubit_start_second + qubit_points)
        qubit_waveform.wave[qubit_start_first:qubit_end_first] = pi_half_amp
        qubit_waveform.wave[qubit_start_second:qubit_end_second] = pi_half_amp
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(readout_waveform)
        ramsey_sequence.add_element(element)
    ramsey_sequence.check()
    return ramsey_sequence


def make_ramsey_gaussian_sequence(pi_duration, sigma_cuttoff,
                                  pi_half_amp=1 / 2, start=0, stop=200e-9,
                                  step=2e-9, channels=[1, 4], pulse_mod=False):
    """
    Two Gaussian pulses on qubit of pi duration of half pi amplitude separated
    by varying duration. Square readout with markers (1 for readout start,
    2 for seq start)
    """
    ramsey_sequence = Sequence(name='ramsey',
                               variable='drive_drive_delay',
                               variable_label='Delay',
                               variable_unit='s',
                               step=step,
                               start=start,
                               stop=stop)

    p_dict = get_calibration_dict()
    pi_half_amp = p_dict['pi_half_pulse_amp']
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
    readout_template.add_marker(1, readout_marker_start_points, marker_points)

    delay_array_points = np.round(
        ramsey_sequence.variable_array / resolution).astype(np.int)

    pi_half_pulse = gaussian(pi_duration, sigma_cuttoff, pi_half_amp,
                             p_dict['sample_rate'])

    pulse_mod_points = int((stop + 4 * pi_duration * sigma_cuttoff) /
                           resolution)

    for i, delay_points in enumerate(delay_array_points):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        if i == 0:
            readout_waveform.add_marker(2, 0, marker_points)
        qubit_start_first = int(pulse_end_points - delay_points - 2 *
                                len(pi_half_pulse))
        qubit_end_first = int(qubit_start_first + len(pi_half_pulse))
        qubit_start_second = int(pulse_end_points - len(pi_half_pulse))
        qubit_end_second = int(qubit_start_second + len(pi_half_pulse))
        qubit_waveform.wave[qubit_start_first:qubit_end_first] = pi_half_pulse
        qubit_waveform.wave[
            qubit_start_second:qubit_end_second] = pi_half_pulse
        if pulse_mod:
            qubit_waveform.add_marker(1, pulse_end_points - pulse_mod_points,
                                      pulse_mod_points)
        element.add_waveform(qubit_waveform)
        element.add_waveform(readout_waveform)
        ramsey_sequence.add_element(element)
    ramsey_sequence.check()
    return ramsey_sequence


##################################################################
# AllXY
##################################################################


def make_allxy_seq(pi_duration, sigma_cuttoff, channels=[1, 2, 4],
                   pulse_mod=False):
    """
    Oh dear.
    """
    p_dict = get_calibration_dict()
    pi_amp = p_dict['pi_pulse_amp']
    pi_half_amp = p_dict['pi_half_pulse_amp']
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    pi_pulse = gaussian(pi_duration, sigma_cuttoff, pi_amp,
                        p_dict['sample_rate'])
    pi_half_pulse = gaussian(pi_duration, sigma_cuttoff, pi_half_amp,
                             p_dict['sample_rate'])

    readout_waveform = Waveform(length=total_points, channel=channels[2])
    readout_waveform.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_waveform.add_marker(1, readout_marker_start_points, marker_points)

    pulse_points = len(pi_pulse)
    pulse_mod_points = int(pi_duration * sigma_cuttoff * 2 / resolution)

    seq = Sequence(name='allxy', variable='operation_combination',
                   variable_label='Operation Combination Id')

    elem_0 = Element()
    x_waveform_0 = Waveform(length=total_points, channel=1)
    y_waveform_0 = Waveform(length=total_points, channel=2)
    readout_first = readout_waveform.copy()
    readout_first.add_marker(2, 0, marker_points)
    if pulse_mod:
        x_waveform_0.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_0.add_waveform(x_waveform_0)
    elem_0.add_waveform(y_waveform_0)
    elem_0.add_waveform(readout_first)
    seq.add_element(elem_0)

    elem_1 = Element()
    x_waveform_1 = Waveform(length=total_points, channel=channels[0])
    y_waveform_1 = Waveform(length=total_points, channel=channels[1])
    x_waveform_1.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_pulse
    x_waveform_1.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_1.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_1.add_waveform(x_waveform_1)
    elem_1.add_waveform(y_waveform_1)
    elem_1.add_waveform(readout_first)
    seq.add_element(elem_1)

    elem_2 = Element()
    x_waveform_2 = Waveform(length=total_points, channel=channels[0])
    y_waveform_2 = Waveform(length=total_points, channel=channels[1])
    y_waveform_2.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_pulse
    y_waveform_2.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_2.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_2.add_waveform(x_waveform_2)
    elem_2.add_waveform(y_waveform_2)
    elem_2.add_waveform(readout_waveform)
    seq.add_element(elem_2)

    elem_3 = Element()
    x_waveform_3 = Waveform(length=total_points, channel=channels[0])
    y_waveform_3 = Waveform(length=total_points, channel=channels[1])
    x_waveform_3.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_pulse
    y_waveform_3.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_3.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_3.add_waveform(x_waveform_3)
    elem_3.add_waveform(y_waveform_3)
    elem_3.add_waveform(readout_waveform)
    seq.add_element(elem_3)

    elem_4 = Element()
    x_waveform_4 = Waveform(length=total_points, channel=channels[0])
    y_waveform_4 = Waveform(length=total_points, channel=channels[1])
    y_waveform_4.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_pulse
    x_waveform_4.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_4.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_4.add_waveform(x_waveform_4)
    elem_4.add_waveform(y_waveform_4)
    elem_4.add_waveform(readout_waveform)
    seq.add_element(elem_4)

    elem_5 = Element()
    x_waveform_5 = Waveform(length=total_points, channel=channels[0])
    y_waveform_5 = Waveform(length=total_points, channel=channels[1])
    x_waveform_5.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_5.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_5.add_waveform(x_waveform_5)
    elem_5.add_waveform(y_waveform_5)
    elem_5.add_waveform(readout_waveform)
    seq.add_element(elem_5)

    elem_6 = Element()
    x_waveform_6 = Waveform(length=total_points, channel=channels[0])
    y_waveform_6 = Waveform(length=total_points, channel=channels[1])
    y_waveform_6.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_6.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_6.add_waveform(x_waveform_6)
    elem_6.add_waveform(y_waveform_6)
    elem_6.add_waveform(readout_waveform)
    seq.add_element(elem_6)

    elem_7 = Element()
    x_waveform_7 = Waveform(length=total_points, channel=channels[0])
    y_waveform_7 = Waveform(length=total_points, channel=channels[1])
    x_waveform_7.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_half_pulse
    y_waveform_7.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_7.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_7.add_waveform(x_waveform_7)
    elem_7.add_waveform(y_waveform_7)
    elem_7.add_waveform(readout_waveform)
    seq.add_element(elem_7)

    elem_8 = Element()
    x_waveform_8 = Waveform(length=total_points, channel=channels[0])
    y_waveform_8 = Waveform(length=total_points, channel=channels[1])
    y_waveform_8.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_half_pulse
    x_waveform_8.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_8.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_8.add_waveform(x_waveform_8)
    elem_8.add_waveform(y_waveform_8)
    elem_8.add_waveform(readout_waveform)
    seq.add_element(elem_8)

    elem_9 = Element()
    x_waveform_9 = Waveform(length=total_points, channel=channels[0])
    y_waveform_9 = Waveform(length=total_points, channel=channels[1])
    x_waveform_9.wave[pulse_end_points - 2 * pulse_points:
                      pulse_end_points - pulse_points] = pi_half_pulse
    y_waveform_9.wave[pulse_end_points - pulse_points:
                      pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_9.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_9.add_waveform(x_waveform_9)
    elem_9.add_waveform(y_waveform_9)
    elem_9.add_waveform(readout_waveform)
    seq.add_element(elem_9)

    elem_10 = Element()
    x_waveform_10 = Waveform(length=total_points, channel=channels[0])
    y_waveform_10 = Waveform(length=total_points, channel=channels[1])
    y_waveform_10.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulse
    x_waveform_10.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_10.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_10.add_waveform(x_waveform_10)
    elem_10.add_waveform(y_waveform_10)
    elem_10.add_waveform(readout_waveform)
    seq.add_element(elem_10)

    elem_11 = Element()
    x_waveform_11 = Waveform(length=total_points, channel=channels[0])
    y_waveform_11 = Waveform(length=total_points, channel=channels[1])
    x_waveform_11.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    y_waveform_11.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_11.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_11.add_waveform(x_waveform_11)
    elem_11.add_waveform(y_waveform_11)
    elem_11.add_waveform(readout_waveform)
    seq.add_element(elem_11)

    elem_12 = Element()
    x_waveform_12 = Waveform(length=total_points, channel=channels[0])
    y_waveform_12 = Waveform(length=total_points, channel=channels[1])
    y_waveform_12.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    x_waveform_12.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_12.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_12.add_waveform(x_waveform_12)
    elem_12.add_waveform(y_waveform_12)
    elem_12.add_waveform(readout_waveform)
    seq.add_element(elem_12)

    elem_13 = Element()
    x_waveform_13 = Waveform(length=total_points, channel=channels[0])
    y_waveform_13 = Waveform(length=total_points, channel=channels[1])
    x_waveform_13.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulse
    x_waveform_13.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_13.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_13.add_waveform(x_waveform_13)
    elem_13.add_waveform(y_waveform_13)
    elem_13.add_waveform(readout_waveform)
    seq.add_element(elem_13)

    elem_14 = Element()
    x_waveform_14 = Waveform(length=total_points, channel=channels[0])
    y_waveform_14 = Waveform(length=total_points, channel=channels[1])
    x_waveform_14.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    x_waveform_14.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_14.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_14.add_waveform(x_waveform_14)
    elem_14.add_waveform(y_waveform_14)
    elem_14.add_waveform(readout_waveform)
    seq.add_element(elem_14)

    elem_15 = Element()
    x_waveform_15 = Waveform(length=total_points, channel=channels[0])
    y_waveform_15 = Waveform(length=total_points, channel=channels[1])
    y_waveform_15.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulse
    y_waveform_15.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_pulse
    if pulse_mod:
        x_waveform_15.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_15.add_waveform(x_waveform_15)
    elem_15.add_waveform(y_waveform_15)
    elem_15.add_waveform(readout_waveform)
    seq.add_element(elem_15)

    elem_16 = Element()
    x_waveform_16 = Waveform(length=total_points, channel=channels[0])
    y_waveform_16 = Waveform(length=total_points, channel=channels[1])
    y_waveform_16.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    y_waveform_16.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_16.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_16.add_waveform(x_waveform_16)
    elem_16.add_waveform(y_waveform_16)
    elem_16.add_waveform(readout_waveform)
    seq.add_element(elem_16)

    elem_17 = Element()
    x_waveform_17 = Waveform(length=total_points, channel=channels[0])
    y_waveform_17 = Waveform(length=total_points, channel=channels[1])
    x_waveform_17.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    if pulse_mod:
        x_waveform_17.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_17.add_waveform(x_waveform_17)
    elem_17.add_waveform(y_waveform_17)
    elem_17.add_waveform(readout_waveform)
    seq.add_element(elem_17)

    elem_18 = Element()
    x_waveform_18 = Waveform(length=total_points, channel=channels[0])
    y_waveform_18 = Waveform(length=total_points, channel=channels[1])
    y_waveform_18.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_pulse
    if pulse_mod:
        x_waveform_18.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_18.add_waveform(x_waveform_18)
    elem_18.add_waveform(y_waveform_18)
    elem_18.add_waveform(readout_waveform)
    seq.add_element(elem_18)

    elem_19 = Element()
    x_waveform_19 = Waveform(length=total_points, channel=channels[0])
    y_waveform_19 = Waveform(length=total_points, channel=channels[1])
    x_waveform_19.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulse
    x_waveform_19.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_19.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_19.add_waveform(x_waveform_19)
    elem_19.add_waveform(y_waveform_19)
    elem_19.add_waveform(readout_waveform)
    seq.add_element(elem_19)

    elem_20 = Element()
    x_waveform_20 = Waveform(length=total_points, channel=channels[0])
    y_waveform_20 = Waveform(length=total_points, channel=channels[1])
    y_waveform_20.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulse
    y_waveform_20.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulse
    if pulse_mod:
        x_waveform_20.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_20.add_waveform(x_waveform_20)
    elem_20.add_waveform(y_waveform_20)
    elem_20.add_waveform(readout_waveform)
    seq.add_element(elem_20)

    return seq


def make_allxySSB_seq(pi_duration, pi_amp, SSBfreq, sigma_cuttoff,
                      channels=[1, 2, 4], pulse_mod=False):
    """
    Oh dear. Part 2
    """
    p_dict = get_calibration_dict()

    pi_amp = p_dict['pi_pulse_amp']
    pi_half_amp = p_dict['pi_half_pulse_amp']
    resolution = 1 / p_dict['sample_rate']
    readout_start = p_dict['pulse_end'] + p_dict['pulse_readout_delay']
    readout_marker_start = readout_start - p_dict['marker_readout_delay']
    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = round(p_dict['readout_time'] / resolution)
    pulse_end_points = round(p_dict['pulse_end'] / resolution)
    marker_points = round(p_dict['marker_time'] / resolution)
    total_points = round(p_dict['cycle_duration'] / resolution)

    pi_pulseI_x = SSB_I_gaussian(pi_duration, sigma_cuttoff, SSBfreq,
                                 pi_amp, p_dict['sample_rate'])
    pi_pulseQ_x = SSB_Q_gaussian(pi_duration, sigma_cuttoff, SSBfreq,
                                 pi_amp, p_dict['sample_rate'])

    pi_half_pulseI_x = SSB_I_gaussian(pi_duration, sigma_cuttoff, SSBfreq,
                                      pi_half_amp, p_dict['sample_rate'])
    pi_half_pulseQ_x = SSB_Q_gaussian(pi_duration, sigma_cuttoff, SSBfreq,
                                      pi_half_amp, p_dict['sample_rate'])

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
        x_waveform_0.add_marker(1, pulse_end_points - pulse_mod_points,
                                pulse_mod_points)
    elem_0.add_waveform(x_waveform_0)
    elem_0.add_waveform(y_waveform_0)
    elem_0.add_waveform(readout_first)
    seq.add_element(elem_0)

    # elem_1 = Element()
    # x_waveform_1 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_1 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_1.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    # y_waveform_1.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    # x_waveform_1.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    # y_waveform_1.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_x
    # if pulse_mod:
    #     x_waveform_1.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_1.add_waveform(x_waveform_1)
    # elem_1.add_waveform(y_waveform_1)
    # elem_1.add_waveform(readout_first)
    # seq.add_element(elem_1)

    # elem_2 = Element()
    # x_waveform_2 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_2 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_2.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    # y_waveform_2.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    # x_waveform_2.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # y_waveform_2.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_y
    # if pulse_mod:
    #     x_waveform_2.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_2.add_waveform(x_waveform_2)
    # elem_2.add_waveform(y_waveform_2)
    # elem_2.add_waveform(readout_waveform)
    # seq.add_element(elem_2)

    # elem_3 = Element()
    # x_waveform_3 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_3 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_3.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    # y_waveform_3.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    # x_waveform_3.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # y_waveform_3.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_y
    # if pulse_mod:
    #     x_waveform_3.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_3.add_waveform(x_waveform_3)
    # elem_3.add_waveform(y_waveform_3)
    # elem_3.add_waveform(readout_waveform)
    # seq.add_element(elem_3)

    # elem_4 = Element()
    # x_waveform_4 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_4 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_4.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    # y_waveform_4.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    # x_waveform_4.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    # y_waveform_4.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseQ_x
    # if pulse_mod:
    #     x_waveform_4.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_4.add_waveform(x_waveform_4)
    # elem_4.add_waveform(y_waveform_4)
    # elem_4.add_waveform(readout_waveform)
    # seq.add_element(elem_4)

    # elem_5 = Element()
    # x_waveform_5 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_5 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_5.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # y_waveform_5.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    # if pulse_mod:
    #     x_waveform_5.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_5.add_waveform(x_waveform_5)
    # elem_5.add_waveform(y_waveform_5)
    # elem_5.add_waveform(readout_waveform)
    # seq.add_element(elem_5)

    # elem_6 = Element()
    # x_waveform_6 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_6 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_6.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    # y_waveform_6.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    # if pulse_mod:
    #     x_waveform_6.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_6.add_waveform(x_waveform_6)
    # elem_6.add_waveform(y_waveform_6)
    # elem_6.add_waveform(readout_waveform)
    # seq.add_element(elem_6)

    # elem_7 = Element()
    # x_waveform_7 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_7 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_7.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # y_waveform_7.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    # x_waveform_7.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # y_waveform_7.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # if pulse_mod:
    #     x_waveform_7.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_7.add_waveform(x_waveform_7)
    # elem_7.add_waveform(y_waveform_7)
    # elem_7.add_waveform(readout_waveform)
    # seq.add_element(elem_7)

    # elem_8 = Element()
    # x_waveform_8 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_8 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_8.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    # y_waveform_8.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    # x_waveform_8.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # y_waveform_8.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # if pulse_mod:
    #     x_waveform_8.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_8.add_waveform(x_waveform_8)
    # elem_8.add_waveform(y_waveform_8)
    # elem_8.add_waveform(readout_waveform)
    # seq.add_element(elem_8)

    # elem_9 = Element()
    # x_waveform_9 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_9 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_9.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # y_waveform_9.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    # x_waveform_9.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # y_waveform_9.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # if pulse_mod:
    #     x_waveform_9.add_marker(1, pulse_end_points - pulse_mod_points,
    #                             pulse_mod_points)
    # elem_9.add_waveform(x_waveform_9)
    # elem_9.add_waveform(y_waveform_9)
    # elem_9.add_waveform(readout_waveform)
    # seq.add_element(elem_9)

    # elem_10 = Element()
    # x_waveform_10 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_10 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_10.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    # y_waveform_10.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    # x_waveform_10.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # y_waveform_10.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # if pulse_mod:
    #     x_waveform_10.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_10.add_waveform(x_waveform_10)
    # elem_10.add_waveform(y_waveform_10)
    # elem_10.add_waveform(readout_waveform)
    # seq.add_element(elem_10)

    # elem_11 = Element()
    # x_waveform_11 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_11 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_11.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    # y_waveform_11.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    # x_waveform_11.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # y_waveform_11.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # if pulse_mod:
    #     x_waveform_11.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_11.add_waveform(x_waveform_11)
    # elem_11.add_waveform(y_waveform_11)
    # elem_11.add_waveform(readout_waveform)
    # seq.add_element(elem_11)

    # elem_12 = Element()
    # x_waveform_12 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_12 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_12.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    # y_waveform_12.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    # x_waveform_12.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # y_waveform_12.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # if pulse_mod:
    #     x_waveform_12.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_12.add_waveform(x_waveform_12)
    # elem_12.add_waveform(y_waveform_12)
    # elem_12.add_waveform(readout_waveform)
    # seq.add_element(elem_12)

    # elem_13 = Element()
    # x_waveform_13 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_13 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_13.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # y_waveform_13.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_x
    # x_waveform_13.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    # y_waveform_13.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_x
    # if pulse_mod:
    #     x_waveform_13.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_13.add_waveform(x_waveform_13)
    # elem_13.add_waveform(y_waveform_13)
    # elem_13.add_waveform(readout_waveform)
    # seq.add_element(elem_13)

    # elem_14 = Element()
    # x_waveform_14 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_14 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_14.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    # y_waveform_14.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    # x_waveform_14.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # y_waveform_14.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # if pulse_mod:
    #     x_waveform_14.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_14.add_waveform(x_waveform_14)
    # elem_14.add_waveform(y_waveform_14)
    # elem_14.add_waveform(readout_waveform)
    # seq.add_element(elem_14)

    # elem_15 = Element()
    # x_waveform_15 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_15 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_15.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_y
    # y_waveform_15.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseQ_y
    # x_waveform_15.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # y_waveform_15.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_pulseI_y
    # if pulse_mod:
    #     x_waveform_15.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_15.add_waveform(x_waveform_15)
    # elem_15.add_waveform(y_waveform_15)
    # elem_15.add_waveform(readout_waveform)
    # seq.add_element(elem_15)

    # elem_16 = Element()
    # x_waveform_16 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_16 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_16.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    # y_waveform_16.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    # x_waveform_16.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # y_waveform_16.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_y
    # if pulse_mod:
    #     x_waveform_16.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_16.add_waveform(x_waveform_16)
    # elem_16.add_waveform(y_waveform_16)
    # elem_16.add_waveform(readout_waveform)
    # seq.add_element(elem_16)

    # elem_17 = Element()
    # x_waveform_17 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_17 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_17.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_x
    # y_waveform_17.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_x
    # if pulse_mod:
    #     x_waveform_17.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_17.add_waveform(x_waveform_17)
    # elem_17.add_waveform(y_waveform_17)
    # elem_17.add_waveform(readout_waveform)
    # seq.add_element(elem_17)

    # elem_18 = Element()
    # x_waveform_18 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_18 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_18.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseI_y
    # y_waveform_18.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_pulseQ_y
    # if pulse_mod:
    #     x_waveform_18.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_18.add_waveform(x_waveform_18)
    # elem_18.add_waveform(y_waveform_18)
    # elem_18.add_waveform(readout_waveform)
    # seq.add_element(elem_18)

    # elem_19 = Element()
    # x_waveform_19 = Waveform(length=total_points, channel=channels[0])
    # y_waveform_19 = Waveform(length=total_points, channel=channels[1])
    # x_waveform_19.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # y_waveform_19.wave[pulse_end_points-2*pulse_points:pulse_end_points-pulse_points] = pi_half_pulseI_x
    # x_waveform_19.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # y_waveform_19.wave[pulse_end_points-pulse_points:pulse_end_points] = pi_half_pulseI_x
    # if pulse_mod:
    #     x_waveform_19.add_marker(1, pulse_end_points - pulse_mod_points,
    #                              pulse_mod_points)
    # elem_19.add_waveform(x_waveform_19)
    # elem_19.add_waveform(y_waveform_19)
    # elem_19.add_waveform(readout_waveform)
    # seq.add_element(elem_19)

    elem_20 = Element()
    x_waveform_20 = Waveform(length=total_points, channel=channels[0])
    y_waveform_20 = Waveform(length=total_points, channel=channels[1])
    x_waveform_20.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulseI_y
    y_waveform_20.wave[pulse_end_points - 2 * pulse_points:
                       pulse_end_points - pulse_points] = pi_half_pulseI_y
    x_waveform_20.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulseI_y
    y_waveform_20.wave[pulse_end_points - pulse_points:
                       pulse_end_points] = pi_half_pulseI_y
    if pulse_mod:
        x_waveform_20.add_marker(1, pulse_end_points - pulse_mod_points,
                                 pulse_mod_points)
    elem_20.add_waveform(x_waveform_20)
    elem_20.add_waveform(y_waveform_20)
    elem_20.add_waveform(readout_waveform)
    seq.add_element(elem_20)

    return seq
