from . import get_calibration_dict
from . import Waveform, Element, Sequence
import numpy as np


def make_floquet_dur_sequence(floquet_freq, start, stop, step,
                              channels=[1, 4]):
    """
    Cosine modulated qubit drive at floquet frequency for varying durations at
    floquet amplitude. Square readout pulse with markers (1 for readout start,
    2 for seq start)
    """
    # validate_pulse_dictionary()
    duration_floquet_sequence = Sequence(name='floquet',
                                         variable='floquet_dur',
                                         variable_label='time',
                                         variable_unit='S',
                                         step=step,
                                         start=start,
                                         stop=stop)
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    resolution = 1 / calib_dict['sample_rate'][qubit]
    readout_start = calib_dict['pulse_end'][qubit] + calib_dict['pulse_readout_delay'][qubit]
    readout_marker_start = readout_start - calib_dict['marker_readout_delay'][qubit]

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(calib_dict['readout_time'][qubit] / resolution))
    pulse_end_points = int(np.ceil(calib_dict['pulse_end'][qubit] / resolution))
    marker_points = int(np.ceil(calib_dict['marker_time'][qubit] / resolution))
    total_points = int(np.ceil(calib_dict['cycle_time'][qubit] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    # readout_template.marker_1[
    #     readout_marker_start_points:
    #     readout_marker_start_points + marker_points] = 1

    floquet_durations = duration_floquet_sequence.variable_array
    for i, floquet_duration in enumerate(floquet_durations):
        element = Element()
        qubit_waveform = Waveform(length=total_points, channel=channels[0])
        readout_waveform = readout_template.copy()
        # if i == 0:
        #     readout_waveform.marker_2[10:10 + marker_points] = 1

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
    # validate_pulse_dictionary()
    frequency_floquet_sequence = Sequence(name='floquet',
                                          variable='floquet_freq',
                                          variable_label='floquet_freq',
                                          variable_unit='Hz',
                                          step=step,
                                          start=start,
                                          stop=stop)
    calib_dict = get_calibration_dict()
    resolution = 1 / calib_dict['sample_rate']
    readout_start = calib_dict['pulse_end'] + calib_dict['pulse_readout_delay']
    readout_marker_start = readout_start - calib_dict['marker_readout_delay']

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(calib_dict['readout_time'] / resolution))
    pulse_end_points = int(np.ceil(calib_dict['pulse_end'] / resolution))
    marker_points = int(np.ceil(calib_dict['marker_time'] / resolution))
    floquet_points = int(np.ceil(floquet_dur / resolution))
    total_points = int(np.ceil(calib_dict['cycle_duration'] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:
        readout_marker_start_points + marker_points] = 1

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
    amplitude_floquet_sequence = Sequence(name='floquet',
                                          variable='floquet_amp',
                                          variable_label='floquet_amp',
                                          variable_unit='V',
                                          step=step,
                                          start=start,
                                          stop=stop)
    calib_dict = get_calibration_dict()

    resolution = 1 / calib_dict['sample_rate']
    readout_start = calib_dict['pulse_end'] + calib_dict['pulse_readout_delay']
    readout_marker_start = readout_start - calib_dict['marker_readout_delay']

    readout_start_points = round(readout_start / resolution)
    readout_marker_start_points = round(readout_marker_start / resolution)
    readout_points = int(np.ceil(calib_dict['readout_time'] / resolution))
    pulse_end_points = int(np.ceil(calib_dict['pulse_end'] / resolution))
    marker_points = int(np.ceil(calib_dict['marker_time'] / resolution))
    floquet_points = int(np.ceil(floquet_dur / resolution))
    total_points = int(np.ceil(calib_dict['cycle_duration'] / resolution))

    readout_template = Waveform(length=total_points, channel=channels[1])

    readout_template.wave[
        readout_start_points:readout_start_points + readout_points] = 1
    readout_template.marker_1[
        readout_marker_start_points:
        readout_marker_start_points + marker_points] = 1

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
