from . import make_readout_wf, get_calibration_dict, \
    make_time_varying_sequence, make_multi_varying_sequence, \
    make_time_multi_varying_sequence, \
    cos_wave, sin_wave, flat, gaussian, SSB_Q_gaussian, \
    SSB_I_gaussian
from . import Segment, Waveform, Element, Sequence


# TODO: checks

################################################################
# Spectroscopy
################################################################

def make_calib_SSB_sequence(freq, amp=1, dur=None, channels=[1, 2]):
    p_dict = get_calibration_dict()
    dur = dur or p_dict['cycle_time']

    seq = Sequence(name='ssb_calib_seq')
    element = Element()
    waveform_i = Waveform(channel=channels[0])
    waveform_q = Waveform(channel=channels[1])
    waveform_i.wave = cos_wave(freq, amp, dur, p_dict['sample_rate'])
    waveform_q.wave = sin_wave(freq, amp, dur, p_dict['sample_rate'])
    element.add_waveform(waveform_i)
    element.add_waveform(waveform_q)
    seq.add_element(element)
    return seq


def make_spectrscopy_SSB_sequence(start, stop, step, channels=[1, 2, 4],
                                  pulse_mod=False):
    p_dict = get_calibration_dict()

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [0],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    time_before_qubit = p_dict['pulse_end'] - p_dict['qubit_spec_time']
    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']
    before_qubit_wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_before_qubit})
    variable_qubit_drive_I_segment = Segment(
        name='SSB_drive_I', gen_func=cos_wave,
        func_args={'amp': 1, 'dur': p_dict['qubit_spec_time']},
        time_markers=pulse_mod_markers)
    variable_qubit_drive_Q_segment = Segment(
        name='SSB_drive_Q', gen_func=sin_wave,
        func_args={'amp': 1, 'dur': p_dict['qubit_spec_time'],
                   'positive': False})
    after_qubit_wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    ssb_I_wf = Waveform(
        channel=channels[0],
        segment_list=[before_qubit_wait_segment,
                      variable_qubit_drive_I_segment,
                      after_qubit_wait_segment])
    ssb_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[before_qubit_wait_segment,
                      variable_qubit_drive_Q_segment,
                      after_qubit_wait_segment])
    readout_wf = make_readout_wf(channel=channels[2])

    ssb_element = Element()
    ssb_element.add_waveform(ssb_I_wf)
    ssb_element.add_waveform(ssb_Q_wf)
    ssb_element.add_waveform(readout_wf)

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    ssb_seq = make_multi_varying_sequence(
        ssb_element, channels[0], 1, 'freq', start, stop, step,
        channels[1], 'freq', start, stop, step,
        readout_ch=channels[2], marker_points=marker_points)
    return ssb_seq

################################################################
# Rabi and T1
################################################################


def make_rabi_carrier_sequence(start, stop, step, pi_amp=None,
                               channels=[1, 4], pulse_mod=False,
                               gaussian=True):
    p_dict = get_calibration_dict()
    pi_amp = pi_amp or p_dict['pi_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})

    if gaussian:
        variable_pi_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp})
        variable_arg = 'sigma'
    else:
        variable_pi_segment = Segment(name='square_pi_pulse', gen_func=flat,
                                      func_args={'amp': pi_amp})
        variable_arg = 'dur'

    wait_segment = Segment(name='wait', gen_func=flat,
                           func_args={'amp': 0, 'dur': time_after_qubit},
                           time_markers=pulse_mod_markers)
    rabi_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, variable_pi_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[1])

    rabi_element = Element()
    rabi_element.add_waveform(rabi_wf)
    rabi_element.add_waveform(readout_wf)
    rabi_element.sample_rate = p_dict['sample_rate']

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    rabi_sequence = make_time_varying_sequence(
        rabi_element, channels[0], 1, variable_arg, start, stop, step, 0,
        p_dict['cycle_time'], readout_ch=channels[1],
        marker_points=marker_points)

    return rabi_sequence


def make_rabi_SSB_sequence(start, stop, step, SSBfreq, channels=[1, 2, 4],
                           gaussian=True, pulse_mod=False,
                           pi_amp=None):
    p_dict = get_calibration_dict()

    pi_amp = pi_amp or p_dict['pi_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})

    if gaussian:
        variable_pi_I_segment = Segment(
            name='gaussian_SSB_pi_I_pulse', gen_func=SSB_I_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp,
                       'SSBfreq': SSBfreq})
        variable_pi_Q_segment = Segment(
            name='gaussian_SSB_pi_Q_pulse', gen_func=SSB_Q_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp,
                       'SSBfreq': SSBfreq})
        variable_arg = 'sigma'
    else:
        variable_pi_I_segment = Segment(
            name='square_SSB_pi_I_pulse', gen_func=cos_wave,
            func_args={'amp': pi_amp, 'freq': SSBfreq})
        variable_pi_Q_segment = Segment(
            name='square__SSB_pi_Q_pulse', gen_func=sin_wave,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'positive': False})
        variable_arg = 'dur'

    wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    rabi_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, variable_pi_I_segment,
                      wait_segment])
    rabi_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment, variable_pi_Q_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[2])

    rabi_element = Element()
    rabi_element.add_waveform(rabi_I_wf)
    rabi_element.add_waveform(rabi_Q_wf)
    rabi_element.add_waveform(readout_wf)

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    rabi_sequence = make_time_multi_varying_sequence(
        rabi_element, channels[0], 1, 'sigma', start, stop, step,
        channels[1], 1, variable_arg, start, stop, step,
        0, 0, p_dict['cycle_time'], marker_ch=channels[2],
        marker_points=marker_points)

    return rabi_sequence


def make_rabi_sequence(start, stop, step, channels=[1, 2, 4], pulse_mod=False,
                       SSBfreq=None, pi_amp=None, gaussian=True):
    if SSBfreq is not None:
        seq = make_rabi_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_amp=pi_amp, gaussian=gaussian)
    else:
        seq = make_rabi_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[2]],
            pulse_mod=pulse_mod, pi_amp=pi_amp, gaussian=gaussian)
    return seq


def make_t1_carrier_sequence(start, stop, step, pi_dur=None, pi_amp=None,
                             channels=[1, 4], gaussain=True, pulse_mod=False):
    p_dict = get_calibration_dict()

    pi_amp = pi_amp or p_dict['pi_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})
    if gaussian:
        pi_sigma = pi_dur or p_dict['pi_pulse_sigma']
        pi_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp,
                       'sigma': pi_sigma})
    else:
        pi_dur = pi_dur or p_dict['pi_pulse_dur']
        pi_segment = Segment(
            name='square_pi_pulse', gen_func=flat,
            func_args={'amp': pi_amp, 'dur': pi_dur})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat,
        func_args={'amp': 0},
        time_markers=pulse_mod_markers)

    wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    t1_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_segment,
                      variable_wait_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[1])

    t1_element = Element()
    t1_element.add_waveform(t1_wf)
    t1_element.add_waveform(readout_wf)

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    t1_sequence = make_time_varying_sequence(
        t1_element, channels[0], 2, 'dur', start, stop, step, 0,
        p_dict['cycle_time'], readout_ch=channels[1],
        marker_points=marker_points)
    return t1_sequence


def make_t1_SSB_sequence(start, stop, step, SSBfreq, pi_dur=None,
                         pi_amp=None, channels=[1, 2, 4], gaussain=True,
                         pulse_mod=False):
    p_dict = get_calibration_dict()

    pi_amp = pi_amp or p_dict['pi_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})

    if gaussian:
        pi_sigma = pi_dur or p_dict['pi_pulse_sigma']
        pi_I_segment = Segment(
            name='gaussian_SSB_pi_I_pulse', gen_func=SSB_I_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp,
                       'SSBfreq': SSBfreq, 'sigma': pi_sigma})
        pi_Q_segment = Segment(
            name='gaussian_SSB_pi_Q_pulse', gen_func=SSB_Q_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'], 'amp': pi_amp,
                       'SSBfreq': SSBfreq, 'sigma': pi_sigma})
    else:
        pi_dur = pi_dur or p_dict['pi_pulse_dur']
        pi_I_segment = Segment(
            name='square_SSB_pi_I_pulse', gen_func=cos_wave,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'dur': pi_dur})
        pi_Q_segment = Segment(
            name='square_SSB_pi_Q_pulse', gen_func=sin_wave,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'dur': pi_dur,
                       'positive': False})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat,
        func_args={'amp': 0},
        time_markers=pulse_mod_markers)

    wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    t1_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_I_segment,
                      variable_wait_segment, wait_segment])
    t1_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment, pi_Q_segment,
                      variable_wait_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[2])

    t1_element = Element()
    t1_element.add_waveform(t1_I_wf)
    t1_element.add_waveform(t1_Q_wf)
    t1_element.add_waveform(readout_wf)

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    t1_sequence = make_time_multi_varying_sequence(
        t1_element, channels[0], 2, 'dur', start, stop, step,
        channels[1], 2, 'dur', start, stop, step,
        0, 0, p_dict['cycle_time'], marker_ch=channels[2],
        marker_points=marker_points)
    return t1_sequence


def make_t1_seq(start, stop, step, SSBfreq=None, pi_dur=None,
                pi_amp=None, channels=[1, 4], gaussain=True,
                pulse_mod=False):
    if SSBfreq is not None:
        seq = make_t1_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_amp=pi_amp, pi_dur=pi_dur, gaussian=gaussian)
    else:
        seq = make_t1_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[2]],
            pulse_mod=pulse_mod, pi_amp=pi_amp, pi_dur=pi_dur,
            gaussian=gaussian)
    return seq


################################################################
# Ramsey and T2
################################################################


def make_ramsey_carrier_sequence(start, stop, step, pi_half_amp=None,
                                 pi_half_dur=None, channels=[1, 4],
                                 pulse_mod=False, gaussian=True):
    p_dict = get_calibration_dict()
    pi_half_amp = pi_half_amp or p_dict['pi_half_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})

    if gaussian:
        pi_half_sigma = pi_half_dur or p_dict['pi_half_pulse_sigma']
        pi_half_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'],
                       'amp': pi_half_amp, 'sigma': pi_half_sigma})
    else:
        pi_half_dur = pi_half_dur or p_dict['pi_half_pulse_dur']
        pi_half_segment = Segment(
            name='square_pi_pulse', gen_func=flat,
            func_args={'amp': pi_half_amp, 'dur': pi_half_dur})

    variable_wait_segment = Segment(
        name='pulse_pulse_delay', gen_func=flat,
        func_args={'amp': 0},
        time_markers=pulse_mod_markers)

    wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    ramsey_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_half_segment,
                      variable_wait_segment, pi_half_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[1])

    ramsey_element = Element()
    ramsey_element.add_waveform(ramsey_wf)
    ramsey_element.add_waveform(readout_wf)
    ramsey_element.sample_rate = p_dict['sample_rate']

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    ramsey_sequence = make_time_varying_sequence(
        ramsey_element, channels[0], 2, 'dur', start, stop, step, 0,
        p_dict['cycle_time'], readout_ch=channels[1],
        marker_points=marker_points)

    return ramsey_sequence


def make_ramsey_SSB_sequence(start, stop, step, SSBfreq, pi_half_amp=None,
                             pi_half_dur=None, channels=[1, 2, 4],
                             pulse_mod=False, gaussian=True):
    p_dict = get_calibration_dict()

    pi_half_amp = pi_half_amp or p_dict['pi_half_pulse_amp']

    time_after_qubit = p_dict['cycle_time'] - p_dict['pulse_end']

    if pulse_mod:
        pulse_mod_markers = {1: {'delay_time': [-p_dict['pulse_mod_time']],
                                 'duration_time': [p_dict['pulse_mod_time']]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat, func_args={'amp': 0})

    if gaussian:
        pi_half_sigma = pi_half_dur or p_dict['pi_half_pulse_sigma']
        pi_half_I_segment = Segment(
            name='gaussian_SSB_pi_half_I_pulse', gen_func=SSB_I_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'],
                       'amp': pi_half_amp, 'SSBfreq': SSBfreq,
                       'sigma': pi_half_sigma})
        pi_half_Q_segment = Segment(
            name='gaussian_SSB_pi_half_Q_pulse', gen_func=SSB_Q_gaussian,
            func_args={'sigma_cutoff': p_dict['sigma_cutoff'],
                       'amp': pi_half_amp, 'SSBfreq': SSBfreq,
                       'sigma': pi_half_sigma})
    else:
        pi_half_dur = pi_half_dur or p_dict['pi_half_pulse_sigma']
        pi_half_I_segment = Segment(
            name='square_SSB_pi_half_I_pulse', gen_func=cos_wave,
            func_args={'amp': pi_half_amp, 'freq': SSBfreq,
                       'dur': pi_half_dur})
        pi_half_Q_segment = Segment(
            name='square_SSB_pi_half_Q_pulse', gen_func=sin_wave,
            func_args={'amp': pi_half_amp, 'freq': SSBfreq,
                       'dur': pi_half_dur, 'positive': False})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat,
        func_args={'amp': 0},
        time_markers=pulse_mod_markers)

    wait_segment = Segment(
        name='wait', gen_func=flat,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    ramsey_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_half_I_segment,
                      variable_wait_segment, pi_half_I_segment, wait_segment])
    ramsey_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment, pi_half_Q_segment,
                      variable_wait_segment, pi_half_Q_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[2])

    ramsey_element = Element()
    ramsey_element.add_waveform(ramsey_I_wf)
    ramsey_element.add_waveform(ramsey_Q_wf)
    ramsey_element.add_waveform(readout_wf)

    marker_points = int(p_dict['marker_time'] * p_dict['sample_rate'])
    ramsey_sequence = make_time_multi_varying_sequence(
        ramsey_element, channels[0], 2, 'dur', start, stop, step,
        channels[1], 2, 'dur', start, stop, step,
        0, 0, p_dict['cycle_time'], marker_ch=channels[2],
        marker_points=marker_points)
    return ramsey_sequence


def make_ramsey_sequence(start, stop, step, SSBfreq=None, pi_half_amp=None,
                         pi_half_dur=None, channels=[1, 2, 4],
                         pulse_mod=False, gaussian=True):
    if SSBfreq is not None:
        seq = make_ramsey_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_half_amp=pi_half_amp, pi_half_dur=pi_half_dur,
            gaussian=gaussian)
    else:
        seq = make_ramsey_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[2]],
            pulse_mod=pulse_mod, pi_half_amp=pi_half_amp,
            pi_half_dur=pi_half_dur, gaussian=gaussian)
    return seq
