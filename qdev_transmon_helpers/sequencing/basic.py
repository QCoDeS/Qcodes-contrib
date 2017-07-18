from . import make_readout_wf, get_calibration_dict, \
    make_time_varying_sequence, make_multi_varying_sequence, \
    make_time_multi_varying_sequence, \
    cos_array, sin_array, flat_array, gaussian_array, cos_gaussian_array, \
    sin_gaussian_array
from . import Segment, Waveform, Element, Sequence

# TODO: T2 echo
# TODO: docstrings
# TODO: checks
# TODO: drag
# TODO: Exception types
# TODO: add multi qubit reaout

###############################################################
# Single Element Sequences
###############################################################


def make_readout_single_sequence(channels=[4]):
    seq = Sequence(name='readout_seq')
    element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    readout_wf = make_readout_wf(first_in_sequence=True, channel=channels[0])
    element.add_waveform(readout_wf)
    seq.add_element(element)
    seq.labels = {'seq_type': 'readout'}
    return seq


def make_readout_SSB_single_sequence(freq_list, channels=[3, 4]):
    seq = Sequence(name='readout_SSB_seq')
    element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    readout_wf_I = make_readout_ssb_wf_I(freq_list, first_in_sequence=True,
                                         channels=channels[0])
    readout_wf_Q = make_readout_ssb_wf_Q(freq_list, channels=channels[1])
    element.add_waveform(readout_wf_I)
    element.add_waveform(readout_wf_Q)
    seq.add_element(element)
    seq.labels = {'seq_type': 'readout_SSB', 'freq_list': freq_list}
    return seq

def make_calib_SSB_single_sequence(freq, amp=1, dur=None, channels=[1, 2]):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    dur = dur or calib_dict['cycle_time'][qubit]

    seq = Sequence(name='ssb_calib_seq')
    element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    waveform_i = Waveform(channel=channels[0])
    waveform_q = Waveform(channel=channels[1])
    waveform_i.wave = cos_array(
        freq, amp, dur, calib_dict['sample_rate'][qubit])
    waveform_q.wave = sin_array(
        freq, amp, dur, calib_dict['sample_rate'][qubit])
    element.add_waveform(waveform_i)
    element.add_waveform(waveform_q)
    seq.add_element(element)
    return seq

################################################################
# Readout SSB sweep
################################################################

def make_readout_SSB_sequence(start, stop, step, channels=[3, 4]):
    if len(channels) < 2:
        raise Exception('2 channels needed for single sideband sequence for'
                        ' I and Q')
    element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    readout_wf_I = make_readout_ssb_wf_I([freq], channels=channels[0])
    readout_wf_Q = make_readout_ssb_wf_Q([freq], channels=channels[1])
    element.add_waveform(readout_wf_I)
    element.add_waveform(readout_wf_Q)
    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    seq = make_multi_varying_sequence(    
        element, channels[0], 1, 'freq', start, stop, step,
        channels[1], 1, 'freq', start, stop, step, name="reabout_SSB_sequence",
        variable_name='LSB_drive_detuning', variable_unit='Hz',
        readout_ch=channels[1], marker_points=marker_points)
    seq.labels = {'seq_type': 'readout_SSB'}
    return ssb_seq

################################################################
# Spectroscopy
################################################################

def make_spectroscopy_SSB_sequence(start, stop, step, channels=[1, 2, 4],
                                   pulse_mod=False):
    if len(channels) < 3:
        raise Exception('3 channels needed for single sideband sequence for'
                        ' I, Q and readout')
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    time_before_qubit = (calib_dict['pulse_end'][qubit] -
                         calib_dict['qubit_spec_time'][qubit])
    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])
    before_qubit_wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_before_qubit})
    variable_qubit_drive_I_segment = Segment(
        name='SSB_drive_I', gen_func=cos_array,
        func_args={'amp': 1, 'dur': calib_dict['qubit_spec_time'][qubit]})
    variable_qubit_drive_Q_segment = Segment(
        name='SSB_drive_Q', gen_func=sin_array,
        func_args={'amp': 1, 'dur': calib_dict['qubit_spec_time'][qubit],
                   'positive': False})
    after_qubit_wait_segment = Segment(
        name='wait', gen_func=flat_array,
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
    readout_wf = make_readout_wf(channel=channels[-1])

    ssb_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    ssb_element.add_waveform(ssb_I_wf)
    ssb_element.add_waveform(ssb_Q_wf)
    ssb_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    ssb_seq = make_multi_varying_sequence(
        ssb_element, channels[0], 1, 'freq', start, stop, step,
        channels[1], 1, 'freq', start, stop, step, name="ssb_sequence",
        variable_name='LSB_drive_detuning', variable_unit='Hz',
        readout_ch=channels[-1], marker_points=marker_points)
    ssb_seq.labels = {'seq_type': 'spectroscopy', 'pulse_mod': pulse_mod}
    return ssb_seq


################################################################
# Rabi and T1
################################################################


def _make_rabi_carrier_sequence(start, stop, step, pi_amp=None,
                                channels=[1, 4], pulse_mod=False,
                                gaussian=True):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_amp = pi_amp or calib_dict['pi_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})

    if gaussian:
        variable_pi_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_amp})
        variable_arg = 'sigma'
    else:
        variable_pi_segment = Segment(
            name='square_pi_pulse', gen_func=flat_array,
            func_args={'amp': pi_amp})
        variable_arg = 'dur'

    wait_segment = Segment(name='wait', gen_func=flat_array,
                           func_args={'amp': 0, 'dur': time_after_qubit},
                           time_markers=pulse_mod_markers)
    rabi_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, variable_pi_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    rabi_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    rabi_element.add_waveform(rabi_wf)
    rabi_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    rabi_sequence = make_time_varying_sequence(
        rabi_element, channels[0], 1, variable_arg, start, stop, step, 0,
        calib_dict['cycle_time'][qubit], name='rabi_seq',
        variable_name='pi_pulse_' + variable_arg, variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    rabi_sequence.labels = {'SSBfreq': None, 'seq_type': 'rabi',
                            'gaussian': gaussian, 'drag': False,
                            'pulse_mod': pulse_mod}
    return rabi_sequence


def _make_rabi_SSB_sequence(start, stop, step, SSBfreq, channels=[1, 2, 4],
                            gaussian=True, pulse_mod=False,
                            pi_amp=None):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_amp = pi_amp or calib_dict['pi_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment_I = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    compensating_wait_segment_Q = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})

    if gaussian:
        variable_pi_I_segment = Segment(
            name='gaussian_SSB_pi_I_pulse', gen_func=cos_gaussian_array,
            func_args={
                'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                'amp': pi_amp, 'SSBfreq': SSBfreq})
        variable_pi_Q_segment = Segment(
            name='gaussian_SSB_pi_Q_pulse', gen_func=sin_gaussian_array,
            func_args={
                'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                'amp': pi_amp, 'SSBfreq': SSBfreq, 'positive': False})
        variable_arg = 'sigma'
    else:
        variable_pi_I_segment = Segment(
            name='square_SSB_pi_I_pulse', gen_func=cos_array,
            func_args={'amp': pi_amp, 'freq': SSBfreq})
        variable_pi_Q_segment = Segment(
            name='square__SSB_pi_Q_pulse', gen_func=sin_array,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'positive': False})
        variable_arg = 'dur'

    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    rabi_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment_I, variable_pi_I_segment,
                      wait_segment])
    rabi_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment_Q, variable_pi_Q_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    rabi_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    rabi_element.add_waveform(rabi_I_wf)
    rabi_element.add_waveform(rabi_Q_wf)
    rabi_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    rabi_sequence = make_time_multi_varying_sequence(
        rabi_element, channels[0], 1, variable_arg, start, stop, step,
        channels[1], 1, variable_arg, start, stop, step,
        0, 0, calib_dict['cycle_time'][qubit], name='rabi_ssb_seq',
        variable_name='pi_pulse_' + variable_arg, variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    rabi_sequence.labels = {'SSBfreq': SSBfreq, 'seq_type': 'rabi',
                            'gaussian': gaussian, 'drag': False,
                            'pulse_mod': pulse_mod}
    return rabi_sequence


def make_rabi_sequence(start, stop, step, SSBfreq=None, channels=[1, 2, 4],
                       pulse_mod=False, pi_amp=None, gaussian=True):
    if SSBfreq is not None:
        if len(channels) < 3:
            raise Exception('at least 3 channels needed for single sideband '
                            'sequence for I, Q and readout')
        seq = _make_rabi_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_amp=pi_amp, gaussian=gaussian)
    else:
        if len(channels) < 2:
            raise Exception('at least 2 channels needed for drive and readout')
        seq = _make_rabi_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[-1]],
            pulse_mod=pulse_mod, pi_amp=pi_amp, gaussian=gaussian)
    return seq


def _make_t1_carrier_sequence(start, stop, step, pi_dur=None, pi_amp=None,
                              channels=[1, 4], gaussian=True, pulse_mod=False):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_amp = pi_amp or calib_dict['pi_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    if gaussian:
        pi_sigma = pi_dur or calib_dict['pi_pulse_sigma'][qubit]
        pi_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian_array,
            func_args={
                'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                'amp': pi_amp, 'sigma': pi_sigma})
    else:
        pi_dur = pi_dur or calib_dict['pi_pulse_dur'][qubit]
        pi_segment = Segment(
            name='square_pi_pulse', gen_func=flat_array,
            func_args={'amp': pi_amp, 'dur': pi_dur})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat_array,
        func_args={'amp': 0})

    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    t1_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_segment,
                      variable_wait_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    t1_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    t1_element.add_waveform(t1_wf)
    t1_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    t1_sequence = make_time_varying_sequence(
        t1_element, channels[0], 2, 'dur', start, stop, step, 0,
        calib_dict['cycle_time'][qubit], name='t1_seq',
        variable_name='pi_pulse_readout_delay', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    t1_sequence.labels = {'SSBfreq': None, 'seq_type': 't1',
                          'gaussian': gaussian, 'drag': False,
                          'pulse_mod': pulse_mod}
    return t1_sequence


def _make_t1_SSB_sequence(start, stop, step, SSBfreq, pi_dur=None,
                          pi_amp=None, channels=[1, 2, 4], gaussian=True,
                          pulse_mod=False):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_amp = pi_amp or calib_dict['pi_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment_I = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    compensating_wait_segment_Q = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})

    if gaussian:
        pi_sigma = pi_dur or calib_dict['pi_pulse_sigma'][qubit]
        pi_I_segment = Segment(
            name='gaussian_SSB_pi_I_pulse', gen_func=cos_gaussian_array,
            func_args={
                'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                'amp': pi_amp, 'SSBfreq': SSBfreq, 'sigma': pi_sigma})
        pi_Q_segment = Segment(
            name='gaussian_SSB_pi_Q_pulse', gen_func=sin_gaussian_array,
            func_args={
                'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                'amp': pi_amp,
                'SSBfreq': SSBfreq, 'sigma': pi_sigma, 'positive': False})
    else:
        pi_dur = pi_dur or calib_dict['pi_pulse_dur'][qubit]
        pi_I_segment = Segment(
            name='square_SSB_pi_I_pulse', gen_func=cos_array,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'dur': pi_dur})
        pi_Q_segment = Segment(
            name='square_SSB_pi_Q_pulse', gen_func=sin_array,
            func_args={'amp': pi_amp, 'freq': SSBfreq, 'dur': pi_dur,
                       'positive': False})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat_array,
        func_args={'amp': 0})

    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    t1_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment_I, pi_I_segment,
                      variable_wait_segment, wait_segment])
    t1_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment_Q, pi_Q_segment,
                      variable_wait_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    t1_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    t1_element.add_waveform(t1_I_wf)
    t1_element.add_waveform(t1_Q_wf)
    t1_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    t1_sequence = make_time_multi_varying_sequence(
        t1_element, channels[0], 2, 'dur', start, stop, step,
        channels[1], 2, 'dur', start, stop, step,
        0, 0, calib_dict['cycle_time'][qubit], name='t1_ssb_seq',
        variable_name='pi_pulse_readout_delay', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    t1_sequence.labels = {'SSBfreq': SSBfreq, 'seq_type': 't1',
                          'gaussian': gaussian, 'drag': False,
                          'pulse_mod': pulse_mod}
    return t1_sequence


def make_t1_sequence(start, stop, step, SSBfreq=None, pi_dur=None,
                     pi_amp=None, channels=[1, 2, 4], gaussian=True,
                     pulse_mod=False):
    if SSBfreq is not None:
        if len(channels) < 3:
            raise Exception('at least 3 channels needed for single sideband '
                            'sequence for I, Q and readout')
        seq = _make_t1_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_amp=pi_amp, pi_dur=pi_dur, gaussian=gaussian)
    else:
        if len(channels) < 2:
            raise Exception('at least 2 channels needed for drive and readout')
        seq = _make_t1_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[-1]],
            pulse_mod=pulse_mod, pi_amp=pi_amp, pi_dur=pi_dur,
            gaussian=gaussian)
    return seq


################################################################
# Ramsey and T2
################################################################


def _make_ramsey_carrier_sequence(start, stop, step, pi_half_amp=None,
                                  pi_half_dur=None, channels=[1, 4],
                                  pulse_mod=False, gaussian=True):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_half_amp = pi_half_amp or calib_dict['pi_half_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})

    if gaussian:
        pi_half_sigma = pi_half_dur or calib_dict['pi_half_pulse_sigma'][qubit]
        pi_half_segment = Segment(
            name='gaussian_pi_pulse', gen_func=gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_half_amp, 'sigma': pi_half_sigma})
    else:
        pi_half_dur = pi_half_dur or calib_dict['pi_half_pulse_dur'][qubit]
        pi_half_segment = Segment(
            name='square_pi_pulse', gen_func=flat_array,
            func_args={'amp': pi_half_amp, 'dur': pi_half_dur})

    variable_wait_segment = Segment(
        name='pulse_pulse_delay', gen_func=flat_array,
        func_args={'amp': 0})

    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    ramsey_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment, pi_half_segment,
                      variable_wait_segment, pi_half_segment,
                      wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    ramsey_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    ramsey_element.add_waveform(ramsey_wf)
    ramsey_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    ramsey_sequence = make_time_varying_sequence(
        ramsey_element, channels[0], 2, 'dur', start, stop, step, 0,
        calib_dict['cycle_time'][qubit], name='ramsey_seq',
        variable_name='pi_half_pulse_pi_half_pulse_delay', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    ramsey_sequence.labels = {'SSBfreq': None, 'seq_type': 'ramsey',
                              'gaussian': gaussian, 'drag': False,
                              'pulse_mod': pulse_mod}
    return ramsey_sequence


def _make_ramsey_SSB_sequence(start, stop, step, SSBfreq, pi_half_amp=None,
                              pi_half_dur=None, channels=[1, 2, 4],
                              pulse_mod=False, gaussian=True):
    calib_dict = get_calibration_dict()
    qubit = calib_dict['current_qubit']
    pi_half_amp = pi_half_amp or calib_dict['pi_half_pulse_amp'][qubit]

    time_after_qubit = (calib_dict['cycle_time'][qubit] -
                        calib_dict['pulse_end'][qubit])

    if pulse_mod:
        pulse_mod_markers = {
            1: {'delay_time': [-1 * calib_dict['pulse_mod_time'][qubit]],
                'duration_time': [calib_dict['pulse_mod_time'][qubit]]}}
    else:
        pulse_mod_markers = None

    compensating_wait_segment_I = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})
    compensating_wait_segment_Q = Segment(
        name='compensating_wait', gen_func=flat_array, func_args={'amp': 0})

    if gaussian:
        pi_half_sigma = pi_half_dur or calib_dict['pi_half_pulse_sigma'][qubit]
        pi_half_I_segment = Segment(
            name='gaussian_SSB_pi_half_I_pulse', gen_func=cos_gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_half_amp, 'SSBfreq': SSBfreq,
                       'sigma': pi_half_sigma})
        pi_half_Q_segment = Segment(
            name='gaussian_SSB_pi_half_Q_pulse', gen_func=sin_gaussian_array,
            func_args={'sigma_cutoff': calib_dict['sigma_cutoff'][qubit],
                       'amp': pi_half_amp, 'SSBfreq': SSBfreq,
                       'sigma': pi_half_sigma, 'positive': False})
    else:
        pi_half_dur = pi_half_dur or calib_dict['pi_half_pulse_sigma'][qubit]
        pi_half_I_segment = Segment(
            name='square_SSB_pi_half_I_pulse', gen_func=cos_array,
            func_args={'amp': pi_half_amp, 'freq': SSBfreq,
                       'dur': pi_half_dur})
        pi_half_Q_segment = Segment(
            name='square_SSB_pi_half_Q_pulse', gen_func=sin_array,
            func_args={'amp': pi_half_amp, 'freq': SSBfreq,
                       'dur': pi_half_dur, 'positive': False})

    variable_wait_segment = Segment(
        name='pulse_readout_delay', gen_func=flat_array,
        func_args={'amp': 0})

    wait_segment = Segment(
        name='wait', gen_func=flat_array,
        func_args={'amp': 0, 'dur': time_after_qubit},
        time_markers=pulse_mod_markers)

    ramsey_I_wf = Waveform(
        channel=channels[0],
        segment_list=[compensating_wait_segment_I, pi_half_I_segment,
                      variable_wait_segment, pi_half_I_segment, wait_segment])
    ramsey_Q_wf = Waveform(
        channel=channels[1],
        segment_list=[compensating_wait_segment_Q, pi_half_Q_segment,
                      variable_wait_segment, pi_half_Q_segment, wait_segment])
    readout_wf = make_readout_wf(channel=channels[-1])

    ramsey_element = Element(sample_rate=calib_dict['sample_rate'][qubit])
    ramsey_element.add_waveform(ramsey_I_wf)
    ramsey_element.add_waveform(ramsey_Q_wf)
    ramsey_element.add_waveform(readout_wf)

    marker_points = int(calib_dict['marker_time'][qubit] *
                        calib_dict['sample_rate'][qubit])
    ramsey_sequence = make_time_multi_varying_sequence(
        ramsey_element, channels[0], 2, 'dur', start, stop, step,
        channels[1], 2, 'dur', start, stop, step,
        0, 0, calib_dict['cycle_time'][qubit], name='ramsey_ssb_seq',
        variable_name='pi_half_pulse_pi_half_pulse_delay', variable_unit='s',
        readout_ch=channels[-1], marker_points=marker_points)
    ramsey_sequence.labels = {'SSBfreq': SSBfreq, 'seq_type': 'ramsey',
                              'gaussian': gaussian, 'drag': False,
                              'pulse_mod': pulse_mod}
    return ramsey_sequence


def make_ramsey_sequence(start, stop, step, SSBfreq=None, pi_half_amp=None,
                         pi_half_dur=None, channels=[1, 2, 4],
                         pulse_mod=False, gaussian=True):
    if SSBfreq is not None:
        if len(channels) < 3:
            raise Exception('at least 3 channels needed for single sideband '
                            'sequence for I, Q and readout')
        seq = _make_ramsey_SSB_sequence(
            start, stop, step, SSBfreq, channels=channels, pulse_mod=pulse_mod,
            pi_half_amp=pi_half_amp, pi_half_dur=pi_half_dur,
            gaussian=gaussian)
    else:
        if len(channels) < 2:
            raise Exception('at least 2 channels needed for drive and readout')
        seq = _make_ramsey_carrier_sequence(
            start, stop, step, channels=[channels[0], channels[-1]],
            pulse_mod=pulse_mod, pi_half_amp=pi_half_amp,
            pi_half_dur=pi_half_dur, gaussian=gaussian)
    return seq
