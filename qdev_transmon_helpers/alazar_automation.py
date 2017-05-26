import numpy as np

from . import get_demod_freq, do_cavity_freq_sweep, find_extreme, \
    set_calibration_val, set_single_demod_freq, make_ssb_qubit_seq, \
    get_calibration_val, measure_ssb, set_up_sequence, \
    make_rabi_gaussian_sequence, make_rabi_sequence, sweep1d

# TODO: docstrings
# TODO: find qubit to compare to average
# TODO: complete do_tracking_ssb_gate_sweep
# TODO: add time tracking automation


def calibrate_cavity(cavity, localos, acq_ctrl, alazar, centre_freq=None,
                     demod_freq=None, calib_update=True, cavity_power=-35,
                     lo_power=15, live_plot=True):
    if centre_freq is None:
        centre_freq = cavity.frequency()
    if demod_freq is None:
        demod_freq = get_demod_freq(cavity, localos, acq_ctrl)
    alazar_mode = alazar.seq_mode()
    alazar.seq_mode(0)
    cavity.status('on')
    cavity.power(cavity_power)
    localos.power(lo_power)
    data, plot = do_cavity_freq_sweep(cavity, localos, centre_freq, acq_ctrl,
                                      cavity_pm=3e6, freq_step=0.1e6,
                                      demod_freq=demod_freq,
                                      live_plot=True, key="mag", save=True)
    cavity_res, mag = find_extreme(
        data, x_key="frequency_set", y_key="mag", extr="min")
    good_cavity_freq = cavity_res + 3e5
    if calib_update:
        set_calibration_val('cavity_freqs', good_cavity_freq)
        set_calibration_val('cavity_pows', cavity_power)
        set_calibration_val('demod_freqs', demod_freq)
    set_single_demod_freq(cavity, localos, [acq_ctrl], demod_freq,
                          cav_freq=good_cavity_freq)
    alazar.seq_mode(alazar_mode)
    print('cavity_freq set to {}, mag = '.format(good_cavity_freq, mag))


def find_qubit(awg, alazar, acq_ctrl, qubit, start_freq=4e9, stop_freq=6e9,
               qubit_power=None, calib_update=True):
    if 'ssb_qubit' not in acq_ctrl.acquisition.setpoint_names[0][0]:
        ssb_seq = make_ssb_qubit_seq()
        set_up_sequence(awg, alazar, [acq_ctrl], ssb_seq, seq_mode=1)
    else:
        alazar.seq_mode(1)
    old_power = qubit.power()
    if qubit_power is None:
        qubit_power = get_calibration_val('spec_powers')
    qubit.status('on')
    qubit.power(qubit_power)
    qubit_freq = None
    qubit_mag = 0
    for centre in np.linspace(start_freq + 100e6, stop_freq - 100e6,
                              num=(stop_freq - start_freq) / 200e6):
        ssb_centre = centre
        data, pl = measure_ssb(qubit, acq_ctrl, ssb_centre, key="mag")
        freq, maximum = find_extreme(data, x_key="set", extr="max")
        if maximum > qubit_mag:
            qubit_freq = freq
            qubit_mag = maximum
    if calib_update:
        set_calibration_val('actual_qubit_positions', qubit_freq)
        set_calibration_val('spec_powers', qubit_power)
    qubit.power(old_power)
    print('qubit found at {}, mag {}'.format(qubit_freq, qubit_mag))
    return qubit_freq, qubit_mag


def do_tracking_ssb_gate_sweep(qubit, cavity, localos, rec_acq_ctrl,
                               ave_acq_ctrl, gate,
                               initial_qubit_freq, initial_cavity_freq,
                               demod_freq, gate_start, gate_stop,
                               gate_step=0.01, live_plot=True):
    raise NotImplementedError


def do_rabis(awg, alazar, acq_ctrl, qubit, start_dur=0, stop_dur=200e-9,
             step_dur=1e-9, pi_pulse_amp=None, qubit_power=None,
             freq_centre=None, freq_pm=10e6, freq_step=2e6, live_plot=True,
             calib_update=True, gaussian=False, pulse_sigmas_number=2,
             pulse_mod=False):
    if gaussian:
        rabis_uploaded = ('gaussian_drive_duration' is
                          acq_ctrl.acquisition.setpoint_names[0][0])
    else:
        rabis_uploaded = ('drive_duration' is
                          acq_ctrl.acquisition.setpoint_names[0][0])
    if (pi_pulse_amp is not None) or not rabis_uploaded:
        pi_pulse_amp = pi_pulse_amp or get_calibration_val(
            'pi_pulse_amplitudes')
        if gaussian:
            rabi_seq = make_rabi_gaussian_sequence(pi_pulse_amp,
                                                   pulse_sigmas_number,
                                                   start=start_dur,
                                                   stop=stop_dur,
                                                   step=step_dur,
                                                   pulse_mod=pulse_mod)
        else:
            rabi_seq = make_rabi_sequence(pi_pulse_amp,
                                          start=start_dur, stop=stop_dur,
                                          step=step_dur, pulse_mod=pulse_mod)
        set_up_sequence(awg, alazar, [acq_ctrl], rabi_seq, seq_mode=1)
    else:
        alazar.seq_mode(1)
    old_power = qubit.power()
    old_frequency = qubit.frequency()
    qubit_power = qubit_power or get_calibration_val('pi_pulse_powers')
    qubit.power(qubit_power)
    qubit.status('on')

    centre = freq_centre or get_calibration_val('actual_qubit_positions')
    freq_start = centre - freq_pm
    freq_stop = centre + freq_pm

    data, plot = sweep1d(acq_ctrl.acquisition, qubit.frequency, freq_start,
                         freq_stop, freq_step,
                         live_plot=live_plot, key="mag", save=True)

    if calib_update:
        set_calibration_val('pi_pulse_powers', qubit_power)
    qubit.power(old_power)
    qubit.frequency(old_frequency)
