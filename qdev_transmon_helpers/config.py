################################################
# instr import 
################################################

alazar_acq_types = ['samp', 'ave', 'rec']

################################################
# calib dict
################################################

vna_dict_keys = ['expected_qubit_freq', 'g_value', 'bare_res_freq',
                 'pushed_res_freq', 'gate_volt']

alazar_dict_keys = ['int_time', 'int_delay', 'cavity_freq',
                    'cavity_pow', 'localos_pow', 'demod_freq',
                    'qubit_freq', 'pi_pulse_pow', 'spec_pow']

default_pulse_dict = {'cycle_time': 20e-6, 'sample_rate': 1e9,
                      'pulse_end': 10e-6, 'pulse_readout_delay': 30e-9,
                      'readout_time': 4e-6, 'readout_amp': 1,
                      'marker_time': 500e-9,
                      'marker_readout_delay': 0, 'qubit_spec_time': 1e-6,
                      'pulse_mod_time': 1.5e-6, 'pi_pulse_amp': 1,
                      'pi_half_pulse_amp': 0.5, 'pi_pulse_sigma': None,
                      'pi_pulse_dur': None, 'sigma_cutoff': 4,
                      'z_pulse_amp': None, 'z_pulse_dur': None,
                      'z_half_pulse_amp': None, 'drag_coef': 0.5}