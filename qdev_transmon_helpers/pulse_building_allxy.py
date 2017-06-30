from . import make_sequence_from_gate_lists

# TODO: pulse_mod

################################################################
# ALLXY
################################################################


allxy_gates = [['id', 'id'],
               ['x_pi', 'x_pi'],
               ['y_pi', 'y_pi'],
               ['x_pi', 'y_pi'],
               ['y_pi', 'x_pi'],
               ['x_pi_half', 'id'],
               ['y_pi_half', 'id'],
               ['x_pi_half', 'y_pi_half'],
               ['y_pi_half', 'x_pi_half'],
               ['x_pi_half', 'y_pi'],
               ['y_pi_half', 'x_pi'],
               ['x_pi', 'y_pi_half'],
               ['y_pi', 'x_pi_half'],
               ['x_pi_half', 'x_pi'],
               ['x_pi', 'x_pi_half'],
               ['y_pi_half', 'y_pi'],
               ['y_pi', 'y_pi_half'],
               ['x_pi', 'id'],
               ['y_pi', 'id'],
               ['x_pi_half', 'x_pi_half'],
               ['y_pi_half', 'y_pi_half']]


def make_allxy_sequence(SSBfreq=None, drag=False, channels=[1, 2, 3, 4],
                        spacing=None):
    seq = make_sequence_from_gate_lists(allxy_gates, SSBfreq=SSBfreq,
                                        drag=drag, variable_label=None,
                                        spacing=None, name='allxy_seq')
    seq.labels = {'seq_type': 'allxy', 'pulse_mod': False}
    return seq
