import qcodes as qc
import pickle
from . import plot_data_single_window, plot_data, \
    get_calibration_val, get_pulse_location

# TODO: rount -> int/ciel
# TODO: test segments
# TODO: make_save_send_load_awg_file -> not doing same thing twice!
# TODO: docstrings
# TODO: checks
# TODO: floquet

#################################################################
# General AWG and Sequence functions
#################################################################


def get_current_seq(awg):
    seq_name = awg.current_seq()
    path = get_pulse_location()
    seq = pickle.load(open(path + seq_name, "rb"))
    return seq


def make_save_send_load_awg_file(awg, sequence, file_name):
    """
    WYSIYWYG

    Args:
        awg instrument for upload
        sequence to be uploaded
    """
    awg.make_and_save_awg_file(*sequence.unwrap(), filename=file_name)
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


def check_seq_uploaded(awg, seq_type, dict_to_check,
                       start=None, stop=None, step=None):
    uploaded_seq = get_current_seq(awg)
    if uploaded_seq.labels['seq_type'] is not seq_type:
        return False
    for k in dict_to_check:
        if k not in uploaded_seq.labels:
            return False
        elif uploaded_seq.labels[k] != dict_to_check[k]:
            return False
    try:
        if ([start, stop, step] !=
                [uploaded_seq.start, uploaded_seq.stop, uploaded_seq.step]):
            return False
    except AttributeError:
        return False
    return True


def sweep_awg_chans(meas_param, awg_ch_1, awg_ch_2, start, stop, step,
                    delay=0.01, live_plot=True, key=None, save=True):
    loop = qc.Loop(awg_ch_1.sweep(start, stop, step)).each(
        qc.Task(awg_ch_2.set, awg_ch_1.get),
        meas_param)
    if live_plot:
        dataset = loop.get_data_set()
        dataset.data_num = dataset.location_provider.counter
        plot = plot_data_single_window(dataset, meas_param, key=key)
        try:
            if save:
                _ = loop.with_bg_task(plot.update, plot.save).run()
            else:
                _ = loop.with_bg_task(plot.update).run()
                print('warning: plots not saved, if you want to save this'
                      'plot run plot.save()')
        except KeyboardInterrupt:
            print("Measurement Interrupted")
        return dataset, plot
    else:
        data = loop.run()
        data.data_num = data.location_provider.counter
        plots = plot_data(data, key=key)
        if (key is not None) and save:
            plots.save()
        else:
            print('warning: plots not saved. To save one choose '
                  'and run plots[i].save()')
        return data, plots
