from . import check_sample_rate, make_save_send_load_awg_file


def set_up_sequence(awg, alazar, acq_controllers, sequence, seq_mode=1):
    """
    Function which checks sample rate compatability between sequence and awg
    setting, uploads sequence to awg, sets the alazar instrument to the
    relevant sequence mode and sets the acquisition controller aquisiion
    parameter to have setpoints based on the sequence variable.

    Args:
        awg instrument (AWG5014)
        alazar instrument
        acq_controllers list
        sequence for upload
        seq_mode (default 1)
    """
    check_sample_rate(awg)
    make_save_send_load_awg_file(awg, sequence)
    alazar.seq_mode(seq_mode)
    record_num = len(sequence)
    if seq_mode == 1:
        for ctrl in acq_controllers:
            try:
                ctrl.records_per_buffer(record_num)
                try:
                    start = sequence.variable_array[0]
                    stop = sequence.variable_array[-1]
                except Exception:
                    start = 0
                    stop = len(sequence) - 1

                ctrl.acquisition.set_base_setpoints(base_name=sequence.variable,
                                                    base_label=sequence.variable_label,
                                                    base_unit=sequence.variable_unit,
                                                    setpoints_start=start,
                                                    setpoints_stop=stop)
            except NotImplementedError:
                pass
    alazar.seq_mode(seq_mode)
    awg.all_channels_on()
    awg.run()
