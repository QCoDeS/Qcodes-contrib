from . import check_sample_rate


def set_up_sequence(awg, alazar, acq_controllers, sequence, seq_mode=1):
    check_sample_rate(awg)
    awg.make_save_send_load_awg_file(awg, sequence)
    alazar.seq_mode(seq_mode)
    record_num = len(sequence)
    if seq_mode == 1:
        for ctrl in acq_controllers:
            try:
                ctrl.records_per_buffer(record_num)
                if sequence.variable is not None:
                    ctrl.acquisition.set_base_setpoints(base_name=sequence.variable,
                                                        base_label=sequence.variable,
                                                        base_unit=sequence.variable_unit,
                                                        setpoints_start=sequence.start,
                                                        setpoints_stop=sequence.stop)
            except NotImplementedError:
                pass
