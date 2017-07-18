def init():
    global files_setup
    files_setup= {'analysis_loc': False,
               'data_loc_fmt': False,
               'python_log_loc': False,
               'ipython_log_loc': False,
               'temp_dict_loc': False,
               'pulse_loc': False}
    global sample_name
    sample_name = None
    global qubit_count
    qubit_count = None
    global current_qubit
    current_qubit = None
