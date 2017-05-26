import qcodes as qc

from . import get_metadata, plot_data, get_sample_name, get_data_duration


def load(counter, plot=True, metadata=True, matplot=False):
    """
    Function like shownum which loads dataset, (optionally)
    QtPlot from default location and (optionally) prints metadata for that
    measurement given current EXPERIMENT_VARS['metadata_list'] entries

    Args:
        counter (int)
        plot(default True): do you want to return plots as well as dataset?
        metadata (default True): do you want to print the metadata?
        matplot (bool) (default False): default is to QtPlot the data
    Returns:
        dataset (qcodes DataSet)
        plot (QtPlot): optional
    """
    str_counter = '{0:03d}'.format(counter)
    data = qc.load_data(
        qc.DataSet.location_provider.fmt.format(sample_name=get_sample_name(),
                                                counter=str_counter))
    data.data_num = counter

    if metadata:
        get_metadata(data, display=True)
        get_data_duration(data)
    if plot:
        plots = plot_data(data, matplot=matplot)
        return data, plots
    else:
        return data


def get_metadata(dataset, display=True, specific_list=None):
    """
    Function which gets the metadata dictionary for a dataset.

    Args:
        dataset(qcodes dataset): dataset for which we want metadata
        display (default True): should metadata be printed as well as returned?
        specific_list (default None):

    Returns:
        meta_dict dictionary of the form:
            {'instr': {'param':
                            {'value': val,
                            'unit': u}}}
            for a parameter with name 'param' belonging to instrument 'instr'
            with value 'val' and unit 'u'
    """
    missing_keys = []
    if isinstance(dataset, int):
        dataset = load(dataset, plot=False, metadata=False)
    snapshot = dataset.snapshot()
    meta_dict = defaultdict(dict)
    for instr, param in specific_list or get_metadata_list():
        try:
            unit = getFromDict(snapshot, ["station",
                                          "instruments",
                                          instr,
                                          "parameters",
                                          param,
                                          "unit"])
            value = getFromDict(snapshot, ["station",
                                           "instruments",
                                           instr,
                                           "parameters",
                                           param,
                                           "value"])
            meta_dict[instr][param] = {}
            meta_dict[instr][param]['value'] = value
            meta_dict[instr][param]['unit'] = unit
        except KeyError:
            missing_keys.append([instr, param])
    if display:
        print_metadata(meta_dict)
    if len(missing_keys) > 0:
        print('\nSome of the specified parameters were not found in the '
              'snapshot metadata for this dataset:')
        for missing in missing_keys:
            print(missing[0] + ' ' + missing[1])
    return meta_dict