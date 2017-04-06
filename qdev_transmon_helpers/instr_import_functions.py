import numpy as np
from functools import partial
import logging
from . import get_qubit_count, config_alazar, get_alazar_seq_mode, \
    set_alazar_seq_mode

alazar_acq_types = ['samp', 'ave', 'rec']

# TODO: AWG external ref source
# TODO: alazar seq mode param default val


def import_decadac(port=5, station=None):
    from qcodes.instrument_drivers.Harvard.Decadac import Decadac
    qubit_count = get_qubit_count()
    slot_count = int(np.ceil(qubit_count / 4))
    dec_slots = []
    dec_chans = []

    for s in range(slot_count):
        dec_slots.append(Decadac('dec_slot_{}'.format(s), port=5, slot=s))
        if station is not None:
            station.add_component(dec_slots[s])
        dec_slots[s].mode(2)

    for c in range(qubit_count):
        slot_i = int(np.floor(c / 4))
        chan = c % 4
        dec_chans.append(getattr(dec_slots[slot_i],
                                 'ch{}_voltage'.format(chan)))
        dec_chans[c].set_step(0.1)
        dec_chans[c].set_delay(0.01)

    logging.info('imported {} decadac slots \'dec_slots\': sorted into {} '
                 'decadac channels \'dec_chans\''.format(len(dec_slots),
                                                         len(dec_chans)))
    print('imported {} decadac slots \'dec_slots\': sorted into {} '
          'decadac channels \'dec_chans\''.format(len(dec_slots),
                                                  len(dec_chans)))
    print('-------------------------')
    return dec_slots, dec_chans


def import_vna(visa_address, name='vna', station=None):
    from qcodes.instrument_drivers.rohde_schwarz.ZNB20 import ZNB20
    vna = ZNB20(name, visa_address)
    if station is not None:
        station.add_component(vna)
    logging.info('imported VNA ZNB20: \'{}\''.format(name))
    print('imported VNA ZNB20: \'vna\''.format(name))
    print('-------------------------')
    return vna


def import_manual_param(name='dummy_time', station=None):
    import qcodes.instrument.parameter as parameter
    dummy_time = parameter.ManualParameter(name=name)
    if station is not None:
        station.add_component(dummy_time)
    logging.info('imported manual parameter: \'{}\''.format(name))
    print('imported manual parameter: \'{}\''.format(name))
    print('-------------------------')
    return dummy_time


def import_alazar(name='alazar', station=None, clock_source='EXTERNAL_CLOCK_10MHz_REF'):
    from qcodes.instrument_drivers.AlazarTech.ATS9360 import AlazarTech_ATS9360
    alazar = AlazarTech_ATS9360(name=name)
    config_alazar(alazar, seq_mode=False, clock_source=clock_source)
    alazar.add_parameter(name='seq_mode',
                         get_cmd=partial(get_alazar_seq_mode, alazar),
                         set_cmd=partial(set_alazar_seq_mode, alazar))
    if station is not None:
        station.add_component(alazar)
    logging.info('imported Alazar ATS9360: \'{}\''.format(name))
    print('imported Alazar ATS9360: \'{}\''.format(name))
    print('-------------------------')
    return alazar


def import_acq_controller(alazar, name=None, ctrl_type='ave', station=None):
    from qcodes.instrument_drivers.AlazarTech.acq_controllers import ATS9360Controller
    ctrl_name = name or (ctrl_type + '_ctrl')
    if ctrl_type is alazar_acq_types[0]:
        ctrl = ATS9360Controller(name=ctrl_name, alazar_name=alazar.name,
                                 integrate_samples=False, average_records=True)
    elif ctrl_type is alazar_acq_types[1]:
        ctrl = ATS9360Controller(name=ctrl_name, alazar_name=alazar.name,
                                 integrate_samples=True, average_records=True)
    elif ctrl_type is alazar_acq_types[2]:
        ctrl = ATS9360Controller(name=ctrl_name, alazar_name=alazar.name,
                                 integrate_samples=True, average_records=False)
    else:
        raise Exception('acquisition controller type must be in {}, received: '
                        '{}'.format(alazar_acq_types, ctrl_type))
    if station is not None:
        station.add_component(ctrl)
    logging.info('imported acquisition controller: \'{}\''.format(ctrl_name))
    print('imported acquisition controller: \'{}\''.format(ctrl_name))
    print('-------------------------')
    return ctrl


def import_rs(visa_address, name='rs_source', station=None):
    from qcodes.instrument_drivers.rohde_schwarz.SGS100A import RohdeSchwarz_SGS100A
    rs = RohdeSchwarz_SGS100A(name, visa_address)
    if station is not None:
        station.add_component(rs)
    logging.info('imported r&s microwave source SGS100A: \'{}\''.format(name))
    print('imported r&s microwave source SGS100A: \'{}\''.format(name))
    print('-------------------------')
    return rs


def import_awg(visa_address, name='awg', timeout=40, station=None):
    from qcodes.instrument_drivers.tektronix.AWG5014 import Tektronix_AWG5014
    awg = Tektronix_AWG5014(name, visa_address, timeout=timeout)
    # TODO
    if station is not None:
        station.add_component(awg)
    logging.info('imported tektronix awg5014c: \'{}\''.format(name))
    print('imported tektronix awg5014c: \'{}\''.format(name))
    print('-------------------------')
    return awg
