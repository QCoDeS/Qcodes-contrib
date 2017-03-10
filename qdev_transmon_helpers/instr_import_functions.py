from qcodes.instrument_drivers.Harvard.Decadac import Decadac
import numpy as np
import qcodes as qc
import logging
from . import get_qubit_count

alazar_acq_types = ['samp', 'rec', 'ave']

# TODO: import alazar function
# TODO: import alazar acq controller function
# TODO: import rohde schwarz function (more?)
# TODO: how should station be dealt with?

def import_decadac(port=5):
    qubit_count = get_qubit_count()
    slot_count = int(np.ceil(qubit_count / 4))
    dec_slots = []
    dec_chans = []

    for s in range(slot_count):
        dec_slots.append(Decadac('dec_slot_{}'.format(s), port=5, slot=s))
        qc.Station.add_component(dec_slots[s])

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


def import_vna(visa_address='TCPIP0::172.20.3.94::inst0::INSTR'):
    from qcodes.instrument_drivers.rohde_schwarz.ZNB20 import ZNB20
    vna = ZNB20('vna', visa_address)
    # qc.Station.add_component(vna)
    logging.info('imported VNA: \'vna\'')
    print('imported VNA: \'vna\'')
    print('-------------------------')
    return vna


def import_dummy_time():
    import qcodes.instrument.parameter as parameter
    dummy_time = parameter.ManualParameter(name="dummy_time")
    print('imported manual parameter: \'dummy_time\'')
    print('-------------------------')
    return dummy_time

def import_alazar():
    raise NotImplementedError

def import_acq_controller(alazar_name, type='ave'):
    return NotImplementedError

