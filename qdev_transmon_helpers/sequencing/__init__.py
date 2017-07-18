import sys
pulse_building_folder = r'A:\\PulseBuilding' #r'C:\\Users\\qcodes-natalie\\PulseBuilding' r'/Users/Natalie/Documents/PhD/Qdev/QcodesRelated/PulseBuilding'
if pulse_building_folder not in sys.path:
    sys.path.insert(0, pulse_building_folder)
from pulse_building import Segment, Waveform, Element, Sequence

from .. import get_calibration_dict, get_allowed_keys

from ..math_functions import *
from .waveform_makers import *
from .helpers import *
from .basic import *
from .allxy import *
from .benchmarking import *
# from .floquet import *
# from .legacy import *