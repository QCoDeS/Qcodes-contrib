import sys
pulse_building_folder = r'/Users/Natalie/Documents/PhD/Qdev/QcodesRelated/PulseBuilding' #r'A:\\PulseBuilding' #r'C:\\Users\\qcodes-natalie\\PulseBuilding' 
if pulse_building_folder not in sys.path:
    sys.path.insert(0, pulse_building_folder)
from pulse_building import Segment, Waveform, Element, Sequence, segment_functions

from .. import get_calibration_dict, get_allowed_keys

from ..math_functions import *
from .waveform_makers import *
from .helpers import *
from .basic import *
from .allxy import *
from .benchmarking import *
from .floquet_new import *
from .majorana import *
# from .legacy import *