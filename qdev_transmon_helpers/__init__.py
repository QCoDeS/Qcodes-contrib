
import sys
pulse_building_folder = r'A:\\PulseBuilding' # r'/Users/Natalie/Documents/PhD/Qdev/QcodesRelated/PulseBuilding' #r'C:\\Users\\qcodes-natalie\\PulseBuilding'
if pulse_building_folder not in sys.path:
    sys.path.insert(0, pulse_building_folder)
from pulse_building import Waveform, Element, Sequence

from .math_functions import *
from .file_functions import *
from .sweep_functions import *
from .calib_dict_functions import *
from .vna_helper_functions import *
from .awg_helper_functions import *
from .alazar_rs_helper_functions import *
from .instr_import_functions import *
from .alazar_awg_wrapper import *
