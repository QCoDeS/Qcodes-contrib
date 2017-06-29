
import sys
pulse_building_folder = r'/Users/Natalie/Documents/PhD/Qdev/QcodesRelated/PulseBuilding' #r'A:\\PulseBuilding' #r'C:\\Users\\qcodes-natalie\\PulseBuilding'
if pulse_building_folder not in sys.path:
    sys.path.insert(0, pulse_building_folder)
from pulse_building import Waveform, Element, Sequence

from .math_functions import *
from .file_functions import *
from .temp_dict_funtions import *
from .measurement_plot_functions import *
from .analysis_plot_functions import *
from .analysis_helpers import *
from .loading_data import *
from .sweep_functions import *
from .vna_helper_functions import *
from .awg_helper_functions import *
from .pulse_building_functions import *
from .pulse_building_basic import *
from .pulse_building_allxy import *
from .alazar_rs_helper_functions import *
from .alazar_awg_wrapper import *
from .instr_import_functions import *
from .alazar_automation import *
from .temp_helpers import *
