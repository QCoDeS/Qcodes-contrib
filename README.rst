QDev Transmon Helpers and Examples
===================================

An experiment layer to go on top of QCoDeS for the transmon team.

Motivation
------
The basic motivation is to provide the layer between QCoDeS and the experimentalist so that the UI is useable and not getting in the way, this includes things like standardising logging procedure, datasaving procedure and the folder structure of an experiment. It also includes wrappers to make loading in data or recognising which analysis plot came from which dataset easier. Totally unfinished.


Modules
-------

math_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module containing basic math functions used for fitting.

file_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module containing experiment level helper functions such as setting the sample_name and qubit_count at the start of an experiment, setting up datasaving locations and logfiles etc, load data as well as plotting functions. 

sweep_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module containing wrapper functions which perform basic do1d and do2d type loops.

calib_dict_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module containing functions for interacting with the pickled dictionary used to save calibration values.

vna_helper_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains quite specific functions for the part of the experiment which mainly uses the R&S Vector Network Analyser. Mainly these functions are for setting instrument setting on the vna or doing specific sweeps but also analysis functions for the data we expect to obtain. Should possibly be split into instrument related functions and analysis related functions.

alazar_rs_helper_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains again specific functions for the time domain part of the experiment, specificaly the microwave drives and lock in
amplifier readout as well as specific measurement setups and fit functions for the resulting data.

instr_import_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module which has wrapper functions for instrument imports. Could be extended to include update initial setting on instrument if needed.

awg_helper_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module which has some wrapper functions for interacting with the awg but mainly consists of pulse building functions for 
running specific experiments. It also includes interaction with the 'pulse dictionary' we use for storing pulse calibration values.

alazar_awg_wrapper.py
^^^^^^^^^^^^^^^^^^^^^^^^^
This is the single function used to link the sequence for upload onto awg with the awg and the alazar readout settings.

Outstanding issues
------------------
- Need a better way to store, access, validate calibration data

- Plotting things (eg would be nice to be able to add markers, take cuts, change axes of live plot etc)

-	Folder structure is currently set as 

	- root/sample_name
	
			- /analysis
	
			- /data
	
			- /log_files
	
			- /pulse_lib
	which is not possible with the https://github.com/qdev-dk/Qcodes/blob/master/qcodes/utils/wrappers.py extension but also might be too experiment specific, need a way to make something more general so that the experimentalist can choose file setup at start.

- The way dataset identification is done is currently by the name of the folder in the specified directory which is pretty horrible. This should improve when QCoDeS has an id for a dataset but will it be readable/searchable etc? Unsolved issue of how to find a dataset and link it to an analysis plot or similar.

- A wrapper funtion for readable metadata and a readable snapshot of the station (or specified instrumente or parameter values from it) is required. This could be used to check the status of instruments easily and repeatably or to check the settings of instruments when a particular dataset was taken. This is a precursur to a 'monitor' and links with the idea of being able to add some instrument data to the pngs that get saved on demand. 
