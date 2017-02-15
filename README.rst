HelperFunctions
===================================

An experiment layer to go on top of QCoDeS for the transmon team.

Needs
------
The basic motivation is to provide the layer between QCoDeS and the experimentalist so that the UI is useable and npt getting in the way, this includes things like standardising loggind, datasaving and folder structure of an experiment and making loading in data or recognising which analysis plot came from which dataset easier. Lots of it is unpolished and unfinishes, its a hacky first round.


Modules
-------

general_helper_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module containing experiment level helper functions such as setting the smaple_name and qubit_count at the start of an experiment, setting up datasaving locations and logfiles etc, load data, etc. 

vna_helper_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains quite specific functions for the part of the experiment which mainly uses the R&S Vector Network Analyser. Mainly these functions are for setting instrument setting on the vna or doing specific sweeps but also analysis functions for the data we expect to obtain. Should possibly be split into instrument related functions and analysis related functions.

instr_import_functions.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This is a module which has wrapper functions for instrument imports. Could be extended to include update initial setting on instrument if needed.


Outstanding issues
------------------
- Currently datasets, metadata and png of plots are saved under only a counter number, this will become problematic if there is more than one plot per dataset, should also at least be possible to have a more explicit name containing some information about the plot

- Should be possible/easier to specify what you want to plot from a dataset on what axes and have this labeled etc

-	Folder structure is currently set as 

	- root/sample_name
	
			- /analysis
	
			- /data
	
			- /log_files
	
			- /pulse_lib
	which is not possible with the https://github.com/qdev-dk/Qcodes/blob/master/qcodes/utils/wrappers.py extension but also might be too experiment specific, need a way to make something more general so that the experimentalist can choose file setup at start.

- The way dataset identification is done is currently by the name of the folder in the specified directory which is pretty horrible. This should improve when QCoDeS has an id for a dataset but will it be readable/searchable etc? Unsolved issue of how to find a dataset and link it to an analysis plot or similar.

- A wrapper funtion for readable metadata and a readable snapshot of the station (or specified instrumente or parameter values from it) is required. This could be used to check the status of instruments easily and repeatably or to check the settings of instruments when a particular dataset was taken. This is a precursur to a 'monitor' and links with the idea of being able to add some instrument data to the pngs that get saved on demand. 

- Further outstanding issues include better data analysis package and plotting defaults and ease of use but are more at home as QCoDeS requirements rather than wrapper functions

