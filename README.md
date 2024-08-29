The code in this respository is designed to process and analyze files that have been generated by Suite2P following imaging of GCaMP signals. The files in the 'library' directory are functions upon which the rest of the files depend and should be added locally somewhere within the current search path in Matlab on your computer. The files in the top level directory can be downloaded and stored wherever you want on your computer.

Currently, the readme file only provides a general overview of how to run the accompanying code. With time, I hope to provide more information about the data structures created and output by each of the Matlab files.

The general workflow from experiment to analysis looks something like this:

1) Image responses of GCaMP-expressing cells to various visual stimuli.
2) Convert the raw scanbox files to hdf5 files.
3) Run Suite2P on all the hdf5 files from an imaging session to generate fluorescence traces for each roi across all experiments. This way, responses from each roi/cell to all visual stimuli and under all conditions (i.e. control vs CNO) can be analyzed and compared.
4) Process the mat file generated by Suite2P by running Suite2P_Analysis_Batch.m. This code creates data structures for all rois that are split according to the different experiments (i.e. different visual stimuli presented and/or imaging conditions).
5) Compile tuning properties for all rois from specific subsets of experiments by running Tuning_Analysis.m on the data structures generated in step 4. This code creates data structures containing preferred SF and TF as well as DSI, OSI, and preferred direction for each roi/cell. It is necessary to compile multiple experiments when doing this analysis because the full range of parameters tested are often collected over multiple experiments.
6) Compare rois/cells in control versus CNO by running Control_CNO_Compare_Compile.m.

The code here runs the analyzis steps outlined on steps 4-6.

Suite2P_Analysis_Batch.m

This code expects several different types of files contained within specific directories. To run this code locally, you will need to copy the files to your computer or to an external hard drive from the Two Photon PC. The files are as follows:
 1) The mat file generated by Suite2P containing fluorescence traces for all rois in a directory named: Suite2P
 2) The mat files generated by scanbox from each experiment in a directory named: Scanbox_Files
 3) The analyzer files for each experiment in a directory named: Analyzer_Files
 4) The text files generated by Spike2 for each experiment in a directory named: Spike2_Files

The code expects the directories containing these files to be located within the following directory structure where Animal_ID is the ID of the animal (typically its ear tag number) and Experiment_Day is the day of the recording for the animal (eg. day1, day2, etc.):

Animal_ID/Experiment_Day/Suite2P

Animal_ID/Experiment_Day/Scanbox_Files

Animal_ID/Experiment_Day/Analyzer_Files

Animal_ID/Experiment_Day/Spike2_Files

Variable Names:

When you run Suite2P_Analysis_Batch.m, you manually specify the location of the files by providing the names of the directories where they are located using the following variables:
1) animal_id: the name of the animal which is typically the mouse ear tag number.
2) experiment_day: typically day1, day2, day3, etc.
3) suite2pDir: the path where the Suite2P mat file is located
4) scanboxDir: the path where the Scanbox_Files are located

There are also four 'flags' representing analysis steps you can 'turn on' by setting to a value of 1. These are:
1) analyze_all_files - when set to 1, the code will use every mat file in the Scanbox_Files directory and use it to analyze the relevant part of the fluorescence traces. When set to 0, the code will ask you to select a specific scanbox file to use for the analysis. Typically, you will want this set to 1.
2) smoothTraces - when set to 1, this will smooth the fluoresence traces
3) subtractNeuropil - when set to 1, this will subtract the neuropil signal from the fluorescence traces
4) photoBleachSubtract - when set to 1, this will use the interleaved blank trials to fit the background fluorescence decay with a polynomial and subtract it from the fluorescence traces to try to linearize the dfof trace.

Output:

The code saves data structures in mat files for each experiment in the following location:

Experiment_Name/Recording_Day/Analysis

The output files are named after the scanbox file that was used to generate them. A typical scanbox name is: day1_000_000.mat. The corresponding mat file in the Analysis directory will be named: day1_000_000_Data.mat

Tuning_Analysis.m

This code takes the output files from Suite2P_Analysis_Batch.m and calculates numerous parameters for each ROI. It also allows you to analyze more than one output file at a time. This is required for stimulus sets that were collected over the course of two or more experiments. For example, when analyzing SF and TF tuning, two different experiments are typically run, each of which covers half of the SFs, TFs, and directions required for this analysis. In order to calculate the SF and TF tuning of each ROI, two files must be loaded and analyzed sequentially. The code should be robust enough to determine how many SFs, TFs, and directions were used in the experiments and analyze them accordingly. However, it is important that you only analyze data from one 'set' of experiments at a time. For example, don't try to analyze an experiment with multiple SFs, TFs, and directions (a SF and TF tuning experiment) with an experiment with one SF, one TF, and multiple directions (a DS and OS tuning experiment).

When you run Tuning_Analysis.m, you manually specify the location and names of the files to be anlyzed using the following variables:
1) animal_id: the name of the animal which is typically the animal's ear tag number.
2) experiment_day: typically day1, day2, day3, etc.
3) files_to_analyze: a list of the matlab files to open and analyze. These are the output files from Suite2P_Analysis_Batch.m and will have names like: day1_000_000_Data.mat, day1_000_001_Data.mat; day1_000_002_Data.mat, etc.
4) analysisDir: the path where the files specified in files_to_analyze are located. This is typically animal_id/experiment_day/Analysis

There are also five flags you can set that impact how the analysis is run and whether output mat files are saved.
1) stationary_trials - when set to 1, the code will only analyze trails when the animal was not running.
2) responsive_thresh - p-value below which ROIs will be considered responsive.
3) plot_cell_traces - when set to 1, the code will plot the preferred response mean dfof trace of each responsive and reliable roi. Currently, this is plotted in a tiled plot that can display up to 80 rois. If you end up having more rois than this, you will need to manually change the number plots in the subplot call in the code. This can be a good sanity check of the rois the code is defining as responsive and reliable.
4) save_data - when set to 1, the code will save out data structures in a mat file in a directory created and located within the directory specified in analysisDir.
5) cno_data - when set to 1, the output mat files will be saved in a directory named 'CNO' and have '_cno' appended to the file name, otherwise they will be saved in a directory named 'Control' and have '_control' appended to the file name.

Output:

Depending on the experiments being analyzed, the code will produce different plots. If you are analyzing a SF and TF tuning set of experiments, a heatmap of the average response of all rois to the range of SFs and TFs wil be generated. If you are analyzing a DS and OS experiment, a plot of the mean +/- sem DSI and OSI, as well as a histogram of the preferred direction for all rois will be generated.

When save_data = 1, the code saves data structures in mat files. The names of these output files depend on how many files were included in files_to_analyze. For each file analyzed, the code extracts the second set of three numbers in the file name (these correspond to the experiment number) and includes it in the file name. If only one file is analyzed, the output file will have a name like:

day1_000_control.mat

If two files are analyzed, the output file will have a name like:

day1_000_002_control.mat

Control_CNO_Compare_Compile.m

This code is still a work in progress, but is fully functional for comparing one control dataset against one CNO dataset. This code will load the output files from Tuning_Analysis.m and compare various parameters between control and CNO conditions. It is robust enough to determine whether the underlying experiments were probing SF and TF tuning or DS and OS of the ROIs. Currently, you can only input one control file and one CNO file. It should go without saying, the control and CNO files should be the same type of experiment (i.e. don't compare a control SF and TF tuning experiment with a CNO DS and OS tuning experiment).

When you run Control_CNO_Compare_Compile.m there are four variables you must manually set that provide the location and names of the files to be analyzed:
1) animal_ids: a list of animal names which are typically the mouse ear tag number. Currently, these have to be the same.
2) experiment_days: a list of experiment days (i.e. day1). Currently, these have to be the same.
3) files_to_compare: a list of files to compare. The first file should be the control data (eg. day1_000_001_control.mat) and the second should be the CNO data (eg. day1_002_003_cno.mat).
4) data_dir: the path to where the mat files are located. This will just be the directory that contains the 'Control' and 'CNO' directories created by Tuning_Analysis.m. The code will load the the data files from the Control or CNO directories automatically.

There is also one flag (compare) that should not be changed at the moment. Leave it set to 1.

Output:

Depending on the type of the experiments being analyzed, the code will produce different plots. If you are anlyzing a SF and TF tuning set of experiments, various plots of the preferred SF, TF, and TF/SF ratio in control and CNO conditions will be generated. If you are analyzing a DS and OS tuning experiment, plots of the DSI and OSI in control and CNO conditions will be generated.
