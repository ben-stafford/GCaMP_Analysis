The code in the respository is designed to process and analyze files that have been generated by Suite2P following imaging of GCaMP signals. Currently, the readme file only provides a general overview of how to run the accompanying code. With time, I hope to provide more information about the data structures created and output by each of the Matlab files. The general workflow looks something like this:

1) Image responses of GCaMP-expressing cells to various visual stimuli.
2) Convert the raw scanbox files to hdf5 files.
3) Run Suite2P on all the hdf5 files from an imaging session to generate fluorescence traces for each roi across all experiments. This way, responses from each roi/cell to all visual stimuli and under all conditions (i.e. control vs CNO) can be analyzed and compared.
4) Process the mat file generated by Suite2P by running Suite2P_Analysis_Batch.m. This code creates data structures for all rois that are split according to the different experiments (i.e. different visual stimuli presented and/or imaging conditions).
5) Compile tuning properties for all rois from specific subsets of experiments by running Tuning_Analysis.m on the data structures generated in step 4. This code creates data structures containing preferred SF and TF as well as DSI, OSI, and preferred direction for each roi/cell. It is necessary to compile multiple experiments when doing this analysis because the full range of parameters tested are often collected over multiple experiments.
6) Compare rois/cells in control versus CNO by running Control_CNO_Compare_Compile.m.

The code here is designed to analyze data from steps 4-6.

Suite2P_Analysis_Batch.m

This code expects several different types of files contained within a specific directories. To run this code locally, you will need to copy the files to your computer or to an external hard drive from the Two Photon PC. The files are as follows:
 1) The mat file generated by Suite2P containing fluorescence traces for all rois in a directory named: Suite2P
 2) The mat files generated by scanbox from each experiment in a directory named: Scanbox_Files
 3) The analyzer files for each experiment in a directory named: Analyzer_Files
 4) The text files generated by Spike2 for each experiment in a directory named: Spike2_Files

The code expects the directories containing these files to be located in the following directory structure:

Experiment_Name/Recording_Day/Suite2P
Experiment_Name/Recording_Day/Scanbox_Files
Experiment_Name/Recording_Day/Analyzer_Files
Experiment_Name/Recording_Day/Spike2_Files

Variable Names:

When you run Suite2P_Analysis_Batch.m, you specify the location of the files by providing the names of the directories where they are located using the following variables:
1) animal_id: the name of the animal which is typically the mouse ear tag number.
2) experiment_day: typically day1, day2, day3, etc.
3) suite2pDir: the path where the Suite2P mat file is located
4) scanboxDir: the path where the Scanbox_Files are located

There are also four 'flags' representing analysis steps you can 'turn on' by setting to a value of 1. These are:
1) analyze_all_files - when set to 1, the code will use every mat file in the Scanbox_Files directory and use it to analyze the relevant part of the fluorescence traces.
2) smoothTraces - when set to 1, this will smooth the fluoresence traces
3) subtractNeuropil - when set to 1, this will subtract the neuropil signal from the fluorescence traces
4) photoBleachSubtract - when set to 1, this will use the interleaved blank trials to fit the background fluorescnce decay with a polynomial and subtract from the fluorescence traces to try to linearize the dfof trace.

Output:

The code saves data structures in mat files for each experiment in the following location:

Experiment_Name/Recording_Day/Analysis

The output files are named after the scanbox file that was used to generate them. A typical scanbox name is: day1_000_000.mat. The corresponding mat file in the Analysis directory will be named: day1_000_000_Data.mat

Tuning_Analysis.m

This code takes the output files from Suite2P_Analysis_Batch.m and calculates numerous parameters for each ROI. It also allows you to analyze more than one output file at a time. This is required for stimulus sets that were collected over the course of two or more experiments. For example, when analyzing SF and TF tuning, two different experiments are typically run, each of which covers half of the SFs, TFs, and directions required for this analysis. In order to calculate the SF and TF tuning of each ROI, two files must loaded and analyzed sequentially.

When you run Tuning_Analysis.m, you specify the location of the files by providing the names and paths to them using the following variables:
1) animal_id: the name of the animal which is typically the animal's ear tag number.
2) experiment_day: typically day1, day2, day3, etc.
3) files_to_analyze: a list of the matlab files to open and analyze. These are the output files from Suite2P_Analysis_Batch.m and will have names like: day1_000_000_Data.mat, day1_000_001_Data.mat; day1_000_002_Data.mat, etc.
4) analysisDir: the path where the files specified in files_to_analyze are located. This is typically animal_id/experiment_day/Analysis

There are also three flags you can set that impact how the analysis is run and whether output mat files are saved.
1) stationary_trials - when set to 1, the code will only analyze trails when the animal was not running.
2) save_data - when set to 1, the code will save out data structures in a mat file in a directory created and located within the dircetory specified in analysisDir.
3) cno_data - when set to 1, the output mat files will be saved in a directory named 'CNO', otherwise they will be saved in a directory named 'Control'.

   Output

Control_CNO_Compare_Compile.m
