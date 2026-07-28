clear

if ismac
    data_path = '~/Dropbox/Salk/GCaMP/Data/';
elseif ispc
    data_path = 'C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data';
end

batch_process = 0; save_data = 1;

if ismac
    cd ~/Dropbox/Salk/GCaMP/Data/
elseif ispc
    cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
end

genotype = 'Tlx3';
experiment_type = '_No_DREADDs'; % Normal Experiment: '_CNO'; Control No DREADDs Experiment: '_No_DREADDs'
stimulus_type = 'SF_TF_Direction'; % 'Speed', 'SF_TF_Direction', or 'Direction' (this uses one DS OS dataset from a speed tuning imaging session)
analysis_type = [genotype, experiment_type];
disp(['Loading ',genotype,experiment_type,'_',stimulus_type,' Data']);
% load(['Control',experiment_type,'_',stimulus_type,'_',genotype,'_List']);
load(['Control_No_DREADDs_',stimulus_type,'_',genotype,'_List'])

% I think I want to actually calculate means for each dataset and then
% report standard error of the mean, rather than lumping all cells from all
% animals together.

dataset_sizes = [];

for a = 1:length(data.animal_ids)
    cd(data.animal_ids{a}); % cd to the animal id directory
    tmp_files = data.file_names{a}; % cell array of control (1) and CNO (2) files
    % Load the first set of control and cno data
    control_file = tmp_files{1};
    cd 'Control';
    load(control_file);
    cont_tuning = tuning;
    clear tuning;
    cd ..
    cno_file = tmp_files{2};
    cd 'CNO'
    load(cno_file);
    cno_tuning = tuning;
    cd ..
    clear tuning

    if a == 1
        all_cont_tuning = cont_tuning;
        all_cno_tuning = cno_tuning;
    else
        all_cont_tuning = [all_cont_tuning,cont_tuning];
        all_cno_tuning = [all_cno_tuning,cno_tuning];
    end

    dataset_sizes = [dataset_sizes, size(cont_tuning,2)];
    
    clear cont_tuning cno_tuning

    cd(data_path);
end
%%
% Add animal ids as a separate field to the data structures
start_idx = 1;
end_idx = 0;
for j = 1:size(dataset_sizes,2)
    end_idx = end_idx + dataset_sizes(j);
    [all_cont_tuning(start_idx:end_idx).animal_ids] = deal(data.animal_ids(j));
    [all_cno_tuning(start_idx:end_idx).animal_ids] = deal(data.animal_ids(j));
    start_idx = start_idx+dataset_sizes(j);
end

% Save the control and cno data
disp('Saving combined control and cno data...')
if ~isfolder([data_path,'/',analysis_type,'_Analysis'])
    mkdir([data_path,'/',analysis_type,'_Analysis']);
end
cd([data_path,'/',analysis_type,'_Analysis'])
save(['All_',analysis_type,'_',stimulus_type,'_Data.mat'],'all_cont_tuning','all_cno_tuning');
