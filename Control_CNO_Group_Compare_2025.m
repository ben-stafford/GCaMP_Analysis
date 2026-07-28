%%
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

genotype = 'Nr5a1';
experiment_type = '_No_DREADDs'; % Empty for normal experiment, '_No_DREADDs', or '_No_CNO'
stimulus_type = 'Direction'; % Speed or Direction
if isempty(experiment_type)
    analysis_type = [genotype];
    load(['Control_CNO_',stimulus_type,'_',genotype,'_List_New']);
else
    analysis_type = [genotype, experiment_type];
    load(['Control',experiment_type,'_',stimulus_type,'_',genotype,'_List_New']);
end
disp(['Analyzing ',genotype,' Data']);

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
%%
% Loop through all_cont_tuning and all_cno_tuning and analyze individual
% datasets separately
% Create data structure to hold the responsive and well fit cells from both
% control and CNO conditions separated out by each dataset analyzed
all_subset_tuning(size(dataset_sizes,2)) = struct();
start_in = 1;
last_in = 0;
for j = 1:size(dataset_sizes,2)
    if j==1
        last_in = dataset_sizes(j);
    else
        last_in = last_in+dataset_sizes(j);
    end
    data_in = start_in:last_in;
    sub_cont_tuning = all_cont_tuning(data_in);
    sub_cno_tuning = all_cno_tuning(data_in);

    cont_isresp_state = get_isresp_by_state(sub_cont_tuning, 'stat_trials', 'max_mean_resp', 0.05);
    cont_wellfit_state = get_well_fit_by_state(sub_cont_tuning,'stat_trials','max_mean_resp');
    cno_isresp_state = get_isresp_by_state(sub_cno_tuning, 'stat_trials', 'max_mean_resp', 0.05);
    control_cno_resp_in = get_logical_comparison_mat(cont_isresp_state,cno_isresp_state);

    n_cont_resp = size(find(control_cno_resp_in(:,1)),1); % number of control responsive cells
    n_both_resp = size(find(control_cno_resp_in(:,3)==2),1); % number of control responsive cells that stayed responsive in cno
    
    resp_both_in = control_cno_resp_in(:,3)==2; % logical index of cells responsive in both conditions
    
    cont_data = [sub_cont_tuning(:).stat_trials];
    cont_both_data = cont_data(resp_both_in);
    cont_both_sf_tf = [cont_both_data(:).sf_tf_mat];
    cont_both_fits = [cont_both_sf_tf(:).max_mean_resp_fits];
    % logical index of the control cells that were responsive in both
    % conditions and well fit
    cont_iswellfit_both = [cont_both_fits(:).well_fit]; % currently unused
    % All of the control cells there were responsive in both conditions and
    % well fit in control
    cont_both_wellfit = cont_both_fits([cont_both_fits(:).well_fit]==1);
        
    cno_data = [sub_cno_tuning(:).stat_trials];
    cno_both_data = cno_data(resp_both_in);
    cno_both_sf_tf = [cno_both_data(:).sf_tf_mat];
    cno_both_fits = [cno_both_sf_tf(:).max_mean_resp_fits];
    % logical index of the CNO cells that were responsive in both
    % conditions and well fit
    cno_iswellfit_both = [cno_both_fits(:).well_fit]; % currently unused
    % All of the CNO cells there were responsive in both conditions and
    % well fit in CNO
    cno_both_wellfit = cno_both_fits([cno_both_fits(:).well_fit]==1);

    % All of the cells that were responsive and well fit in both control
    % and CNO
    cont_cno_resp_wellfit_both_in = get_logical_comparison_mat([cont_both_fits(:).well_fit],[cno_both_fits(:).well_fit]);
    n_cont_resp_wellfit = size(find(cont_cno_resp_wellfit_both_in(:,1)),1); % number of control cells 
    n_cno_resp_wellfit = size(find(cont_cno_resp_wellfit_both_in(:,2)),1); % number of cno cells
    n_both_resp_wellfit = size(find(cont_cno_resp_wellfit_both_in(:,3)==2),1);
    resp_wellfit_both_in = cont_cno_resp_wellfit_both_in(:,3)==2; % logical index of cells responsive and well fit in both conditions

    % Create data structure of cells that were responsive and well fit in
    % both conditions from the current dataset being analyzed
    tmp_cont_data = [cont_both_fits(resp_wellfit_both_in)];
    tmp_cno_data = [cno_both_fits(resp_wellfit_both_in)];

    % Place in data structure containing data from each dataset analyzed
    all_subset_tuning(j).cont_resp_wellfit = tmp_cont_data;
    all_subset_tuning(j).cno_resp_wellfit = tmp_cno_data;

    start_in = last_in+1; % at the end of the loop do this
end
%%
% This calculates means for each animal, and then compares means of means
params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma', 'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)', 'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
cont_mean_mat = nan(size(all_subset_tuning,2),size(params_to_plot,2));
cno_mean_mat = nan(size(cont_mean_mat));
% for k = 1:size(params_to_plot,2)
%     for j = 1:size(all_subset_tuning,2)
%         tmp_cont_data = all_subset_tuning(j).cont_resp_wellfit;
%         tmp_cont_param = [tmp_cont_data(:).(params_to_plot{k})];
%         mean_cont_param = mean(tmp_cont_param);
%         cont_mean_mat(j,k) = mean_cont_param;
%         tmp_cno_data = all_subset_tuning(j).cno_resp_wellfit;
%         tmp_cno_param = [tmp_cno_data(:).(params_to_plot{k})];
%         mean_cno_param = mean(tmp_cno_param);
%         cno_mean_mat(j,k) = mean_cno_param;
%     end
% 
%     subplot(1,3,j); hold on;
%     two_condition_bar_plot([cont_both_fits(:).(params_to_plot{j})],[cno_both_fits(:).(params_to_plot{j})],1.0,[0.75 0.75 0.75]);
%     xl(-0.5,1.5);
%     xticks([0 1]);
%     xticklabels({'Control','CNO'});
%     ylabel(ylabels{j});
% end


% This combines all cells across all animals and performs analysis
% logical index of responsive cells in control and cno
cont_isresp_state = get_isresp_by_state(all_cont_tuning, 'stat_trials', 'max_mean_resp', 0.05);
cno_isresp_state = get_isresp_by_state(all_cno_tuning, 'stat_trials', 'max_mean_resp', 0.05);
% Compare logical values
control_cno_resp_in = [cont_isresp_state', cno_isresp_state'];
% third column will determine if responsive only control (1), only cno
% (-1), both (2), or none (nan)
control_cno_resp_in = [control_cno_resp_in, (nan(size(control_cno_resp_in,1),1))];
for i=1:size(control_cno_resp_in,1)
    if control_cno_resp_in(i,1)==0 && control_cno_resp_in(i,2)==0 % none responsive
        control_cno_resp_in(i,3) = nan;
    elseif control_cno_resp_in(i,1)==1 && control_cno_resp_in(i,2)==0 % only control responsive
        control_cno_resp_in(i,3) = 1;
    elseif control_cno_resp_in(i,1)==1 && control_cno_resp_in(i,2)==1 % both responsive
        control_cno_resp_in(i,3) = 2;
    elseif control_cno_resp_in(i,1)==0 && control_cno_resp_in(i,2)==1 % only cno responsive
        control_cno_resp_in(i,3) = -1;
    end
end

n_cont_resp = size(find(control_cno_resp_in(:,1)),1); % number of control responsive cells
n_both_resp = size(find(control_cno_resp_in(:,3)==2),1); % number of control responsive cells that stayed responsive in cno
n_cno_resp = size(find(control_cno_resp_in(:,2)),1); % number of cno responsive cells
n_only_cont_resp = size(find(control_cno_resp_in(:,3)==1),1); % number of cells responsive only in control
n_only_cno_resp = size(find(control_cno_resp_in(:,3)==-1),1); % number of cells responsive only in cno

resp_both_in = control_cno_resp_in(:,3)==2; % logical index of cells responsive in both conditions
resp_cont_only_in = control_cno_resp_in(:,3)==1; % logical index of cells responsive only in control
resp_cno_only_in = control_cno_resp_in(:,3)==-1; % logical index of cells responsive only in cno

cont_data = [all_cont_tuning(:).stat_trials];
cont_both_data = cont_data(resp_both_in);
cont_only_data = cont_data(resp_cont_only_in);
cont_both_sf_tf = [cont_both_data(:).sf_tf_mat];
cont_both_fits = [cont_both_sf_tf(:).max_mean_resp_fits];
cont_only_sf_tf = [cont_only_data(:).sf_tf_mat];
cont_only_fits = [cont_only_sf_tf(:).max_mean_resp_fits];
cont_only_wellfit = [cont_only_fits([cont_only_fits.well_fit]==1)];
cont_all_data = cont_data(cont_isresp_state);
cont_all_sf_tf = [cont_all_data(:).sf_tf_mat];
cont_all_fits = [cont_all_sf_tf(:).max_mean_resp_fits];
cont_all_wellfit = [cont_all_fits([cont_all_fits.well_fit]==1)];

cont_peak_sf_tf = get_peak_sf_tf(cont_both_data, 'max_mean_resp');

cno_data = [all_cno_tuning(:).stat_trials];
cno_both_data = cno_data(resp_both_in);
cno_only_data = cno_data(resp_cno_only_in);
cno_both_sf_tf = [cno_both_data(:).sf_tf_mat];
cno_both_fits = [cno_both_sf_tf(:).max_mean_resp_fits];
cno_only_sf_tf = [cno_only_data(:).sf_tf_mat];
cno_only_fits = [cno_only_sf_tf(:).max_mean_resp_fits];
cno_only_wellfit = [cno_only_fits([cno_only_fits.well_fit]==1)];
cno_all_data = cno_data(cno_isresp_state);
cno_all_sf_tf = [cno_all_data(:).sf_tf_mat];
cno_all_fits = [cno_all_sf_tf(:).max_mean_resp_fits];
cno_all_wellfit = [cno_all_fits([cno_all_fits.well_fit]==1)];

cno_peak_sf_tf = get_peak_sf_tf(cno_both_data, 'max_mean_resp');

params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma', 'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)', 'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,size(ylabels,2),j); hold on;
    two_independent_condition_bar_plot([cont_all_fits(:).(params_to_plot{j})],[cno_all_fits(:).(params_to_plot{j})],4.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'All Control','All CNO'});
    ylabel(ylabels{j});
end

params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma', 'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)', 'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,size(ylabels,2),j); hold on;
    two_independent_condition_bar_plot([cont_only_fits(:).(params_to_plot{j})],[cno_only_fits(:).(params_to_plot{j})],4.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Only Control','Only CNO'});
    ylabel(ylabels{j});
end

% Make mutliple plots
% Line Plots
params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    two_condition_plot([cont_both_fits(:).(params_to_plot{j})],[cno_both_fits(:).(params_to_plot{j})],1.0,[0.75 0.75 0.75]);
    plot(zeros(length([cont_both_fits(:).(params_to_plot{j})])),[cont_both_fits(:).(params_to_plot{j})],'ko','MarkerSize',8.0);
    plot(ones(length([cno_both_fits(:).(params_to_plot{j})])),[cno_both_fits(:).(params_to_plot{j})],'ro','MarkerSize',8.0);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end

params_to_plot = {'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    two_condition_plot([cont_both_fits(:).(params_to_plot{j})],[cno_both_fits(:).(params_to_plot{j})],1.0,[0.75 0.75 0.75]);
    plot(zeros(length([cont_both_fits(:).(params_to_plot{j})])),[cont_both_fits(:).(params_to_plot{j})],'ko','MarkerSize',8.0);
    plot(ones(length([cno_both_fits(:).(params_to_plot{j})])),[cno_both_fits(:).(params_to_plot{j})],'ro','MarkerSize',8.0);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end

% Bar Plots
params_to_plot = {'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    two_condition_bar_plot([cont_both_fits(:).(params_to_plot{j})],[cno_both_fits(:).(params_to_plot{j})],1.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end

params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    two_condition_bar_plot([cont_both_fits(:).(params_to_plot{j})],[cno_both_fits(:).(params_to_plot{j})],1.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end

% Plotting peak SF and TF for all cells responsive in control and CNO
% High SF tuned
figure(); clf; hold on;
all_cont_data = cont_peak_sf_tf((cont_peak_sf_tf(:,1) >= 0.04),1)';
all_cno_data = cno_peak_sf_tf((cont_peak_sf_tf(:,1) >= 0.04),1)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,0.1);
yticks([0.01 0.02 0.04 0.08]);
xticklabels({'Control','CNO'});
ylabel('Peak SF (cpd)');
% Low SF tuned
figure(); clf; hold on;
all_cont_data = cont_peak_sf_tf((cont_peak_sf_tf(:,1) <= 0.02),1)';
all_cno_data = cno_peak_sf_tf((cont_peak_sf_tf(:,1) <= 0.02),1)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,0.1);
yticks([0.01 0.02 0.04 0.08]);
xticklabels({'Control','CNO'});
ylabel('Peak SF (cpd)');
% All SF tuned
figure(); clf; hold on;
all_cont_data = cont_peak_sf_tf(:,1)';
all_cno_data = cno_peak_sf_tf(:,1)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,0.1);
yticks([0.01 0.02 0.04 0.08]);
xticklabels({'Control','CNO'});
ylabel('Peak SF (cpd)');

% Low TF tuned
figure(); clf; hold on;
all_cont_data = 60./cont_peak_sf_tf((cont_peak_sf_tf(:,2) >= 60),2)';
all_cno_data = 60./cno_peak_sf_tf((cont_peak_sf_tf(:,2) >= 60),2)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,4.25);
yticks([0.5 1 2 4]);
xticklabels({'Control','CNO'});
ylabel('Peak TF (hz)');
% High TF tuned
figure(); clf; hold on;
all_cont_data = 60./cont_peak_sf_tf((cont_peak_sf_tf(:,2) <= 30),2)';
all_cno_data = 60./cno_peak_sf_tf((cont_peak_sf_tf(:,2) <= 30),2)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,4.25);
yticks([0.5 1 2 4]);
xticklabels({'Control','CNO'});
ylabel('Peak TF (hz)');
% All TF tuned
figure(); clf; hold on;
all_cont_data = 60./cont_peak_sf_tf(:,2)';
all_cno_data = 60./cno_peak_sf_tf(:,2)';
two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
xl(-0.5,1.5);
xticks([0 1]);
yl(0,4.25);
yticks([0.5 1 2 4]);
xticklabels({'Control','CNO'});
ylabel('Peak TF (hz)');

params_to_plot = {'tf_peak', 'tf_pref', 'tf_sigma'};
ylabels = {'Peak TF (hz)', 'Preferred TF (hz)', 'TF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    % control data
    all_cont_data = [];
    all_cno_data = [];
    for k=1:size(all_subset_tuning,2)
        all_cont_data = [all_cont_data,[all_subset_tuning(k).cont_resp_wellfit.(params_to_plot{j})]];
        all_cno_data = [all_cno_data,[all_subset_tuning(k).cno_resp_wellfit.(params_to_plot{j})]];
    end
    two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end
sgtitle('Responsive and Well-Fit');

params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma'};
ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)'};
figure(); clf;
for j=1:size(params_to_plot,2)
    subplot(1,3,j); hold on;
    % control data
    all_cont_data = [];
    all_cno_data = [];
    for k=1:size(all_subset_tuning,2)
        all_cont_data = [all_cont_data,[all_subset_tuning(k).cont_resp_wellfit.(params_to_plot{j})]];
        all_cno_data = [all_cno_data,[all_subset_tuning(k).cno_resp_wellfit.(params_to_plot{j})]];
    end
    two_condition_bar_plot(all_cont_data,all_cno_data,1.0,[0.75 0.75 0.75]);
    xl(-0.5,1.5);
    xticks([0 1]);
    xticklabels({'Control','CNO'});
    ylabel(ylabels{j});
end
sgtitle('Responsive and Well-Fit');
