clear

batch_process = 0; save_data = 0; plot_on = 0;

if batch_process == 1
    if ismac
        cd ~/Dropbox/Salk/GCaMP/Data/
    elseif ispc
        cd C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data
    end    
    load('Tlx3_Control_CNO_SF_TF_Direction_Batch_List')
else
    % data = struct();
    % data.animal_ids = {'129';'129';'7823';'7823';'69';'69';'71';'71';'66';'66'};
    % data.experiment_days = {'day2';'day2';'day3';'day3';'day2';'day2';'day2';'day2';'day2';'day2'};
    % % data.analysis_files = {{'_000_003_Data.mat', '_000_004_Data.mat'}}; % speed uses two files
    % % data.analysis_files = {{'_000_002_Data.mat'}}; % direction uses one or three files
    % data.analysis_files = {
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     {'_000_004_Data.mat', '_000_005_Data.mat', '_000_006_Data.mat'};
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     {'_000_004_Data.mat', '_000_005_Data.mat', '_000_006_Data.mat'};
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     {'_000_004_Data.mat', '_000_005_Data.mat', '_000_006_Data.mat'};
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     {'_000_004_Data.mat', '_000_005_Data.mat', '_000_006_Data.mat'};
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     {'_000_004_Data.mat', '_000_005_Data.mat', '_000_006_Data.mat'};
    %     }; % direction SF x TF uses three files
    % % data.stim_type = {'direction'}; % speed or direction
    % data.stim_type = {'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction'}; % speed or direction
    % data.condition = {'control';'cno';'control';'cno';'control';'cno';'control';'cno';'control';'cno'}; % 'cno' or 'control'
    % % data.condition = {'cno'}; % 'cno' or 'control'
    % % data.stim_type = {'sf_tf_direction'}; % uses three files

    data = struct();
    data.animal_ids = {'712'};
    data.experiment_days = {'day1'};
    data.analysis_files = {{'_000_004_Data.mat', '_000_005_Data.mat'}}; % speed uses two files
    % data.analysis_files = {{'_000_002_Data.mat'}}; % direction uses one or three files
    % data.analysis_files = {
    %     {'_000_000_Data.mat', '_000_001_Data.mat', '_000_002_Data.mat'};
    %     }; % direction SF x TF uses three files
    data.stim_type = {'speed'}; % speed, direction, or sf_tf_direction
    % data.stim_type = {'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction';'sf_tf_direction'}; % speed or direction
    % data.condition = {'control';'cno';'control';'cno';'control';'cno';'control';'cno';'control';'cno'}; % 'cno' or 'control'
    data.condition = {'cno'}; % 'cno' or 'control'
    % data.stim_type = {'sf_tf_direction'}; % uses three files
    
end

% Parameters to adjust how analysis is performed
responsive_thresh = 0.05; % p-value for responsiveness ANOVA
dsi_thresh = 0.4;
osi_thresh = 0.2;
%%
for d = 1:length(data.animal_ids)
    %%
    % Get parameters from data structure created above
    animal_id = data.animal_ids{d};
    experiment_day = data.experiment_days{d};
    analysis_files = data.analysis_files{d};
    stim_type = data.stim_type{d};
    condition = data.condition{d};
    
    disp(['Analyzing ', animal_id, ' ', stim_type, ' tuning']);

    files_to_analyze = {};
    for i = 1:length(analysis_files)
        files_to_analyze{i} = [experiment_day,analysis_files{i}];
        disp(files_to_analyze{i});
    end

    % Dropbox
    if ismac
        analysis_dir = ['~/Dropbox/Salk/GCaMP/Data/', animal_id, '/', experiment_day,'/Analysis'];
        save_dir = ['~/Dropbox/Salk/GCaMP/Data/', animal_id, '/'];
    elseif ispc
        analysis_dir = ['C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data/', animal_id, '/', experiment_day,'/Analysis'];
        save_dir = ['C:/Users/SNLC_BS_PC/Dropbox/Salk/GCaMP/Data/', animal_id, '/'];
    end

    % Local to PC
    if ispc
        save_dir_local = ['D:/GCaMP_Data/', animal_id, '/', experiment_day,'/Analysis'];
    end

    cd (analysis_dir)
%%
    % Load the first file to find out how many cells/ROIs there were so that
    % the correct sized data structure can be created
    mat_obj = matfile(files_to_analyze{1});

    % Create a structure to contain the merged data from all the files
    % being analyzed.
    tuning(length(mat_obj.gcamp_data)) = struct();

    % Create variable to keep track of the number of stimuli presented in
    % each file that is appended so subsets of data can be extracted later
    % if needed.
    in_stim_start = 0;
    in_stim_end = 0;
    in_stim_mat = nan(length(files_to_analyze),2);

    clear mat_obj; % clear the matfile object to save memory

    disp('Extracting GCaMP data for all files');

    for f = 1:length(files_to_analyze)
        file = files_to_analyze{f};

        % Load the GCaMP Data
        load (file)

        if f == 1
            % Put the cell ids into the tuning data structure
            [tuning.cell_ids] = deal(gcamp_data.cell_ids);
        end

        % Get the number of stimuli presented by checking the size of the
        % first cell in gcamp_data. Since all cells 'saw' the same stimuli,
        % only the first one needs to be checked. Depending on whether this
        % is the first file, adjust the index values accordingly
        if f == 1
            in_stim_start = 1;
            in_stim_end = size(gcamp_data(1).dfof_traces,1);
        else
            in_stim_start = in_stim_start+in_stim_end;
            in_stim_end = in_stim_end+size(gcamp_data(1).dfof_traces,1)+1;
        end

        in_stim_mat(f,:) = [in_stim_start in_stim_end];

        tuning = extract_gcamp_data(f,gcamp_data,tuning,in_stim_mat,stimulus_info);
        % for g = 1:length(gcamp_data)
        %     tuning(g).stimulus_index_mat = in_stim_mat;
        %     dfof_trace = gcamp_data(g).dfof_traces;
        %     dfof_trace = [stimulus_info dfof_trace];
        % 
        %     % Get the mean, max mean, and max response
        %     max_resp_trials = nan(size(dfof_trace,1),5);
        %     max_mean_resp_trials = nan(size(dfof_trace,1),5);
        %     mean_resp_trials = nan(size(dfof_trace,1),5);
        %     % First get the max response from each trial in dfof_trace
        %     stim_start = dfof_trace(1,2)-dfof_trace(1,1);
        %     for t = 1:size(dfof_trace,1)
        %         tmp_trace = dfof_trace(t,9:end-7);
        %         if dfof_trace(t,6) ~= 500
        %             tmp_max_resp = max(tmp_trace(1,stim_start:end-7)); % max response
        %             tmp_mean_resp = mean(tmp_trace(1,stim_start:end-7)); % mean response
        %             in_max_resp = find(tmp_trace == max(tmp_trace(1,stim_start:end-7))); % minus 7 so max can't be too close to end of trace for averaging
        %             tmp_max_mean_resp = mean(tmp_trace(in_max_resp-7:in_max_resp+7)); % max mean response
        %         else
        %             tmp_max_mean_resp = mean(tmp_trace); % if blank just average the entire 'response'
        %         end
        %         max_resp_trials(t,1:3) = dfof_trace(t,6:8);
        %         max_mean_resp_trials(t,1:3) = dfof_trace(t,6:8);
        %         mean_resp_trials(t,1:3) = dfof_trace(t,6:8);
        %         max_resp_trials(t,4) = dfof_trace(t,5); % make moving signal the fourth column
        %         max_mean_resp_trials(t,4) = dfof_trace(t,5);
        %         mean_resp_trials(t,4) = dfof_trace(t,5);
        %         max_resp_trials(t,5) = tmp_max_resp;
        %         max_mean_resp_trials(t,5) = tmp_max_mean_resp;
        %         mean_resp_trials(t,5) = tmp_mean_resp;
        %     end
        %     if f == 1 % if first file put dfof data in structure
        %         tuning(g).traces = dfof_trace;
        %         tuning(g).all_trials.max_resp = max_resp_trials;
        %         tuning(g).all_trials.max_mean_resp = max_mean_resp_trials;
        %         tuning(g).all_trials.mean_resp = mean_resp_trials;
        %     else % otherwise get the existing dfof data and append the new dfof data to it
        %         if ~(size(dfof_trace) == size(tuning(g).traces))
        %             unequal = find(~(size(dfof_trace) == size(tuning(g).traces)));
        %             for i = 1:length(unequal)
        %                 if unequal(i) == 1 % all the existing data is smaller so needs to be padded
        %                     if size(dfof_trace,1) < size(tuning(g).traces,1)
        %                         dfof_trace(end+1:size(tuning(g).traces,1),:) = missing;
        %                     end
        %                     % create function to recreate the tuning data
        %                     % structure
        %                 elseif unequal(i) == 2 % dfof_trace needs to be adjusted
        %                     if size(dfof_trace,2) < size(tuning(g).traces,2)
        %                         dfof_trace(:,end+1:end+1) = missing;
        %                     elseif size(dfof_trace,2) > size(tuning(g).traces,2)
        %                         dfof_trace = dfof_trace(:,1:end-1);
        %                     end
        %                 end
        %             end
        %         end
        %         tuning(g).traces = [tuning(g).traces; dfof_trace];
        %         tuning(g).all_trials.max_resp = [tuning(g).all_trials.max_resp; max_resp_trials];
        %         tuning(g).all_trials.max_mean_resp = [tuning(g).all_trials.max_mean_resp; max_mean_resp_trials];
        %         tuning(g).all_trials.mean_resp = [tuning(g).all_trials.mean_resp; mean_resp_trials];
        %     end
        % end
    
    end
%%
    clear  in_stim_mat in_stim_start in_stim_end

    disp('Getting stationary and running trials')
    % Get the stationary and running trials and add that data to the data
    % structure
    for t = 1:size(tuning,2)
        tuning = get_stat_run_trials(t,tuning);
        % for n = 1:size(state_names,1) % each dfof field
        %     resp_names = fieldnames(tuning(t).(state_names{n}));
        %     for f = 1:length(resp_names) % each response field
        %         tmp_resp_data = tuning(t).(state_names{n}).(resp_names{f});
        %         stims = unique(tmp_resp_data(:,1:3,1), 'rows');
        %         stim_resp_mat = nan(max_stims,size(stims,1));
        %         for s = 1:size(stims,1)
        %             in_stim = ismember(tmp_resp_data(:,1:3),stims(s,:),'rows');
        %             % subset of dfof traces corresponding to all repeats of
        %             stim_resp = tmp_resp_data(in_stim,:);
        %             stim_resp_mat(1:size(stim_resp,1),s) = stim_resp(:,5);
        %         end
        %         % calculate statistics
        %         [p_anova, anovatab, stats] = anova1('kruskalwallis', stim_resp_mat, [],'off');
        %         % tuning(t).(state_names{n}).([resp_names{f}, '_p_val']) = p_anova;
        %         tuning(t).(state_names{n}).p_vals.(resp_names{f}) = p_anova;
        %     end
        % end
    end
    
    % Get fieldnames for use later
    names = fieldnames(tuning);
    state_names = names(contains(names, 'trials'));
    resp_names = fieldnames(tuning(1).(state_names{1}));
    resp_names = resp_names(contains(resp_names,'resp'));
    resp_names = resp_names(~contains(resp_names,'by_stim'));

    % Check how many trials there are in stationary and running states
    % First get how many total stimulus combinations were presented
    all_stat_trials_by_stim = get_trials_by_stim(tuning, 'stat_trials');
    all_run_trials_by_stim = get_trials_by_stim(tuning, 'run_trials');

    % Now, loop through all the states (all, stat, run) and
    % calculate responsive values for each cell. This is the coarse filter
    % where all trials from all conditions are compared

    % Get max number of repeats of each stimulus condition including
    % blanks
    stim_data = tuning(1).all_trials.max_resp;
    stims = unique(stim_data(:,1:3,1), 'rows');
    max_stims = 0;
    for s = 1:size(stims,1)
        in_stim = ismember(stim_data(:,1:3),stims(s,:),'rows');
        max_stims = max(max_stims, size(stim_data(in_stim,:),1));
    end

    disp('Finding responsive cells')
    for t = 1:size(tuning,2) % each cell
        tuning = get_resp_p_val(t,tuning,state_names,max_stims,responsive_thresh);
    %     for n = 1:size(state_names,1) % each dfof field
    %         resp_names = fieldnames(tuning(t).(state_names{n}));
    %         resp_names = resp_names(~contains(resp_names, 'p_vals')); % resp names
    %         for f = 1:length(resp_names) % each response field
    %             tmp_resp_p_val = tuning(t).(state_names{n}).p_vals.(resp_names{f});
    %             tmp_resp_data = tuning(t).(state_names{n}).(resp_names{f});
    %             stims = unique(tmp_resp_data(:,1:3,1), 'rows');
    %             stim_resp_mat = nan(max_stims,size(stims,1));
    %             for s = 1:size(stims,1)
    %                 in_stim = ismember(tmp_resp_data(:,1:3),stims(s,:),'rows'); % logical index of stimulus conditions
    %                 stim_resp = tmp_resp_data(in_stim,:); % all responses to one stimulus condition
    %                 stim_resp_mat(1:size(stim_resp,1),s) = stim_resp(:,5);
    %             end
    %             tuning(t).(state_names{n}).([resp_names{f}, '_by_stim']) = [stims, nanmean(stim_resp_mat)']; % mean response to each stimulus condition
    %             if tmp_resp_p_val < responsive_thresh % only do second responsive test if initially deemed responsive
    %                 peak_resp_cond = max(nanmean(stim_resp_mat(:,1:end))); % last column are blanks so don't include
    %                 in_peak = find(nanmean(stim_resp_mat(:,1:end-1)) == peak_resp_cond);
    %                 peak_cond_trials = stim_resp_mat(:,in_peak);
    %                 peak_cond_trial_mat = [peak_cond_trials,stim_resp_mat(:,end)];
    %                 p_val = kruskalwallis(peak_cond_trial_mat,[],'off');
    %                 tuning(t).(state_names{n}).p_vals.([resp_names{f}, '_peak']) = p_val;
    %             else
    %                 tuning(t).(state_names{n}).p_vals.([resp_names{f}, '_peak']) = 1; % not significant because cell did not pass initial responsive test
    %             end
    %         end
    %     end  
    end

    % Loop through each field and calculate mean response to each stimulus condition and p-value of 'peak' response
    for t = 1:size(tuning,2) % each cell
        tuning = get_mean_by_stim_and_peak_resp_p_val(t,tuning,state_names,max_stims,responsive_thresh);
    end

    % Print out number of responsive cells for stationary trials
    isresp_stationary = get_isresp_by_state(tuning, 'stat_trials', 'max_mean_resp', 0.05);
    disp(['Number of responsive cells (stationary trials): ',num2str(size(find(isresp_stationary),2))])
    

    % To access the stringent responsiveness p-values for stationary states
    % calculated from the max response:
    % all_stat_data = [tuning(:).stat_trials]; % all the stationary data
    % all_p_vals = [all_stat_data(:).max_resp_peak_p_val]; % all the p-values
    % in_resp_stat = [all_stat_data(:).max_resp_peak_p_val]' < 0.05; % logical index of the significant p-values

    % Pull out responsive cells and calculate dsi and osi. If less than 8
    % directions, use traditional dsi and osi calculation.
   %%
    % Generate SF x TF matrices for each responsive cell. If not
    % responsive, pass an empty matrix to pad the data structure.
    if strcmp(data.stim_type(d),'speed')
        disp('Generating SF x TF matrix')

        % Need to check if any of the stationary or running trials are
        % missing any stimulus combinations (i.e. 0 trials of that state
        % for a specific combination). If a stimulus combination has zero
        % trials, the analysis below will crash.
        state_names_to_analyze = {'all_trials'};
        stat_by_stim = get_trials_by_stim(tuning, 'stat_trials');
        if isempty(find(stat_by_stim(:,4)==0))
            state_names_to_analyze{end+1} = 'stat_trials';
        else
            disp('No stationary trials for some stimulus combinations!')
            disp('No SF x TF fitting will be done for stationary trials')
        end
        run_by_stim = get_trials_by_stim(tuning, 'run_trials');
        if isempty(find(run_by_stim(:,4)==0))
            state_names_to_analyze{end+1} = 'run_trials';
        else
            disp('No running trials for some stimulus combinations!')
            disp('No SF x TF fitting will be done for running trials')
        end
        state_names_to_analyze = state_names_to_analyze';
        for t = 1:size(tuning,2)
            for n = 1:size(state_names_to_analyze,1)
                for f = 1:size(resp_names,1)
                    tmp_resp_data = tuning(t).(state_names_to_analyze{n}).([resp_names{f},'_by_stim']);
                    sf_tf_combs = unique(tmp_resp_data(:,2:3,1), 'rows');
                    dir_combs = unique(tmp_resp_data(:,1,1), 'rows');
                    sf_tf_resp_mat = nan(size(sf_tf_combs,1)-1,3); % one row for each sf x tf combination, and one column each for sf, tf, mean response
                    for s = 1:(size(sf_tf_combs,1)-1) % each sf x tf combination excluding blanks
                        dir_resp_mat = nan((size(dir_combs,1)-1),2); % excluding blanks
                        for r = 1:(size(dir_combs,1)-1) % each direction, excluding blanks
                            in_stim = ismember(tmp_resp_data(:,1:3),[dir_combs(r,:) sf_tf_combs(s,:)],'rows'); % logical index of all trials from one dir x sf x tf combination
                            % subset of responses corresponding to all repeats
                            % of a given stimulus
                            dir_resp_mat(r,1) = dir_combs(r,:);
                            dir_resp_mat(r,2) = tmp_resp_data(in_stim,4);
                        end
                        if size(dir_resp_mat,1) < 8
                            [dsi, osi, pref_dir, pref_resp] = get_trad_dsi_osi(dir_resp_mat(:,2),dir_resp_mat(:,1));
                        else
                            [osi, dsi, pref_dir, L_ori, L_dir] = get_CV_DirCV(dir_resp_mat(:,2),dir_resp_mat(:,1));
                        end
                        if dsi > dsi_thresh || osi > osi_thresh % if a cell was tuned at a given sf x tf, combination, use it's preferred direction
                            sf_tf_dir_resp = pref_resp;
                        else
                            sf_tf_dir_resp = nanmean(dir_resp_mat(:,2)); % if a cell was NOT tuned at a given sf x tf, combination, average all its responses
                        end
                        sf_tf_resp_mat(s,1:2) = [sf_tf_combs(s,:)];
                        sf_tf_resp_mat(s,3) = sf_tf_dir_resp;
                    end
                    tuning(t).(state_names_to_analyze{n}).sf_tf_mat.([resp_names{f},'_sf_tf_matrix']) = sf_tf_resp_mat;
                end
            end
            if ~isequal(state_names,state_names_to_analyze)
                % Check if enough stationary trials
                enough_stat = contains(state_names_to_analyze,state_names{2});
                enough_run = contains(state_names_to_analyze,state_names{3});
                if isempty(find(enough_stat)) % not enough stationary trials so insert empty matrix in structure
                    tuning(t).('stat_trials').sf_tf_mat.([resp_names{f},'_dir_tuning']) = nan; % if cell wasn't responsive pass empty matrix
                elseif isempty(find(enough_run)) % not enough running trials so insert empty matrix in structure
                    tuning(t).('run_trials').sf_tf_mat.([resp_names{f},'_dir_tuning']) = nan; % if cell wasn't responsive pass empty matrix
                end
            end
        end
        dsi_pd = []; % not a direction stimulus set so pass empty matrices for dsi and osi probability dictributions
        osi_pd = [];
    elseif strcmp(data.stim_type(d),'direction') || strcmp(data.stim_type(d),'sf_tf_direction') % if direction tuning experiment, there is only one SF or TF tested, so no SF x TF matrix
        % Need to check if any of the stationary or running trials are
        % missing any stimulus combinations (i.e. 0 trials of that state
        % for a specific combination). If a stimulus combination has zero
        % trials, the analysis below will crash.
        %%
        state_names_to_analyze = {'all_trials'};
        stat_by_stim = get_trials_by_stim(tuning, 'stat_trials');
        if isempty(find(stat_by_stim(:,4)==0))
            state_names_to_analyze{end+1} = 'stat_trials';
        else
            disp('No stationary trials for some stimulus combinations!')
        end
        run_by_stim = get_trials_by_stim(tuning, 'run_trials');
        if isempty(find(run_by_stim(:,4)==0))
            state_names_to_analyze{end+1} = 'run_trials';
        else
            disp('No running trials for some stimulus combinations!')
        end
        state_names_to_analyze = state_names_to_analyze';

        % Get permutation dsi and osi probability distributions.
        pd_dsi_osi_struct = get_dsi_osi_prob_dist(tuning, state_names_to_analyze);

        textprogressbar_2025('Sum of Von Mises Fitting: ');
        pause(1);
        for t = 1:size(tuning,2)
            textprogressbar_2025(t,size(tuning,2));
            for n = 1:size(state_names_to_analyze,1)
                for f = 1:size(resp_names,1)
                    % if tuning(t).(state_names_to_analyze{n}).p_vals.([resp_names{f},'_peak']) < responsive_thresh % if cell is responsive do calculation
                        tmp_resp_data = tuning(t).(state_names_to_analyze{n}).([resp_names{f},'_by_stim']);
                        sf_tf_combs = unique(tmp_resp_data(:,2:3,1), 'rows');
                        dir_combs = unique(tmp_resp_data(:,1,1), 'rows');
                        sf_tf_resp_mat = nan(size(sf_tf_combs,1)-1,9); % one row for each sf x tf combination, and one column each for sf, tf, dsi, osi, preferred direction, preferred response, mean response
                        %%
                        for s = 1:(size(sf_tf_combs,1)-1) % each sf x tf combination excluding blanks
                            dir_resp_mat = nan((size(dir_combs,1)-1),2); % excluding blanks
                            for r = 1:(size(dir_combs,1)-1) % each direction, excluding blanks
                                in_stim = ismember(tmp_resp_data(:,1:3),[dir_combs(r,:) sf_tf_combs(s,:)],'rows'); % logical index of all trials from one dir x sf x tf combination
                                % subset of responses corresponding to all repeats
                                % of a given stimulus
                                dir_resp_mat(r,1) = deg2rad(dir_combs(r,:));
                                dir_resp_mat(r,2) = tmp_resp_data(in_stim,4);
                            end
                            if min(dir_resp_mat(:,2)) < 0
                                dir_resp_mat(:,2) = dir_resp_mat(:,2)-(min(dir_resp_mat(:,2))); % eliminate negative values
                            end
                            % traditional dsi and osi calculation
                            [dsi, osi, pref_dir, pref_resp] = get_trad_dsi_osi(dir_resp_mat(:,2),dir_resp_mat(:,1));
                            % CV dsi osi calculation
                            [cv_osi, cv_dsi, cv_pref_dir, L_ori, L_dir] = get_CV_DirCV(dir_resp_mat(:,2),dir_resp_mat(:,1));
                            [x_data_interp, y_data_interp, x_fit, y_fit,...
                                        vm_pref_dir, vm_fwhm, vm_osi_fit, vm_dsi_fit, resnorm] = fit_sum_von_mises(dir_resp_mat(:,1),dir_resp_mat(:,2),360,0);
                            sf_tf_resp_mat(s,1:2) = [sf_tf_combs(s,:)];
                            sf_tf_resp_mat(s,3:5) = [dsi osi pref_dir];
                            sf_tf_resp_mat(s,6:7) = [vm_pref_dir vm_fwhm];
                            % closest real value to interpolated preferred
                            % direction
                            [~,closest_in] = min(abs(vm_pref_dir-rad2deg(dir_resp_mat(:,1))));
                            pref_resp = dir_resp_mat(closest_in,2);
                            [pref_amp, pref_idx] = max(y_fit);
                            sf_tf_resp_mat(s,8:9) = [pref_resp mean(dir_resp_mat(:,2))];
                            sf_tf_resp_mat(s,10:12) = [cv_osi cv_dsi cv_pref_dir];
                        end
                        %%
                        dir_struct = [];
                        dir_struct.sf_tf = sf_tf_resp_mat(:,1:2);
                        dir_struct.dsi_osi_pdir = sf_tf_resp_mat(:,3:5);
                        dir_struct.cv_dsi_cv_osi_cv_pdir = sf_tf_resp_mat(:,10:12);
                        dir_struct.vm_pdir_fwhm = sf_tf_resp_mat(:,6:7);
                        dir_struct.pref_resp_mean_resp = sf_tf_resp_mat(:,8:9);
                        dir_struct.dsi_pd = pd_dsi_osi_struct.(state_names_to_analyze{n}).dsi_pd;
                        dir_struct.osi_pd = pd_dsi_osi_struct.(state_names_to_analyze{n}).osi_pd;
                        tuning(t).(state_names_to_analyze{n}).dir_mat.([resp_names{f},'_dir_tuning']) = dir_struct;
                end
            end
            if ~isequal(state_names,state_names_to_analyze)
                % Check if enough stationary trials
                enough_stat = contains(state_names_to_analyze,state_names{2});
                enough_run = contains(state_names_to_analyze,state_names{3});
                if isempty(find(enough_stat)) % not enough stationary trials so insert empty matrix in structure
                    dir_struct = [];
                    tuning(t).('stat_trials').dir_mat.([resp_names{f},'_dir_tuning']) = dir_struct; % if cell wasn't responsive pass empty matrix
                elseif isempty(find(enough_run)) % not enough running trials so insert empty matrix in structure
                    dir_struct = [];
                    tuning(t).('run_trials').dir_mat.([resp_names{f},'_dir_tuning']) = dir_struct; % if cell wasn't responsive pass empty matrix
                end
            end
        end
        textprogressbar_2025(' done');
    end

    % Create file name for saving tuning data structure
    file_name = experiment_day;
    for i = 1:length(analysis_files)
        tmp_name = analysis_files{i};
        split_name = split(tmp_name,'_');
        file_name = [file_name,'_',split_name{3}];
    end
    fig_plot_dir = file_name;
    file_name = [file_name,'_',stim_type,'_',condition];
    if strcmp(condition, 'control')
        if ~isfolder([save_dir,'Control/'])
            mkdir([save_dir,'Control/'])
        end
        % save([save_dir,'Control/',file_name],'tuning');
        % disp(['Saving ',file_name,'.mat to ',save_dir,'Control']);
        if ispc
            if ~isfolder([save_dir_local,'/Control/'])
                mkdir([save_dir_local,'/Control/'])
            end
            % save([save_dir_local,'/Control/',file_name],'tuning');
            % disp(['Saving ',file_name,'.mat locally to PC at ',save_dir_local,'/Control']);
        end
    elseif strcmp(condition, 'cno')
        if ~isfolder([save_dir,'CNO/'])
            mkdir([save_dir,'CNO/'])
        end
        % save([save_dir,'CNO/',file_name],'tuning');
        % disp(['Saving ',file_name,'.mat to ',save_dir,'CNO']);
        if ispc
            if ~isfolder([save_dir_local,'/CNO/'])
                mkdir([save_dir_local,'/CNO/'])
            end
            % save([save_dir_local,'/CNO/', file_name],'tuning');
            % disp(['Saving ',file_name,'.mat locally to PC at ',save_dir_local,'/CNO']);
        end
    end
    
    if strcmp(data.stim_type(d),'speed')
        % Go through each state and response type and do sf x tf fitting only
        % for responsive cells
        if ispc % only make directories for saving plots if on PC so they can be saved to big HDD
            % plot_on = 1; % only save figs if analyzing on PC
            if strcmp(condition, 'control')
                path_to_save = [save_dir_local,'/Control/Speed_Tuning_Plots/',fig_plot_dir];
                if ~isfolder([save_dir_local,'/Control/Speed_Tuning_Plots/',fig_plot_dir])
                    mkdir([save_dir_local,'/Control/Speed_Tuning_Plots/',fig_plot_dir]);
                end
            elseif strcmp(condition, 'cno')
                path_to_save = [save_dir_local,'/CNO/Speed_Tuning_Plots/',fig_plot_dir];
                if ~isfolder([save_dir_local,'/CNO/Speed_Tuning_Plots/',fig_plot_dir])
                    mkdir([save_dir_local,'/CNO/Speed_Tuning_Plots/',fig_plot_dir]);
                end
            end
        else
            plot_on = 0; % just in case set to zero so no plots are saved to laptop
            path_to_save = 'null'; % empty string to pass to function below so it works
        end
        % Do the fitting now that we know whether to save figs and where to
        % save them
        for i = 1:size(state_names_to_analyze,1)
            disp(['Fitting SF x TF matrices for ',state_names_to_analyze{i}]);
            for j = 1:size(resp_names,1)
                % isresp_state = get_isresp_by_state(tuning, state_names_to_analyze{i}, resp_names{j}, responsive_thresh);
                % [all_speed_fits, tuning] = get_speed_fits_by_state(tuning, isresp_state, state_names_to_analyze{i}, resp_names{j}, plot_on, 0, path_to_save);
                % Might be better idea to run fitting on all cells, since
                % responsive thresholding may not work that well?
                all_cells = ones(1,size(tuning,2)); % logical index of 1 for each cell to force 2D fitting of all cells
                [all_speed_fits, tuning] = get_speed_fits_by_state(tuning, all_cells, state_names_to_analyze{i}, resp_names{j}, plot_on, 0, path_to_save);
            end
        end
    end

    % Save the tuning data structure. Location dependent on condition type
    % and whether analysis was run on PC.
    if strcmp(condition, 'control')
        save([save_dir,'Control/',file_name],'tuning');
        disp(['Saving ',file_name,'.mat to ',save_dir,'Control']);
        if ispc
            save([save_dir_local,'/Control/',file_name],'tuning');
            disp(['Saving ',file_name,'.mat locally to PC at ',save_dir_local,'/Control']);
        end
    elseif strcmp(condition, 'cno')
        save([save_dir,'CNO/',file_name],'tuning');
        disp(['Saving ',file_name,'.mat to ',save_dir,'CNO']);
        if ispc
            save([save_dir_local,'/CNO/', file_name],'tuning');
            disp(['Saving ',file_name,'.mat locally to PC at ',save_dir_local,'/CNO']);
        end
    end

    % Example plots depending on stimulus type
    if strcmp(data.stim_type(d),'speed')
        %%
        resp_cell_idx = find(get_isresp_by_state(tuning,'stat_trials','max_resp',0.05));
        % Make a heatmap of the first responsive cell
        tmp = tuning(resp_cell_idx(3)).stat_trials.sf_tf_mat.max_resp_sf_tf_matrix;
        sfs = unique(tmp(:,1));
        tmp_mat = nan(4,4);
        for i = 1:size(sfs,1)
            in_sf = ismember(tmp(:,1),sfs(i));
            tmp_mat(5-i,:) = tmp(in_sf,3); % low sf is bottom of y-axis this way
        end
        % tfs = tmp(in_sf,2);
        % figure(1); clf; heatmap(tfs,flipud(sfs),tmp_mat);
        % figure(2); clf; heatmap(flipud(60./tfs),flipud(sfs),fliplr(tmp_mat));
    
        tfs = unique(tmp(:,2));
        tmp_mat_alt = nan(4,4);
        for i = 1:size(tfs,1)
            in_tf = ismember(tmp(:,2),tfs(i));
            tmp_mat_alt(5-i,:) = tmp(in_tf,3); % low tf is bottom of y-axis this way
        end
        sfs = tmp(in_tf,1);
        % figure(3); clf; heatmap(sfs,flipud(tfs),tmp_mat_alt); % frames
        figure(4); clf; heatmap(sfs,flipud(60./tfs),tmp_mat_alt); % hertz
    
        % Need tfs in hertz for fitting to work
        tfs_hz = flipud(60./tfs);
        % And need tfs to increase from top to bottom in matrix
        [X,Y] = meshgrid(sfs, tfs_hz);
        x_data = zeros(size(X,1),size(Y,2),2);
        x_data(:,:,1) = X;
        x_data(:,:,2) = Y;
    
        speed_fit = get_speed_fits_2025(tmp_mat_alt,sfs,tfs_hz,1,1); % data structure with fitted parameters from 2D gaussian fit
        [well_fit, well_fit_eps] = find_well_fit_2025(speed_fit); % logical value indicating whether cell was well fit by 2D gaussian
        
        % Make a bar plot of all the responsive and well fit cells
        well_fit_state = get_well_fit_by_state(tuning(resp_cell_idx), 'stat_trials', 'max_resp');
        resp_tuning = tuning(resp_cell_idx);
        tmp_data = [resp_tuning((well_fit_state)).stat_trials];
        tmp_data_sf_tf = [tmp_data(:).sf_tf_mat];
        tmp_data_fits = [tmp_data_sf_tf(:).max_resp_fits];
        params_to_plot = {'sf_peak', 'sf_pref', 'sf_sigma'};
        ylabels = {'Peak SF (cpd)', 'Preferred SF (cpd)', 'SF Width (octaves)'};
        figure(); clf;
        for j=1:size(params_to_plot,2)
            subplot(1,3,j); hold on;
            one_condition_bar_plot([tmp_data_fits(:).(params_to_plot{j})],8.0,[0.75 0.75 0.75]);
            xl(-0.5,0.5);
            % xticks([0 1]);
            % xticklabels({'Control','CNO'});
            ylabel(ylabels{j});
        end
        %%
    elseif contains(data.stim_type(d),'direction')
        % Plot DSI and OSI as a function of SF or TF, depending on stimulus
        % that was run
        all_stat_data = [tuning(:).stat_trials]; % all the stationary data
        all_p_vals = nan(1,size(all_stat_data,2));
        for i = 1:size(all_stat_data,2)
            all_p_vals(1,i) = all_stat_data(i).p_vals.max_mean_resp_peak;
        end
        in_resp_stat = all_p_vals < 0.05; % logical index of the significant p-values
        in_resp_stat = find(in_resp_stat);
        sf_tf = tuning(in_resp_stat(1)).stat_trials.dir_mat.max_mean_resp_dir_tuning.sf_tf;
        % Find the SF and TF that occur most frequently
        [sf_cnt,sf_vals] = groupcounts(sf_tf(:,1));
        sf_most = sf_vals(find(sf_cnt==max(sf_cnt)));
        [tf_cnt,tf_vals] = groupcounts(sf_tf(:,2));
        tf_most = tf_vals(find(tf_cnt==max(tf_cnt)));
        % Get SF vals from the most frequent TF value, and the TF vals from
        % the most frequent SF value
        sf_in = find(sf_tf(:,2)==tf_most); % positions in dsi_osi and vm_pdir_fwhm matrices corresponding to the sf direction analysis
        tf_in = find(sf_tf(:,1)==sf_most); % positions in dsi_osi and vm_pdir_fwhm matrices corresponding to the tf direction analysis
        % Matrices holding DSI and OSI values as a function of SF
        sf_dsi_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        sf_osi_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        sf_fwhm_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        sf_pdir_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        for i = 1:size(in_resp_stat,2)
            sf_dsi_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.dsi_osi_pdir(sf_in,1)';
            sf_osi_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.dsi_osi_pdir(sf_in,2)';
            sf_fwhm_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.vm_pdir_fwhm(sf_in,2);
            sf_pdir_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.vm_pdir_fwhm(sf_in,1);
        end
        % Matrices holding DSI and OSI values as a function of TF
        tf_dsi_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        tf_osi_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        tf_fwhm_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        tf_pdir_mat = nan(size(in_resp_stat,2),size(sf_in,1));
        for i = 1:size(in_resp_stat,2)
            tf_dsi_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.dsi_osi_pdir(tf_in,1)';
            tf_osi_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.dsi_osi_pdir(tf_in,2)';
            tf_fwhm_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.vm_pdir_fwhm(tf_in,2);
            tf_pdir_mat(i,:) = tuning(in_resp_stat(i)).stat_trials.dir_mat.max_mean_resp_dir_tuning.vm_pdir_fwhm(tf_in,1);
        end

        % DSI plots
        sf_dsi_in = (sf_dsi_mat > 0.1);
        sf_dsi_vals = sf_dsi_mat;
        sf_dsi_vals(~sf_dsi_in) = nan;
        figure(); clf; hold on;
        multiple_independent_condition_bar_plot(sf_dsi_vals,4,[0.75 0.75 0.75]);
        x_vals =(1:1:size(sf_dsi_vals,2));
        xticks(x_vals);
        xticklabels(unique(sf_tf(:,1)));
        ylabel('DSI'); xlabel('SF (cpd)');

        tf_dsi_in = (tf_dsi_mat > 0.1);
        tf_dsi_vals = tf_dsi_mat;
        tf_dsi_vals(~tf_dsi_in) = nan;
        figure(); clf; hold on;
        multiple_independent_condition_bar_plot(tf_dsi_vals,4,[0.75 0.75 0.75]);
        x_vals =(1:1:size(tf_dsi_vals,2));;
        xticks(x_vals);
        xticklabels(unique(sf_tf(:,2)));
        ylabel('DSI'); xlabel('TF (hz)');

        % OSI plots
        sf_osi_in = (sf_osi_mat > 0.1);
        sf_osi_vals = sf_osi_mat;
        sf_osi_vals(~sf_osi_in) = nan;
        figure(); clf; hold on;
        multiple_independent_condition_bar_plot(sf_osi_vals,4,[0.75 0.75 0.75]);
        x_vals =(1:1:size(sf_osi_vals,2));
        xticks(x_vals);
        xticklabels(unique(sf_tf(:,1)));
        ylabel('OSI'); xlabel('SF (cpd)');

        tf_osi_in = (tf_osi_mat > 0.1);
        tf_osi_vals = tf_osi_mat;
        tf_osi_vals(~tf_osi_in) = nan;
        figure(); clf; hold on;
        multiple_independent_condition_bar_plot(tf_osi_vals,4,[0.75 0.75 0.75]);
        x_vals =(1:1:size(tf_osi_vals,2));
        xticks(x_vals);
        xticklabels(unique(sf_tf(:,2)));
        ylabel('OSI'); xlabel('TF (hz)');

    end

    clear tuning dir_struct gcamp_data pd_dsi_osi_struct all_stat_data
end