% cd /Volumes/EPHYS_000/GCaMP_Data/6510/day1/day1_000/Analysis/

clear all

animal_id = '6948';
experiment_day = 'day2';

files_to_analyze = {    
    [experiment_day,'_000_000_Data.mat'];
    [experiment_day,'_000_002_Data.mat'];

    
%     'day2_000_003_Data.mat';
};

% Flag to indicate whether to save the data and whether the data is from control or CNO conditions
save_data = 0;
cno_data = 0;

% Flag to indicate whether you want to analyze just stationary trials or
% all trials
stationary_trials = 1;
responsive_thresh = 0.05;
plot_cell_traces = 0;

% Path to data

% Local Directory
% analysisDir = ['~/Dropbox/Salk/GCaMP/GCaMP_Data/', animal_id, '/', experiment_day,...
%     '/Analysis'];

% External Hard Drive
analysisDir = ['/Volumes/EPHYS_000/GCaMP_Data/', animal_id, '/', experiment_day,...
    '/Analysis'];

% Flag to turn on and off making heatmaps of the mean response of all cells
% to each direction
plot_heatmap_all_directions = 0;

cd (analysisDir)

% Load the first file to find out how many cells/ROIs there were so that
% the correct sized data structure can be created
load (files_to_analyze{1})

% Create data structure to hold the output of this analysis.
tuning(length(gcamp_data)) = struct();

% Create columns to hold the data from the analysis
tuning(1).cell_ids = [];

% Create cell array to hold the names of the files being analyzed
file_names = {};

all_stim_cond_cell = {};
file_stim_cond_struct = struct();

for f = 1:length(files_to_analyze)
    
    file = files_to_analyze{f};
    
    % Load the GCaMP Data
    load (file)
    
    short_name = [file(1:4),file(9:12)];
    
    % Data structure field names
    sf_tf_dir_mean_traces = [file(1:4), file(9:12), '_mean_traces'];
    sf_tf_dir_all_traces = [file(1:4), file(9:12), '_sf_tf_dir_all_traces'];
    sf_tf_tuning = [file(1:4), file(9:12), '_sf_tf_means'];
    sf_tf_all_means = [file(1:4), file(9:12), '_sf_tf_all_mean_resp'];
    blank_field = [file(1:4), file(9:12), '_blank_means'];
    max_dfof_field = [file(1:4), file(9:12), '_max_dfof_trials'];
    
    file_names{f} = [file(1:4), file(9:12)];
    
    % Add field to the data structure
    tuning(1).(sf_tf_dir_mean_traces) = [];
    tuning(1).(sf_tf_dir_all_traces) = [];
    tuning(1).(sf_tf_tuning) = [];
    tuning(1).(sf_tf_all_means) = [];
    tuning(1).(blank_field) = [];
    tuning(1).(max_dfof_field) = [];
    
    % Flags to set what type of stimulus set was run
    sf_analysis = 0;
    tf_analysis = 0;

    % Get unique stimulus values and exclude blanks
    sfs = unique(stimulus_info(:,7));
    sfs = sfs(sfs < 500);
    tfs = unique(stimulus_info(:,8));
    tfs = tfs(tfs < 500);
    dirs = unique(stimulus_info(:,6));
    blank = dirs(dirs == 500);
    dirs = dirs(dirs < 500);
    
    % Determine how many repeats of each stimulus condition
    n_repeats = length(find(stimulus_info(:,7) == sfs(1) & stimulus_info(:,8) == tfs(1) & stimulus_info(:,6) == dirs(1)));
    n_blanks = length(find(stimulus_info(:,7) == blank & stimulus_info(:,8) == blank & stimulus_info(:,6) == blank));
    
    % Create a nSFs, nTFs, nDirs, 3 matrix to hold the stimulus conditions
    % for each cell
    stim_cond_mat = zeros(length(sfs), length(tfs), length(dirs), 3); % matrix of all stimulus combinations
    
    % For each cell in the current file, get the mean response for all
    % repeats of each stimulus condition as well as the mean dfof trace of
    % all repeats and save out into the mean_dfof_mean_mat and mean_dfof_trace_mat
    % matrices respectively
    
    for i = 1:length(gcamp_data)
%     for i = 1:10 % for testing
        
        % Create a nSFs, nTFs, nDirs, length of dfof trace matrix to hold the
        % data from each cell.
        mean_dfof_trace_mat = nan(length(sfs), length(tfs), length(dirs), size(gcamp_data(1).dfof_traces,2));
        
        % Creat a nSFs, nTFs, nDirs, nRepeats, length of dfof matrix to
        % hold all the dfof traces from each stimulus condition for each
        % cell
        all_dfof_trace_cond_mat = nan(length(sfs), length(tfs), length(dirs), n_repeats, size(gcamp_data(1).dfof_traces,2));

        % Create a nSFs, nTFs, nDirs matrix to hold the mean of the mean dfof
        % of the response to each stimlus set from each cell.
        mean_dfof_mean_mat = nan(length(sfs), length(tfs), length(dirs));

        % Create a nSFs, nTFs, nDirs, nRepeats matrix to hold all of the mean
        % dfof response to each stimulus set from each cell
        dfof_mean_mat = nan(length(sfs), length(tfs), length(dirs), n_repeats);

        % Create matrix to hold mean response to blanks
        blank_mean_mat = nan(1,n_blanks);
        
        % Put cell ids in the tuning structure
        tuning(i).cell_ids = gcamp_data(i).cell_ids;
        % Get data from gcamp_data
        all_means = gcamp_data(i).mean_dfof_resp;
        all_traces = gcamp_data(i).dfof_traces;
        all_max_dfofs = gcamp_data(i).max_dfof_resp;
        
        % Get data from blank trials
        blanks_idx = find(stimulus_info(:,6) == blank);
        blank_mean_mat(1,:) = all_means(blanks_idx);
        tuning(i).(blank_field) = blank_mean_mat;
        
        % Get the max dfof from all trials excluding blanks
        all_max_trial_dfof_mat = all_max_dfofs(find(stimulus_info(:,6) ~= blank));
        tuning(i).(max_dfof_field) = all_max_trial_dfof_mat;
        
        for j = 1:length(sfs)
            tmp_sf = sfs(j);
            for k = 1:length(tfs)
                tmp_tf = tfs(k);
                for m = 1:length(dirs)
                    tmp_dir = dirs(m);
                    if stationary_trials == 1
                        tmp_idxs = find(stimulus_info(:,7) == tmp_sf & stimulus_info(:,8) == tmp_tf ...
                            & stimulus_info(:,6) == tmp_dir & stimulus_info(:,5) == 0); % only stationary trials
                    else
                        tmp_idxs = find(stimulus_info(:,7) == tmp_sf & stimulus_info(:,8) == tmp_tf ...
                            & stimulus_info(:,6) == tmp_dir); % running and stationary trials
                    end
                    tmp_means = all_means(tmp_idxs);
                    tmp_traces = all_traces(tmp_idxs,:);
                    mean_tmp_means = nanmean(tmp_means);
                    sem_tmp_means = sem(tmp_means);
                    mean_tmp_trace = nanmean(tmp_traces);
                    
                    % Put all the dfof traces from stationary trials in the
                    % matrix
                    all_dfof_trace_cond_mat(j,k,m,(1:size(tmp_traces,1)),:) = tmp_traces;
                    
                    % Put the mean of all mean dfofs for each stimulus in matrix
                    mean_dfof_mean_mat(j,k,m) = mean_tmp_means; % j is sf increment, k is tf increment, m is direction increment
                    mean_dfof_trace_mat(j,k,m,:) = mean_tmp_trace;
                    dfof_mean_mat(j,k,m,1:size(tmp_means,2)) = tmp_means;
                    
                    % Put the stimulus conditions in the stimulus condition
                    % matrix
                    if i == 1 % only need to do this once since all cells were presented the same stimulus set
                        stim_cond_mat(j,k,m,1) = tmp_sf;
                        stim_cond_mat(j,k,m,2) = tmp_tf;
                        stim_cond_mat(j,k,m,3) = tmp_dir;
                    end
                end
            end
        end
        
        tuning(i).(sf_tf_tuning) = mean_dfof_mean_mat;
        tuning(i).(sf_tf_dir_mean_traces) = mean_dfof_trace_mat;
        tuning(i).(sf_tf_dir_all_traces) = all_dfof_trace_cond_mat;
        tuning(i).(sf_tf_all_means) = dfof_mean_mat;
        
        % If we want to plot the SF x TF heatmap at one direction:
%         sf_tf_mat_one = squeeze(mean_dfof_mean_mat(:,:,1)); % first direction analyzed
%         figure(2); heatmap(sf_tf_mat_one);
%         sf_tf_mat_two = squeeze(mean_dfof_mean_mat(:,:,2)); % second direction analyzed
%         figure(3); heatmap(sf_tf_mat_two);
        
        % If we want to plot the dfof traces from a range of SFs at one TF
        % and both directions
%         figure(2);
%         subplot(1,2,1); plot(squeeze(mean_dfof_trace_mat(:,4,2,:))');
%         % Plot the same but at the other direction
%         subplot(1,2,2); plot(squeeze(mean_dfof_trace_mat(:,4,1,:))');
        
        % If we want to plot the dfof traces from a range of SFs and TFs
        % at one direction
%         figure(5);
%         for n = 1:4
%             plot_num = 0;
%             for p = 1:4
%                 subplot(4,4,n+plot_num); plot(squeeze(mean_dfof_trace_mat(n,p,1,:))');
%                 plot_num = plot_num + 4;
%                 % stimulus conditions for each plot
%                 stim_conds = squeeze(stim_cond_mat(n,p,1,:));
%                 title(num2str(stim_conds'));
%             end
%         end

    end

    if plot_heatmap_all_directions == 1
        % Get responsive and reliable index values for the current file
        resp_rel_idx = find([tuning(1,:).(resp_rel_field)] == 1);

        sf_tf_mat = [];

        for i = 1:length(resp_rel_idx)
            sf_tf_mat = cat(3,sf_tf_mat,tuning(resp_rel_idx(i)).(sf_tf_tuning));
        end

        % Get mean sf x tf matrix for one direction
        sf_tf_mat_one_dir = sf_tf_mat(:,:,1:length(dirs):end);
        sf_tf_mat_one_dir_mean = nanmean(sf_tf_mat_one_dir,3);

        % Plot heatmap of mean responses from all neurons from the current file
        % for one direction
        figure(6); subplot(2,2,f); heatmap(sf_tf_mat_one_dir_mean);

        % Get mean sf x tf matrix for the other direction
        sf_tf_mat_two_dir = sf_tf_mat(:,:,2:length(dirs):end);
        sf_tf_mat_two_dir_mean = nanmean(sf_tf_mat_two_dir,3);

        % Plot heatmap of mean responses from all neurons from the current file
        % for the other direction
        figure(6); subplot(2,2,f+2); heatmap(sf_tf_mat_two_dir_mean);
    end
    % Need to save out stimulus conditions matrix for each file so they can
    % be compared to determine which dimension to concatenate the data. The
    % data will be concatenated along the 'incomplete' dimension. For
    % example, if each data set only covered 2 out of 4 directions, we need
    % to concatenate along that dimension. If each data set only covered 2
    % out of 4 SFs, then we need to concatenate along that dimension
%     file_stim_cond_mat = [file_stim_cond_mat,stim_cond_mat];
    file_stim_cond_struct(f).sfs = squeeze(stim_cond_mat(:,1,1,1))';
    file_stim_cond_struct(f).tfs = squeeze(stim_cond_mat(1,:,1,2));
    file_stim_cond_struct(f).dirs = squeeze(stim_cond_mat(1,1,:,3))';

    file_sfs_mat = squeeze(stim_cond_mat(:,1,1,1));
    file_tfs_mat = squeeze(stim_cond_mat(1,:,1,2));
    file_dirs_mat = squeeze(stim_cond_mat(1,1,:,3));
%     all_stim_cond_mat = cat(cat_dim,all_stim_cond_mat,stim_cond_mat);
    all_stim_cond_cell{f} = stim_cond_mat;

end

%%
% Check the stimulus conditions mat and set the concatenation dimension
if length(files_to_analyze) > 1 % only matters which dimension if there's more than one file being analyzed
    dim_check_mat = ones(1,3);
    if length(unique([file_stim_cond_struct(1:end).dirs])) ~= length([file_stim_cond_struct(1:end).dirs])...
            || length(unique([file_stim_cond_struct(1:end).dirs])) == 1
        dim_check_mat(1,3) = 0;
    end
    if length(unique([file_stim_cond_struct(1:end).tfs])) ~= length([file_stim_cond_struct(1:end).tfs])...
            || length(unique([file_stim_cond_struct(1:end).tfs])) == 1
        dim_check_mat(1,2) = 0;
    end
    if length(unique([file_stim_cond_struct(1:end).sfs])) ~= length([file_stim_cond_struct(1:end).sfs])...
            || length(unique([file_stim_cond_struct(1:end).sfs])) == 1
        dim_check_mat(1,1) = 0;
    end
    cat_dim = find(dim_check_mat == 1);
else
    cat_dim = 1;
end
all_stim_cond_mat = [];

for f = 1:length(all_stim_cond_cell)
    if f == 1
        all_stim_cond_mat = all_stim_cond_cell{f};
    else
        all_stim_cond_mat = cat(cat_dim,all_stim_cond_mat,all_stim_cond_cell{f});
    end
end

clear gcamp_data % clear to save memory

%%
% Need to go through and get mean of all dfof_means for each cell at all
% directions for each SF/TF combination to check for DS or OS. This is so a
% true average response to each SF/TF can be calculated and one SF x TF
% heatmap can be generated for each neuron (i.e. not one for each
% direction). Only do this for responsive and reliable cells?


% Matrix to hold SF/TF matrices for all cells to quickly calculate mean and
% plot heatmap when finished
all_sf_tf_matrix = [];

% Get the responsive and reliable cells from both datasets
all_resp_rel_idx = [];

% Create structure array to hold responsive and reliable metrics
% all_resp_rel(length(files_to_analyze)) = struct();
% all_resp_rel(1).resp_rel_struct = [];

all_resp_rel_idx = [];
good_rois = 1;

for f = 1:length(files_to_analyze)
    file = files_to_analyze{f};
    resp_rel_struct = getResponsiveReliableMetrics_2024(file,tuning);
%     all_resp_rel(f).resp_rel_struct = resp_rel_struct;
    resp_idx = find([resp_rel_struct(1,:).responsive] < responsive_thresh); % responsive
    rel_idx = find([resp_rel_struct(1,:).reliable] == 1); % reliable
    good_roi_idx = intersect(resp_idx,rel_idx); % responsive and reliable
    high_idx = find([resp_rel_struct(1,:).high_fluorescence] == 0); % did not have high fluorscence
    good_roi_idx = intersect(good_roi_idx,high_idx); % responsive, reliable, and not too high fluorscence
    
    if isempty(good_roi_idx)
        error('Analysis did not finish. No cells have been deemed responsive and reliable. Consider adjusting the responsive_thresh variable.')
    end
    
    % Now take the subset of good rois and use them to recalculate the low
    % fluorescence rois
    tmp_resp_rel = getResponsiveReliableMetrics_2024(file,tuning,good_roi_idx);
    low_idx = find([tmp_resp_rel(1,:).low_fluorescence] == 0);
    % Remove the low fluorescence rois
    resp_rel_roi_idx = good_roi_idx(low_idx);
    
    all_resp_rel_idx = [all_resp_rel_idx,resp_rel_roi_idx];
end

% Now get just the unique responsive and reliable cells since they can show
% up in both datasets
all_resp_rel_idx = unique(all_resp_rel_idx);

all_tuning(length(all_resp_rel_idx)) = struct();
all_tuning(1).cell_ids = [];
all_tuning(1).sf_tf_dirs_mean_resp = [];
all_tuning(1).sf_tf_dirs_all_mean_resp = [];
all_tuning(1).sf_tf_dirs_all_traces = [];

all_tuning(1).sf_tf_mean_resp = [];
all_tuning(1).sf_tf_cv_mean_resp = [];
all_tuning(1).sf_tf_dsi_osi = [];
all_tuning(1).sf_tf_cv_dsi_osi = [];

%%
% Matrix to hold SF/TF matrices for all cells to quickly calculate mean and
% plot heatmap when finished

all_sf_tf_matrix = [];

all_dirs_tested = squeeze(all_stim_cond_mat(1,1,:,3));
all_rads_tested = circ_ang2rad(all_dirs_tested);

if length(all_dirs_tested) <= 2
    disp('Only two directions, so no OSI is calculated');
end

if length(all_dirs_tested) <= 4
    disp('Only four directions, so no Von Mises fitting is performed');
end

% The code below is calculating the mean responses for each SF/TF
% combination.

for c = 1:length(tuning)
%     disp(['Cell: ',num2str(c),'/',num2str(length(tuning))]);
% for c = 1 % for testing
    
%     tmp_idx = all_resp_rel_idx(c);
    tmp_idx = c;
    
    all_tuning(c).cell_ids = tuning(tmp_idx).cell_ids;
    
    % Get the sf x tf matrix for all directions for the cell being analyzed
    sf_tf_all_dirs_mat = []; % matrix to hold the all direction data
    sf_tf_all_dirs_mean_resp_mat = []; % matrix to hold mean responses from all trials at all directions
    sf_tf_all_dirs_traces_mat = []; % matrix to hold traces from all stimulus conditions
    for f = 1:length(files_to_analyze)
       file = files_to_analyze{f};
       sf_tf_tuning = [file(1:4), file(9:12), '_sf_tf_means'];
       sf_tf_all_means = [file(1:4), file(9:12), '_sf_tf_all_mean_resp'];
       sf_tf_dir_all_traces = [file(1:4), file(9:12), '_sf_tf_dir_all_traces'];
       
       sf_tf_tmp = tuning(tmp_idx).(sf_tf_tuning);
       sf_tf_dir_all_means = tuning(tmp_idx).(sf_tf_all_means);
       sf_tf_dir_traces = tuning(tmp_idx).(sf_tf_dir_all_traces);
       
       sf_tf_all_dirs_mat = cat(cat_dim,sf_tf_all_dirs_mat,sf_tf_tmp);
       sf_tf_all_dirs_mean_resp_mat = cat(cat_dim,sf_tf_all_dirs_mean_resp_mat,sf_tf_dir_all_means);
       sf_tf_all_dirs_traces_mat = cat(cat_dim,sf_tf_all_dirs_traces_mat,sf_tf_dir_traces);
    end
    
    % Save the all direction matrix for use later
    all_tuning(c).sf_tf_dirs_mean_resp = sf_tf_all_dirs_mat;
    all_tuning(c).sf_tf_dirs_all_mean_resp = sf_tf_all_dirs_mean_resp_mat;
    all_tuning(c).sf_tf_dirs_all_traces = sf_tf_all_dirs_traces_mat;
end
%%

% Group all responses at all sf and tf and all directions for all cells to
% perform permutation analysis for calculation of pDS and pOS.

all_sf_tf_dirs_all_mean_resp = [];

for r = 1:length(all_resp_rel_idx)
    tmp_idx = all_resp_rel_idx(r);
    tmp_sf_tf_dirs_all_mean_resp = all_tuning(tmp_idx).sf_tf_dirs_all_mean_resp;
    all_sf_tf_dirs_all_mean_resp = cat(4,all_sf_tf_dirs_all_mean_resp,tmp_sf_tf_dirs_all_mean_resp);
end

cv_dsi_mat = nan(1,1);
cv_osi_mat = nan(1,1);

for i = 1:1000
    % Select random sf x tf combination each time
    rand_sf = randi([1,size(all_sf_tf_dirs_all_mean_resp,1)]);
    rand_tf = randi([1,size(all_sf_tf_dirs_all_mean_resp,2)]);
    tmp_permute = squeeze(all_sf_tf_dirs_all_mean_resp(rand_sf,rand_tf,:,:));
    rows = randperm(size(tmp_permute,1));
    cols = randperm(size(tmp_permute,2));
    rand_tmp = tmp_permute(:,cols);
    rand_tmp = rand_tmp(rows,:);
    mean_rand_tmp = nanmean(rand_tmp');
    if min(mean_rand_tmp) < 0
        mean_rand_tmp = mean_rand_tmp-(min(mean_rand_tmp));
    end
    [cv_osi_perm, cv_dsi_perm, pref_dir, L_ori, L_dir] = get_CV_DirCV(mean_rand_tmp', all_rads_tested);
    cv_dsi_mat(1,i) = cv_dsi_perm;
    cv_osi_mat(1,i) = cv_osi_perm;
end

dsi_pd = fitdist(cv_dsi_mat','Normal');
osi_pd = fitdist(cv_osi_mat','Normal');


%%
for c = 1:length(all_tuning)
    disp(['Cell: ',num2str(c),'/',num2str(length(all_tuning))]);
% for c = 1 % for testing
    
%     tmp_idx = all_resp_rel_idx(c);
    tmp_idx = c;
    
    sf_tf_all_dirs_mat = all_tuning(tmp_idx).sf_tf_dirs_mean_resp;
    
    % Create matrix to hold the mean response at all SF/TF combinations.
    % This is the sf_tf_all_dirs_mat collapsed along the direction axis
    % because we are either averaging all directions (untuned) or using the
    % preferred response (tuned)    
    sf_tf_mean_mat = nan(size(squeeze(all_tuning(c).sf_tf_dirs_mean_resp(:,:,1))));
    sf_tf_cv_mean_mat = nan(size(sf_tf_mean_mat));
    
    % Create a matrix to hold the dsi, osi, and preferred direction for each
    % sf x tf combination. Third dimension is: 1:dsi, 2:osi, 3:direction.
    % The second matrix is for storing dsi snd osi calculated using CV
    sf_tf_dsi_osi_mat = nan(size(sf_tf_mean_mat,1),size(sf_tf_mean_mat,2),3);
    sf_tf_cv_dsi_osi_mat = nan(size(sf_tf_mean_mat,1),size(sf_tf_mean_mat,2),3);
    sf_tf_fwhm_pref_dir_mat = nan(size(sf_tf_mean_mat,1),size(sf_tf_mean_mat,2),2);

    for j = 1:size(sf_tf_all_dirs_mat,1) % sfs
%     for j = 1 % for testing
        for k = 1:size(sf_tf_all_dirs_mat,2) % tfs
%         for k = 1 % for testing
            tmp_sf_tf_mat = nan(1,size(sf_tf_all_dirs_mat,3));
            for m = 1:size(sf_tf_all_dirs_mat,3) % directions
                tmp_sf_tf_mat(1,m) = sf_tf_all_dirs_mat(j,k,m);
            end
            % Determine how many directions were presented between a preferred
            % direction and the opposite direction based on the total number of
            % directions presented.
            null_shift = length(tmp_sf_tf_mat)/2;
            orth_shift = null_shift/2;

            % Get the preferred response
            pref_resp = max(tmp_sf_tf_mat);
            pref_idx = find(tmp_sf_tf_mat == max(tmp_sf_tf_mat));

            % Get the null response
            null_shift = circshift(tmp_sf_tf_mat,null_shift);
            null_resp = null_shift(pref_idx);

            if length(all_dirs_tested) > 2
                % Get the orthogonal response
                orth_shift_one = circshift(tmp_sf_tf_mat,orth_shift);
                orth_resp_one = orth_shift_one(pref_idx);
                orth_shift_two = circshift(tmp_sf_tf_mat,-orth_shift);
                orth_resp_two = orth_shift_two(pref_idx);

                orth_resp = nanmean([orth_resp_one,orth_resp_two]);
                cell_osi = (pref_resp - orth_resp)/(pref_resp + orth_resp);
                sf_tf_dsi_osi_mat(j,k,2) = cell_osi;
            else
                cell_osi = nan;
            end

            % Calculate DSI and OSI
            cell_dsi = (pref_resp - null_resp)/(pref_resp + null_resp);
            
            sf_tf_dsi_osi_mat(j,k,1) = cell_dsi;
            
            if cell_dsi > 0.4 || cell_osi > 0.2
                sf_tf_mean_mat(j,k) = pref_resp;
                sf_tf_dsi_osi_mat(j,k,3) = all_dirs_tested(pref_idx);           
            else
                sf_tf_mean_mat(j,k) = nanmean(tmp_sf_tf_mat);
                sf_tf_dsi_osi_mat(j,k,3) = 500; % not DS or OS, so set to 500           
            end
            % Eliminate negative dfof values
            if min(tmp_sf_tf_mat) < 0
                tmp_sf_tf_mat = tmp_sf_tf_mat-(min(tmp_sf_tf_mat));
            end
            [cv_osi, cv_dsi, pref_dir, L_ori, L_dir] = get_CV_DirCV(tmp_sf_tf_mat, all_rads_tested');
            p_dsi = 1-normcdf(cv_dsi,dsi_pd.mu,dsi_pd.sigma);
            p_osi = 1-normcdf(cv_osi,osi_pd.mu,osi_pd.sigma);
            if p_dsi < 0.05 || p_osi < 0.05
                sf_tf_cv_mean_mat(j,k) = pref_resp;
            else
                sf_tf_cv_mean_mat(j,k) = nanmean(tmp_sf_tf_mat);
            end
            sf_tf_cv_dsi_osi_mat(j,k,1) = cv_dsi;
            sf_tf_cv_dsi_osi_mat(j,k,2) = cv_osi;
            sf_tf_cv_dsi_osi_mat(j,k,3) = pref_dir;
            
            % Do Von Mises fitting to get fwhm and preferred direction, but
            % only if cell is responsive and reliable to save computational
            % time and also because there's no point in fitting
            % non-responsive cells. Also, do not do von mises fitting if
            % there are fewer than 8 directions because the fwhm from the
            % fit would be meaningless anyway.
            if ismember(c,all_resp_rel_idx) == 1 && length(all_rads_tested) >= 8
                [x_data_interp, y_data_interp, x_fit, y_fit,...
                    pref_dir, fwhm, osi_fit, dsi_fit, resnorm] = fit_sum_von_mises(all_rads_tested',tmp_sf_tf_mat,360,0);
                sf_tf_fwhm_pref_dir_mat(j,k,1) = fwhm;
                sf_tf_fwhm_pref_dir_mat(j,k,2) = pref_dir;
            else
                sf_tf_fwhm_pref_dir_mat(j,k,1) = nan;
                sf_tf_fwhm_pref_dir_mat(j,k,2) = 500;
            end
        end
    end
    
    all_tuning(c).sf_tf_mean_resp = sf_tf_mean_mat;
    all_tuning(c).sf_tf_cv_mean_resp = sf_tf_cv_mean_mat;
    all_tuning(c).sf_tf_dsi_osi = sf_tf_dsi_osi_mat;
    all_tuning(c).sf_tf_cv_dsi_osi = sf_tf_cv_dsi_osi_mat;
    all_tuning(c).sf_tf_fwhm_pref_dir = sf_tf_fwhm_pref_dir_mat;

    all_sf_tf_matrix = cat(cat_dim,all_sf_tf_matrix,sf_tf_mean_mat);
end

%%
all_sfs = all_stim_cond_mat(:,1,1,1);
all_tfs = all_stim_cond_mat(1,:,1,2);
if length(all_sfs) > 1 && length(all_tfs) > 1
    figure(6); clf;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_mean_resp));
    mean_resp = nanmean(all_sf_tf,3);
    heatmap(rot90(mean_resp,2), fliplr(60./all_tfs), flipud(all_sfs), 1, 'Colormap', [], 'Colorbar', true);
    title('Traditional DS or OS Calculation')
    figure(7); clf;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_cv_mean_resp));
    mean_resp = nanmean(all_sf_tf,3);
    heatmap(rot90(mean_resp,2), fliplr(60./all_tfs), flipud(all_sfs), 1, 'Colormap', [], 'Colorbar', true);
    title('CV DS or OS Calculation');
elseif length(all_sfs) == 1 && length(all_tfs) == 1
    % DSI and OSI plot
    all_dsi_osi_dir = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_dsi_osi));
    mean_dsi = nanmean(all_dsi_osi_dir(1:3:end));
    sem_dsi = sem(all_dsi_osi_dir(1:3:end));
    mean_osi = nanmean(all_dsi_osi_dir(2:3:end));
    sem_osi = sem(all_dsi_osi_dir(2:3:end));
    figure(1); clf; hold on
    dsi_osi = bar([mean_dsi; mean_osi]);
    dsi_osi.BarWidth = 0.5;
    errorbar(1,mean_dsi,sem_dsi,'k-');
    errorbar(2,mean_osi,sem_osi,'k-');
    axis tight
    xl(0.5,2.5);
    xticks([1,2]);
    xticklabels(['DSI'; 'OSI']);
    
    % Preferred Direction Plot
    pref_dirs = all_dsi_osi_dir(3:3:end);
    tuned_dirs = pref_dirs;
    % Check for untuned cells
    untuned_idx = find(pref_dirs == 500);
    if ~isempty(untuned_idx)
        tuned_dirs = pref_dirs;
        tuned_dirs(untuned_idx) = []; % remove the untuned directions from matrix
    end
    figure(2); clf; hold on;
    h = histogram(tuned_dirs,length(dirs),'BinEdges',(dirs(1)-22.5:45:dirs(end)+22.5)); % (dirs(1)-22.5:45:dirs(end)+22.5) (-45:90:315)
    xticks(dirs);
    axis tight
    xl(dirs(1)-22.5,dirs(end)+22.5);
    yl(0,max(h.Values)+size(tuned_dirs,1)*0.05);
    xlabel('Preferred Direction');
    ylabel('nCells');
elseif length(all_sfs) == 1 && length(all_tfs) > 1 % if one SF and multiple TFs make a heatmap
    figure(6); clf; hold on;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_mean_resp));
    heatmap((nanmean(all_sf_tf,2)), all_sfs, (60./all_tfs), 1, 'Colormap', [], 'Colorbar', true);
    title('Traditional DS or OS Calculation')
    figure(7); clf; hold on;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_cv_mean_resp));
    heatmap((nanmean(all_sf_tf,2)), all_sfs, (60./all_tfs), 1, 'Colormap', [], 'Colorbar', true);
    title('CV DS or OS Calculation')
elseif length(all_sfs) > 1 && length(all_tfs) == 1  % if one TF and multiple SFs make a heatmap
    figure(6); clf; hold on;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_mean_resp));
    [sort_sfs,sort_idx] = sort(all_sfs);
    mean_sf_tf = nanmean(all_sf_tf,2);
    heatmap((mean_sf_tf(sort_idx)), (60./all_tfs), sort_sfs, 1, 'Colormap', [], 'Colorbar', true);
    title('Traditional DS or OS Calculation')
    figure(7); clf; hold on;
    all_sf_tf = squeeze(cat(3,all_tuning(all_resp_rel_idx).sf_tf_cv_mean_resp));
    [sort_sfs,sort_idx] = sort(all_sfs);
    mean_sf_tf = nanmean(all_sf_tf,2);
    heatmap((mean_sf_tf(sort_idx)), (60./all_tfs), sort_sfs, 1, 'Colormap', [], 'Colorbar', true);
    title('CV DS or OS Calculation')

end

% Plot all traces of the max response at the preferred direction (or
% average of all directions if not DS or OS)
if plot_cell_traces == 1
    figure(3); clf;
    % Determine number of rows and columns for subplot based on number of
    % responsive and reliable cells
    n_rois_ten = ceil(size(all_resp_rel_idx,2)/10)*10;
    for c = 1:length(all_resp_rel_idx)
%     for c = 1:20 % for testing
        tmp_sf_tf_mat = all_tuning(all_resp_rel_idx(c)).sf_tf_mean_resp;
        tmp_traces_all = all_tuning(all_resp_rel_idx(c)).sf_tf_dirs_all_traces;
        tmp_dsi_osi = all_tuning(all_resp_rel_idx(c)).sf_tf_dsi_osi;
        
        % Find the max response from the SF x TF matrux
        tmp_max_sf_tf_idx = find(tmp_sf_tf_mat == max(tmp_sf_tf_mat(:)));
        % Get n-dimensional matrix position of the max value
        tmp_sf_tf_dim = [size(tmp_sf_tf_mat)]; % dimensions of sf_tf matrix
        [A,B] = ind2sub(tmp_sf_tf_dim,tmp_max_sf_tf_idx);
        
        % Get the traces, DSI, OSI, preferred direction, and stimulus
        % conditions of the max response in the matrix
        tmp_max_traces = squeeze(tmp_traces_all(A,B,:,:,:));
        tmp_max_dsi_osi = squeeze(tmp_dsi_osi(A,B,:));
        tmp_stim_cond = squeeze(all_stim_cond_mat(A,B,:,:));
        
        % Check to see if the cell was DS or OS
        if tmp_max_dsi_osi(1) > 0.4 && tmp_max_dsi_osi(2) > 0.2
            dir_idx = find(tmp_stim_cond(:,3) == tmp_max_dsi_osi(3));
        else
            disp('Untuned cell')
            dir_idx = (1:1:4);
        end
        % Get the specific traces (or averaged traces if untuned) for
        % plotting
        tmp_max_dir_traces = squeeze(nanmean(tmp_max_traces(dir_idx,:,:),1));
        subplot(n_rois_ten/10,10,c); hold on; plot(nanmean(tmp_max_dir_traces)); title(['Cell: ', num2str(all_resp_rel_idx(c))]);
    end
end
% Save out data
if save_data == 1
    disp('Saving Data');
    % First get shortened versions of the file analyzed to include in the name
    % of the saved .mat file
    save_file_name = experiment_day;
    for f = 1:length(file_names)
        tmp_name = file_names{f};
        file_num = tmp_name(end-3:end);
        save_file_name = [save_file_name,file_num];
    end

    if cno_data == 0 % assume this is control data
        if ~exist('Control', 'dir')
           mkdir('Control')
        end
        cd 'Control'
        save([save_file_name, '_control'], 'all_tuning', 'all_stim_cond_mat', 'all_resp_rel_idx', 'dsi_pd', 'osi_pd');
    elseif cno_data == 1
        if ~exist('CNO', 'dir')
           mkdir('CNO')
        end
        cd 'CNO'
        save([save_file_name, '_cno'], 'all_tuning', 'all_stim_cond_mat', 'all_resp_rel_idx', 'dsi_pd', 'osi_pd');
    end
end