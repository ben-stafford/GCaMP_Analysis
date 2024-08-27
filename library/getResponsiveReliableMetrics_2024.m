function resp_rel_struct = getResponsiveReliableMetrics_2024(file, tuning, cells)
    try
        if ~exist('cells','var')
            cell_idx = length(tuning);
%             disp('Using tuning variable')
            cells = (1:1:length(tuning));
        else
            cell_idx = length(cells);
%             disp('Using cells variable')
        end
        % Create structure to hold output of the function
        resp_rel_struct(cell_idx) = struct();
        resp_rel_struct(1).d_prime_mat = [];
        resp_rel_struct(1).reliable = [];
        resp_rel_struct(1).responsive = [];
        resp_rel_struct(1).high_fluorescence = [];
        resp_rel_struct(1).low_fluorescence = [];

        % Get the sf x tf matrix for all directions for the cell being analyzed

        sf_tf_field = [file(1:4), file(9:12), '_sf_tf_means'];
        sf_tf_all_mean_resp_field = [file(1:4), file(9:12), '_sf_tf_all_mean_resp'];
        sf_tf_blanks_field = [file(1:4), file(9:12), '_blank_means'];
        sf_tf_max_dfof_field = [file(1:4), file(9:12), '_max_dfof_trials'];

        % Low fluorescence check. This needs to be done after removing the high
        % fluorscence, unreliable, and unresponsive rois. Otherwise, I think
        % too many wonky rois are included and none of the rois pass this test.

        % any roi where the max dfof of all its mean responses to all stimulus
        % combinations was less than the 6% of the maximum mean dfof of all the
        % rois in the imaging session.

        tmp_low_fluor_data = [];
        for c = 1:cell_idx
            idx = cells(c);
            tmp_low_fluor_data = cat(3,tmp_low_fluor_data,[tuning(1,idx).(sf_tf_field)]);
        end
        max_val = max(tmp_low_fluor_data(:));
        low_fluor = max_val * 0.06;
        idx_low_fluor = [];
        for c = 1:cell_idx
            idx = cells(c);
            tmp_thresh_data = [tuning(1,idx).(sf_tf_all_mean_resp_field)];
            max_tmp_data = max(tmp_thresh_data(:));
            if max_tmp_data < low_fluor
                idx_low_fluor = [idx_low_fluor, idx];
                resp_rel_struct(c).low_fluorescence = 1;
%                 disp('Low fluorescence detected'); % will detect a lot if no previous dfof exclusion has been done
            else
                resp_rel_struct(c).low_fluorescence = 0;
            end
        end

        high_fluor_idx = [];

        roi_rel_mat = zeros(cell_idx,1);
        roi_resp_mat = zeros(cell_idx,2);

        for c = 1:cell_idx
    %     for c = 1:5 % for testing
            idx = cells(c);
            sf_tf_all_dirs_resp_mat = tuning(idx).(sf_tf_field);
            sf_tf_all_means_resp = tuning(idx).(sf_tf_all_mean_resp_field);
            all_blanks = tuning(idx).(sf_tf_blanks_field);

            % High dfof removal

            % roi is removed if it had a single mean trial dfof that exceeded
            % the median of the maximum dfof for all trials to all stimulus
            % combinations for the same roi + the 95th percentile of the
            % maximum mean dfof for all trials for the same roi

            % trial dF/F that exceeded the median maximum trial dF/F per ROI + the
            % 95% percentile maximum trial dF/F per ROI

            % Get all max dfof values for the cell
            all_max_dfof = tuning(idx).(sf_tf_max_dfof_field);
            max_mean_dfof = max(tuning(idx).(sf_tf_all_mean_resp_field)(:)); % all responses excluding blanks
            median_dfof = median(all_max_dfof);
            percentile_max_dfof = prctile(all_max_dfof,95);
            high_dfof_thresh = median_dfof + percentile_max_dfof;

            if max_mean_dfof > high_dfof_thresh
                high_fluor_idx = [high_fluor_idx, idx];
                resp_rel_struct(c).high_fluorescence = 1;
%                 disp('High fluorescence detected')
            else
                resp_rel_struct(c).high_fluorescence = 0;
            end

            % Responsive calculation. One way ANOVA using all repeats of each
            % stimulus condition and the blanks.

            % Run one way ANOVA using all SF, TF, and direction combinations

            n_stim_conds = size(sf_tf_all_means_resp,1) * size(sf_tf_all_means_resp,2) * size(sf_tf_all_means_resp,3);

            tmp_anova_mat = nan(size(sf_tf_all_means_resp,4),n_stim_conds + 1);

            stim_cond = 1;
            % Create matrix to run anova such that each column contains all the
            % responses to one stimulus combination
            for j = 1:size(sf_tf_all_dirs_resp_mat,1) % sfs
    %         for j = 1 % for testing
                for k = 1:size(sf_tf_all_dirs_resp_mat,2) % tfs
    %             for k = 1 % for testing
                    for m = 1:size(sf_tf_all_means_resp,3) % directions
                        tmp_sf_tf_dir_all_resp_mat = squeeze(sf_tf_all_means_resp(j,k,m,:));
                        tmp_anova_mat(:,stim_cond) = tmp_sf_tf_dir_all_resp_mat;
                        stim_cond = stim_cond + 1;
                    end
                end
            end

            % Now add blanks
            for b = 1:length(all_blanks)
                tmp_anova_mat(b,end) = all_blanks(b);
            end

            [p_anova, anovatab, stats] = anova1(tmp_anova_mat, [],'off');
            roi_resp_mat(idx,1) = p_anova;
            resp_rel_struct(c).responsive = p_anova;
            if p_anova < 0.05
                roi_resp_mat(idx,2) = 1;
%                 disp('Responsive ROI found')
            else
                roi_resp_mat(idx,2) = 0;
            end

            
            % d' calculation. If d' is less than 0.5, cell is scored as
            % unreliable.

            % Regarding the preferred stimulus, should this just be the maximum
            % mean dfof of responses to all stimulus conditions, or should DSI
            % and OSI be calculated find the preferred stimulus? Currently
            % written to base calculation off of DSI/OSI calculation.

            % Matrix to hold the reliability score for each sf x tf
            % combination
            sf_tf_rel_mat = nan(size(squeeze(sf_tf_all_dirs_resp_mat(:,:,1))));
            

            for j = 1:size(sf_tf_all_dirs_resp_mat,1) % sfs
    %         for j = 1 % for testing
                for k = 1:size(sf_tf_all_dirs_resp_mat,2) % tfs
    %             for k = 1 % for testing
                    tmp_sf_tf_mat = nan(1,size(sf_tf_all_dirs_resp_mat,3));
                    for m = 1:size(sf_tf_all_dirs_resp_mat,3) % directions
                        tmp_sf_tf_mat(1,m) = sf_tf_all_dirs_resp_mat(j,k,m);
                    end     

                    % Determine how many directions were tested
                    all_dirs_tested = size(sf_tf_all_dirs_resp_mat,3);
                    
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

                        orth_resp = mean([orth_resp_one,orth_resp_two]);
                        cell_osi = (pref_resp - orth_resp)/(pref_resp + orth_resp);
                    else
                        cell_osi = nan;
                    end

                    % Calculate DSI and OSI
                    cell_dsi = (pref_resp - null_resp)/(pref_resp + null_resp);

                    if cell_dsi > 0.4 || cell_osi > 0.2                    
                        % cell is tuned, use preferred direction for resp/rel
                        pref_resp_means = squeeze(sf_tf_all_dirs_resp_mat(j,k,pref_idx,:));

                        % calculate d'
                        d_prime = (mean(pref_resp_means) - mean(all_blanks)) / ...
                            (std(pref_resp_means) + std(all_blanks));

                        sf_tf_rel_mat(j,k) = d_prime;

                    else                    
                        % cell is not tuned, take average of all responses
                        all_resp_means = mean(squeeze(sf_tf_all_dirs_resp_mat(j,k,:,:)));

                        % calculate d'
                        d_prime = (mean(all_resp_means) - mean(all_blanks)) / ...
                            (std(all_resp_means) + std(all_blanks));

                        sf_tf_rel_mat(j,k) = d_prime;
                    end
                end
            end
                        
            % If any of the responses have a d' greater than 0.5 the roi is
            % considered reliable. Should this be changed to requiring all of
            % them? Or the mean of all d' scores for all stimulus conditions?
            resp_rel_struct(c).d_prime_mat = sf_tf_rel_mat;
            if ~isempty(find(sf_tf_rel_mat >= 1.0))
                roi_rel_mat(idx) = 1;
                resp_rel_struct(c).reliable = 1;
%                 disp('Reliable ROI found')
            else
                resp_rel_struct(c).reliable = 0;
            end
        end
        
    catch
        
    end
end
