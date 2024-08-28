clear all

% Flag to either compile the data (i.e. open several files from the same
% conditions and save the data in a new file) or compare the data (i.e.
% open files from different conditions (control vs CNO) and compare them.
compare = 1;

    animal_ids = {
        '6510'
        '6510'
        };
        
    experiment_days = {
        'day1'
        'day1'
        };

    files_to_compare = {
        'day1_000_control.mat'
        'day1_005_cno.mat'
        };

if compare % comparing datasets
    % Load the saved data files for control and cno data and store them as
    % separate variables
    for f = 1:length(files_to_compare)
        data_dir = ['/Volumes/EPHYS_000/GCaMP_Data/', animal_ids{f}, '/', experiment_days{f},...
        '/Analysis'];

        cd (data_dir)
        
        file = files_to_compare{f};
        if ~isempty(strfind(file, 'control'))
            load (['Control/', file]);
            control_tuning = all_tuning;
            control_stim_cond = all_stim_cond_mat;
            control_resp_rel_idx = all_resp_rel_idx;
        elseif ~isempty(strfind(file, 'cno'))
            load (['CNO/', file]);
            cno_tuning = all_tuning;
            cno_stim_cond = all_stim_cond_mat;
            cno_resp_rel_idx = all_resp_rel_idx;
        end
    end

    % Clear loaded variables to save memory and reduce variables
    clear all_tuning all_stim_cond_mat all_resp_rel_idx

    % Quick check to make sure control and cno datasets were from the same
    % stimulus conditions
    if isequal(control_stim_cond,cno_stim_cond)
        disp('Datasets have the same stimulus conditions')
    else
        disp('Datasets have different stimulus conditions. Analysis might not work')
    end

    % Get all the responsive and reliable cells in both control and cno
    all_resp_rel_idx = [control_resp_rel_idx cno_resp_rel_idx];
    all_resp_rel_idx = unique(all_resp_rel_idx);
%     all_resp_rel_idx = control_resp_rel_idx;

    % Create matrix to hold preferred SF, TF, direction, and DSI and OSI values
    % for each cell
    cont_cno_sf_tf_struct = struct();
    cont_cno_sf_tf_struct(1).control_sf = [];
    cont_cno_sf_tf_struct(1).control_tf = [];
    cont_cno_sf_tf_struct(1).control_dsi = [];
    cont_cno_sf_tf_struct(1).control_osi = [];
    cont_cno_sf_tf_struct(1).control_dir = [];
    cont_cno_sf_tf_struct(1).cno_sf = [];
    cont_cno_sf_tf_struct(1).cno_tf = [];
    cont_cno_sf_tf_struct(1).cno_dsi = [];
    cont_cno_sf_tf_struct(1).cno_osi = [];
    cont_cno_sf_tf_struct(1).cno_dir = [];


    for c = 1:length(all_resp_rel_idx)
    % for c = 1 % for testing
        cell_idx = all_resp_rel_idx(c);

        % Get the control data from the cell
        cont_sf_tf_mat = control_tuning(cell_idx).sf_tf_mean_resp;
        cont_sf_tf_dirs_mat = control_tuning(cell_idx).sf_tf_dirs_mean_resp;
        cont_dsi_osi_mat = control_tuning(cell_idx).sf_tf_dsi_osi;

        % Find the index of the max value from one cell
        cont_max_sf_tf_idx = find(cont_sf_tf_mat == max(cont_sf_tf_mat(:)));
        cont_max_sf_tf_dirs_idx = find(cont_sf_tf_dirs_mat == max(cont_sf_tf_dirs_mat(:)));

        % Get n-dimensional matrix position of the max value
        cont_tmp_sf_tf_dim = [size(cont_sf_tf_mat)]; % dimensions of sf_tf matrix
        cont_tmp_sf_tf_dirs_dim = [size(cont_sf_tf_dirs_mat)]; % dimensions of sf_tf_dirs matrix

        [A,B] = ind2sub(cont_tmp_sf_tf_dim,cont_max_sf_tf_idx);
        cont_tmp_sf_tf_vals_mat = squeeze(control_stim_cond(A,B,1,:));
        cont_tmp_dsi_osi_dir_mat = squeeze(cont_dsi_osi_mat(A,B,:));
        [A,B,C] = ind2sub(cont_tmp_sf_tf_dirs_dim,cont_max_sf_tf_dirs_idx);
        cont_tmp_sf_tf_dirs_vals_mat = squeeze(control_stim_cond(A,B,C,:));


        % Get the CNO data from the cell
        cno_sf_tf_mat = cno_tuning(cell_idx).sf_tf_mean_resp;
        cno_sf_tf_dirs_mat = cno_tuning(cell_idx).sf_tf_dirs_mean_resp;
        cno_dsi_osi_mat = cno_tuning(cell_idx).sf_tf_dsi_osi;

        % Find the index of the max value from one cell
        cno_max_sf_tf_idx = find(cno_sf_tf_mat == max(cno_sf_tf_mat(:)));
        cno_max_sf_tf_dirs_idx = find(cno_sf_tf_dirs_mat == max(cno_sf_tf_dirs_mat(:)));

        % Get n-dimensional matrix position of the max value
        cno_tmp_sf_tf_dim = [size(cno_sf_tf_mat)]; % dimensions of sf_tf matrix
        cno_tmp_sf_tf_dirs_dim = [size(cno_sf_tf_dirs_mat)]; % dimensions of sf_tf_dirs matrix

        [A,B] = ind2sub(cno_tmp_sf_tf_dim,cno_max_sf_tf_idx);
        cno_tmp_sf_tf_vals_mat = squeeze(cno_stim_cond(A,B,1,:));
        cno_tmp_dsi_osi_dir_mat = squeeze(cno_dsi_osi_mat(A,B,:));
        [A,B,C] = ind2sub(cno_tmp_sf_tf_dirs_dim,cno_max_sf_tf_dirs_idx);
        cno_tmp_sf_tf_dirs_vals_mat = squeeze(cno_stim_cond(A,B,C,:));

        cont_cno_sf_tf_struct(c).control_sf = cont_tmp_sf_tf_vals_mat(1);
        cont_cno_sf_tf_struct(c).control_tf = cont_tmp_sf_tf_vals_mat(2);
        cont_cno_sf_tf_struct(c).control_dsi = cont_tmp_dsi_osi_dir_mat(1);
        cont_cno_sf_tf_struct(c).control_osi = cont_tmp_dsi_osi_dir_mat(2);
        cont_cno_sf_tf_struct(c).control_dir = cont_tmp_dsi_osi_dir_mat(3);
        cont_cno_sf_tf_struct(c).cno_sf = cno_tmp_sf_tf_vals_mat(1);
        cont_cno_sf_tf_struct(c).cno_tf = cno_tmp_sf_tf_vals_mat(2);
        cont_cno_sf_tf_struct(c).cno_dsi = cno_tmp_dsi_osi_dir_mat(1);
        cont_cno_sf_tf_struct(c).cno_osi = cno_tmp_dsi_osi_dir_mat(2);
        cont_cno_sf_tf_struct(c).cno_dir = cno_tmp_dsi_osi_dir_mat(3);

    end
    
    if size(cont_sf_tf_mat,1) > 1 && size(cont_sf_tf_mat,2) > 1
        % Make some plots of control vs cno data
        figure(1); clf;
        % sf data
        all_control_sf = [cont_cno_sf_tf_struct(1,:).control_sf];
        all_cno_sf = [cont_cno_sf_tf_struct(1,:).cno_sf];
        mean_control_sf = mean(all_control_sf);
        sem_control_sf = sem(all_control_sf);
        mean_cno_sf = mean(all_cno_sf);
        sem_cno_sf = sem(all_cno_sf);

        % tf data
        all_control_tf = [cont_cno_sf_tf_struct(1,:).control_tf];
        all_cno_tf = [cont_cno_sf_tf_struct(1,:).cno_tf];
        mean_control_tf = nanmean(60./all_control_tf);
        sem_control_tf = nansem(60./all_control_tf);
        mean_cno_tf = nanmean(60./all_cno_tf);
        sem_cno_tf = nansem(60./all_cno_tf);

        % tf/sf data
        tf_sf_control = (60./all_control_tf)./all_control_sf;
        tf_sf_cno = (60./all_cno_tf)./all_cno_sf;
        mean_tf_sf_control = nanmean(tf_sf_control);
        sem_tf_sf_control = nansem(tf_sf_control);
        mean_tf_sf_cno = nanmean(tf_sf_cno);
        sem_tf_sf_cno = nansem(tf_sf_cno);

        % all tf/sf ratios in case they are needed
        all_sfs = squeeze(control_stim_cond(:,:,1,1));
        all_tfs = 60./squeeze(control_stim_cond(:,:,1,2));
        tf_sf_ratios = unique(all_tfs./all_sfs);

        subplot(1,3,1); hold on;    
        plot(zeros(length(all_control_sf)), all_control_sf(:),'ko','Markersize',5);
        plot(ones(length(all_cno_sf))*0.5, all_cno_sf(:),'ro','Markersize',5);
        for i=1:length(all_control_sf)
            plot([0 0.5],[all_control_sf(i) all_cno_sf(i)], 'k-');
        end
        plot(0,mean_control_sf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_sf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_sf,sem_control_sf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_sf,sem_cno_sf,'r-',  'Linewidth', 1.5);
        yl(0,0.1);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred SF');

        subplot(1,3,2); hold on;    
        plot(zeros(length(all_control_tf)), (60./all_control_tf),'ko','Markersize',5);
        plot(ones(length(all_cno_tf))*0.5, (60./all_cno_tf),'ro','Markersize',5);
        for i=1:length(all_control_tf)
            plot([0 0.5],[(60/all_control_tf(i)) (60/all_cno_tf(i))], 'k-');
        end
        plot(0,mean_control_tf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_tf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_tf,sem_control_tf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_tf,sem_cno_tf,'r-',  'Linewidth', 1.5);
        yl(0,10);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF');

        subplot(1,3,3); hold on;    
        plot(zeros(length(tf_sf_control)), tf_sf_control(:),'ko','Markersize',5);
        plot(ones(length(tf_sf_cno))*0.5, tf_sf_cno(:),'ro','Markersize',5);
        for i=1:length(tf_sf_control)
            plot([0 0.5],[tf_sf_control(i) tf_sf_cno(i)], 'k-');
        end
        plot(0,mean_tf_sf_control,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_tf_sf_cno,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_tf_sf_control,sem_tf_sf_control,'k-','Linewidth',1.5);
        errorbar(0.5,mean_tf_sf_cno,sem_tf_sf_cno,'r-',  'Linewidth', 1.5);
        yl(0,900);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF/SF (º/s)');

        % Make some plots of control vs cno data
        figure(2); clf;
        all_control_dsi = [cont_cno_sf_tf_struct(1,:).control_dsi];
        all_cno_dsi = [cont_cno_sf_tf_struct(1,:).cno_dsi];
        all_control_osi = [cont_cno_sf_tf_struct(1,:).control_osi];
        all_cno_osi = [cont_cno_sf_tf_struct(1,:).cno_osi];
        all_control_dir = [cont_cno_sf_tf_struct(1,:).control_dir];
        all_cno_dir = [cont_cno_sf_tf_struct(1,:).cno_dir];
        
        subplot(1,2,1); hold on;    
        plot(zeros(length(all_control_dsi)), all_control_dsi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_dsi))*0.5, all_cno_dsi(:),'ro','Markersize',5);
        for i=1:length(all_control_dsi)
            plot([0 0.5],[all_control_dsi(i) all_cno_dsi(i)], 'k-');
        end
        plot(0,nanmean(all_control_dsi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_dsi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_dsi),nansem(all_control_dsi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_dsi),nansem(all_cno_dsi),'r-',  'Linewidth', 1.5);
        yl(-1.0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('DSI');
        
        subplot(1,2,2); hold on;    
        plot(zeros(length(all_control_osi)), all_control_osi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_osi))*0.5, all_cno_osi(:),'ro','Markersize',5);
        for i=1:length(all_control_osi)
            plot([0 0.5],[all_control_osi(i) all_cno_osi(i)], 'k-');
        end
        plot(0,nanmean(all_control_osi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_osi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_osi),nansem(all_control_osi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_osi),nansem(all_cno_osi),'r-',  'Linewidth', 1.5);
        yl(0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('OSI');
        

        % Make some histograms of the SF and TF
        % SF
        sorted_sfs = unique(sort(all_sfs));
        figure(3); clf;
        subplot(1,3,1); hold on;
        histogram(log2(all_control_sf),length(sorted_sfs),...
            'BinEdges',[log2(sorted_sfs(1))-0.5:1:log2(sorted_sfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(all_control_sf)), 0.5, 'v', 'MarkerSize', 5, 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(all_cno_sf),length(sorted_sfs),...
            'BinEdges',[log2(sorted_sfs(1))-0.5:1:log2(sorted_sfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(all_cno_sf)), 0.5, 'v', 'MarkerSize', 5, 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred SF');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(sorted_sfs(1))-0.75 log2(sorted_sfs(end))+0.75]);
        xticks(log2([sorted_sfs]));
        xticklabels(([sorted_sfs]));


        % TF
        sorted_tfs = unique(sort(all_tfs));
        subplot(1,3,2); hold on;
        histogram(log2(60./all_control_tf),4,...
            'BinEdges',[log2(sorted_tfs(1))-0.5:1:log2(sorted_tfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(60./all_control_tf)), 0.5, 'v', 'MarkerSize', 5,  'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(60./all_cno_tf),4,...
            'BinEdges',[log2(sorted_tfs(1))-0.5:1:log2(sorted_tfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(60./all_cno_tf)), 0.5, 'v', 'MarkerSize', 5,  'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred TF');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(sorted_tfs(1))-0.75 log2(sorted_tfs(end))+0.75]);
        xticks(log2([sorted_tfs]));
        xticklabels(([sorted_tfs]));

        % TF/SF (approximation of cell's speed tuning)
        tf_over_sf = unique(sort(all_tfs./all_sfs));
        subplot(1,3,3); hold on;
        histogram(log2(tf_sf_control),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(tf_sf_control)), 0.5, 'v', 'MarkerSize', 5,  'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(tf_sf_cno),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(tf_sf_cno)), 0.5, 'v', 'MarkerSize', 5,  'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred TF/SF (º/s)');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(tf_over_sf(1))-0.75 log2(tf_over_sf(end))+0.75]);
        xticks(log2(tf_over_sf(1)):1:log2(tf_over_sf(end)));
        xticklabels(2.^(log2(tf_over_sf(1)):1:log2(tf_over_sf(end))));

        % Make a 2D histogram of cells in control vs CNO
        figure(4); clf; hold on;
        edges = [0 0.015 0.025 0.045 0.085]; % set the bin edges (sfs)
        cell_counts = histcounts2(all_cno_sf,all_control_sf,edges,edges); % matrix of sf combinations
        x = reshape(repmat((1:4)', 4, 1),1,[]); % x values for plotting
        y = reshape(repmat(1:4, 4, 1),1,[]); % y values for plotting
        cell_counts_reshaped = reshape(flip(cell_counts)',1,[]); % convert cell count matrix to vector
        scatter(x,y,cell_counts_reshaped*300,'Facecolor','k');
        plot([0 5],[0 5],'k--');
        xticks(1:1:4); yticks(1:1:4);
        xticklabels([0.01 0.02 0.04 0.08]);
        xlabel('preferred SF (control)');
        yticklabels([0.01 0.02 0.04 0.08]);
        ylabel('preferred SF (CNO)');
    elseif size(cont_sf_tf_mat,1) == 1 && size(cont_sf_tf_mat,2) > 1
        % Make some plots of control vs cno data
        figure(1); clf;
        % sf data
        all_control_sf = [cont_cno_sf_tf_struct(1,:).control_sf];
        all_cno_sf = [cont_cno_sf_tf_struct(1,:).cno_sf];
        mean_control_sf = mean(all_control_sf);
        sem_control_sf = sem(all_control_sf);
        mean_cno_sf = mean(all_cno_sf);
        sem_cno_sf = sem(all_cno_sf);

        % tf data
        all_control_tf = [cont_cno_sf_tf_struct(1,:).control_tf];
        all_cno_tf = [cont_cno_sf_tf_struct(1,:).cno_tf];
        mean_control_tf = nanmean(60./all_control_tf);
        sem_control_tf = nansem(60./all_control_tf);
        mean_cno_tf = nanmean(60./all_cno_tf);
        sem_cno_tf = nansem(60./all_cno_tf);

        % tf/sf data
        tf_sf_control = (60./all_control_tf)./all_control_sf;
        tf_sf_cno = (60./all_cno_tf)./all_cno_sf;
        mean_tf_sf_control = nanmean(tf_sf_control);
        sem_tf_sf_control = nansem(tf_sf_control);
        mean_tf_sf_cno = nanmean(tf_sf_cno);
        sem_tf_sf_cno = nansem(tf_sf_cno);

        % all tf/sf ratios in case they are needed
        all_sfs = squeeze(control_stim_cond(:,:,1,1));
        all_tfs = 60./squeeze(control_stim_cond(:,:,1,2));
        tf_sf_ratios = unique(all_tfs./all_sfs);

        subplot(1,3,1); hold on;    
        plot(zeros(length(all_control_sf)), all_control_sf(:),'ko','Markersize',5);
        plot(ones(length(all_cno_sf))*0.5, all_cno_sf(:),'ro','Markersize',5);
        for i=1:length(all_control_sf)
            plot([0 0.5],[all_control_sf(i) all_cno_sf(i)], 'k-');
        end
        plot(0,mean_control_sf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_sf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_sf,sem_control_sf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_sf,sem_cno_sf,'r-',  'Linewidth', 1.5);
        yl(0,0.1);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred SF');

        subplot(1,3,2); hold on;    
        plot(zeros(length(all_control_tf)), (60./all_control_tf),'ko','Markersize',5);
        plot(ones(length(all_cno_tf))*0.5, (60./all_cno_tf),'ro','Markersize',5);
        for i=1:length(all_control_tf)
            plot([0 0.5],[(60/all_control_tf(i)) (60/all_cno_tf(i))], 'k-');
        end
        plot(0,mean_control_tf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_tf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_tf,sem_control_tf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_tf,sem_cno_tf,'r-',  'Linewidth', 1.5);
        yl(0,10);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF');

        subplot(1,3,3); hold on;    
        plot(zeros(length(tf_sf_control)), tf_sf_control(:),'ko','Markersize',5);
        plot(ones(length(tf_sf_cno))*0.5, tf_sf_cno(:),'ro','Markersize',5);
        for i=1:length(tf_sf_control)
            plot([0 0.5],[tf_sf_control(i) tf_sf_cno(i)], 'k-');
        end
        plot(0,mean_tf_sf_control,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_tf_sf_cno,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_tf_sf_control,sem_tf_sf_control,'k-','Linewidth',1.5);
        errorbar(0.5,mean_tf_sf_cno,sem_tf_sf_cno,'r-',  'Linewidth', 1.5);
        yl(0,900);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF/SF (º/s)');

        % Make some plots of control vs cno data
        figure(2); clf;
        all_control_dsi = [cont_cno_sf_tf_struct(1,:).control_dsi];
        all_cno_dsi = [cont_cno_sf_tf_struct(1,:).cno_dsi];
        all_control_osi = [cont_cno_sf_tf_struct(1,:).control_osi];
        all_cno_osi = [cont_cno_sf_tf_struct(1,:).cno_osi];
        all_control_dir = [cont_cno_sf_tf_struct(1,:).control_dir];
        all_cno_dir = [cont_cno_sf_tf_struct(1,:).cno_dir];
        
        subplot(1,2,1); hold on;    
        plot(zeros(length(all_control_dsi)), all_control_dsi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_dsi))*0.5, all_cno_dsi(:),'ro','Markersize',5);
        for i=1:length(all_control_dsi)
            plot([0 0.5],[all_control_dsi(i) all_cno_dsi(i)], 'k-');
        end
        plot(0,nanmean(all_control_dsi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_dsi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_dsi),nansem(all_control_dsi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_dsi),nansem(all_cno_dsi),'r-',  'Linewidth', 1.5);
        yl(-1.0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('DSI');
        
        subplot(1,2,2); hold on;    
        plot(zeros(length(all_control_osi)), all_control_osi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_osi))*0.5, all_cno_osi(:),'ro','Markersize',5);
        for i=1:length(all_control_osi)
            plot([0 0.5],[all_control_osi(i) all_cno_osi(i)], 'k-');
        end
        plot(0,nanmean(all_control_osi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_osi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_osi),nansem(all_control_osi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_osi),nansem(all_cno_osi),'r-',  'Linewidth', 1.5);
        yl(0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('OSI');
        

        % Make some histograms of the SF and TF
        % No SF histogram, because there is only one SF being tested
        % TF
        sorted_tfs = sort(all_tfs);
        subplot(1,3,2); hold on;
        histogram(log2(60./all_control_tf),length(sorted_tfs),...
            'BinEdges',[log2(sorted_tfs(1))-0.5:1:log2(sorted_tfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(60./all_control_tf)), 0.5, 'v', 'MarkerSize', 5,  'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(60./all_cno_tf),length(sorted_tfs),...
            'BinEdges',[log2(sorted_tfs(1))-0.5:1:log2(sorted_tfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(60./all_cno_tf)), 0.5, 'v', 'MarkerSize', 5,  'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred TF');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(sorted_tfs(1))-0.75 log2(sorted_tfs(end))+0.75]);
        xticks(log2([sorted_tfs]));
        xticklabels(([sorted_tfs]));
        yticks([0:0.1:0.5]);

        % TF/SF (approximation of cell's speed tuning)
        tf_over_sf = unique(sort(all_tfs./all_sfs));
        subplot(1,3,3); hold on;
        histogram(log2(tf_sf_control),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(tf_sf_control)), 0.5, 'v', 'MarkerSize', 5,  'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(tf_sf_cno),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(tf_sf_cno)), 0.5, 'v', 'MarkerSize', 5,  'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred TF/SF (º/s)');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(tf_over_sf(1))-0.75 log2(tf_over_sf(end))+0.75]);
        xticks(log2(tf_over_sf(1)):1:log2(tf_over_sf(end)));
        xticklabels(2.^(log2(tf_over_sf(1)):1:log2(tf_over_sf(end))));

%         % Make a 2D histogram of cells in control vs CNO
%         figure(4); clf; hold on;
%         edges = [0 0.015 0.025 0.045 0.085]; % set the bin edges (sfs)
%         cell_counts = histcounts2(all_cno_sf,all_control_sf,edges,edges); % matrix of sf combinations
%         x = reshape(repmat((1:4)', 4, 1),1,[]); % x values for plotting
%         y = reshape(repmat(1:4, 4, 1),1,[]); % y values for plotting
%         cell_counts_reshaped = reshape(flip(cell_counts)',1,[]); % convert cell count matrix to vector
%         scatter(x,y,cell_counts_reshaped*300,'Facecolor','k');
%         plot([0 5],[0 5],'k--');
%         xticks(1:1:4); yticks(1:1:4);
%         xticklabels([0.01 0.02 0.04 0.08]);
%         xlabel('preferred SF (control)');
%         yticklabels([0.01 0.02 0.04 0.08]);
%         ylabel('preferred SF (CNO)');
    elseif size(cont_sf_tf_mat,1) > 1 && size(cont_sf_tf_mat,2) == 1
        % Make some plots of control vs cno data
        figure(1); clf;
        % sf data
        all_control_sf = [cont_cno_sf_tf_struct(1,:).control_sf];
        all_cno_sf = [cont_cno_sf_tf_struct(1,:).cno_sf];
        mean_control_sf = mean(all_control_sf);
        sem_control_sf = sem(all_control_sf);
        mean_cno_sf = mean(all_cno_sf);
        sem_cno_sf = sem(all_cno_sf);

        % tf data
        all_control_tf = [cont_cno_sf_tf_struct(1,:).control_tf];
        all_cno_tf = [cont_cno_sf_tf_struct(1,:).cno_tf];
        mean_control_tf = nanmean(60./all_control_tf);
        sem_control_tf = nansem(60./all_control_tf);
        mean_cno_tf = nanmean(60./all_cno_tf);
        sem_cno_tf = nansem(60./all_cno_tf);

        % tf/sf data
        tf_sf_control = (60./all_control_tf)./all_control_sf;
        tf_sf_cno = (60./all_cno_tf)./all_cno_sf;
        mean_tf_sf_control = nanmean(tf_sf_control);
        sem_tf_sf_control = nansem(tf_sf_control);
        mean_tf_sf_cno = nanmean(tf_sf_cno);
        sem_tf_sf_cno = nansem(tf_sf_cno);

        % all tf/sf ratios in case they are needed
        all_sfs = squeeze(control_stim_cond(:,:,1,1));
        all_tfs = 60./squeeze(control_stim_cond(:,:,1,2));
        tf_sf_ratios = unique(all_tfs./all_sfs);

        subplot(1,3,1); hold on;    
        plot(zeros(length(all_control_sf)), all_control_sf(:),'ko','Markersize',5);
        plot(ones(length(all_cno_sf))*0.5, all_cno_sf(:),'ro','Markersize',5);
        for i=1:length(all_control_sf)
            plot([0 0.5],[all_control_sf(i) all_cno_sf(i)], 'k-');
        end
        plot(0,mean_control_sf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_sf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_sf,sem_control_sf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_sf,sem_cno_sf,'r-',  'Linewidth', 1.5);
        yl(0,0.1);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred SF');

        subplot(1,3,2); hold on;    
        plot(zeros(length(all_control_tf)), (60./all_control_tf),'ko','Markersize',5);
        plot(ones(length(all_cno_tf))*0.5, (60./all_cno_tf),'ro','Markersize',5);
        for i=1:length(all_control_tf)
            plot([0 0.5],[(60/all_control_tf(i)) (60/all_cno_tf(i))], 'k-');
        end
        plot(0,mean_control_tf,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_cno_tf,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_control_tf,sem_control_tf,'k-','Linewidth',1.5);
        errorbar(0.5,mean_cno_tf,sem_cno_tf,'r-',  'Linewidth', 1.5);
        yl(0,10);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF');

        subplot(1,3,3); hold on;    
        plot(zeros(length(tf_sf_control)), tf_sf_control(:),'ko','Markersize',5);
        plot(ones(length(tf_sf_cno))*0.5, tf_sf_cno(:),'ro','Markersize',5);
        for i=1:length(tf_sf_control)
            plot([0 0.5],[tf_sf_control(i) tf_sf_cno(i)], 'k-');
        end
        plot(0,mean_tf_sf_control,'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,mean_tf_sf_cno,'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,mean_tf_sf_control,sem_tf_sf_control,'k-','Linewidth',1.5);
        errorbar(0.5,mean_tf_sf_cno,sem_tf_sf_cno,'r-',  'Linewidth', 1.5);
        yl(0,900);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('preferred TF/SF (º/s)');

        % Make some plots of control vs cno data
        figure(2); clf;
        all_control_dsi = [cont_cno_sf_tf_struct(1,:).control_dsi];
        all_cno_dsi = [cont_cno_sf_tf_struct(1,:).cno_dsi];
        all_control_osi = [cont_cno_sf_tf_struct(1,:).control_osi];
        all_cno_osi = [cont_cno_sf_tf_struct(1,:).cno_osi];
        all_control_dir = [cont_cno_sf_tf_struct(1,:).control_dir];
        all_cno_dir = [cont_cno_sf_tf_struct(1,:).cno_dir];
        
        subplot(1,2,1); hold on;    
        plot(zeros(length(all_control_dsi)), all_control_dsi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_dsi))*0.5, all_cno_dsi(:),'ro','Markersize',5);
        for i=1:length(all_control_dsi)
            plot([0 0.5],[all_control_dsi(i) all_cno_dsi(i)], 'k-');
        end
        plot(0,nanmean(all_control_dsi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_dsi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_dsi),nansem(all_control_dsi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_dsi),nansem(all_cno_dsi),'r-',  'Linewidth', 1.5);
        yl(-1.0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('DSI');
        
        subplot(1,2,2); hold on;    
        plot(zeros(length(all_control_osi)), all_control_osi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_osi))*0.5, all_cno_osi(:),'ro','Markersize',5);
        for i=1:length(all_control_osi)
            plot([0 0.5],[all_control_osi(i) all_cno_osi(i)], 'k-');
        end
        plot(0,nanmean(all_control_osi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_osi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_osi),nansem(all_control_osi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_osi),nansem(all_cno_osi),'r-',  'Linewidth', 1.5);
        yl(0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('OSI');
        

        % Make some histograms of the SF and TF
        % SF
        sorted_sfs = unique(sort(all_sfs));
        figure(3); clf;
        subplot(1,3,1); hold on;
        histogram(log2(all_control_sf),length(sorted_sfs),...
            'BinEdges',[log2(sorted_sfs(1))-0.5:1:log2(sorted_sfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(all_control_sf)), 0.5, 'v', 'MarkerSize', 5, 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(all_cno_sf),length(sorted_sfs),...
            'BinEdges',[log2(sorted_sfs(1))-0.5:1:log2(sorted_sfs(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(all_cno_sf)), 0.5, 'v', 'MarkerSize', 5, 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred SF');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(sorted_sfs(1))-0.75 log2(sorted_sfs(end))+0.75]);
        xticks(log2([sorted_sfs]));
        xticklabels(([sorted_sfs]));


        % TF
        % No TF plot because only one TF was used

        % TF/SF (approximation of cell's speed tuning)
        tf_over_sf = unique(sort(all_tfs./all_sfs));
        subplot(1,3,3); hold on;
        histogram(log2(tf_sf_control),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [0 0 0]);
        plot(nanmean(log2(tf_sf_control)), 0.5, 'v', 'MarkerSize', 5,  'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        histogram(log2(tf_sf_cno),length(tf_over_sf),...
            'BinEdges',[log2(tf_over_sf(1))-0.5:1:log2(tf_over_sf(end))+0.5],...
            'Normalization','probability', ...
            'DisplayStyle','stairs', 'LineWidth', 1.5, 'EdgeColor', [1 0 0]);
        plot(nanmean(log2(tf_sf_cno)), 0.5, 'v', 'MarkerSize', 5,  'Color', [1 0 0], 'MarkerFaceColor', [1 0 0]);
        xlabel('preferred TF/SF (º/s)');
        ylabel('proportion');
        ylim([0 1.0]);
        xlim([log2(tf_over_sf(1))-0.75 log2(tf_over_sf(end))+0.75]);
        xticks(log2(tf_over_sf(1)):1:log2(tf_over_sf(end)));

    else
        % Make some plots of control vs cno data
        figure(2); clf;
        all_control_dsi = [cont_cno_sf_tf_struct(1,:).control_dsi];
        all_cno_dsi = [cont_cno_sf_tf_struct(1,:).cno_dsi];
        all_control_osi = [cont_cno_sf_tf_struct(1,:).control_osi];
        all_cno_osi = [cont_cno_sf_tf_struct(1,:).cno_osi];
        all_control_dir = [cont_cno_sf_tf_struct(1,:).control_dir];
        all_cno_dir = [cont_cno_sf_tf_struct(1,:).cno_dir];
        
        subplot(1,2,1); hold on;    
        plot(zeros(length(all_control_dsi)), all_control_dsi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_dsi))*0.5, all_cno_dsi(:),'ro','Markersize',5);
        for i=1:length(all_control_dsi)
            plot([0 0.5],[all_control_dsi(i) all_cno_dsi(i)], 'k-');
        end
        plot(0,nanmean(all_control_dsi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_dsi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_dsi),nansem(all_control_dsi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_dsi),nansem(all_cno_dsi),'r-',  'Linewidth', 1.5);
        yl(-1.0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('DSI');
        
        subplot(1,2,2); hold on;    
        plot(zeros(length(all_control_osi)), all_control_osi(:),'ko','Markersize',5);
        plot(ones(length(all_cno_osi))*0.5, all_cno_osi(:),'ro','Markersize',5);
        for i=1:length(all_control_osi)
            plot([0 0.5],[all_control_osi(i) all_cno_osi(i)], 'k-');
        end
        plot(0,nanmean(all_control_osi),'ko','Linewidth',1.5,'Markersize',5);
        plot(0.5,nanmean(all_cno_osi),'ro','Linewidth',1.5,'Markersize',5);
        errorbar(0,nanmean(all_control_osi),nansem(all_control_osi),'k-','Linewidth',1.5);
        errorbar(0.5,nanmean(all_cno_osi),nansem(all_cno_osi),'r-',  'Linewidth', 1.5);
        yl(0,5.0);
        xl(-0.5,1.0);
        xticks([0,0.5]);
        xticklabels({'Control'; 'CNO'});
        ylabel('OSI');
    end
    
else % compiling datasets
    
    animal_ids = {
        '6627'
        '6627'
        };
        
    experiment_days = {
        'day4'
        'day4'
        };

    files_to_compare = {
        'day4_000_001_control.mat'
        'day4_000_001_control.mat'
        };
    
    % Load the saved data files for control and cno data and store them as
    % separate variables
    for f = 1:length(files_to_compare)
        data_dir = ['/Volumes/EPHYS_000/GCaMP_Data/', animal_ids{f}, '/', experiment_days{f},...
        '/Analysis'];

        cd (data_dir)
        
        file = files_to_compare{f};
        
        if ~isempty(strfind(file, 'control'))
            load (['Control/', file]);
        elseif ~isempty(strfind(file, 'cno'))
            load (['CNO/', file]);
        end
        
        if f == 1
            % Just get the responsive and reliable cells and make a new
            % structure
            resp_rel_tuning = all_tuning(all_resp_rel_idx);
            % If the first file, pass the responsive and reliable structure
            % to the compiled structure
            compiled_tuning = resp_rel_tuning;
            stim_cond = all_stim_cond_mat;
            resp_rel_idx = all_resp_rel_idx;
        else
            resp_rel_tuning = all_tuning(all_resp_rel_idx);
            compiled_tuning = [compiled_tuning, resp_rel_tuning];
            stim_cond = cat(4, stim_cond, all_stim_cond_mat);
            resp_rel_idx = [resp_rel_idx, all_resp_rel_idx];
        end
    end

    % Clear loaded variables to save memory and reduce variables
    clear all_tuning all_stim_cond_mat all_resp_rel_idx
end 

%% bubble plot sizes match scatter plot sizes (linear scale)
% x = squeeze(control_stim_cond(:,1,1));
% y = squeeze(control_stim_cond(:,1,1));
% 
% 
% [X, Y] = ndgrid(0:0.1, 0:0.1);
% x = X(:);
% y = Y(:);
% s = 1:numel(x);
% figure;
% p = bubble(x, y, s, 'b', 'EdgeColor', 'none');
% hold on;
% h = scatter(x, y, s, 'r');

% mean_dsi = nanmean(all_dsi_osi_dir(1:3:end));
% sem_dsi = sem(all_dsi_osi_dir(1:3:end));
% mean_osi = nanmean(all_dsi_osi_dir(2:3:end));
% sem_osi = sem(all_dsi_osi_dir(2:3:end));
% figure(1); clf; hold on
% dsi_osi = bar([mean_dsi; mean_osi]);
% dsi_osi.BarWidth = 0.5;
% errorbar(1,mean_dsi,sem_dsi,'k-');
% errorbar(2,mean_osi,sem_osi,'k-');
% axis tight
% xl(0.5,2.5);
% xticks([1,2]);
% xticklabels(['DSI'; 'OSI']);