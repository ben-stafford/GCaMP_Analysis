% Load and analyze data from Suite2P files
clear all

experimentName = '6627';
recordingDay = 'day5';

% Change this to match the directory that contains the Suite2P mat file you
% want to load and analyze

% suite2pDir = ['~/Dropbox/Salk/GCaMP/GCaMP_Data/', experimentName, '/', recordingDay,...
%     '/Suite_2P'];
% 
% scanboxDir = ['~/Dropbox/Salk/GCaMP/GCaMP_Data/', experimentName, '/', recordingDay,...
%     '/Scanbox_Files'];

suite2pDir = ['/Volumes/EPHYS_000/GCaMP_Data/', experimentName, '/', recordingDay,...
    '/Suite_2P'];

scanboxDir = ['/Volumes/EPHYS_000/GCaMP_Data/', experimentName, '/', recordingDay,...
    '/Scanbox_Files'];

% Flag to analyze all files in the Scanbox_Files directory
analyze_all_files = 1;

% Flag to turn on/off 5 frame smoothing of fluorescence traces
smoothTraces = 1;

% Flag to turn on/off neuropil subtraction
subtractNeuropil = 1;

% Flag to turn on/off photobleaching decay
photoBleachSubtract = 0;


cd (suite2pDir)

% Get the mat file located in the directory
matFiles = dir('*.mat');
% Check for weird hidden files if the data has been copied from my Mac to
% the external drive
if size(matFiles,1) > 1
    goodMatFiles = {};
    file_names = extractfield(matFiles, 'name');
    for i = 1:length(file_names)
        if contains(file_names{i},'._') == 0
            goodMatFiles{end + 1} = file_names{i};
        end
    end
    oldMatFileName = goodMatFiles{1};
else
    oldMatFileName = matFiles(1).name;
end    
% oldMatFileName = matFiles(1).name;
dirPath = matFiles(1).folder;
dirPathStr = split('/',dirPath);

disp('Loading Suite 2P Data');

% % Load the mat file after renaming
load(oldMatFileName);

% This uses the mat file generated by the python version of Suite2P. It is
% different from the one generated by the Matlab version of Suite2P in that
% it has no data structure in it. Instead, it just has a handful of
% variables that correspond to the fluorscence traces, neuropil, and the
% rois that were deemed to be cells.
cell_ids = find(iscell(:,1)>0);

% This pulls out signals only from ROIs that have been deemed to
% actually be neurons
tracesAll = F(cell_ids,:);

% Neuropil subtraction
if subtractNeuropil
    disp('Performing neuropil subtraction');
    tracesAll = tracesAll - Fneu(cell_ids,:);
end

% smooth traces by mean of 5 frame window.
if smoothTraces
    disp('Smoothing fluorescence traces')
    [tracesAll, w] = smoothData(tracesAll,2,'movmean',5);
end

% Get list of filenames that were used in the Suite2P analysis

% Check to see if this list is already a variable in the mat file generated
% by Suite2P

if exist('h5fileList','var')
    hd5Files = h5fileList;
else    
    hd5FileDir = extractBefore(dirPath,char(dirPathStr{end}));
    cd (char(hd5FileDir));
    cd ('Scanbox_Files');
    hd5Files = dir('*.mat');
    % Check for hidden files that begin with '._' and exclude
    goodHd5Files = {};
    hd5file_names = extractfield(hd5Files, 'name');
    for i = 1:length(hd5file_names)
        if contains(hd5file_names{i},'._') == 0
            goodHd5Files{end + 1} = hd5file_names{i};
        end
    end
    hd5Files = goodHd5Files;
end

% Get info about the number of frames in the file you are analyzing

files_to_analyze = {};

if ~analyze_all_files
    disp(['Which file are you analyzing?']);
    [file, path] = uigetfile('*.mat');
    files_to_analyze{1} = file;
else
    files_to_analyze = hd5Files;
end

for m = 1:length(files_to_analyze)
    
    % cd to the Scanbox_Files directory at the start of each loop
    cd (char(hd5FileDir));
    cd ('Scanbox_Files');
    
    file_being_analyzed = files_to_analyze{m};
    
    file = file_being_analyzed;
    
    disp(['Analyzing: ', file]);
    
    for i = 1:length(hd5Files)
        fileStr = hd5Files{i};
        if strcmp(fileStr, file)
            nFile = i;
        end
    end

    tFrames = 0;

    % If the file being analyzed isn't the first one in the Suite2P analysis,
    % this code sums the total number of frames that were imaged before the
    % start of the file being analyzed. One frame is subtracted from the total
    % to try to compensate for small timing differences,
    if nFile > 1
        for i = 1:nFile-1
            matFileToLoad = [hd5Files{i}(1:end-4),'.mat'];
            load(matFileToLoad);
    %         disp([hd5Files{i}(1:end-4),'.mat'])
    %         info.totalFrame
            tFrames = tFrames + (info.totalFrame - 1);
        end
    end

    % Now load the mat file that corresponds to the file being analyzed and
    % then get the stimulus times, in frames, generated by scanbox
    load([file(1:end-4), '.mat']);

    % Stim times from TTLs generated by scanbox
    rawStimFrames = info.frame;

    % Stim times downsampled based on the number of imaging planes. This
    % method isn't 100% precise, but is probably OK.
    if ~isempty(info.otparam)
        stimFrames = floor(rawStimFrames/info.otparam(3));
    else
        stimFrames = rawStimFrames;
    end

    % The stimFrames matrix contains the time in frames of each trial.
    % Sometimes, the first value 0, which is incorrect. I don't yet know why
    % this happens, but my workaround is to check to see if the first value is
    % zero and, if it is, to trim the matrix by one value.
    if stimFrames(1) == 0
        stimFrames = stimFrames(2:end);
    end

    % Add the number of frames that have elapsed before the file that is being
    % analyzed occurs so that the stimulus on and off times are shifted to
    % align them with the fluorescence traces.
    stimFramesShift = stimFrames + tFrames;

    % Get the chunk of the raw fluorescence trace generated by Suite2P that
    % corresponds to the file you are analyzing for all the rois
    % tracesExpt = tracesAll(:,tFrames + 1:stimFrames(end) + (tFrames + 1));

    % Modification in 2024 because the length of tracesAll and the last values
    % of stimFrames sometimes aren't the same. Typically, they are off by only
    % one (1 frame), so the current hack is to add or subtract one from the
    % last value of stimFrames so it matches the length of tracesAll.
    if stimFrames(end) + tFrames == length(tracesAll)
        tracesExpt = tracesAll(:,tFrames + 1:stimFrames(end) + tFrames);
    elseif stimFrames(end) + tFrames < length(tracesAll)
        tracesExpt = tracesAll(:,tFrames + 1:stimFrames(end) + (tFrames + 1));
    elseif stimFrames(end) + tFrames > length(tracesAll)
        tracesExpt = tracesAll(:,tFrames + 1:stimFrames(end) + (tFrames - 1));
    end


    fpsRaw = 15;
    % Downsample the frame rate depending on how many planes you are
    % imaging
    if info.volscan~=0
        fps=fpsRaw./info.otparam(3);
    else
        fps = fpsRaw;
    end

    % Load analyzer file automatically
    analyzerFile = [file(1:end-4), '.analyzer'];
    analyzerPath = getPathSubset(pwd,1);
    cd(analyzerPath);
    cd('Analyzer_Files');
    analyzerPath = pwd;
    load([analyzerPath, '/',  analyzerFile], '-mat');
    [trial_num, stimTime] = looper_2024([analyzerPath, '/', analyzerFile]);

    % trial_num is a matrix where each row contains the direction, sf, and
    % tf of each stimulus presented
    
    nStimParams = unique(trial_num, 'rows');

    % Get the parameters of the stimuli, including pre/post-stim times as
    % well as stimulus duration
    preTime = stimTime(1);
    postTime = stimTime(3);
    stimTime = stimTime(2);

    % Check to see how many ttls there are per stimulus presentation. At some
    % point, Peichao modified scanbox so that it generates 4 ttls/stimulus
    % presentation instead of 2.
    % nTtlsPerStim = length(stimFramesShift)/length(trial_num);
    nTtlsPerStim = 4;

    % Load Spike2 text file automatically (save as spreadsheet text file)
    spikeFile = [file(1:end-4), '.txt'];
    spikePath = getPathSubset(pwd,1);
    cd(spikePath);
    cd('Spike2_Files');
    spikePath = pwd;
    stim_info = readtable([spikePath, '/', spikeFile]);
    stim_info = table2array(stim_info);

    % Voltage signal from photodiode.
    galvo = stim_info(:,6);
    % Voltage signal from one movement encoder
    movtA = round(stim_info(:,3)./max(stim_info(:,3))); %normalize to binary 0,1
    % Voltage signal from the other movement encoder
    movtB = round(stim_info(:,2)./max(stim_info(:,2)));
    % Sample number of the voltage measurements
    t_samp = stim_info(:,1);
    % Voltage signal from the scanbox output indicating start/stop time of
    % visual stimulus
    stim = stim_info(:,5);

    % Determine sampling rate of the digitizer by finding the sample number
    % that corresponds to 1 second.
    fs = find(t_samp == 1)-1;
    fs = fs(1);

    % Create array of times using the calculated sample rate and the total
    % number of samples in the encoder file.
    t_int = (1:size(stim,1)) * fs;

    % %plot to sanity check
    % figure; 
    % subplot(3,1,1); hold all; plot(t_int, galvo,'r'); plot(t_int, stim,'b');
    % subplot(3,1,2); plot(t_int, movtA);
    % subplot(3,1,3); plot(t_int, movtB);

    % Set threshold value for determining the start/stop time of each trial
    threshold = 3;
    % Find samples where visual stimulus was above a threshold voltage
    pulses = find(stim>threshold);

    % Find samples where  visual stimulus ended by finding the places where the
    % difference between values in pulses matrix were greater than one.
    downEdges = find(diff(pulses)>1)';
    startEdges = downEdges(1:2:end);
    endEdges = startEdges + 1;

    % Create matrix of start/stop samples for stimulus by sorting the start times
    % and stop times sequentially.
    allEdges = sort([startEdges endEdges]);

    % Times on/off
    on_off = pulses(allEdges);

    % On and off times in samples
    on_times = pulses(startEdges);
    off_times = pulses(endEdges);

    % Duration of stimuli in seconds
    times_s = (off_times-on_times)./fs(1);
    % Duration of inter-trial intervals in seconds
    gap = diff(off_times)./fs(1);

    % Determine which trials were running or stationary. Running if mean speed >0.5 cm/s
    % speed = + = wheel moving forward = backward running, -= wheel moving bw=
    % forward running.

    % Total number of stimuli presented, including blanks
    alltrials = size(trial_num,1);

    % Determine trials during which the mouse was moving at all. Output is a
    % matrix with a binary value for each trial
    [moveTrials, allSpeeds, meanSpeeds] = Intan_digital_movement_HW(alltrials,movtA, movtB, [on_times off_times], 0.5,fs); %size(trial_pos,1)

    % Create a matrix that contains all relevant parameters for each
    % stimulus presentation.
    %
    % Col 1: trial start in frames
    % Col 2: stimulus start in frames
    % Col 3: stimulus end in frames
    % Col 4: trial end in frames
    % Col 5: binary value indicating if animal was moving (1 = moving)
    % Col 6: direction of stimulus
    % Col 7: SF of stimulus
    % Col 8: TF of stimulus
    
    trialInfoMat = [stimFrames(1:nTtlsPerStim:end), stimFrames(2:nTtlsPerStim:end), stimFrames(3:nTtlsPerStim:end),...
        stimFrames(4:nTtlsPerStim:end), moveTrials, trial_num];

    % Subtract off photobleaching decay using polynomial fit to the blank trials
    if photoBleachSubtract == 1

        tracesBlankSub = nan(size(tracesExpt));

        % Get the stimulus times that correspond to blank trials
        blankTrialsIndex = find(trialInfoMat(:,6) == 500);
        blankTrials = trialInfoMat(blankTrialsIndex,1:2);

        % Shift back tFrames because the raw trace won't be shifted.
        blankTrials = blankTrials;

        blankXvals = [];

        for i = 1:length(blankTrials)
            blankXvals = [blankXvals, blankTrials(i,1):1:blankTrials(i,2)];
        end

        for i = 1:length(cell_ids)
            % Create matrix the size of the experiment being analyzed in frames
            decayTrace = nan(1,stimFrames(end,1));

            % Get the raw fluorescence trace for the cell being analyzed
            rawCellTrace = tracesExpt(i,:);

            for j = 1:length(blankTrials)
                blankRawTrace = rawCellTrace(blankTrials(j,1):blankTrials(j,2));
                decayTrace(blankTrials(j,1):blankTrials(j,2)) = blankRawTrace;

            end

            blankYvalInd = find(~isnan(decayTrace));
            blankYvals = decayTrace(blankYvalInd);
            blankFitVals = polyfit(blankXvals,blankYvals,1);
            fitXvals = (1:1:length(decayTrace) + 1);
            blankFit = blankFitVals(1)*fitXvals + blankFitVals(2);

            rawCellTraceSub = rawCellTrace - blankFit;
            tracesBlankSub(i,:) = rawCellTraceSub;

        end

        tracesExpt = tracesBlankSub;
    end


    % Split raw fluorescence traces and calculate dF/F

    % The trial times in frames in the scanbox file don't always produce
    % trials of the exact same number of samples. They seem to alternate by
    % one frame every other trial, so trim all traces to be the smallest
    % trail size.
    min_trial_length = floor(min(trialInfoMat(:,4)-trialInfoMat(:,1)));
    
    % Create NaN array to hold raw fluorescence traces for each trial for each cell
    rawCellTrials = nan(size(cell_ids,1), size(trial_num,1), min_trial_length);
    
    % Create NaN array to hold dF/F traces for each trial for each cell
    dfofCellTrials = nan(size(cell_ids,1), size(trial_num,1), min_trial_length);

    % Create NaN array to hold mean and max dF/F during stimulus presentation for each trial for each
    % cell
    meanDfofCellTrials = nan(size(cell_ids,1), size(trial_num,1), 1);
    maxDfofCellTrials = nan(size(cell_ids,1), size(trial_num,1), 1);


    for j = 1:size(cell_ids,1)
        % Get the full fluorescence trace for the cell being examined
        cellTrace = tracesExpt(j,:);

        % Go through each stimulus presentation and cut out the trace aligned
        % with the stimulus start and stop times in frames
        for i = 1:size(trialInfoMat,1)
            % trial start, end, and stimulus start times in frames
            trial_start = trialInfoMat(i,1);
            trial_end = trial_start + (min_trial_length - 1);
            trial_length = trial_end - trial_start;
            stim_start = trialInfoMat(i,2) - trialInfoMat(i,1);
            
            % Trace for trial being examined
            tmpTrace = cellTrace(trial_start:trial_end);
            % Baseline from trial for dF/F calculation
            tmpBaseline = mean(tmpTrace(1:stim_start));
            % Calculate dF/F
            dfofTrace = (tmpTrace-tmpBaseline)./tmpBaseline;
            % Find max dF/F during stimulus presentation, including post-stim
            maxDfofTrace = max(dfofTrace(stim_start:end));
            % Find mean df/F during stimulus presentation, including post-stim
            meanDfofTrace = mean(dfofTrace(stim_start:end));

            % Place values in matrices for later analysis
            rawCellTrials(j,i,1:size(tmpTrace,2)) = tmpTrace;
            dfofCellTrials(j,i,1:size(dfofTrace,2)) = dfofTrace;
            meanDfofCellTrials(j,i,:) = meanDfofTrace;
            maxDfofCellTrials(j,i,:) = maxDfofTrace;

        end

    end

    % Get the trials that occured when the mouse was running or stationary
    moveTrials = find(trialInfoMat(:,5) == 1);
    statTrials = find(trialInfoMat(:,5) == 0);
    allTrials = find(trialInfoMat(:,5) < 2);

    disp(['Calculating Responsive and Reliable Metrics'])

    % Get matrices of responsivity and reliability parameters for each cell
%     cellRespRelMoveMat = getResponsiveReliableMetrics_2024(cell_ids,meanDfofCellTrials,maxDfofCellTrials,dfofCellTrials,trial_num,preTime,stimTime,fps,moveTrials);
%     cellRespRelStatMat = getResponsiveReliableMetrics_2024(cell_ids,meanDfofCellTrials,maxDfofCellTrials,dfofCellTrials,trial_num,preTime,stimTime,fps,statTrials);
    cellRespRelMat = getResponsiveReliableMetrics_2024(cell_ids,meanDfofCellTrials,maxDfofCellTrials,dfofCellTrials,trial_num,preTime,stimTime,fps,allTrials);

%     disp(['Calculating Tuning'])
% 
%     dsiOsiMoveMat = getAllCellTuning_2024(cellRespRelMoveMat, cell_ids, nStimParams, trial_num, meanDfofCellTrials, moveTrials,0);
%     dsiOsiStatMat = getAllCellTuning_2024(cellRespRelStatMat, cell_ids, nStimParams, trial_num, meanDfofCellTrials, statTrials,0);
%     dsiOsiMat = getAllCellTuning_2024(cellRespRelMat, cell_ids, nStimParams, trial_num, meanDfofCellTrials, allTrials,0);

    % Creating a data structure to store this analysis that should be more
    % modular to account for running stimuli in which multiple parameters vary
    % (i.e. direction and SF, or SF and TF, or SF, TF, and direction). This is
    % organized by cell, which means trial data gets repeated a lot so it's not
    % very memory efficient.
    
    % rename a few variables before saving them
    cellIds = cell_ids;
    trialStim = trial_num;
    trialMoveInfo = trialInfoMat(:,5:6);
    stimInfo = [fps, preTime, stimTime, postTime];
    
    % Create a directory to save the data in
    cd ..
    if ~exist('Analysis', 'dir')
       mkdir('Analysis')
    end
    cd 'Analysis'

    
    disp(['Saving files'])

    
    % Clear the old data structure if it exists
    if exist('gcamp_data','var')
        clear gcamp_data;
    end
    gcamp_data(length(cellIds)) = struct();

    % Start making columns in the data structure to store different parts of
    % the data

    % Cell IDs
    gcamp_data(1).cell_ids = [];
    % Dfof trace for each trial presented for each cell. Each entry (row) will
    % be m x n matrix, where m = number of trials and n = the dfof trace
    % for that trial
    gcamp_data(1).dfof_traces = [];
    % matrix of mean dfofs for all trials for each cell
    gcamp_data(1).mean_dfof_resp = [];
    % matrix of max dfof for all trials for each cell
    gcamp_data(1).max_dfof_resp = [];
    % Responsive p-values for t-test
    gcamp_data(1).ttest_resp = [];
    % Responsive p-values for anova
    gcamp_data(1).anova_resp = [];
    % d-prime scores
    gcamp_data(1).dprime = [];
    % reliability threshold pass (1 = reliable cell/ROI)
    gcamp_data(1).reliable = [];


    % Loop through the list of cells and put data in the structure
    for i = 1:length(cellIds)
        gcamp_data(i).cell_ids = cellIds(i);
        gcamp_data(i).dfof_traces = squeeze(dfofCellTrials(i,:,:));
        gcamp_data(i).mean_dfof_resp = squeeze(meanDfofCellTrials(i,:,:));
        gcamp_data(i).max_dfof_resp = squeeze(maxDfofCellTrials(i,:,:));
        gcamp_data(i).ttest_resp = cellRespRelMat(i,2);
        gcamp_data(i).anova_resp = cellRespRelMat(i,3);
        gcamp_data(i).dprime = cellRespRelMat(i,4);
        gcamp_data(i).reliable = cellRespRelMat(i,9);
    end
    
    % Get maximum mean fluorescence of all ROIs
    all_mean_dfof = [];
    all_mean_dfof = [all_mean_dfof, gcamp_data(1,:).mean_dfof_resp];
    max_mean_dfof = max(all_mean_dfof);
    fluor_thresh = max_mean_dfof * 0.06;
    
    for i = 1:size(gcamp_data,2)
        max_roi_mean_dfof = max(gcamp_data(i).mean_dfof_resp);
        if max_roi_mean_dfof <= fluor_thresh
            gcamp_data(i).good_fluor = 0; % not reliable, too low
        else
            gcamp_data(i).good_fluor = 1; % reliable, above threshold
        end
    end


    % Stimulus conditions for each trial presented. This is a matrix that
    % contains (by column) the stimulus pre, start, stop, and post times (in frames), a
    % binary value indicating whether the animal was moving, then the various
    % stimulus conditions for each trial of the stimulus

    % Also save out the matrix that contains the stimulus conditions for each
    % trial presented. This is a matrix that contains (by column) the (1) pre,
    % (2) start, (3) stop, and (4) post stim times, (5) a binary value indicating
    % whether the animal was moving, then the parameters of the visual stimulus
    % that was presented in the following order: (6) orientation, (7) spatial
    % frequency, (8) temporal frequency.

    stimulus_info = trialInfoMat;

    save([file(1:end-4), '_Data'], 'gcamp_data', 'stimulus_info');

end


