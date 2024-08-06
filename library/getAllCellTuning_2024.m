% Analyzes responses of cells to drifting grating stimuli to determine
% direction and orientation tuning of the cells.
% Col 1 is the cell ID. Col 2 is the DSI.
% Col 3 is the OSI.

function dsiOsiMat = getAllCellTuning_2024(respCellMat, cellIds, nStimParams, trialStim, meanDfofCellTrials, trials, nRespThresh)
    try
        %% Calculate tuning of cells

        stimForTrials = trialStim(trials);
        respRelCells = cellIds;

        % Create matrix to hold output. Col 1 is the cell ID. Col 2 is the DSI.
        % Col 3 is the OSI.
        dsiOsiMat = nan(length(respRelCells),3);
        dsiOsiMat(:,1) = respRelCells(:,1);

        % Take the list of responsive and reliable cells and calculate DSI and OSI
        % indices for them using the preferred directions stored in respCellMat
        for i = 1:length(respRelCells)
            cellInd = find(respCellMat(:,1) == respRelCells(i));
            prefDir = respCellMat(cellInd,6);
            prefResp = respCellMat(cellInd,5);

            % Get directions, but exclude last one, because these are the blanks
            stimDirs = nStimParams(1:end-1,1);
            % Determine how many directions were presented between a preferred
            % direction and the opposite direction based on the total number of
            % directions presented.
            nullShift = length(stimDirs)/2;
            orthShift = nullShift/2;

            % Get the index positions of the preferred direction, null direction,
            % and two orthogonal directions
            if prefDir == 500
                nullDir = 500;
                orthDirOne = 500;
                orthDirTwo = 500;
            else
                prefIndex = find(stimDirs == prefDir);
                nullStimDirs = circshift(stimDirs,nullShift);
                nullDir = nullStimDirs(prefIndex);
                orthStimDirsOne = circshift(stimDirs,orthShift);
                orthDirOne = orthStimDirsOne(prefIndex);
                orthStimDirsTwo = circshift(stimDirs,-orthShift);
                orthDirTwo = orthStimDirsTwo(prefIndex);
            end

            % Get the trials that correspond to stimulus presented at preferred,
            % null, and orthogonal directions.
            prefTrials = find(stimForTrials(:,1) == prefDir);
            nullTrials = find(stimForTrials(:,1) == nullDir);
            orthOneTrials = find(stimForTrials(:,1) == orthDirOne);
            orthTwoTrials = find(stimForTrials(:,1) == orthDirTwo);

            % Get mean dF/F values for the cell being investigated
            dfofAllMeans = squeeze(meanDfofCellTrials(i,trials,:))';

            % Get all mean values for a given stimulus (i.e. direction)
            allDfofPref = dfofAllMeans(prefTrials,1);
            allDfofNull = dfofAllMeans(nullTrials,1);
            allDfofOrthOne = dfofAllMeans(orthOneTrials,1);
            allDfofOrthTwo = dfofAllMeans(orthTwoTrials,1);

            dfofPref = mean(allDfofPref);
            dfofNull = mean(allDfofNull);
            dfofOrth = mean([allDfofOrthOne',allDfofOrthTwo']);

            dsi = (dfofPref - dfofNull)/(dfofPref + dfofNull);
            osi = (dfofPref - dfofOrth)/(dfofPref + dfofOrth);

            if length(allDfofPref) >= nRespThresh && length(allDfofNull) >= nRespThresh
                dsiOsiMat(i,2) = dsi;
            end
            if length(allDfofPref) >= nRespThresh && length(allDfofOrthOne) >= nRespThresh && length(allDfofOrthTwo) >= nRespThresh
                dsiOsiMat(i,3) = osi;
            end
        end
    catch
        
    end
end
