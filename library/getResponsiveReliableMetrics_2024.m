% Calculate responsiveness of cells
%
% This is done by grouping all mean responses of a given stimulus
% (direction) and comparing all directions and blanks using an ANOVA

% Create matrix to hold output. Col 1 is the cell number. Col 2 is the
% p-value of the ttest comparing the max response against the blanks. Col 3
% is is the p-value of the ANOVA comparing all the directions against each
% other. Col 4 is the d' for the cell. Col 5 is the mean response to the
% stimulus that produced the maximum response. Col 6-8 are the parameters
% of the stimulus (6: direction, 7: SF, 8: TF) that produced the maximal
% response

function respCellMat = getResponsiveReliableMetrics_2024(cellIds, meanDfofCellTrials, maxDfofCellTrials,...
    dfofCellTrials, trialStim, preTime, stimTime, fps, trials)
    try
        respCellMat = nan(length(cellIds),8);
        respCellMat(:,1) = cellIds;

        for i = 1:size(cellIds,1)
            % Get all the mean responses for the cell
            dfofAllMeans = squeeze(meanDfofCellTrials(i,trials,:));
            stimForTrials = trialStim(trials,:);

            % Create a matrix to hold the means and the stimulus parameters
            dfofMeanParam = nan(length(dfofAllMeans),4);

            % Populate with mean values and stimulus parameters
            dfofMeanParam(:,1) = dfofAllMeans';
            dfofMeanParam(:,2) = stimForTrials(:,1);
            dfofMeanParam(:,3) = stimForTrials(:,2);
            dfofMeanParam(:,4) = stimForTrials(:,3);

            % Find the stimulus that generated the maximum mean response
            maxRespIndex = find(dfofAllMeans == max(dfofAllMeans));
            maxRespStim = stimForTrials(maxRespIndex,:);
            allMaxResp = dfofAllMeans(find(stimForTrials(:,1) == maxRespStim(1) & stimForTrials(:,2) == maxRespStim(2) & stimForTrials(:,3) == maxRespStim(3)));
            allBlankResp = dfofAllMeans(find(stimForTrials(:,1) == 500));

            [t, p_ttest] = ttest2(allMaxResp,allBlankResp);

            % Sort based on stimulus parameters
            dfofMeanParam = sortrows(dfofMeanParam,2);

            % Now exclude blanks
%             nonBlankTrials = find(dfofMeanParam(:,2) < 500);
%             dfofMeanTrials = dfofMeanParam(nonBlankTrials,:);

            % Now perform one-way ANOVA on the data
        %     [p, anovatab, stats] = anova1(dfofMeanTrials(:,1),dfofMeanTrials(:,2),'off');

            [p_anova, anovatab, stats] = anova1(dfofMeanParam(:,1),dfofMeanParam(:,2),'off');

            % Store the p-values in the matrix
            respCellMat(i,2) = p_ttest;
            respCellMat(i,3) = p_anova;

            nStimParams = unique(trialStim, 'rows');

            meanRespVals = nan(length(nStimParams),4);

            % Get all the mean dfof responses for each stimulus combination
            % and average them together and store the mean value and direction in a
            % matrix
            for j = 1:length(nStimParams)
                tmpVals = dfofMeanParam(find(dfofMeanParam(:,2) == nStimParams(j,1) & dfofMeanParam(:,3) == nStimParams(j,2) & dfofMeanParam(:,4) == nStimParams(j,3)),:);
                tmpMeanResp = mean(tmpVals(:,1));
                meanRespVals(j,1) = tmpMeanResp;
                meanRespVals(j,2) = nStimParams(j,1);
                meanRespVals(j,3) = nStimParams(j,2);
                meanRespVals(j,4) = nStimParams(j,3);
            end

            % Find the direction with the maximum mean dfof value
            maxIndex = find(meanRespVals(:,1) == max(meanRespVals(:,1)));

            % Store the max mean dfof values and the corresponding direction in the
            % matrix
            respCellMat(i,5) = meanRespVals(maxIndex,1);
            respCellMat(i,6) = meanRespVals(maxIndex,2);
            respCellMat(i,7) = meanRespVals(maxIndex,3);
            respCellMat(i,8) = meanRespVals(maxIndex,4);
            
        end
    catch
        warning('Could not compute anovas likely because there were no trials that occured during the specified behavior.');        
        i
    end
      
    try
        % Calculate d' for cells

        % Create NaN matrix to hold d' values for each cell
        dPrimeMat = nan(size(cellIds,1),1);

        for i = 1:size(cellIds,1)
            % Get the max dF/F values for all trials
            cellMaxDfof = squeeze(maxDfofCellTrials(i,trials,:));
            % Find the index of the max dF/F
            prefIndex = find(cellMaxDfof == max(cellMaxDfof));
            % Get the stimulus parameter (direction) for the preferred stimulus
            prefStim = stimForTrials(prefIndex,:);

            % If the preferred direction is a blank trial, then the d' has to be
            % zero and the cell is not reliable (and probably not responsive
            % either)
            if prefStim == 500
                dPrimeMat(i,:) = 0;
            else
                % Mean and STD of responses to preferred direction
                allPrefIndex = find(stimForTrials(:,1) == prefStim(1) & stimForTrials(:,2) == prefStim(2)...
                    & stimForTrials(:,3) == prefStim(3));
                dfofTrials = squeeze(dfofCellTrials(i,trials,:));
                prefRespAll = mean(dfofTrials(allPrefIndex,(fps*preTime+1):(fps*preTime+1) + (fps*stimTime)));
                prefRespMean = mean(prefRespAll);
                prefRespStd = std(prefRespAll);

                % Mean and STD of responses to blanks
                allBlankIndex = find(stimForTrials(:,1) == 500);
                blankRespAll = mean(dfofTrials(allBlankIndex,(fps*preTime+1):(fps*preTime+1) + (fps*stimTime)));
                blankRespMean = mean(blankRespAll);
                blankRespStd = std(blankRespAll);    

                dPrimeMat(i,:) = (prefRespMean - blankRespMean) / sqrt((prefRespStd + blankRespStd)/2);
            end
        end

        respCellMat(:,4) = dPrimeMat;
            
        
    catch
        warning('Could not compute dprime likely because there were no trials that occured during the specified behavior.');
        cellIds(i)
    end
end