function [trial_num, stim_time] = looper_2024(fa,stf_stim)
    %Looper file

    if nargin<2
        stf_stim=1;
    end
    load (fa, '-mat');

    bflag = 0; %are there blank trials?
    if strcmp(Analyzer.loops.conds{length(Analyzer.loops.conds)}.symbol{1},'blank')
        bflag = 1;
    end
    cond = length(Analyzer.loops.conds); %total number of combinations of conditions (last cond is blank condition)
    stim_types = Analyzer.loops.conds{1}.symbol;
    % Create the matrix to hold all the stimulus conditions
    trial_num = ones(3,(length(Analyzer.loops.conds)-1)*length(Analyzer.loops.conds{1}.repeats)+(length(Analyzer.loops.conds{cond}.repeats)))'; %9x20=180 even tho 175 total, bc 20 repeats for stim + 15 blanks
    if ~any(strcmp(stim_types, 'ori'))
        trial_num(:,1) = Analyzer.P.param{1,21}{3};
    end
    if ~any(strcmp(stim_types, 's_freq'))
        trial_num(:,2) = Analyzer.P.param{1,25}{3};
    end
    if ~any(strcmp(stim_types, 't_period'))
        trial_num(:,3) = Analyzer.P.param{1,31}{3};
    end

    % s_freq = Analyzer.P.param{1,25};
    % t_freq = Analyzer.P.param{1,31};

    if bflag == 1 %if so, here's how you deal with it
        for i = 1:length(Analyzer.loops.conds)-1 %this is -1 because there are 1 blanks
    %         % Create the matrix to hold all the stimulus conditions
    %         if i == 1
    %         %         trial_num=ones(length(Analyzer.loops.conds{i}.symbol),(length(Analyzer.loops.conds)-1)*length(Analyzer.loops.conds{i}.repeats)+(length(Analyzer.loops.conds{cond}.repeats)))'; %9x20=180 even tho 175 total, bc 20 repeats for stim + 15 blanks
    %             trial_num = ones(3,(length(Analyzer.loops.conds)-1)*length(Analyzer.loops.conds{i}.repeats)+(length(Analyzer.loops.conds{cond}.repeats)))'; %9x20=180 even tho 175 total, bc 20 repeats for stim + 15 blanks
    %         end

            for k = 1:length(Analyzer.loops.conds{i}.symbol)
                % This gets the value of the stimulus condition
                trial_val = Analyzer.loops.conds{i}.val{k};
                % This gets the string representing the stimulus parameter (i.e.
                % ori, sf, or tf
                trial_type = Analyzer.loops.conds{i}.symbol{k};
                % Find out which type it is and adjust column position accordingly
                if strcmp(trial_type,'ori') == 1
                    col = 1;
                elseif strcmp(trial_type,'s_freq') == 1
                    col = 2;
                elseif strcmp(trial_type,'t_period') == 1
                    col = 3;
                end

                for j = 1:length(Analyzer.loops.conds{i}.repeats)
                    aux_trial = Analyzer.loops.conds{i}.repeats{j}.trialno;
                    trial_num(aux_trial,col) = trial_val;
                end
        %       trial_vals(:,k)= Analyzer.loops.conds{i}.val{k};
            end
        %     for j=1:length(Analyzer.loops.conds{i}.repeats)
        %       aux_trial=Analyzer.loops.conds{i}.repeats{j}.trialno;
        %       trial_num(aux_trial,:)=trial_vals;
        %     end 
        end
    else %if there's no blank trials, use all conds:
        for i=1:length(Analyzer.loops.conds)
           if i==1;trial_num=ones(length(Analyzer.loops.conds{i}.symbol),length(Analyzer.loops.conds)*length(Analyzer.loops.conds{i}.repeats))';end ;
              for k=1:length(Analyzer.loops.conds{i}.symbol)
                  trial_vals(:,k)= Analyzer.loops.conds{i}.val{k};
              end
              for j=1:length(Analyzer.loops.conds{i}.repeats)
                  aux_trial=Analyzer.loops.conds{i}.repeats{j}.trialno;
                  trial_num(aux_trial,:)=trial_vals;
              end 
        end
    end

    predelay = cell2mat(Analyzer.P.param{1}(3));
    postdelay = cell2mat(Analyzer.P.param{2}(3));
    stim_time = cell2mat(Analyzer.P.param{3}(3));
    stim_time = [predelay stim_time postdelay];

    [r,c]=find(trial_num==1); %this takes all blank trials (if there were any) and sets it to 256

    %if analyzer file is from a sf, tf, ori presenting stimulus set, find the
    %blanks and set to 500 bc 1 is too similar to a real value.
    if stf_stim==1
        trial_num(r,:)=500; % sort either symbol by 256 for blank trials 
    end
end
