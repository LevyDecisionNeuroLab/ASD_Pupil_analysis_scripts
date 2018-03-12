% This scripts load and analyze preprocessed choice data 

clearvars
close all

root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata_extracted\'); 

 choiceall = zeros(100*length(subj),5);

 for subjidx = 1:length (subj)
%      subjidx = 1;
    dataname = ['ASD' num2str(subj(subjidx)) '_Initial.mat'];
    load([datafold,dataname])

    %% Create a matrix for writing data into text
    % Column name: 1-subj, 2-al, 3-val,4-choice(lottery 1,reference 0, no response 2), 5-RT
    for j = 1:length(sInitial.Trial)
        choiceall(j+100*(subjidx-1),1) = subj(subjidx);
        choiceall(j+100*(subjidx-1),2) = sInitial.AL(j);
        choiceall(j+100*(subjidx-1),3) = sInitial.Val(j);
        choiceall(j+100*(subjidx-1),4) = sInitial.Choice(j);       
        choiceall(j+100*(subjidx-1),5) = sInitial.Rt(j);
    end
    
 end
 
 
%% Write data, choice for all subject and all trials

% results file
 fid = fopen(['ASDChoiceAll.txt'],'w');
 fprintf(fid,'%s\t %s\t %s\t %s\t %s\n', 'Subj','Al', 'Val', 'Choice', 'Rt');
 fprintf(fid, '%d\t %d\t %d\t %d\t %d\n',choiceall');
 fclose(fid)
%% Write data into text file, prepare the logistic regression, for individual subject




