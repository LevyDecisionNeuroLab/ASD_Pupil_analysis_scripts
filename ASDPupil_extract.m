% This script extract pupil size data from the Tobii Studio exported .txt file

%% History:
% Ruonan written 12.07.2017
% Ruonan updated 03.09.2018

clearvars
close all
%% Load individual data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
datafold = fullfile(root,'Matlab data','pupildata_extracted\'); 
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

%% Extract pupil data from different time period
for subjidx = 1:length(subj)
    
    subjidx = 2;
    
    filename = fullfile(root, ['ASD' num2str(subj(subjidx)) '.txt']);
    s = tdfread(filename);
    
    %% 1. Initial pupil response
    % Extract the 2sITI + 3sLottery presentation + 5s Delay
    sInitial = struct; % sInitial stores pupil size for 8s time line
    sInitial.subjid = subj(subjidx);
    i = 0; % counting index for going over each element in s in sequence
    j = 0; % No. for trial 1~25
    % sInitial.Trial = repmat(char(0),25,26);
    while i < length(s.Timestamp)
        i=i+1;
        if s.EventKey(i) ~= 8 || s.Descriptor(i,1)~= 'a'; continue; end
        if s.EventKey(i) == 8 && s.Descriptor(i,1) == 'a' 
            j=j+1;
            sInitial.Trial(j,:) = s.Descriptor(i,:);
        end
        k = 1; % No for each Pupil measurement in the time line
        for m = i-120:i+480
            sInitial.PupilLeft(j,k) = s.PupilLeft(m);
            sInitial.PupilRight(j,k)= s.PupilRight(m);
            sInitial.Timestamp(j,k) = s.Timestamp(m);
            k = k+1;
        end
        while s.Descriptor(i,3) ~= 'I'
            i=i+1;
        end
    end
    sInitial.PupilLeft(sInitial.PupilLeft == -1 | sInitial.PupilLeft == 0)=NaN;
    sInitial.PupilRight(sInitial.PupilRight == -1 | sInitial.PupilRight == 0)=NaN;
    sInitial.Timestamp(sInitial.Timestamp == 0)=NaN;


    for k = 1:size(sInitial.Timestamp,1)
        sInitial.Timestamp(k,:) = sInitial.Timestamp(k,:)-sInitial.Timestamp(k,1);
    end

    % % plot raw data
    % ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),sInitial.Timestamp(1,:),'name', 'PupilSize_Initial');
    % plot(ts_initialLeft)

    % % plot average with STD shade
    % average = nanmean(sInitial.PupilLeft(:,:));
    % std = nanstd(sInitial.PupilLeft(:,:));
    % % a function from toolbox kakearney boudedline-pkg
    % figure
    % boundedline(sInitial.Timestamp(1,:), average, std,'-b.')

    % ambig level and value for each trial
    sInitial.AL = str2num(sInitial.Trial(:,6:8));
    sInitial.Val = str2num(sInitial.Trial(:,10:11));

    savefilename = fullfile(datafold, ['ASD' num2str(subj(subjidx)) '_Initial.mat']);
    save(savefilename, 'sInitial')

    %% consider discard the following codes:
    % %% 2. ITI Pupil size last 2s
    % sITI = struct;
    % sITI.PupilLeft_int = sInitial.PupilLeft_int(:,1:120);
    % sITI.Timestamp = sInitial.Timestamp(:,1:120);
    % sITI.Trial = sInitial.Trial;
    % for k = 1:size(sITI.Timestamp,1)
    %     sITI.Timestamp(k,:) = sITI.Timestamp(k,:)-sITI.Timestamp(k,1);
    % end
    % for a = 1:size(sITI.PupilLeft_int,1)
    %     sITI.PupilLeftMean(a)= nanmean(sITI.PupilLeft_int(a,:));
    %     sITI.PupilLeftStd(a)= nanstd(sITI.PupilLeft_int(a,:));
    % %     sITI.PupilRightMean(a)= nanmean(sITI.PupilRight(a,:));
    % %     sITI.PupilRightStd(a)= nanstd(sITI.PupilRight(a,:));
    % end
    % 
    % 
    % %% 3. Lottery presentation Pupil size, last 1.8s (3s-1.2s)
    % sLott = struct;
    % sLott.PupilLeft_int = sInitial.PupilLeft_int(:,193:300);
    % sLott.Timestamp = sInitial.Timestamp(:,193:300);
    % sLott.Trial = sInitial.Trial;
    % for k = 1:size(sLott.Timestamp,1)
    %     sLott.Timestamp(k,:) = sLott.Timestamp(k,:)-sLott.Timestamp(k,1);
    % end
    % for a = 1:size(sLott.PupilLeft_int,1)
    %     sLott.PupilLeftMean(a)= nanmean(sLott.PupilLeft_int(a,:));
    %     sLott.PupilLeftStd(a)= nanstd(sLott.PupilLeft_int(a,:));
    % %     sITI.PupilRightMean(a)= nanmean(sITI.PupilRight(a,:));
    % %     sITI.PupilRightStd(a)= nanstd(sITI.PupilRight(a,:));
    % end
    % % normalized pupil size during Lott
    % for i= 1 : length(sITI.PupilLeftMean)
    %     sLott.PupilLeftNorm(i) = sLott.PupilLeftMean(i)/sITI.PupilLeftMean(i);
    % end
    % 
    % sLott.PupilLeftNormmean = nanmean(sLott.PupilLeftNorm);
    % 
    % 
    % %% 4. Delay pupil size,last 3.8 s
    % sDelay = struct;
    % sDelay.PupilLeft_int = sInitial.PupilLeft_int(:,length(sInitial.PupilLeft_int)-227:length(sInitial.PupilLeft_int));
    % sDelay.Timestamp = sInitial.Timestamp(:,length(sInitial.PupilLeft_int)-227:length(sInitial.PupilLeft_int));
    % sDelay.Trial = sInitial.Trial;
    % for k = 1:size(sDelay.Timestamp,1)
    %     sDelay.Timestamp(k,:) = sDelay.Timestamp(k,:)-sDelay.Timestamp(k,1);
    % end
    % for a = 1:size(sDelay.PupilLeft_int,1)
    %     sDelay.PupilLeftMean(a)= nanmean(sDelay.PupilLeft_int(a,:));
    %     sDelay.PupilLeftStd(a)= nanstd(sDelay.PupilLeft_int(a,:));
    % end
    % % normalized pupil size during delay
    % for i= 1 : length(sITI.PupilLeftMean)
    %     sDelay.PupilLeftNorm(i) = sDelay.PupilLeftMean(i)/sITI.PupilLeftMean(i);
    % end
end