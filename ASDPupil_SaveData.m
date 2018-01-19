%% Load individual data & specify preprocessing filter type
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% Pupil data preprocessing filter type
filter = struct;
filter.filterType = 'sgolay';
filter.order = 3; % order of polynomial for sgolay filter?
filter.framelen = 21; % length of window? must be odd number
filter.clearWin = 0; % delete the n surrounding data points of a blink
filter.velThreshold = 2; % de-blinking relative velocity threshold
graph = false;


for subjidx = 1:length(subj)
    
    filename = fullfile(root, ['ASD' num2str(subj(subjidx)) '.txt']);
    s = tdfread(filename);

    %% Extract pupil data from different time period

    %% 1. Initial pupil response
    % Extract the 2sITI + 3sLottery presentation + 5s Delay
    sInitial = struct; % sInitial stores pupil size for 8s time line
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


    % dblink, interpolate and smooth data(preprocessing). Filtered data signal, and filted z-scored data signal
    sInitial.PupilLeft_filt=zeros(size(sInitial.PupilLeft));
    for i=1:size(sInitial.PupilLeft,1)
        [sInitial.PupilLeft_filt(i,:),sInitial.PupilLeft_filtz(i,:),sInitial.Velleft_filt(i,:)] =...
            pupilPrepro(sInitial.PupilLeft(i,:),sInitial.PupilLeft(i,:), sInitial.Timestamp(i,:), filter, graph);
    end


    % ambig level and value for each trial
    sInitial.AL = str2num(sInitial.Trial(:,6:8));
    sInitial.Val = str2num(sInitial.Trial(:,10:11));
    
    save(fullfile(root,'Matlab data','pupildata',['ASD' num2str(subj(subjidx)) '_Initial.mat']),'sInitial')
    
    clear sInitial
    clear s
end