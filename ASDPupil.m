%% History:
% Ruonan written 12.07.2017

% clearvars
% clc
%% Load individual data & specify preprocessing filter type
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
filename = fullfile(root, 'ASD25.txt');
% cd 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
s = tdfread(filename);

% Pupil data preprocessing filter type
filter = struct;
filter.filterType = 'sgolay';
filter.order = 3; % order of polynomial for sgolay filter?
filter.framelen = 21; % length of window? must be odd number
filter.clearWin = 0; % delete the n surrounding data points of a blink
filter.velThreshold = 2; % de-blinking relative velocity threshold


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

% plot raw data
ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),sInitial.Timestamp(1,:),'name', 'PupilSize_Initial');
plot(ts_initialLeft)

average = nanmean(sInitial.PupilLeft(:,:));
std = nanstd(sInitial.PupilLeft(:,:));
% a function from toolbox kakearney boudedline-pkg
figure
boundedline(sInitial.Timestamp(1,:), average, std,'-b.')


% dblink, interpolate and smooth data(preprocessing). Filtered data signal, and filted z-scored data signal
sInitial.PupilLeft_filt=zeros(size(sInitial.PupilLeft));
for i=1:size(sInitial.PupilLeft,1)
    [sInitial.PupilLeft_filt(i,:),sInitial.PupilLeft_filtz(i,:)] = pupilPrepro(sInitial.PupilLeft(i,:), sInitial.Timestamp(i,:), filter);
end

% timeseries and figure
ts_initialLeft_int = timeseries(sInitial.PupilLeft_int(:,:),sInitial.Timestamp(1,:),'name', 'PupilSize_Initial_int');
figure('Name','InitialResponseLeft','NumberTitle','off')
plot(ts_initialLeft_int)

% plot averaged time course with shaded error bar
average = nanmean(sInitial.PupilLeft_filtz(:,:));
std = nanstd(sInitial.PupilLeft_filtz(:,:));
se = std ./ sqrt(size(sInitial.PupilLeft_filtz,1));
% a function from toolbox kakearney boudedline-pkg
figure
boundedline(sInitial.Timestamp(1,:), average, se,'-b.')
figure
boundedline(sInitial.Timestamp(1,:), average, std,'-r.')

% plot(sInitial.Timestamp(1,:), average,'LineStyle', '-', 'Marker', '.');
% hold on
% errorbar(sInitial.Timestamp(1,:), average,std,'CapSize',0.1)

%% 2. ITI Pupil size last 2s
sITI = struct;
sITI.PupilLeft_int = sInitial.PupilLeft_int(:,1:120);
sITI.Timestamp = sInitial.Timestamp(:,1:120);
sITI.Trial = sInitial.Trial;
for k = 1:size(sITI.Timestamp,1)
    sITI.Timestamp(k,:) = sITI.Timestamp(k,:)-sITI.Timestamp(k,1);
end
for a = 1:size(sITI.PupilLeft_int,1)
    sITI.PupilLeftMean(a)= nanmean(sITI.PupilLeft_int(a,:));
    sITI.PupilLeftStd(a)= nanstd(sITI.PupilLeft_int(a,:));
%     sITI.PupilRightMean(a)= nanmean(sITI.PupilRight(a,:));
%     sITI.PupilRightStd(a)= nanstd(sITI.PupilRight(a,:));
end


% 3. Lottery presentation Pupil size, last 1.8s (3s-1.2s)
sLott = struct;
sLott.PupilLeft_int = sInitial.PupilLeft_int(:,193:300);
sLott.Timestamp = sInitial.Timestamp(:,193:300);
sLott.Trial = sInitial.Trial;
for k = 1:size(sLott.Timestamp,1)
    sLott.Timestamp(k,:) = sLott.Timestamp(k,:)-sLott.Timestamp(k,1);
end
for a = 1:size(sLott.PupilLeft_int,1)
    sLott.PupilLeftMean(a)= nanmean(sLott.PupilLeft_int(a,:));
    sLott.PupilLeftStd(a)= nanstd(sLott.PupilLeft_int(a,:));
%     sITI.PupilRightMean(a)= nanmean(sITI.PupilRight(a,:));
%     sITI.PupilRightStd(a)= nanstd(sITI.PupilRight(a,:));
end
% normalized pupil size during Lott
for i= 1 : length(sITI.PupilLeftMean)
    sLott.PupilLeftNorm(i) = sLott.PupilLeftMean(i)/sITI.PupilLeftMean(i);
end

sLott.PupilLeftNormmean = nanmean(sLott.PupilLeftNorm);




% 4. Delay pupil size,last 3.8 s
sDelay = struct;
sDelay.PupilLeft_int = sInitial.PupilLeft_int(:,length(sInitial.PupilLeft_int)-227:length(sInitial.PupilLeft_int));
sDelay.Timestamp = sInitial.Timestamp(:,length(sInitial.PupilLeft_int)-227:length(sInitial.PupilLeft_int));
sDelay.Trial = sInitial.Trial;
for k = 1:size(sDelay.Timestamp,1)
    sDelay.Timestamp(k,:) = sDelay.Timestamp(k,:)-sDelay.Timestamp(k,1);
end
for a = 1:size(sDelay.PupilLeft_int,1)
    sDelay.PupilLeftMean(a)= nanmean(sDelay.PupilLeft_int(a,:));
    sDelay.PupilLeftStd(a)= nanstd(sDelay.PupilLeft_int(a,:));
end
% normalized pupil size during delay
for i= 1 : length(sITI.PupilLeftMean)
    sDelay.PupilLeftNorm(i) = sDelay.PupilLeftMean(i)/sITI.PupilLeftMean(i);
end


% Time series and figures
ts_ITILeft_int = timeseries(sITI.PupilLeft_int,sITI.Timestamp(1,:),'name', 'PupilSize_ITI_int');
figure('Name','ITILeft','NumberTitle','off')
% hold on
plot (ts_ITILeft_int)

ts_LottLeft_int = timeseries(sLott.PupilLeft_int,sLott.Timestamp(1,:),'name', 'PupilSize_Lott_int');
figure('Name','LottLeft','NumberTitle','off')
% hold on
plot (ts_LottLeft_int)

ts_DelayLeft_int = timeseries(sDelay.PupilLeft_int,sDelay.Timestamp(1,:),'name', 'PupilSize_Delay_int');
figure('Name','DelayLeft','NumberTitle','off')
% hold on
plot (ts_DelayLeft_int)


%% Trial information, summary pupil size into ambig and val levels

% ambig level and value for each trial
sTrial = struct;
sTrial.AL = str2num(sInitial.Trial(:,6:8));
sTrial.Val = str2num(sInitial.Trial(:,10:11));

% Pupil size by ambig level
for i = 1: size(sTrial.AL,1)
    if sTrial.AL(i) == 0; sTrial.a(i) = 1;
    elseif sTrial.AL(i) == 24; sTrial.a(i) = 2;
    elseif sTrial.AL(i) == 50; sTrial.a(i) = 3;
    elseif sTrial.AL(i) == 74; sTrial.a(i) = 4;
    elseif sTrial.AL(i) == 100; sTrial.a(i) = 5;
    end
    
    if sTrial.Val(i) == 4; sTrial.b(i) = 1;
    elseif sTrial.Val(i) == 5; sTrial.b(i) = 2;
    elseif sTrial.Val(i) == 6; sTrial.b(i) = 3;
    elseif sTrial.Val(i) == 7; sTrial.b(i) = 4;
    elseif sTrial.Val(i) == 8; sTrial.b(i) = 5;
    elseif sTrial.Val(i) == 9; sTrial.b(i) = 6;
    elseif sTrial.Val(i) == 10; sTrial.b(i) = 7;
    elseif sTrial.Val(i) == 11; sTrial.b(i) = 8;
    elseif sTrial.Val(i) == 12; sTrial.b(i) = 9;
    elseif sTrial.Val(i) == 13; sTrial.b(i) = 10;
    elseif sTrial.Val(i) == 14; sTrial.b(i) = 11;
    elseif sTrial.Val(i) == 15; sTrial.b(i) = 12;
    elseif sTrial.Val(i) == 16; sTrial.b(i) = 13;
    elseif sTrial.Val(i) == 17; sTrial.b(i) = 14;
    elseif sTrial.Val(i) == 18; sTrial.b(i) = 15;
    elseif sTrial.Val(i) == 23; sTrial.b(i) = 16;
    elseif sTrial.Val(i) == 34; sTrial.b(i) = 17;
    elseif sTrial.Val(i) == 39; sTrial.b(i) = 18;
    elseif sTrial.Val(i) == 57; sTrial.b(i) = 19;
    elseif sTrial.Val(i) == 68; sTrial.b(i) = 20;
    end
end
pupilITI = zeros(5,20);
pupilLott = zeros(5,20);
pupilDelay = zeros(5,20);
for i = 1: size(sTrial.AL,1)
    pupilITI(sTrial.a(i),sTrial.b(i)) = sITI.PupilLeftMean(i);
    pupilLott(sTrial.a(i),sTrial.b(i)) = sLott.PupilLeftMean(i);
    pupilDelay(sTrial.a(i),sTrial.b(i)) = sDelay.PupilLeftMean(i);
    RT(sTrial.a(i),sTrial.b(i)) = sChoice.Rt(i);
end


%%  Devide the lott and delay period pupil size into ambiguity levels, compute the pupil mean at each time point, and the std at each time point

% Lottery presentation period
for i = 1:5
    for j = 1:20
        pupilLottNorm (i,j) = pupilLott(i,j)/pupilITI(i,j);
    end
end

for i = 1: size (sLott.PupilLeft_int,1)
    sLott.PupilLeft_intNorm(i,:) = sLott.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
end

temp1 = sLott.PupilLeft_intNorm;
a1 = find(sTrial.a == 1);
a2 = find(sTrial.a == 2);
a3 = find(sTrial.a == 3);
a4 = find(sTrial.a == 4);
a5 = find(sTrial.a == 5);

pupilLottAL(1,:) = nanmean(temp1(a1,:));
pupilLottAL(2,:) = nanmean(temp1(a2,:));
pupilLottAL(3,:) = nanmean(temp1(a3,:));
pupilLottAL(4,:) = nanmean(temp1(a4,:));
pupilLottAL(5,:) = nanmean(temp1(a5,:));

pupilLottALstd(1,:) = nanstd(temp1(a1,:));
pupilLottALstd(2,:) = nanstd(temp1(a2,:));
pupilLottALstd(3,:) = nanstd(temp1(a3,:));
pupilLottALstd(4,:) = nanstd(temp1(a4,:));
pupilLottALstd(5,:) = nanstd(temp1(a5,:));

ts_pupilLottAL = timeseries(pupilLottAL,sLott.Timestamp(1,:),'name', 'PupilSize_Lott_byAL');
figure('Name','Pupil_Lott_byAL','NumberTitle','off')
plot (ts_pupilLottAL)

% delay period
for i = 1:5
    for j = 1:20
        pupilDelayNorm (i,j) = pupilDelay(i,j)/pupilITI(i,j);
    end
end

for i = 1: size (sDelay.PupilLeft_int,1)
    sDelay.PupilLeft_intNorm(i,:) = sDelay.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
end

temp2 = sDelay.PupilLeft_intNorm;
a1 = find(sTrial.a == 1);
a2 = find(sTrial.a == 2);
a3 = find(sTrial.a == 3);
a4 = find(sTrial.a == 4);
a5 = find(sTrial.a == 5);

pupilDelayAL(1,:) = nanmean(temp2(a1,:));
pupilDelayAL(2,:) = nanmean(temp2(a2,:));
pupilDelayAL(3,:) = nanmean(temp2(a3,:));
pupilDelayAL(4,:) = nanmean(temp2(a4,:));
pupilDelayAL(5,:) = nanmean(temp2(a5,:));

pupilDelayALstd(1,:) = nanstd(temp2(a1,:));
pupilDelayALstd(2,:) = nanstd(temp2(a2,:));
pupilDelayALstd(3,:) = nanstd(temp2(a3,:));
pupilDelayALstd(4,:) = nanstd(temp2(a4,:));
pupilDelayALstd(5,:) = nanstd(temp2(a5,:));


ts_pupilDelayAL = timeseries(pupilDelayAL,sDelay.Timestamp(1,:),'name', 'PupilSize_Delay_byAL');
figure('Name','Pupil_Delay_byAL','NumberTitle','off')
plot (ts_pupilDelayAL)


