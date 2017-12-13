%% this script uses a sliding window, and a linear model pupil ~ val+al, for each subject

 clear
% clc
% cd 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
% s = tdfread('ASD22.txt');

%% Extract pupil data from different time period

%1. Initial pupil response
% Extract the 2sITI + 3sLottery presentation + 5s Delay
i=12;
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];
datafold = 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData\Matlab data\Matlab extracted original data';
dataname = ['\ASD',num2str(subj(i)), '_origData.mat']
load([datafold,dataname])
sInitial = struct; % sInitial stores pupil size for 8s time line
i = 0; % counting index for going over each element in s in sequence
j = 0; % No. for trial 1~25
% sInitial.Trial = repmat(char(0),25,26);
while i < length(s.Timestamp)
    i=i+1;
    if s.EventKey(i) ~= 8 | s.Descriptor(i,1)~= 'a'; continue; end
    if s.EventKey(i) == 8 & s.Descriptor(i,1) == 'a' 
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


% dblink, smooth and interpolate data(preprocessing)
sInitial.PupilLeft_int=zeros(size(sInitial.PupilLeft));
for i=1:100
    sInitial.PupilLeft_int(i,:) = blinkinterp(sInitial.PupilLeft(i,:), sInitial.Timestamp(i,:));
end
% timeseries and figure
ts_initialLeft_int = timeseries(sInitial.PupilLeft_int(:,:),sInitial.Timestamp(1,:),'name', 'PupilSize_Initial_int');
figure('Name','InitialResponseLeft','NumberTitle','off')
% hold on
plot (ts_initialLeft_int)

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

% normalize the whole time series by ITI mean
for i = 1: size (sInitial.PupilLeft_int,1)
    sInitial.PupilLeft_intNorm(i,:) = sInitial.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
end

%% Trial information, summary pupil size into ambig and val levels

% ambig level and value for each trial
sTrial = struct;
sTrial.AL = str2num(sInitial.Trial(:,6:8));
sTrial.Val = str2num(sInitial.Trial(:,10:11));

%% Linear model with a sliding window 
x = [sTrial.AL sTrial.Val];
windl = 6; % actual length = 6+1+6 = 13 = 218 ms
signal = sInitial.PupilLeft_intNorm;
regcoeff = zeros(length(signal),2);
pval = zeros(length(signal),2);
for i = windl+1 : length(signal)-windl
    y = nanmean(signal(:,i-windl:i+windl),2); % average of the time window
    mdl = LinearModel.fit(x,y);
    coeff = table2array(mdl.Coefficients);
    regcoeff(i,1) = coeff(2,1); % regression coefficient for AL
    regcoeff(i,2) = coeff(3,1); % regression coefficient for Val
    pval(i,1) = coeff(2,4); % p value for AL
    pval(i,2) = coeff(3,4); % p value for Val
end

% plot 
valsignif = find(pval(:,2)<0.05 & pval (:,2)>0);
alsignif = find(pval(:,1)<0.05 & pval (:,1)>0);

figure
plot(1:length(regcoeff),regcoeff)
hold on
yaxis = ylim;
if isempty(valsignif) ==0
    color = [1 0 0];
    aval=area(valsignif,repmat(yaxis(1),length(valsignif),1),'basevalue',yaxis(2),...
        'FaceColor',color,'EdgeColor','none','ShowBaseLine','off');
    % make area color transparent
    drawnow; pause(0.05);  % This needs to be done for transparency to work
    aval.Face.ColorType = 'truecoloralpha';
    aval.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
else
    disp('No siginificant time window for val')
end
hold on

if isempty(alsignif)==0
    color = [0 0 1];
    aal=area(alsignif,repmat(yaxis(1),length(alsignif),1),'basevalue',yaxis(2),...
        'FaceColor',color,'EdgeColor','none','ShowBaseLine','off');
    % make area color transparent
    drawnow; pause(0.05);  % This needs to be done for transparency to work
    aal.Face.ColorType = 'truecoloralpha';
    aal.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
else
    disp('No siginificant time window for al')
end




