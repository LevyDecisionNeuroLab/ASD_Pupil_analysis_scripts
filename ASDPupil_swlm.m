%% this script uses a sliding window, and a linear model pupil ~ val+al, for each subject

 clearvars
 close all
% clc

%% load data
subjidx=14;
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];
datafold = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData\Matlab data';
dataname = ['\ASD',num2str(subj(subjidx)), '_origData.mat']
load([datafold,dataname])

%% Pupil data preprocessing filter type
filter = struct;
filter.filterType = 'sgolay';
filter.order = 3; % order of polynomial for sgolay filter?
filter.framelen = 21; % length of window? must be odd number
filter.clearWin = 1; % delete the n surrounding data points of a blink
filter.velThreshold = 3; % de-blinking relative velocity threshold

screensize = get(groot, 'Screensize');
%% Extract pupil data from different time period
%% 1. Initial pupil response
% Extract the 2sITI + 3sLottery presentation + 5s Delay
sInitial = struct; % sInitial stores pupil size for 10s time line
sInitial.subj = subj(subjidx);
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

% reset timestamp to start from 0
for k = 1:size(sInitial.Timestamp,1)
    sInitial.Timestamp(k,:) = sInitial.Timestamp(k,:)-sInitial.Timestamp(k,1);
end

% dblink, interpolate and smooth data(preprocessing). Filtered data signal, and filted z-scored data signal
sInitial.PupilLeft_filt=zeros(size(sInitial.PupilLeft));
for i=1:size(sInitial.PupilLeft,1)
    [sInitial.PupilLeft_filt(i,:),sInitial.PupilLeft_filtz(i,:)] = pupilPrepro(sInitial.PupilLeft(i,:), sInitial.Timestamp(i,:), filter);
end

% check the data quality: after preproc, how many trials remain
trial2analyze = size(sInitial.PupilLeft_filtz,1);
for i = 1:size(sInitial.PupilLeft_filtz,1)
    if sum(isnan(sInitial.PupilLeft_filtz(i,:)))== length(sInitial.PupilLeft_filtz(i,:))
        trial2analyze = trial2analyze - 1;
    end
end

%% plot time course
% plot raw data
ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),sInitial.Timestamp(1,:),'name', ['ASD ' num2str(subj(subjidx)) ' PupilSize_Initial']);
% ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),mean(sInitial.Timestamp(:,:),1));

figraw = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
plot(ts_initialLeft);
figrawname = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_raw.fig' ]);
saveas(figraw, figrawname)

average = nanmean(sInitial.PupilLeft(:,:));
std = nanstd(sInitial.PupilLeft(:,:));
% a function from toolbox kakearney boudedline-pkg
figrawaverage = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
boundedline(sInitial.Timestamp(1,:), average, std,'-b.')
title(['ASD ' num2str(subj(subjidx)) ' averaged raw timecourse with std'])
figrawaveragename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_rawAverageStd.fig' ]);
saveas(figrawaverage, figrawaveragename)

% timeseries and figure
ts_initialLeft_filt = timeseries(sInitial.PupilLeft_filt(:,:),sInitial.Timestamp(1,:),'name', ['ASD ' num2str(subj(subjidx)) ' PupilSize_Initial_filt']);
figfiltaverage = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
plot(ts_initialLeft_filt)
figfiltaveragename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtered.fig' ]);
saveas(figfiltaverage, figfiltaveragename)


% plot averaged time course with shaded error bar
average = nanmean(sInitial.PupilLeft_filtz(:,:));
std = nanstd(sInitial.PupilLeft_filtz(:,:));
se = std ./ sqrt(size(sInitial.PupilLeft_filtz,1));
% a function from toolbox kakearney boudedline-pkg
figfiltzaveragese = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
boundedline(sInitial.Timestamp(1,:), average, se,'-b.')
% title(['ASD ' num2str(subj(subjidx)) ' averaged filtered z-scored timecourse with se'])
% figfiltzaveragesename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtzAverageSe.fig' ]);
% saveas(figfiltzaveragese, figfiltzaveragesename)

figfiltzaveragestd = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
boundedline(sInitial.Timestamp(1,:), average, std,'-r.')
title(['ASD ' num2str(subj(subjidx)) ' averaged filtered z-scored timecourse with std'])
figfiltzaveragestdname = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtzAverageStd.fig' ]);
saveas(figfiltzaveragestd, figfiltzaveragestdname)


%% Trial information, summary pupil size into ambig and val levels

% ambig level and value for each trial
sInitial.al = str2num(sInitial.Trial(:,6:8));
sInitial.val = str2num(sInitial.Trial(:,10:11));

%% Linear model with a sliding window 
x = [sInitial.al sInitial.val];
windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
signal = sInitial.PupilLeft_filtz;
timestamp = sInitial.Timestamp(1,:);
% regression coefficient
regcoeff = zeros(length(signal),2);
% p-value
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

%% plot regression
valsignif = find(pval(:,2)<0.05 & pval (:,2)>0);
alsignif = find(pval(:,1)<0.05 & pval (:,1)>0);

figRegress = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
% plot regression coefficient
% ambiguity level, blue
plot([1:length(regcoeff)]*1000/60,regcoeff(:,1),'Color','b')
hold on
% value, red
plot([1:length(regcoeff)]*1000/60,regcoeff(:,2),'Color','r')
yLimit = ylim;
xLimit = xlim;

% paint the windows with significance (uncorrected)
% ambiguity level, blue
if isempty(alsignif)==0
    color = [0 0 1];
    drawy = linspace(yLimit(1),yLimit(2));
    for k = 1:length(alsignif)
        drawx = repmat(alsignif(k)*1000/60,length(drawy),2);
        plot(drawx,drawy,'Color',[0 0 1 0.2])
    end

%     aal=area(alsignif,repmat(yaxis(1),length(alsignif),1),'basevalue',yaxis(2),...
%         'FaceColor',color,'EdgeColor','none','ShowBaseLine','off');
%     % make area color transparent
%     drawnow; pause(0.05);  % This needs to be done for transparency to work
%     aal.Face.ColorType = 'truecoloralpha';
%     aal.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
else
    disp('No siginificant time window for al')
end

% value, red
if isempty(valsignif) ==0
    color = [1 0 0];
    drawy = linspace(yLimit(1),yLimit(2));
    for k = 1:length(valsignif)
        drawx = repmat(valsignif(k)*1000/60,length(drawy),2);
        plot(drawx,drawy,'Color',[1 0 0 0.2])
    end
%     aval=area(valsignif,repmat(yaxis(1),length(valsignif),1),'basevalue',yaxis(2),...
%         'FaceColor',color,'EdgeColor','none','ShowBaseLine','off');
%     % make area color transparent
%     drawnow; pause(0.05);  % This needs to be done for transparency to work
%     aval.Face.ColorType = 'truecoloralpha';
%     aval.Face.ColorData(4) = 255 * 0.3; % Your alpha value is the 0.3
else
    disp('No siginificant time window for val')
end

txtPar1 = ['window length = ' num2str(windl*2+1) '; included trials = ' num2str(trial2analyze)];
text(xLimit(2)*5/9,yLimit(1)+(yLimit(2)-yLimit(1))/8,txtPar1)

title(['ASD ' num2str(subj(subjidx)), ' regression Pupil ~ val(red) +al(blue)'])

figname = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_regression.fig' ]);
saveas (figRegress, figname);

%% Save the sInitial data structure
sInitialName = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_procData.mat' ]);
save(sInitialName, 'sInitial')

%% 2. Average ITI Pupil size last 2s
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
% % normalize the whole time series by ITI mean
% for i = 1: size (sInitial.PupilLeft_int,1)
%     sInitial.PupilLeft_intNorm(i,:) = sInitial.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
% end







