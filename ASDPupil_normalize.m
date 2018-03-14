% This script normalize the pupil size data by subtracing mean, after
% initial preprocessing

 clearvars
 close all

%% load data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata_preprocessed\'); 
figfold = fullfile(root, 'Matlab data', 'pupil average trace');
outputfold = fullfile(root,'Matlab data','pupildata_normalized\');

% get screen size for the convenience of plotting
screensize = get(groot, 'Screensize');

for subjidx = 1:length(subj)
%     subjidx = 1; % for testing, should be commented out
    dataname = ['ASD' num2str(subj(subjidx)) '_Initial_prepro.mat'];
    load([datafold,dataname])

    % baseline as the mean of the last 1000ms before stimulus
    baseline = nanmean(sInitial.filtered.Pupil(:,60:119),2);
    % if baseline of a trial is NaN, that means it has too many missing
    % values, so the trial should be excluded. And it will be done by the
    % next step, because if you subtract a NaN, the result will be NaN as
    % well
    sInitial.filtered.Pupil_norm = sInitial.filtered.Pupil - baseline;
    
%     plot(sInitial.Timestamp(1,:),sInitial.filtered.Pupil_norm);
%     plot(sInitial.Timestamp(1,:),sInitial.filtered.Pupil);
% 
%     plot(sInitial.Timestamp(1,:),nanmean(sInitial.filtered.Pupil_norm,1));
%     plot(sInitial.Timestamp(1,:),nanmean(sInitial.filtered.Pupil,1));
    
    % data qualitiy: proportion of missing data for each trial
    % first column-left, second column-right, third column-averaged
    signalQual = zeros(size(sInitial.PupilLeft,1),3);
    for i = 1:size(sInitial.PupilLeft,1)
        signalQual(i,1) = sum(isnan(sInitial.PupilLeft(i,:)))/length(sInitial.PupilLeft(i,:));
        signalQual(i,2) = sum(isnan(sInitial.PupilRight(i,:)))/length(sInitial.PupilRight(i,:));
        signalQual(i,3) = mean(signalQual(i,1:2));
    end
    
    sInitial.DataQual = signalQual;
    
    savefilename = fullfile(outputfold, ['ASD' num2str(subj(subjidx)) '_Initial_norm.mat']);
    save(savefilename, 'sInitial') 
    
    %% plot time course
    % % plot raw data
    % ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),sInitial.Timestamp(1,:),'name', ['ASD ' num2str(subj(subjidx)) ' PupilSize_Initial']);
    % % ts_initialLeft = timeseries(sInitial.PupilLeft(:,:),mean(sInitial.Timestamp(:,:),1));
    % 
    % figraw = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % plot(ts_initialLeft);
    % figrawname = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_raw.fig' ]);
    % % saveas(figraw, figrawname)
    % 
    % average = nanmean(sInitial.PupilLeft(:,:));
    % std = nanstd(sInitial.PupilLeft(:,:));
    % % a function from toolbox kakearney boudedline-pkg
    % figrawaverage = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % boundedline(sInitial.Timestamp(1,:), average, std,'-b.')
    % title(['ASD ' num2str(subj(subjidx)) ' averaged raw timecourse with std'])
    % figrawaveragename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_rawAverageStd.fig' ]);
    % % saveas(figrawaverage, figrawaveragename)
    % 
    % % timeseries and figure
    % ts_initialLeft_filt = timeseries(sInitial.PupilLeft_filt(:,:),sInitial.Timestamp(1,:),'name', ['ASD ' num2str(subj(subjidx)) ' PupilSize_Initial_filt']);
    % figfiltaverage = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % plot(ts_initialLeft_filt)
    % figfiltaveragename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtered.fig' ]);
    % % saveas(figfiltaverage, figfiltaveragename)
    % 
    % 
    % % plot averaged time course with shaded error bar
    signal2plot = sInitial.filtered.Pupil_norm;
    average = nanmean(signal2plot(:,:));
    std = nanstd(signal2plot(:,:));
    se = std ./ sqrt(size(signal2plot,1));
    % a function from toolbox kakearney boudedline-pkg
    figaveragese = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    boundedline(sInitial.Timestamp(1,:), average, se,'-b.')
    title(['ASD ' num2str(subj(subjidx)) ' averaged normazlied timecourse with se'])
    figaveragesename = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_NormazliedAverageSe.bmp' ]);
    saveas(figaveragese, figaveragesename)
    
    figaveragestd = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    boundedline(sInitial.Timestamp(1,:), average, std,'-r.')
    title(['ASD ' num2str(subj(subjidx)) ' averaged normalized timecourse with std'])
    figaveragestdname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_NormazliedAverageStd.bmp' ]);
    saveas(figaveragestd, figaveragestdname)

end