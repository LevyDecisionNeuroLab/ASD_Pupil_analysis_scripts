% This script preprocess ASD pupil data, and then save it to the
% 'pupil_preprocessed' folder

%% History:
% Ruonan written 12.07.2017
% Ruonan updated 03.09.2018

clearvars
close all
%% Load individual data & specify preprocessing filter type
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
inputfold = fullfile(root,'Matlab data','pupildata_extracted\'); 
outputfold = fullfile(root,'Matlab data','pupildata_preprocessed\'); 
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

% Pupil data preprocessing filter type
filter = struct;
filter.filterType = 'sgolay';
filter.order = 3; % order of polynomial for sgolay filter?
filter.framelen = 21; % length of window? must be odd number
filter.clearWin = 2; % delete the n surrounding data points of a blink
filter.velThreshold = 2; % de-blinking relative velocity threshold
graph = false;

%% preprocessing
for subjidx = 1:length(subj)
    subjidx = 14;
    inputdataname = ['ASD' num2str(subj(subjidx)) '_Initial.mat'];
    load([inputfold,inputdataname])

    % dblink, interpolate and smooth data(preprocessing). 
    % Output: Filtered combine, filtered left, filtered right, filtered
    % combined z-scored, filtered combined velocity
    sInitial.filtered = struct;
    sInitial.filtered.Pupil=zeros(size(sInitial.PupilLeft));
    sInitial.filtered.PupilLeft=zeros(size(sInitial.PupilLeft));
    sInitial.filtered.PupilRight=zeros(size(sInitial.PupilLeft));
    sInitial.filtered.Pupilz=zeros(size(sInitial.PupilLeft));
    sInitial.filtered.vel=zeros(size(sInitial.PupilLeft));

    for i=1:size(sInitial.PupilLeft,1)
        [sInitial.filtered.Pupil(i,:),sInitial.filtered.PupilLeft(i,:),sInitial.filtered.PupilRight(i,:),sInitial.filtered.Pupilz(i,:),sInitial.filtered.vel(i,:)]=...
            combineLeftRight(sInitial.PupilLeft(i,:),sInitial.PupilRight(i,:),sInitial.Timestamp(i,:),filter,graph);
    end

    savefilename = fullfile(outputfold, ['ASD' num2str(subj(subjidx)) '_Initial_prepro.mat']);
    save(savefilename, 'sInitial')

    % % timeseries and figure
    % ts_initialLeft_int = timeseries(sInitial.PupilLeft_filt(:,:),sInitial.Timestamp(1,:),'name', 'PupilSize_Initial_int');
    % figure('Name','InitialResponseLeft','NumberTitle','off')
    % plot(ts_initialLeft_int)

    % % plot averaged time course with shaded error bar
    % average = nanmean(sInitial.PupilLeft_filt(:,:));
    % std = nanstd(sInitial.PupilLeft_filt(:,:));
    % se = std ./ sqrt(size(sInitial.PupilLeft_filt,1));
    % % a function from toolbox kakearney boudedline-pkg
    % figure
    % boundedline(sInitial.Timestamp(1,:), average, se,'-b.')
    % figure
    % boundedline(sInitial.Timestamp(1,:), average, std,'-r.')

    % plot time course, binned by values


    % plot(sInitial.Timestamp(1,:), average,'LineStyle', '-', 'Marker', '.');
    % hold on
    % errorbar(sInitial.Timestamp(1,:), average,std,'CapSize',0.1)
    end
