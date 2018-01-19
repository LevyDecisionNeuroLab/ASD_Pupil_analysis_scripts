clearvars
close all
%% Load individual data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];
datafold = fullfile(root,'Matlab data','pupildata\');

for subj2plot = 1:length(subj)
    pupilLeft_filt = [];
    al = [];
    val = [];

    for subjidx = 1:length(subj)
        dataname = ['ASD' num2str(subj(subjidx)) '_Initial.mat']; 
        load([datafold,dataname])

        pupilLeft_filt = [pupilLeft_filt; sInitial.PupilLeft_filt];
        al = [al; sInitial.AL];
        val = [val; sInitial.Val];  
    end
    
    % normalize by subtractin the mean of the 2s ITI period
    pupilLeft_base = nanmean(pupilLeft_filt(:,1:120),2);
    pupilLeft_filt = pupilLeft_filt - pupilLeft_base;

    pupilLeft_filt = pupilLeft_filt(100*(subj2plot-1)+1:100*subj2plot,:);
    al = al(100*(subj2plot-1)+1:100*subj2plot,:);
    val = val(100*(subj2plot-1)+1:100*subj2plot,:);

    % pupil size by value levels
    uniqueval = unique(val);
    % exclude $4 and $5
    uniqueval = uniqueval(3:20);
    
    bins = 6; % how many groups to draw
    for i = 1:bins
        averageByVal(i,:) = nanmean(pupilLeft_filt(val <= uniqueval(i*18/bins),:));
    %     stdByVal(i,:) = nanstd(pupilLeft_filt(val <= uniqueval(i*4),:));
    %     seByVal(i,:) = std ./ sqrt(size(pupilLeft_filt,1));
    end

    % plot
    colors = {};
    for i = 1:bins
        colors{i} = [i/bins, 0, 0];
    end
    
    f = figure;
    for i = 1:bins
        plot(sInitial.Timestamp(1,:), averageByVal(i,:), 'LineStyle', '-', 'Marker', '.', 'Color',colors{bins+1-i})
        hold on
    end
    title(['ASD' num2str(subj(subj2plot)) ' pupil by value(18 values) ITI normalized']);


    saveas(f,fullfile(root, 'Matlab data', 'pupil trace by value',['ASD' num2str(subj(subj2plot)) ' pupil by value18_normalized.png']));
end
