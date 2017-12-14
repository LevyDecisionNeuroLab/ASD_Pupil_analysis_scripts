%% this script opens saved figures

 clearvars
 close all

 %% load data

% subjidx = 1;
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];
datafold = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData\Matlab data';
fig2open = '_regression.fig'; % name of the figures to open
% fig2open = '_raw.fig'; % name of the figures to open

for subjidx = 1:length(subj)
    figname = fullfile(datafold,['ASD',num2str(subj(subjidx)), fig2open]);
    openfig(figname);
end
