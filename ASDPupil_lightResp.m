%% this script looks at pupil light response

 clearvars
 close all

%% load data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata_normalized\'); 
figfold = fullfile(root, 'Matlab data', 'pupil average trace');

% load behavioral model fitting result
behavFilename = fullfile(root,'ASDBehavModelFit.txt');
modelFit = tdfread(behavFilename);

% get screen size for the convenience of plotting
screensize = get(groot, 'Screensize');

outputdata = struct;
outputdata.subjId = zeros(100*length(subj),1);
outputdata.pupilPeakDil = zeros(100*length(subj),1);
outputdata.al = zeros(100*length(subj),1);
outputdata.val = zeros(100*length(subj),1);
outputdata.subjValLott = zeros(100*length(subj),1);
outputdata.subjValChosen = zeros(100*length(subj),1);
outputdata.alpha = zeros(100*length(subj),1);
outputdata.beta = zeros(100*length(subj),1);

for subjidx = 1:length(subj)
    %% load data
    
%     subjidx = 1; % for testing, should be commented out

    dataname = ['ASD' num2str(subj(subjidx)) '_Initial_norm.mat'];
    load([datafold,dataname])
   
    %% behavioral model fitting result
    alpha = modelFit.alpha(modelFit.SubjID == subj(subjidx));
    beta = modelFit.beta(modelFit.SubjID == subj(subjidx));
    gamma = modelFit.slopeP(modelFit.SubjID == subj(subjidx));
    
    %% calculate subjective value
    % calculate subjective value of lottery
    subjValLott = ambig_utility(0,sInitial.Val,repmat(0.5,length(sInitial.Val),1),sInitial.AL,...
        alpha, beta,'ambigNrisk');

    % calculate subjective value of lottery    
    subjValRef = ambig_utility(0,repmat(5,length(sInitial.Val),1),repmat(0.5,length(sInitial.Val),1),...
        zeros(length(sInitial.Val),1),alpha, beta,'ambigNrisk');
    
    % subjective value of the chose option
    subjValChosen = subjValLott .* sInitial.Choice + subjValRef .* (1-sInitial.Choice);    

    %% clean signal and excluding bad trials    
    % Choose the signal to look at
    pupil = sInitial.filtered.Pupil_norm;
    vel = sInitial.filtered.vel;

    % data qualitiy: proportion of missing data for each trial
    % first column-left, second column-right, third column-averaged
    signalQual = sInitial.DataQual;
    
    % check the data quality: after preproc, how many trials remain
    % standard for excluding a trial: in the chosen time window, if over 18% data is missing 
    threshold = 0.3;
    trial2analyze = sum(signalQual(:,3) < threshold);
    pupil(signalQual(:,3) > threshold,:) = nan;
    vel(signalQual(:,3) > threshold,:) = nan;    

    %% plot averaged velocity and pupil with shaded error bar
    
%     figaveragestd = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
%     
%     signal2plot = pupil;
%     average = nanmean(signal2plot(:,:));
%     std = nanstd(signal2plot(:,:));
%     se = std ./ sqrt(size(signal2plot,1));
%     ax1 = subplot(2,1,1);
%     boundedline(ax1,sInitial.Timestamp(1,:), average, std,'-r.')
%     title(['ASD ' num2str(subj(subjidx)) ' averaged pupil size with std'])
%     
%     signal2plot = vel;
%     average = nanmean(signal2plot(:,:));
%     std = nanstd(signal2plot(:,:));
%     se = std ./ sqrt(size(signal2plot,1));
%     ax2 = subplot(2,1,2);
%     boundedline(ax2,sInitial.Timestamp(1,:), average, std,'-b.')
%     title(['ASD ' num2str(subj(subjidx)) ' averaged velocity with std'])
%     
%     figname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_PupilNVel.bmp' ]);
%     saveas(figaveragestd, figname)
    
    %% calculate peak change in pupil size's light response
    % the peak of the 1st light resposne(seeing lottery) should be a little
    % before 3000ms according to the velocity graph.
    
    % find in the averaged trace, the lowest point which should be the peak
    pupilAver = nanmean(pupil);
    peak = 120 + find(pupilAver(121:210) == min(pupilAver(121:210))); % lowest point pupil trace during 2-3.5s
%     peaktime = peak/60
    
    pupilPeakDil = nanmean(pupil(:,peak-5:peak+5),2)-nanmean(pupil(:,95:119),2); % peak - 0.4s baseline
    
    %% organize data into the outputdata struct, standardize within individual
    outputdata.subjId((subjidx-1)*100+1:subjidx*100) = repmat(subj(subjidx),100,1);
    
    pupilPeakDilmu = nanmean(pupilPeakDil);
    pupilPeakDilstd = nanstd(pupilPeakDil);
    pupilPeakDilz = (pupilPeakDil - pupilPeakDilmu)/pupilPeakDilstd;
    outputdata.pupilPeakDil((subjidx-1)*100+1:subjidx*100) = pupilPeakDilz;
    
    outputdata.al((subjidx-1)*100+1:subjidx*100) = zscore(sInitial.AL);
    outputdata.val((subjidx-1)*100+1:subjidx*100) = zscore(sInitial.Val);
    outputdata.subjValLott((subjidx-1)*100+1:subjidx*100) = zscore(subjValLott);
    outputdata.subjValChosen((subjidx-1)*100+1:subjidx*100) = zscore(subjValChosen);
    outputdata.alpha((subjidx-1)*100+1:subjidx*100) = repmat(alpha,100,1);
    outputdata.beta((subjidx-1)*100+1:subjidx*100) = repmat(beta,100,1);
end

%% linear regression

% create table
mydata = table(outputdata.subjId,outputdata.al,outputdata.val,outputdata.subjValLott,...
    outputdata.subjValChosen,outputdata.alpha,outputdata.beta,outputdata.pupilPeakDil,...
    'VariableNames',{'subjId','al','val','subjValLott','subjValChosen','alpha','beta','pupilPeakDil'});

% mixed effect model
lme1 = fitlme(mydata, 'pupilPeakDil~al+val+(1|subjId)')

lme2 = fitlme(mydata, 'pupilPeakDil~subjValChosen+(1|subjId)')
lme3 = fitlme(mydata, 'pupilPeakDil~subjValLott+(1|subjId)')


%% local function to calculate subjective value
function y = ambig_utility(base,v,p,AL,alpha,beta,model);

if (strcmp(model,'ambiguity') || strcmp(model,'ambigNrisk')) || strcmp(model,'ambigNriskFixSlope')
    % the model we are using
    y = (p - beta .* (AL./2)) .* v .^alpha + (1-p - beta .* (AL./2)) .* base .^alpha;
elseif strcmp(model,'ambigPower')
    y = p .^ (1+beta.*AL) .* v .^alpha; % change that
elseif strcmp(model,'discounting')
    %y = v ./ (1 + alpha.*log(1+(1-p+beta.*AL./2)./(p-beta.*AL./2)));
    y = v ./ (1 + alpha.*(1-p+beta.*AL./2)./(p-beta.*AL./2));
    %y = v ./ (1 + alpha.*(1-p)./p);
end
end