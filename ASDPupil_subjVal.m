%% this script analyze pupil size and subjective value

 clearvars
 close all

%% load data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata_normalized\'); 
figfold = fullfile(root, 'Matlab data', 'pupil regression fig','SubjVal_normalized');

% load behavioral model fitting result
behavFilename = fullfile(root,'ASDBehavModelFit.txt');
modelFit = tdfread(behavFilename);

% get screen size for the convenience of plotting
screensize = get(groot, 'Screensize');

for subjidx = 1:length(subj)
    %% load data
    
%     subjidx = 1; % for testing, should be commented out

    dataname = ['ASD' num2str(subj(subjidx)) '_Initial_norm.mat'];
    load([datafold,dataname])
    
    % behavioral model fitting result
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

    %% 
    
    
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

    %% sliding window regression, pupil ~ subjVal, Regression uses standardized data
  
    x = [zscore(subjValChosen)];
    windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
    timestamp = sInitial.Timestamp(1,:);
    % regression coefficient
    regcoeff = zeros(length(pupil),1);
    % p-value
    pval = zeros(length(pupil),1);

    for i = windl+1 : length(pupil)-windl
        y = nanmean(pupil(:,i-windl:i+windl),2); % average of the time window
        % z-score y
        ymu = nanmean(y);
        ysigma = nanstd(y);
        yz=(y-repmat(ymu,length(y),1))./repmat(ysigma,length(y),1);
        % regression using standardized data
        mdl = LinearModel.fit(x,yz);
        regcoeff(i,1) = mdl.Coefficients.Estimate(2); % regression coefficient for subjective value
        pval(i,1) = mdl.Coefficients.pValue(2); % p value for subjective value
    end

    % plot regression
    subjValLottsignif = find(pval(:,1)<0.05 & pval (:,1)>0);    
    
    figRegress = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % plot regression coefficient
    % subjective value, blue
    plot([1:length(regcoeff)]*1000/60,regcoeff(:,1),'Color','b')
    hold on

    yLimit = ylim;
    xLimit = xlim;
    % paint the windows with significance (uncorrected)
    % ambiguity level, blue
    if isempty(subjValLottsignif)==0
        color = [0 0 1];
        drawy = linspace(yLimit(1),yLimit(2));
        for k = 1:length(subjValLottsignif)
            drawx = repmat(subjValLottsignif(k)*1000/60,length(drawy),2);
            plot(drawx,drawy,'Color',[0 0 1 0.2])
        end
    else
        disp('No siginificant time window for subjValChosen')
    end


    txtPar1 = ['window length = ' num2str(windl*2+1) '; included trials = ' num2str(trial2analyze)];
    text(xLimit(2)*5/9,yLimit(1)+(yLimit(2)-yLimit(1))/8,txtPar1)

    title(['ASD ' num2str(subj(subjidx)), ' regression Pupil ~ subjective value chosen (blue)'])

%     figname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_pupilsize.fig' ]);
%     saveas (figRegress, figname);    
    bmpname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_SubjValChosen.bmp' ]);
    saveas (figRegress, bmpname);        
end

% local function to calculate subjective value
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