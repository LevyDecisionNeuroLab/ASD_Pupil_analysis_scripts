%% this script uses a sliding window, and a linear model, for each subject

 clearvars
 close all

%% load data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata_normalized\'); 
figfold = fullfile(root, 'Matlab data', 'pupil regression fig','Left right eye average_normalized');

% get screen size for the convenience of plotting
screensize = get(groot, 'Screensize');

for subjidx = 1:length(subj)
    %% load data and exclude bad trials
    
%     subjidx = 1; % for testing, should be commented out

    dataname = ['ASD' num2str(subj(subjidx)) '_Initial_norm.mat'];
    load([datafold,dataname])

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


    %% Linear model pupil ~ al + val, with a sliding window 
    x = [sInitial.AL sInitial.Val];
    windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
    timestamp = sInitial.Timestamp(1,:);
    % regression coefficient
    regcoeff = zeros(length(pupil),2);
    % p-value
    pval = zeros(length(pupil),2);

    for i = windl+1 : length(pupil)-windl
        y = nanmean(pupil(:,i-windl:i+windl),2); % average of the time window
        mdl = LinearModel.fit(x,y);
        coeff = table2array(mdl.Coefficients);
        regcoeff(i,1) = coeff(2,1); % regression coefficient for AL
        regcoeff(i,2) = coeff(3,1); % regression coefficient for Val
        pval(i,1) = coeff(2,4); % p value for AL
        pval(i,2) = coeff(3,4); % p value for Val
    end

    % plot regression
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

    figname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_pupilsize.fig' ]);
    saveas (figRegress, figname);    
    bmpname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_pupilsize.bmp' ]);
    saveas (figRegress, bmpname);    

    %% Linear model pupil_velocity ~ al + val, with a sliding window 
    x = [sInitial.AL sInitial.Val];
    windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
    timestamp = sInitial.Timestamp(1,:);
    % pre-allocate regression coefficient
    regcoeff = zeros(length(vel),2);
    % pre-allocate p-value
    pval = zeros(length(vel),2);

    for i = windl+1 : length(vel)-windl
        y = nanmean(vel(:,i-windl:i+windl),2); % average of the time window
        mdl = LinearModel.fit(x,y);
        coeff = table2array(mdl.Coefficients);
        regcoeff(i,1) = coeff(2,1); % regression coefficient for AL
        regcoeff(i,2) = coeff(3,1); % regression coefficient for Val
        pval(i,1) = coeff(2,4); % p value for AL
        pval(i,2) = coeff(3,4); % p value for Val
    end

    % plot regression
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

    title(['ASD ' num2str(subj(subjidx)), ' regression PupilVelocity ~ val(red) +al(blue)'])

    figname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_velocity.fig' ]);
    saveas (figRegress, figname);
    bmpname = fullfile(figfold,['ASD' num2str(subj(subjidx)), '_regression_velocity.bmp' ]);
    saveas (figRegress, bmpname);
    
    
end





