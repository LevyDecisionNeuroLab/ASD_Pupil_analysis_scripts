%% this script uses a sliding window, and a linear model, for each subject

 clearvars
 close all

%% load data
root = 'D:\Ruonan\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';
subj = [21 22 23 24 25 26 27 28 29 30 31 32 33 34];

datafold = fullfile(root,'Matlab data','pupildata\'); 
figfold = fullfile(root, 'Matlab data', 'pupil regression fig');

% get screen size for the convenience of plotting
screensize = get(groot, 'Screensize');

for subjidx = 1:length(subj)
    dataname = ['ASD' num2str(subj(subjidx)) '_Initial.mat'];
    load([datafold,dataname])

    % check the data quality: after preproc, how many trials remain
    trial2analyze = size(sInitial.PupilLeft_filtz,1);
    for i = 1:size(sInitial.PupilLeft_filtz,1)
        if sum(isnan(sInitial.PupilLeft_filtz(i,:)))== length(sInitial.PupilLeft_filtz(i,:))
            trial2analyze = trial2analyze - 1;
        end
    end

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
    % average = nanmean(sInitial.PupilLeft_filtz(:,:));
    % std = nanstd(sInitial.PupilLeft_filtz(:,:));
    % se = std ./ sqrt(size(sInitial.PupilLeft_filtz,1));
    % % a function from toolbox kakearney boudedline-pkg
    % figfiltzaveragese = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % boundedline(sInitial.Timestamp(1,:), average, se,'-b.')
    % % title(['ASD ' num2str(subj(subjidx)) ' averaged filtered z-scored timecourse with se'])
    % % figfiltzaveragesename = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtzAverageSe.fig' ]);
    % % saveas(figfiltzaveragese, figfiltzaveragesename)
    % 
    % figfiltzaveragestd = figure('Position', [screensize(3)/4 screensize(4)/4 screensize(3)/2 screensize(4)/2]);
    % boundedline(sInitial.Timestamp(1,:), average, std,'-r.')
    % title(['ASD ' num2str(subj(subjidx)) ' averaged filtered z-scored timecourse with std'])
    % figfiltzaveragestdname = fullfile(datafold,['ASD' num2str(subj(subjidx)), '_filtzAverageStd.fig' ]);
    % % saveas(figfiltzaveragestd, figfiltzaveragestdname)

    %% Linear model pupil ~ al + val, with a sliding window 
    x = [sInitial.AL sInitial.Val];
    windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
    signal = sInitial.PupilLeft_filt;
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

    %% Linear model pupil_velocity ~ al + val, with a sliding window 
    x = [sInitial.AL sInitial.Val];
    windl = 0; % actual length = (windl+1+windl) * 1000/60 ms
    signal = sInitial.Velleft_filt;
    timestamp = sInitial.Timestamp(1,:);
    % pre-allocate regression coefficient
    regcoeff = zeros(length(signal),2);
    % pre-allocate p-value
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
end





