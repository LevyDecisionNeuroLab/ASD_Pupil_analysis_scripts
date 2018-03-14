function [sampfilt,sampfiltleft,sampfiltright,sampfiltz,velfilt] = combineLeftRight(sampleft,sampright,sampt,filter,graph)

% This function combines left and right eyes data, after the initial
% preprocessing of both eyes

%% For tsting the function, giving input to test the function. 
% Comment out this section when using the function
% 
% % close all
% 
% trialn = 1;
% sampleft = sInitial.PupilLeft(trialn,:);
% sampright= sInitial.PupilRight(trialn,:);
% sampt = sInitial.Timestamp(trialn,:);
% 
% filter.order = 3; % order of polynomial for sgolay filter?
% filter.framelen = 21; % length of windew? must be odd number
% filter.clearWin = 0; % delete the n surrounding data points of a blink
% filter.velThreshold = 2; % de-blinking velocity threshold
% filter.filterType = 'sgolay';
% % filter.filterType = 'hannWindow';
% graph = true;

%% Preprocessing of both eyes' pupil size
[sampfiltleft,sampfiltzleft,sampinterpleft,sampdbleft,velfiltleft,veldbleft,velleft] = pupilPrepro(sampleft,sampt,filter);
[sampfiltright,sampfiltzright,sampinterpright,sampdbright,velfiltright,veldbright,velright] = pupilPrepro(sampright,sampt,filter);

%% Combine left and right data and compute velocity profile
% averaging
sampfilt = zeros(size(sampfiltleft));
sampfilt = mean([sampfiltleft;sampfiltright]);

% z-score
% it is difficult to deal with NaN values with the function zscore
% sampz = zscore(sampfilt,0,2);
if any(isnan(sampfilt(:)))
    sampfiltmu=nanmean(sampfilt);
    sampfiltsigma=nanstd(sampfilt);
    sampfiltz=(sampfilt-repmat(sampfiltmu,1,length(sampfilt)))./repmat(sampfiltsigma,1,length(sampfilt));
else
    sampfiltz=zscore(sampfilt,0,2);
end

% velocity profile of the combined filtered timecourse
velfilt = zeros(size(sampfilt)); % velocity profile, mm/s
% calculate the velocity
for i=2:length(sampfilt)
    velfilt(i)= (sampfilt(i)-sampfilt(i-1))*1000/(sampt(i)-sampt(i-1));
end
velfilt(1)=velfilt(2);

%% Plot individual pupil time series results
if graph
    screensize = get(groot, 'Screensize');

    %% velocity profile
    figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12])
    % left side raw and deblinked
    ax1 = subplot(3,1,1);
    plot(ax1,sampt,velleft, 'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,veldbleft,'LineStyle', '-.', 'Marker', 'o')
    title(ax1,'left velocity')
    % right side raw and deblinked
    ax2 = subplot(3,1,2);
    plot(ax2,sampt,velright, 'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,veldbright,'LineStyle', '-.', 'Marker', 'o')
    title('right velocity')
    % filtered
    ax3 = subplot(3,1,3);
    plot(ax3,sampt,velfilt, 'LineStyle', '-.', 'Marker', 'o')
    title('combined filtered velocity')

    %% Plot pupil size time series
    figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12])
    % left, raw, deblinked, interpolated, filtered
    ax1 = subplot(3,1,1);
    plot(ax1,sampt,sampleft,'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,sampdbleft,'LineStyle', '-.', 'Marker', 'o')
    plot(sampt,sampinterpleft,'LineStyle', '-', 'Marker', 'o')
    plot(sampt,sampfiltleft,'LineStyle', '-', 'Marker', 'x')
    title(ax1,'left signal')  
    % ylim([3.2,4.4]);
    
    % right, raw, deblinked, interpolated, filtered
    ax2 = subplot(3,1,2);
    plot(ax2,sampt,sampright,'LineStyle', '-.', 'Marker', 'o')    
    hold on
    plot(sampt,sampdbright,'LineStyle', '-.', 'Marker', 'o')
    plot(sampt,sampinterpright,'LineStyle', '-', 'Marker', 'o')
    plot(sampt,sampfiltright,'LineStyle', '-', 'Marker', 'x')
    title(ax2,'right signal')
    xLimit = xlim;
    yLimit = ylim;
    if strcmp(filter.filterType,'sgolay')
        txtPar = ['sgolay: order = ', num2str(filter.order), '; framlen = ', num2str(filter.framelen)];
    elseif strcmp(filter.filterType, 'hannWindow')
        txtPar = ['hannWindow: framlen = ', num2str(filter.framelen)];
    end
    text(xLimit(2)*5/9,yLimit(1)+(yLimit(2)-yLimit(1))/8,txtPar)

    % combined filtered 
    ax3 = subplot(3,1,3);
    plot(ax3,sampt, sampfilt,'LineStyle', '-', 'Marker', '.')
    hold on
    title(ax3,'left-right combined filtered signal')
    
    %% plot and check how different the left and right pupil traces are
    figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12])
    % left filtered
    ax1 = subplot(3,1,1);
    plot(ax1,sampt,sampfiltleft, 'LineStyle', '-.', 'Marker', 'o')
    title(ax1,'left filtered')
    % right filtered
    ax2 = subplot(3,1,2);
    plot(ax2,sampt,sampfiltright, 'LineStyle', '-.', 'Marker', 'o')
    title('right filtered')
    % left and right filtered
    ax3 = subplot(3,1,3);
    plot(ax3,sampt,sampfiltleft, 'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,sampfiltright, 'LineStyle', '-.', 'Marker', 'o')
    plot(sampt, sampfilt,'LineStyle', '-', 'Marker', '.')
    title('left, right and combined filtered')

end
