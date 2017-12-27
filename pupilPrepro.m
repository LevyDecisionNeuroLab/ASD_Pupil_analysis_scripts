function [sampfilt,sampfiltz] = pupilPrepro(samp,sampt,filter)

%% This function preprocesses the raw pupil data
%
% Input:
%   - samp: pupil size data
%   - sampt: pupil size data time point
%   - filter: filter fype and parameters
%       - filterType: string
%       - order: polynomial order for sgolay filter
%       - framelen: window size
%       - clearWin: number of surrounding data points to delete on both
%         left and right sides. actually deleting (2*clearWin + 1) data points
%       - velThreshold: how many std is compared with the velocity difference from mean
% Output:
%   - sampfilt: filtered pupil size data
%   - sampfiltz: z-scored filtered pupil size data
%   - velfilt: velocity of z-scored filtered pupil size data. TO BE ADDED

%% Preprocessing steps
% 
% 1. Check data quality:
%       If over 40% data points are missing, delete the trial completely. 
% 2. Deblink: 
%       First calculate the velocity of each time point. If
%       the difference between velocity at time t and the mean velocity is larger
%       than the filter.velThreshold*velocity standard deviation, delete this
%       data point at time t plus the surrounding filter.clearWin data points on both left
%       and right side of time t.
% 3. Interpolate: 
%       linear interpolation of missing data points
% 4. Filter: 
%       filter the interpolated time course depending on the filter type
%       and filter set-up parameters, which typically includes window length
% 5. Z-score: 
%       z-score the filtered time course
% 6. Plot: 
%       plot the velocity profile. 
%       plot the time course in three panels:
%       upper: raw time course; 
%       middle: raw time course, interpolated time course, filtered time course overlaid; 
%       lower: z-scored filtered time course.

%% Histroy
% Ruonan Jia 12.07.2017 written

%% For tsting the function, giving input to test the function
% samp = sInitial.PupilLeft(1,:);
% sampt = sInitial.Timestamp(1,:);
% % filter.order = 3; % order of polynomial for sgolay filter?
% filter.framelen = 21; % length of windew? must be odd number
% filter.clearWin = 2; % delete the n surrounding data points of a blink
% filter.velThreshold = 2; % de-blinking velocity threshold
% filter.filterType = 'sgolay';
% % filter.filterType = 'hannWindow';

% decide if missing too much data
dataQual = 1-sum(isnan(samp))/length(samp);
if dataQual >0.6
    %% De-blink
    % Detect blink
    vel = zeros(size(samp)); % velocity profile, mm/s
    % calculate the velocity
    for i=2:length(samp)
        vel(i)= (samp(i)-samp(i-1))*1000/(sampt(i)-sampt(i-1));
    end
    vel(1)=vel(2);

    veldb=vel;
    sampdb = samp;
    
    velStd = nanstd(vel);
    velMean = nanmean(vel);
    
    % detect fast speed, and delete the surrounding data points
    for i = 1:length(vel)
        if abs(vel(i)-velMean)/velStd > filter.velThreshold || isnan(vel(i))
            if i-filter.clearWin > 1 && i+filter.clearWin < length(vel)
                sampdb(i-filter.clearWin:i+filter.clearWin)=NaN; 
                veldb(i-filter.clearWin:i+filter.clearWin)=NaN; 
            elseif i-filter.clearWin < 1
                sampdb(1:i+filter.clearWin)=NaN; 
                veldb(1:i+filter.clearWin)=NaN; 
            elseif i+filter.clearWin > length(vel)
                sampdb(i-filter.clearWin:length(vel))=NaN; 
                veldb(i-filter.clearWin:length(vel))=NaN; 
            end
        end
    end
    
    %% Interpolate missing data
    miss = find (isnan(sampdb)); % missing data
    exist = find (isnan(sampdb)==0); % existing data
    % whole = [exist(1):length(sampfil)]; % the whole time series,excluding the first NaNs

    if length(exist) > 0
        % % Cubic Spline interpolation
        % interp = spline(sampt(exist),sampfil(exist),sampt(miss));
        % sampinterp=sampfil;
        % for i=1:length(miss)
        %     if miss(i)<exist(1); continue; end
        %     if miss(i)>exist(length(exist)); break; end
        %     sampinterp(miss(i))=interp(i);
        % end

        % Linear interpolation
        interp = interp1(sampt(exist),sampdb(exist),sampt(miss));
        sampinterp=sampdb;
        for i=1:length(miss)
            if miss(i)<exist(1); continue; end
            if miss(i)>exist(length(exist)); break; end
            sampinterp(miss(i))=interp(i);
        end
        % for i=1:length(miss)
        %     sampinterp(miss(i))=interp(i);
        % end
    else % if after filtering, all data are missing
         sampinterp=sampdb;
    end

    %% Filtering
    if strcmp(filter.filterType,'sgolay')
        % savitzky-golay smoothing
        % what order and framelen to use?
        sampfilt = sgolayfilt(sampinterp, filter.order, filter.framelen);
    elseif strcmp(filter.filterType, 'hannWindow')
        % filtering using a hanning window average
        hwin = hann(filter.framelen);
        hwin = hwin/sum(hwin);
        sampfilt = conv(sampinterp,hwin,'same');
    end

    sampfilt(1)=sampfilt(2);
    sampfilt(length(sampfilt))=sampfilt(length(sampfilt)-1);

    %% z-score
    % it is difficult to deal with NaN values with the function zscore
    % sampz = zscore(sampfilt,0,2);

    if any(isnan(sampfilt(:)))
        sampfiltmu=nanmean(sampfilt);
        sampfiltsigma=nanstd(sampfilt);
        sampfiltz=(sampfilt-repmat(sampfiltmu,1,length(sampfilt)))./repmat(sampfiltsigma,1,length(sampfilt));
    else
        sampfiltz=zscore(sampfilt,0,2);
    end

    %% Plot individual pupil time series results
    % velocity profile
    figure
    plot(sampt,vel, 'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,veldb,'LineStyle', '-.', 'Marker', 'o')
    
    % pupil size time series
    screensize = get(groot, 'Screensize');
    figure('Position', [screensize(3)/4 screensize(4)/22 screensize(3)/2 screensize(4)*10/12])
    ax1 = subplot(3,1,1);
    plot(ax1,sampt,samp,'LineStyle', '-.', 'Marker', 'o')
    title(ax1,'raw signal')
    
    % ylim([3.2,4.4]);
    ax2 = subplot(3,1,2);
    plot(ax2,sampt,samp,'LineStyle', '-.', 'Marker', 'o')
    hold on
    plot(sampt,sampdb,'LineStyle', '-.', 'Marker', 'o')
    plot(sampt,sampinterp,'LineStyle', '-', 'Marker', 'o')
    plot(sampt,sampfilt,'LineStyle', '-', 'Marker', 'x')
    title(ax2,'filtered signal')
    xLimit = xlim;
    yLimit = ylim;
    if strcmp(filter.filterType,'sgolay')
        txtPar = ['sgolay: order = ', num2str(filter.order), '; framlen = ', num2str(filter.framelen)];
    elseif strcmp(filter.filterType, 'hannWindow')
        txtPar = ['hannWindow: framlen = ', num2str(filter.framelen)];
    end
    text(xLimit(2)*5/9,yLimit(1)+(yLimit(2)-yLimit(1))/8,txtPar)
    
    % print z-scored data
    ax3 = subplot(3,1,3);
    plot(ax3,sampt, sampfiltz,'LineStyle', '-', 'Marker', '.')
    hold on
    title(ax3,'z-scored filtered signal')
    xLimit = xlim;
    yLimit = ylim;
    text(xLimit(2)*5/9,yLimit(1)+(yLimit(2)-yLimit(1))/8,'z-scored')
else
    % if data is bad, discard
    sampfilt = ones(size(samp))*NaN;
    sampfiltz = ones(size(samp))*NaN;
end
end