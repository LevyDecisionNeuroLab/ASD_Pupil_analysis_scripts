function [sampfilt,sampfiltz,sampinterp,sampdb,velfilt,veldb,vel] = pupilPrepro(samp,sampt,filter)

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
%       - graph: true or fales, whether to output the figures of individual trial
% Output:
%   - sampfilt: filtered pupil size data
%   - sampfiltz: z-scored filtered pupil size data
%   - velfilt: velocity of z-scored filtered pupil size data.

%% Preprocessing steps
% 
% 1. Check data quality:
%       ??If over 40% data points are missing, delete the trial completely. 
% 2. Deblink: 
%       First calculate the velocity of each time point. If
%       the difference between velocity at time t and the mean velocity is larger
%       than the filter.velThreshold*velocity standard deviation, delete this
%       data point at time t plus the surrounding filter.clearWin data points on both left
%       and right side of time t.
% 
%    The following step was tried but was not adopeted eventually. Because
%    the left and right pupil sizes are actually not the same:
%           Combining binocular data. After de-blinking, average the left and right data if both data 
%           exist; or use one side of data if the other side is missing; or leave as blank if both 
%           sides data are missing. Assuming that pupil size are same of both sides of the eyes 
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
% 
% * Not done in the function:
% 1. Baseline correction by subtraction
% 2. convert to percent signal change


%% Histroy
% Ruonan Jia 12.07.2017 written


% decide if missing too much data
% it depends on the analysis window
% dataQual = 1-sum(isnan(samp))/length(samp);
% if dataQual >0.6

    %% De-blink and interpolate missing data
    % Detect blink both left and right pupil data
    [vel,veldb,sampdb,sampinterp] = pupilDeblink(samp,sampt,filter);   

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
    
    % fix beginning and end
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
    
    %% velocity profile of the filtered timecourse
    velfilt = zeros(size(sampfilt)); % velocity profile, mm/s
    % calculate the velocity
    for i=2:length(sampfilt)
        velfilt(i)= (sampfilt(i)-sampfilt(i-1))*1000/(sampt(i)-sampt(i-1));
    end
    velfilt(1)=velfilt(2);
    

   
% else
%     % if data is bad, discard
%     sampfilt = ones(size(samp))*NaN;
%     sampfiltz = ones(size(samp))*NaN;
% end