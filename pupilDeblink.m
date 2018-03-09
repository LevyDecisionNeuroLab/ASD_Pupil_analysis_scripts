function [vel,veldb,sampdb,sampinterp] = pupilDeblink(samp,sampt,filter)

%% This function de-blink and interpolate the raw pupil data
%
% Input:
%   - samp: pupil size data
%   - sampt: pupil size data time point
%   - filter: filter fype and parameters
%       - velThreshold: how many std is compared with the velocity difference from mean
% Output:
%   - vel: velocity profile before deblinking
%   - veldb: velocity profile after deblinking
%   - sampdb: deblinked pupil size data

%% De-blink
% Detect blink both left and right pupil data
 
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

if ~isempty(exist)
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
else % if after de-blink, all data are missing
     sampinterp=sampdb;
end
