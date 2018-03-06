function [vel,veldb,sampdb] = pupilDeblink(samp,sampt,filter)

%% This function de-blink the raw pupil data
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
