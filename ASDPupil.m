%% This script analyze the already-preprocessed pupil size data

% clearvars
% clc

%% 
% pupil by value
uniqueval = unique(sTrial.Val);
for i = 1:5
    averageByVal(i,:) = nanmean(sInitial.PupilLeft_filt(sTrial.Val <= uniqueval(i*4),:));
    stdByVal(i,:) = nanstd(sInitial.PupilLeft_filt(sTrial.Val <= uniqueval(i*4),:));
    seByVal(i,:) = std ./ sqrt(size(sInitial.PupilLeft_filt,1));
end

colors = {[0.2 0 0],[0.4 0 0],[0.6 0 0],[0.8 0 0],[1 0 0]};
figure
for i = 1:5
    plot(sInitial.Timestamp(1,:), averageByVal(i,:), 'LineStyle', '-', 'Marker', '.', 'Color',colors{6-i})
    hold on
end

% Time series and figures
ts_ITILeft_int = timeseries(sITI.PupilLeft_int,sITI.Timestamp(1,:),'name', 'PupilSize_ITI_int');
figure('Name','ITILeft','NumberTitle','off')
% hold on
plot (ts_ITILeft_int)

ts_LottLeft_int = timeseries(sLott.PupilLeft_int,sLott.Timestamp(1,:),'name', 'PupilSize_Lott_int');
figure('Name','LottLeft','NumberTitle','off')
% hold on
plot (ts_LottLeft_int)

ts_DelayLeft_int = timeseries(sDelay.PupilLeft_int,sDelay.Timestamp(1,:),'name', 'PupilSize_Delay_int');
figure('Name','DelayLeft','NumberTitle','off')
% hold on
plot (ts_DelayLeft_int)

%% Trial information, summary pupil size into ambig and val levels

% Pupil size by ambig level
for i = 1: size(sTrial.AL,1)
    if sTrial.AL(i) == 0; sTrial.a(i) = 1;
    elseif sTrial.AL(i) == 24; sTrial.a(i) = 2;
    elseif sTrial.AL(i) == 50; sTrial.a(i) = 3;
    elseif sTrial.AL(i) == 74; sTrial.a(i) = 4;
    elseif sTrial.AL(i) == 100; sTrial.a(i) = 5;
    end
    
    if sTrial.Val(i) == 4; sTrial.b(i) = 1;
    elseif sTrial.Val(i) == 5; sTrial.b(i) = 2;
    elseif sTrial.Val(i) == 6; sTrial.b(i) = 3;
    elseif sTrial.Val(i) == 7; sTrial.b(i) = 4;
    elseif sTrial.Val(i) == 8; sTrial.b(i) = 5;
    elseif sTrial.Val(i) == 9; sTrial.b(i) = 6;
    elseif sTrial.Val(i) == 10; sTrial.b(i) = 7;
    elseif sTrial.Val(i) == 11; sTrial.b(i) = 8;
    elseif sTrial.Val(i) == 12; sTrial.b(i) = 9;
    elseif sTrial.Val(i) == 13; sTrial.b(i) = 10;
    elseif sTrial.Val(i) == 14; sTrial.b(i) = 11;
    elseif sTrial.Val(i) == 15; sTrial.b(i) = 12;
    elseif sTrial.Val(i) == 16; sTrial.b(i) = 13;
    elseif sTrial.Val(i) == 17; sTrial.b(i) = 14;
    elseif sTrial.Val(i) == 18; sTrial.b(i) = 15;
    elseif sTrial.Val(i) == 23; sTrial.b(i) = 16;
    elseif sTrial.Val(i) == 34; sTrial.b(i) = 17;
    elseif sTrial.Val(i) == 39; sTrial.b(i) = 18;
    elseif sTrial.Val(i) == 57; sTrial.b(i) = 19;
    elseif sTrial.Val(i) == 68; sTrial.b(i) = 20;
    end
end
pupilITI = zeros(5,20);
pupilLott = zeros(5,20);
pupilDelay = zeros(5,20);
for i = 1: size(sTrial.AL,1)
    pupilITI(sTrial.a(i),sTrial.b(i)) = sITI.PupilLeftMean(i);
    pupilLott(sTrial.a(i),sTrial.b(i)) = sLott.PupilLeftMean(i);
    pupilDelay(sTrial.a(i),sTrial.b(i)) = sDelay.PupilLeftMean(i);
    RT(sTrial.a(i),sTrial.b(i)) = sChoice.Rt(i);
end

%%  Devide the lott and delay period pupil size into ambiguity levels, compute the pupil mean at each time point, and the std at each time point

% Lottery presentation period
for i = 1:5
    for j = 1:20
        pupilLottNorm (i,j) = pupilLott(i,j)/pupilITI(i,j);
    end
end

for i = 1: size (sLott.PupilLeft_int,1)
    sLott.PupilLeft_intNorm(i,:) = sLott.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
end

temp1 = sLott.PupilLeft_intNorm;
a1 = find(sTrial.a == 1);
a2 = find(sTrial.a == 2);
a3 = find(sTrial.a == 3);
a4 = find(sTrial.a == 4);
a5 = find(sTrial.a == 5);

pupilLottAL(1,:) = nanmean(temp1(a1,:));
pupilLottAL(2,:) = nanmean(temp1(a2,:));
pupilLottAL(3,:) = nanmean(temp1(a3,:));
pupilLottAL(4,:) = nanmean(temp1(a4,:));
pupilLottAL(5,:) = nanmean(temp1(a5,:));

pupilLottALstd(1,:) = nanstd(temp1(a1,:));
pupilLottALstd(2,:) = nanstd(temp1(a2,:));
pupilLottALstd(3,:) = nanstd(temp1(a3,:));
pupilLottALstd(4,:) = nanstd(temp1(a4,:));
pupilLottALstd(5,:) = nanstd(temp1(a5,:));

ts_pupilLottAL = timeseries(pupilLottAL,sLott.Timestamp(1,:),'name', 'PupilSize_Lott_byAL');
figure('Name','Pupil_Lott_byAL','NumberTitle','off')
plot (ts_pupilLottAL)

% delay period
for i = 1:5
    for j = 1:20
        pupilDelayNorm (i,j) = pupilDelay(i,j)/pupilITI(i,j);
    end
end

for i = 1: size (sDelay.PupilLeft_int,1)
    sDelay.PupilLeft_intNorm(i,:) = sDelay.PupilLeft_int(i,:)/sITI.PupilLeftMean(i);
end

temp2 = sDelay.PupilLeft_intNorm;
a1 = find(sTrial.a == 1);
a2 = find(sTrial.a == 2);
a3 = find(sTrial.a == 3);
a4 = find(sTrial.a == 4);
a5 = find(sTrial.a == 5);

pupilDelayAL(1,:) = nanmean(temp2(a1,:));
pupilDelayAL(2,:) = nanmean(temp2(a2,:));
pupilDelayAL(3,:) = nanmean(temp2(a3,:));
pupilDelayAL(4,:) = nanmean(temp2(a4,:));
pupilDelayAL(5,:) = nanmean(temp2(a5,:));

pupilDelayALstd(1,:) = nanstd(temp2(a1,:));
pupilDelayALstd(2,:) = nanstd(temp2(a2,:));
pupilDelayALstd(3,:) = nanstd(temp2(a3,:));
pupilDelayALstd(4,:) = nanstd(temp2(a4,:));
pupilDelayALstd(5,:) = nanstd(temp2(a5,:));


ts_pupilDelayAL = timeseries(pupilDelayAL,sDelay.Timestamp(1,:),'name', 'PupilSize_Delay_byAL');
figure('Name','Pupil_Delay_byAL','NumberTitle','off')
plot (ts_pupilDelayAL)


