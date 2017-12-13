% clear
% cd 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData';

 subj = [21 22 23 24 25 26 28 29 30 31 32 33 34];
 choiceall = zeros(100*length(subj),5);

 for n = 1:length (subj)
     ['ASD' num2str(subj(n)) '.txt']
    s = tdfread(['ASD' num2str(subj(n)) '.txt']);

    %% Response and RT for choice

    % Extract data for the key press for answering question
    sChoice = struct; % sChoice stores information during choice
    % sAnswer.Acc = zeros(1,25);
    i = 0; % counting index for going over each element in s in sequence
    j = 0; % No. for trial 1~100
    t0= 0;
    t1= 0; 
    sChoice.Rt = zeros(1,100);

    while i < length(s.Timestamp)
        i=i+1;
        if s.EventKey ~= 8 | s.Descriptor(i,1) ~= 'a'; continue; end
        j=j+1;
        sChoice.Trial(j,:) = s.Descriptor(i,:);
        while s.EventKey(i) ~= 8 | s.Descriptor(i,1) ~= 'c' 
            i=i+1;
        end
        t0 = s.Timestamp(i);
        while s.EventKey(i)~=4
            i=i+1;
        end
        sChoice.Press(j,:) = s.Descriptor(i,:);
        while s.EventKey(i)~= 16 | s.Descriptor(i,1) ~= 'c'
            i=i+1;
        end
        t1=s.Timestamp(i-1);
        sChoice.Rt(j) = t1-t0;
    end

    % Change response string to number
    for j=1:100
        sChoice.Resp(j) = str2num(sChoice.Press(j,7));
    end

    %ambig and val
    for j=1:100
        sChoice.ar(j,:) = sChoice.Trial(j,1:11);
    end

    %lottery side
    sChoice.Side = [];
    for i = 1:100
        sChoice.Side = [sChoice.Side; 'Right'];
    end

    %% Create a matrix for writing data into text
    % Column name: 1-subj, 2-al, 3-val,4-choice(lottery 1,reference 0, no response 2), 5-RT
    for j = 1:length(sChoice.Trial)
        choiceall(j+100*(n-1),1) = subj(n);
        choiceall(j+100*(n-1),2) = str2num(sChoice.Trial(j,6:8));
        choiceall(j+100*(n-1),3) = str2num(sChoice.Trial(j,10:11));
        
        if (strcmp(sChoice.Side(j,:),'Left') && sChoice.Resp(j)==1) || (strcmp(sChoice.Side(j,:),'Right') && sChoice.Resp(j)==2)
        % lottery chosen
        choiceall(j+100*(n-1),4) = 1;
        elseif sChoice.Resp(j) == 0;
        % no response
        choiceall(j+100*(n-1),4) = 2;
        else
        % reference chosen
        choiceall(j+100*(n-1),4) = 0;
        end
        
        choiceall(j+100*(n-1),5) = sChoice.Rt(j);
    end
 end

%% Write data, choice for all subject and all trials

% results file
 fid = fopen(['ASDChoiceAll.txt'],'w');
 fprintf(fid,'%s\t %s\t %s\t %s\t %s\n', 'Subj','Al', 'Val', 'Choice', 'Rt');
 fprintf(fid, '%d\t %d\t %d\t %d\t %d\n',choiceall');
 fclose(fid)
%% Write data into text file, prepare the logistic regression, for individual subject

% choice_result = [sChoice.Resp; sChoice.ar'; sChoice.Side'; sChoice.Rt];
% choice_result = [sChoice.Resp; sChoice.Rt];
% 
% choice_result = [sChoice.ar'; sChoice.Side'];
% 
% 
% 
% 
% path = 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\AS_DecisionTobiiData\Decision Tobii Data Matlab analysis\';
% choice_file = 'ASD21.txt';
% 
% % results file
% fid = fopen([path choice_file],'w');
% fprintf(fid,'%s\n%s\t%s\t%s\t%s\n', 'Decision','Resp', 'RiskAmbigLevel', 'LotterySide', 'RT')
% fprintf(fid, '%d\t%s\t%s\t%d\n',choice_result)
% fclose(fid)


