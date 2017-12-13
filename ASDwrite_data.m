%% Print the data into text file for each subject

result = zeros(7,100);
result(1,:) = repmat(str2num(subject),1,100);
result(6,:) = repmat(aP,1,100);
result(7,:) = repmat(bP,1,100);

result(2,:) = sDelay.PupilLeftNorm;
result(3,:) = sTrial.AL;
result(4,:) = sTrial.Val;

for i=1:size(result,2)
    if result(3,i) == 0
        result(5,i) = 0;
    elseif result(3,i) == 24
        result(5,i) = ambigatt(1);
    elseif result(3,i) == 50
        result(5,i) = ambigatt(2);
    elseif result(3,i) == 74
        result(5,i) = ambigatt(3);
    elseif result(3,i) == 100
        result(5,i) = ambigatt(4);
    end
end

path = 'C:\Users\rj299\Documents\Projects in the lab\Ambiguity-as-stressor Project\Tobii script\AS_PatternPilotData\Statistical analysis\Decision behav_pupil\';
results_file = 'ASD21.txt';

% results file
fid = fopen([path results_file],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\n', 'subj', 'pupil', 'al', 'val', 'ambigatt', 'alpha', 'beta')
fprintf(fid, '%d\t%f\t%d\t%d\t%f\t%f\t%f\n',result)
fclose(fid)