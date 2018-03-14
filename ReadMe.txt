Functions:
pupilDeblink: de-blink and interpolate pupil trace
pupilPrepro: preprocess pupil size from singe eye
combineLeftRight: preprocess pupil size from both eyes and combine them


Pupil size analysis steps:
Step1. Run ASDPupil_extract, which extracts pupil data and choice data from the Tobii exported files
Step2. Run ASDPupil_prepro, which preprocesses pupil data and combine left and right eyes
Step3. Run ASDPupil_normalize, which normalizes the pupil trace for each trial by subtract the mean of the 1000ms pre-stimulus pupil data
After the previous three steps, pupil size data are ready for analysis. Analyses scripts are:
1. ASDPupil_swlm: sliding window regression
2. ASDPupil_chocice: pupil size analyzed by the choice of subjects.