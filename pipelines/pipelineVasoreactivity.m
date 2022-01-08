%pipelineVasoreactivity - calculates features from Vasoreactivity data
%
% Syntax:  output1 = function_name(input1,input2)
%
% Inputs:
%    path               - path to the .rls file
%    fileName           - name of the .rls file
%    Fs                 - sample frequency of the data
%
% Outputs:
%   ttp                 - time to peak from the start of the stimulation
%   auc_beforePeak      - Area Under the Curve from time of stimulation 
%                       to the peak
%   slope               - the slope of the "sigmoid" curve
%   rCBF                - relative cerebral blood flow map between an average 
%                       frame at the start of the simulation with an
%                       average frame at the peak.
%
% Other m-files required: processVasoreactivity.m
% Subfunctions: subfunctions are called inside other m-files
% MAT-files required: none
%
% See also: Other m-files required.

% Authors: Alberto Gonzalez Olmos, DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 17-November-2021

%------------- BEGIN CODE --------------

clear all, clc, close all,
% %% Add LSCI library path
lsciLibPath='./../../'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

path = 'C:\Users\Alberto\Desktop\github\data\CAIAD batch 2 samples\A54\Month1\';
name = '2021831_1109_vasoreactivity_10fps_20min.rls';

fileName = [path,name];

Fs = 10;
durationRecordingMinutes = 20 + 2; % adding 2 extra minutes
durationRecordingSeconds = durationRecordingMinutes * 60;


%See read RLS file for detailed description of arguments
batchSize = 1000; 
procType = 'fastcpu';
skipFrames = 0;
movMeanKernel = 1;
downSampling = 1; % set to 1 for no downsampling
ROI = []; 
maxFrames = durationRecordingSeconds * Fs;
%   ROI             - Region of Interest:
%                       [firstRow,      lastRow ;
%                       firstColumn,   lastColumn]

[sLSCI,sampling,timeStampsOut,sizeT] = runBatchLSCI(fileName, 'spatial', 5, batchSize, procType, ...
                                                   skipFrames, movMeanKernel, downSampling, ROI, maxFrames);

% Loading time stamps
load([path,'Month1_A54_VR.mat']);

% Get time stamps from the stimulation
stimStartTime = datestr(blocktimes,'HH:MM:SS.FFF');
stimStartTimeSeconds = str2num(stimStartTime(1:2))*360 + ...
                       str2num(stimStartTime(4:5))*60 + ...
                       str2num(stimStartTime(7:end));
stimSampleFrequency = samplerate(1);

dataStimEnd = durationRecordingSeconds * stimSampleFrequency;

stimData = data(datastart(1):dataStimEnd);

% %Get time in seconds from the start of LSCI recording
bfiTimeZero=(double(timeStampsOut)-double(timeStampsOut(1)))./1000;
bfiTime = datetime(timeStampsOut,'ConvertFrom','epochtime','Epoch',datetime(1970,1,1),'TicksPerSecond',1e3);
bfiTime = datestr(bfiTime,'HH:MM:SS.FFF');
bfiTimeSeconds = zeros(1,length(bfiTime));
for iTimeSeconds = 1:1:length(bfiTime)
    date = bfiTime(iTimeSeconds,:);
    bfiTimeSeconds(iTimeSeconds) = str2num(date(1:2))*360 + ...
                                   str2num(date(4:5))*60 + ...
                                   str2num(date(7:end));
end

% Frame at the start of the stimulus
diffTimes = abs(bfiTimeSeconds - stimStartTimeSeconds);
bfiStimStartFrame = find(diffTimes == min(diffTimes));

% Considering only the frames from the start of the stimulus
bfiTimeSecondsStartStim = (bfiTimeSeconds(bfiStimStartFrame:end));
sLsciStartStim = sLSCI(:,:,bfiStimStartFrame:end);
sLsciBaseline = sLSCI(:,:,1:bfiStimStartFrame);

clearvars -except ...
    bfiTimeSecondsStartStim ...
    sLsciStartStim ...
    sLsciBaseline ...
    fileName path

% calculating:
% ttp               - time to peak from the start of the stimulation
% auc_beforePeak    - Area Under the Curve from time of stimulation 
%                       to the peak
% slope             - the slope of the "sigmoid" curve
% rCBF              - relative cerebral blood flow map between an average 
%                       frame at the start of the simulation with an
%                       average frame at the peak.

[ttp, auc_beforePeak, slope, rCBF] = processVasoreactivity (bfiTimeSecondsStartStim, sLsciStartStim);






