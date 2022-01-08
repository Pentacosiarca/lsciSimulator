%pipelinePulsatility - calculates Pulsatility Index and Resistivity Index
%after finding the local minima of normal pulses
%
% Syntax:  output1 = function_name(input1,input2)
%
% Setup parameters:
%    path           - path to the .rls file
%    fileName       - name of the .rls file
%    Fs             - sample frequency of the data
%
% Outputs:
%   pulsatility Index
%   resistivity Index 
%
% Other m-files required: loadSpatialContrastBfi.m, findSignalMinima.m, 
% calculatePulsatilityIndex.m, calculateResistivityIndex.m
% Subfunctions: subfunctions are called inside other m-files
% MAT-files required: none
%
% See also: The other m-files required.

% Authors: Alberto Gonzalez Olmos, DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 17-November-2021

%------------- BEGIN CODE --------------

clear all, clc, close all,
% %% Add LSCI library path
lsciLibPath='./'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

path = 'C:\Users\Alberto\Desktop\github\data\dataFromNumberFive\20211213\baslerPulsatility\';
name = '20211213_4.rls';
fileName = [path,name];

Fs = 194;

%See read RLS file for detailed description of arguments
batchSize = 100; 
procType = 'fastcpu';
skipFrames = 0;
movMeanKernel = 1;
downSampling = 1; % set to 1 for no downsampling
ROI = [];
%   ROI             - Region of Interest:
%                       [firstRow,      lastRow ;
%                       firstColumn,   lastColumn]

data = runBatchLSCI(fileName, 'spatial', 5, batchSize, procType, ...
                    skipFrames, movMeanKernel, downSampling, ROI);



bfiSP=nanmean(data,[1,2]);
bfiSP=squeeze(1./(bfiSP.^2));
bfiSP(bfiSP==Inf)=[];
clear data

locMinima = findSignalMinima (bfiSP, Fs);

% save('locMinPulsatility','locMinima')
% load('./../dataTemporalSpatialContrast')
% load('locMinPulsatility')

figure,hold on,
plot(bfiSP),
scatter(locMinima(1,:),bfiSP(locMinima(1,:)),'or'),
scatter(locMinima(2,:),bfiSP(locMinima(2,:)),'.b'),
legend('signal','Start Pulse, local minima','End Pulse, local minima')
title('Score of Minima based on difference of Max, Min and Min, Min'),
hold off,

%% Pulsatility Index and Resistivitiy Index

pulsatilityIndex = calculatePulsatilityIndex (bfiSP, locMinima);
resistivityIndex = calculateResistivityIndex (bfiSP, locMinima);












