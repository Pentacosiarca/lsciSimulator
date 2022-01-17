clear all, clc, close all,

% camera parameters
cameraNoise = 0; % to implement
exposureTime = 0.005; % in seconds
Fs = 1000; % Sampling frequency (samples per second) 
nPeriods = 10; 

%% processing the simulation data

savePath = './../data/simulation/';
load([savePath,'periodic_signal_50_by_50.mat'])

lsciLibPath='./'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

framesToAverage = exposureTime *Fs;

dataLsci = nan(size(timeIntensityAutocorrelationFunction,1),...
               size(timeIntensityAutocorrelationFunction,2),...
               floor(size(timeIntensityAutocorrelationFunction,3)/framesToAverage));
for iAverage = 1:1:floor(size(timeIntensityAutocorrelationFunction,3)/framesToAverage)
    batchAverage = (iAverage-1)*framesToAverage + 1 : iAverage*framesToAverage;
    dataLsci(:,:,iAverage) = mean(timeIntensityAutocorrelationFunction(:,:,batchAverage),3);
end

procType = 'fastcpu';
dataSLSCI=getSLSCI(dataLsci,5,procType);

bfiSpatial = squeeze(1./mean(dataSLSCI,[1 2]).^2);
figure,plot(bfiSpatial)
title('bfi of spatial contrast')

dataTLSCI=getTLSCI(dataLsci,25,procType);

bfiTemporal = squeeze(1./mean(dataTLSCI,[1 2]).^2);
figure,plot(bfiTemporal)
title('bfi of temporal contrast')

figure,plot(bfiSpatial)
hold on
plot(bfiTemporal)
