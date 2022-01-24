clear all, clc, close all,

% camera parameters
cameraNoise = 0; % to implement
exposureTime = 0.005; % in seconds
Fs = 1000; % Sampling frequency (samples per second) 
nPeriods = 10; 
mainHeartFreq = 10;

%% processing the simulation data

savePath = './../data/simulation/';
load([savePath,'periodic_signal_500_by_500_fs_200_per_5.mat'])

lsciLibPath='./'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

recordingTime = nPeriods / mainHeartFreq;
recordingSamples = recordingTime * Fs;
framesToAverage = floor(size(timeIntensityAutocorrelationFunction,3) / recordingSamples);

dataLsci = nan(size(timeIntensityAutocorrelationFunction,1),...
               size(timeIntensityAutocorrelationFunction,2),...
               floor(size(timeIntensityAutocorrelationFunction,3)/framesToAverage));
for iAverage = 1:1:floor(size(timeIntensityAutocorrelationFunction,3)/framesToAverage)
    batchAverage = (iAverage-1)*framesToAverage + 1 : iAverage*framesToAverage;
    dataLsci(:,:,iAverage) = gpuArray(mean(timeIntensityAutocorrelationFunction(:,:,batchAverage),3));
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


%% create gif

isCreateGif = 0;
if isCreateGif
        kindOfMotion = 'periodic';
        fileName = [savePath,kindOfMotion,'_lscigif.gif'];
        h = figure('Visible','off');
        for iFrame = 1:size(dataLsci,3)
            sensorImage = dataLsci(:,:,iFrame);
            imagesc(sensorImage)
            title(['frame: ',num2str(iFrame)])
            frame = getframe(h); 
            im = frame2im(frame); 
            [imind,cm] = rgb2ind(im,256); 
            % Write to the GIF File 
            if iFrame == 1 
              imwrite(imind,cm,fileName,'gif', 'Loopcount',inf); 
            else 
              imwrite(imind,cm,fileName,'gif','WriteMode','append'); 
            end 
        end
end
