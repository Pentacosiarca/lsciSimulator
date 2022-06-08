%% pipeline for segmentation


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

data = runBatchLSCI(fileName, 'temporal', 25, batchSize, procType, ...
                    skipFrames, movMeanKernel, downSampling, ROI);

[dataLSCI,~,~,~]=readRLS(fileName);
dataLSCI = single(dataLSCI(:,:,1:1000));
dataSTD = std(dataLSCI,[],3);
dataMEAN = mean(dataLSCI,3);
dataLSCI_temporalContrast = dataSTD./dataMEAN;

[bloodFlow,parenchyma]=selectAreasToCompare(dataLSCI_temporalContrast);

[dataBloodFlow,~,~,~]=readRLS(fileName,0,[],bloodFlow);
[dataParenchyma,~,~,~]=readRLS(fileName,0,[],parenchyma);

h = figure;
for iter = 1:size(dataParenchyma,3)
    figure(h)
    subplot(2,1,1)
    histogram(dataBloodFlow(:,:,iter),'BinLimits',[40,150])
    subplot(2,1,2)
    histogram(dataParenchyma(:,:,iter),'BinLimits',[40,150])
end

dataSLSCI = getSLSCI(dataBloodFlow,7,procType);

figure,imagesc(dataSLSCI(:,:,1))

bfi = double(squeeze(mean(dataSLSCI,[1 2])));
bfi = 1./bfi.^2;
figure,plot(bfi)

    figure
    subplot(2,1,1)
    histogram(dataBloodFlow,'BinLimits',[40,150])
    subplot(2,1,2)
    histogram(dataParenchyma,'BinLimits',[40,150])

figure,histogram(dataBloodFlow)
hold on
histogram(dataParenchyma)

seedClusters

