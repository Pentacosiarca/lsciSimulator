%% pipeline for segmentation


clear all, clc, close all,
% %% Add LSCI library path
lsciLibPath='./'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

path = 'C:\Users\Alberto\Desktop\github\data\dataFromNumberFive\20211213\baslerPulsatility\';
name = '20211213_4.rls';
fileName = [path,name];

% visualizeData;

filePathHist = 'C:\Users\Alberto\Desktop\github\data\segmentation\';
if ~exist([filePathHist,'dataHist.mat'],'file')
    [data,~,~,~]=readRLS(fileName);
    dataHist = pixelToHistogram(data);
else
    load([filePathHist,'dataHist.mat'])
end

nCols = size(dataHist,1); % columns
nRows = size(dataHist,2); % rows
nBins = size(dataHist,3);

centroids = seedCentroids(sizeX,sizeY);

vectorizedIndex = vectorizedCentroids (centroids,nRows);

converge = 1;
while (converge)
    distanceP2C = distancePixelToCentroid (dataHist,vectorizedIndex);

end

