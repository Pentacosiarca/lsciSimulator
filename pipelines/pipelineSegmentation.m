%% pipeline for segmentation


clear all, clc, close all,
% %% Add LSCI library path
lsciLibPath='./../'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

path = 'C:\Users\Alberto\Desktop\github\data\dataFromNumberFive\20211213\baslerPulsatility\';
name = '20211213_4.rls';
fileName = [path,name];

% visualizeData;

% filePathHist = 'C:\Users\Alberto\Desktop\github\data\segmentation\';
% if ~exist([filePathHist,'dataHist.mat'],'file')
%     [data,~,~,~]=readRLS(fileName);
%     dataHist = pixelToHistogram(data);
% else
%     load([filePathHist,'dataHist.mat'])
% end
% 
% nCols = size(dataHist,1); % columns
% nRows = size(dataHist,2); % rows
% nBins = size(dataHist,3);



[data,~,~,~]=readRLS(fileName);
data = single(data(:,:,1:1000));
dataSTD = std(data,[],3);
dataMEAN = mean(data,3);
dataLSCI_temporalContrast = dataSTD./dataMEAN;

[bloodFlow,~]=selectAreasToCompare(dataLSCI_temporalContrast);

[dataSegment,~,~,~]=readRLS(fileName,0,[],bloodFlow);

dataTLSCI = getTLSCI(dataSegment,size(dataSegment,3),'fastcpu');

[dataSegment,~,~,~]=readRLS(fileName);

dataHistSegment = pixelToHistogram(dataSegment);

maxHist = max(dataHistSegment,[],3);
minHist = min(dataHistSegment,[],3);
normHist = bsxfun(@minus,dataHistSegment,minHist) ./ (maxHist - minHist);

nCols = size(normHist,1); % columns
nRows = size(normHist,2); % rows
nBins = size(normHist,3);

% centroids = seedCentroids(nRows,nCols);
% nCentroids = size(centroids,1);

% vectorizedIndex = vectorizedCentroids (centroids,nRows);

% converge = 1;
% while (converge)
%     distanceP2C = distancePixelToCentroid (dataHist,vectorizedIndex);
%     distanceP2C = distancePixelToCentroid_loop (normHist,centroids);

clusterHist = clusterBasedHist (normHist);

    figure,subplot(1,2,1)
    imagesc(clusterHist)
    subplot(1,2,2)
    imagesc(dataTLSCI)

    figure,imagesc(dataLSCI_temporalContrast)

%     centroids = relocateCentroids(distanceP2C,nCentroids);

% end



%% other
isOther = 0;
if isOther
figure,plot(squeeze(dataHistBloodFlow(1,1,:)))

dataSLSCI = getTLSCI(dataSegment,size(dataSegment,3),'fastcpu');
figure,imagesc(dataSLSCI)

dataHistBloodFlow = pixelToHistogram(dataSegment);
dataHistBloodFlowReshape = reshape(dataHistBloodFlow,[size(dataHistBloodFlow,1)*size(dataHistBloodFlow,2),size(dataHistBloodFlow,3)]);

dataHistParenchyma = pixelToHistogram(dataParenchyma);
dataHistParenchymaReshape = reshape(dataHistParenchyma,[size(dataHistParenchyma,1)*size(dataHistParenchyma,2),size(dataHistParenchyma,3)]);

net = selforgmap([2 size(dataHist,3)]);
net = train(net,[dataHistBloodFlowReshape;dataHistParenchymaReshape]);

y = net([dataHistBloodFlowReshape;dataHistParenchymaReshape]);
classes = vec2ind(y);

load iris_dataset

figure;
imagesc(dataLSCI_temporalContrast(idx==1),'r.','MarkerSize',12)
hold on
imagesc(dataLSCI_temporalContrast(idx==2),'b.','MarkerSize',12)
hold off

figure,imagesc(dataLSCI_temporalContrast)

plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off
end


