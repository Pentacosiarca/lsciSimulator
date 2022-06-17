%% analyze and save typical histogram for blood flow and parenchyma


clear all, clc, close all,
% %% Add LSCI library path
lsciLibPath='./../'; %path to the root folder with LSCI processing scripts
addpath(genpath(lsciLibPath));

path = 'C:\Users\Alberto\Desktop\github\data\dataFromNumberFive\20211213\baslerPulsatility\';
name = '20211213_4.rls';
fileName = [path,name];

[data,~,~,~]=readRLS(fileName);
data = single(data(:,:,1:1000));
dataSTD = std(data,[],3);
dataMEAN = mean(data,3);
dataLSCI_temporalContrast = dataSTD./dataMEAN;

[bloodFlow,parenchyma]=selectAreasToCompare(dataLSCI_temporalContrast);

[dataBloodFlow,~,~,~]=readRLS(fileName,0,[],bloodFlow);
[dataParenchyma,~,~,~]=readRLS(fileName,0,[],parenchyma);

dataHistBloodFlow = pixelToHistogram(dataBloodFlow);

maxBF = max(dataHistBloodFlow,[],3);
minBF = min(dataHistBloodFlow,[],3);
normBF = bsxfun(@minus,dataHistBloodFlow,minBF) ./ (maxBF - minBF);

averageHistBF = squeeze(mean(normBF,[1 2]));

dataHistParenchyma = pixelToHistogram(dataParenchyma);

maxPR = max(dataHistParenchyma,[],3);
minPR = min(dataHistParenchyma,[],3);
normPR = bsxfun(@minus,dataHistParenchyma,minPR) ./ (maxPR - minPR);

averageHistPR = squeeze(mean(normPR,[1 2]));


figure,hold on
plot(averageHistBF,'b')
plot(averageHistPR,'r')
hold off
legend('region with blood flow, average histogram','parenchyma, average histogram')

save('averageHistograms.mat',"averageHistPR","averageHistBF",'-mat')




