clear all, clc, close all
addpath(genpath('./../'))
% load ("speckleVScontrastK")
% 
% figure,scatter(speckleVScontrastK(2,:),speckleVScontrastK(1,:)),
% xlabel('Speckle Size')
% ylabel('K')

loadPath = 'C:\Users\Alberto\Desktop\github\data\simulation\';
% load([loadPath,'periodic_signal_10_by_10_fs_200_per_50'])

load([loadPath,'periodic_signal_10_by_10_fs_200_per_20_02'])

% contrastKernel = 5;
% batchSize = 100; 
% procType = 'fastgpu';
% 
% data = getSLSCI(timeIntensityAutocorrelationFunction,contrastKernel,procType,batchSize);

% % blood flow index signal from spatial contrast
% % find infinite values and discard them
% bfiTS=squeeze(mean(data,[1,2]));
% bfiTS=1./(bfiTS.^2);
% bfiTsInf = find(bfiTS == inf);
% if ~isempty(bfiTsInf)
%     data(:,:,bfiTsInf) = [];
%     bfiTimeSecondsStartStim(bfiTsInf) = [];
%     
%     % recalculate bfi without infinite values
%     bfiTS=squeeze(mean(data,[1,2]));
%     bfiTS=1./(bfiTS.^2);
% end
% 
% figure,plot(bfiTS)

dataLsci = timeIntensityAutocorrelationFunction;
clear timeIntensityAutocorrelationFunction

% dataLsci = dataLsci(:,:,1:3:end);
downsampleRatio = floor(length(periodicMovement)/size(dataLsci,3));
periodicDownSampled = downsample(periodicMovement,downsampleRatio);
periodicDownSampled(size(dataLsci,3)+1:end) = [];

x = 1:length(periodicDownSampled);
[~,loc] = findpeaks(-periodicDownSampled);

% locMinMedian = findSignalMinima (periodicDownSampled, 200);

locMinima = islocalmax(-periodicDownSampled);
locMinima(1) = 1;
locMinima(end) = 1;

isException = 1;
if isException
deleteLocMin = 11:20:loc(end);
locMinima(deleteLocMin) = 0;
end

figure,plot(x,periodicDownSampled,x(locMinima),periodicDownSampled(locMinima),'r*')

positionLocMin = find(locMinima>0);



%% --------------------------------
% mean values of Temporal Contrast
% --------------------------------
isTempContrast = 1;
if isTempContrast

typeOfAveraging = 'meanContrast';
typeOfContrast = 'temporal';

typeOfAveraging = cell({typeOfAveraging, typeOfContrast});
averagePulse = calculateAveragePulse (dataLsci,positionLocMin,typeOfAveraging);

% bfiTSave=squeeze(mean(averagePulse,[1,2]));
bfiTSave=squeeze(averagePulse(end/2,end/2,:));

bfiTSave=1./(bfiTSave.^2);
figure,hold on,
yyaxis left
plot(periodicDownSampled(positionLocMin(1):positionLocMin(2))),
yyaxis right
plot(bfiTSave)
legend('ideal pulse','calculated average pulse')
hold off
title(['average of ',typeOfContrast,' Contrast'])
end

%% --------------------------------
% mean values of Spatial Contrast
% --------------------------------
isSpatialContrast = 1;
if isSpatialContrast
typeOfAveraging = 'meanContrast';
typeOfContrast = 'spatial';

typeOfAveraging = cell({typeOfAveraging, typeOfContrast});
averagePulse = calculateAveragePulse (dataLsci,positionLocMin,typeOfAveraging);

% bfiTSave=squeeze(mean(averagePulse,[1,2]));
bfiTSave=squeeze(averagePulse(end/2,end/2,:));
bfiTSave=1./(bfiTSave.^2);
figure,hold on,
yyaxis left
plot(periodicDownSampled(positionLocMin(1):positionLocMin(2))),
yyaxis right
plot(bfiTSave)
legend('ideal pulse','calculated average pulse')
hold off
title(['average of ',typeOfContrast,' Contrast'])
end

%% --------------------------------
% Temporal Contrast average
% --------------------------------
isTempContrastAverage = 1;
if isTempContrastAverage
typeOfAveraging = cell({'temporalContrastMean'});
averagePulse = calculateAveragePulse (dataLsci,positionLocMin,typeOfAveraging);

% bfiTSave=squeeze(mean(averagePulse,[1,2]));
bfiTSave=squeeze(averagePulse(end/2,end/2,:));
bfiTSave=1./(bfiTSave.^2);
figure,hold on,
yyaxis left
plot(periodicDownSampled(positionLocMin(1):positionLocMin(2))),
yyaxis right
plot(bfiTSave)
hold off
legend('ideal pulse','calculated average pulse')
title('Lossless Temporal contrast average')
end

%% --------------------------------
% LSCI signal
% --------------------------------
isLSCIsignal = 1;
if isLSCIsignal
% bfiTSave=squeeze(mean(dataLsci,[1,2]));
bfiTSave=squeeze(dataLsci(end/2,end/2,:));
bfiTSave=1./(bfiTSave.^2);
figure,hold on,
yyaxis left
plot(periodicDownSampled(positionLocMin(1):positionLocMin(2))),
yyaxis right
plot(bfiTSave(positionLocMin(1):positionLocMin(2)))
hold off
legend('ideal pulse','calculated average pulse')
title('Raw LSCI signal')
end




