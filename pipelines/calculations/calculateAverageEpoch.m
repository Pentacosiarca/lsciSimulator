function [averageEpochBaseline, averageEpochStimulation] = calculateAverageEpoch (epochBaseline, epochStimulation, dataLsci, dataBfi, bfiTimeSeconds, Fs)

baselineFrames = epochBaseline*Fs;
stimulationFrames = epochStimulation*Fs;
epochLength = baselineFrames + stimulationFrames;


%% detecting motion-like artifacts in blood flow index
nEpoch = floor(length(dataBfi)/(epochLength));
lengthSignalEpoch = nEpoch * epochLength;
dataBfiEpoch = dataBfi (1:lengthSignalEpoch);
medfiltBfi = medfilt1(dataBfiEpoch,100);
medianMedFiltBfi = median(medfiltBfi);
medianThreshold = 1.1;
findAboveThreshold = find(medfiltBfi > medianMedFiltBfi*medianThreshold);
epochZeros = zeros(1,length(medfiltBfi));
epochZeros(findAboveThreshold) = 1;
epochSegmentZeros = zeros(nEpoch,epochLength);
for iSegment = 1:1:nEpoch
    epochBatch = ((iSegment-1)*epochLength) + 1: ((iSegment-1)*epochLength) + epochLength;
    epochSegmentZeros(iSegment,:)  = epochZeros(epochBatch);
end
epochSegmentZerosSum = sum(epochSegmentZeros,2);
idxAbnormalEpoch = find(epochSegmentZerosSum > 2*median(epochSegmentZerosSum));


%% average spatial contrast epoch
averageEpoch = nan(size(dataLsci,1),size(dataLsci,1),epochLength);
epochNumFrame = (epochOrder * epochLength) + 1;
tic;
for iFrame = 1:1:epochLength
    idxFrame = epochNumFrame + (iFrame - 1);
    frame=gpuArray(single(squeeze(mean(dataLsci(:,:,idxFrame),3))));
    averageEpoch(:,:,iFrame) = gather(frame);
end
averageFrameTime = toc;
disp(['Calculating the average epoch took: ',num2str(averageFrameTime),' seconds.'])

% clear dataLsci
averageEpochBaseline = averageEpoch(:,:,1:baselineFrames);
averageEpochStimulation = averageEpoch(:,:,baselineFrames+1:epochLength);





