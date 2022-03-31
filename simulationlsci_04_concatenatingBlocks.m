function simulationlsci_04_concatenatingBlocks (blocksFolder)

blockFiles = dir([blocksFolder,'block*']);

load([blocksFolder,blockFiles(1).name])

data = timeIntensityAutocorrelationFunction;

blockSize = size(data,3);
masterBlock = nan(size(data,1),...
                  size(data,2),...
                  blockSize*length(blockFiles));

masterBlock(:,:,1:blockSize) = data;

for iBlock = 2:length(blockFiles)
    load([blocksFolder,blockFiles(iBlock).name])
    masterBlock(:,:,(iBlock-1) * blockSize+1 : iBlock * blockSize) = data;

end

fileName = 'BlockOfSpecklePatterns.mat';
save([blocksFolder,fileName],'masterBlock','-mat','-v7.3');