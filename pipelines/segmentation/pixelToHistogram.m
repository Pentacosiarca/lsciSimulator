function dataHist = pixelToHistogram(data)

limBinHist = [40,100];
nBins = 48;

dataHist = zeros(size(data,1),size(data,2),nBins);

tic;

for iRow = 1:size(data,1)
    for iCol = 1:size(data,2)
         hist = histogram(data(iRow,iCol,:),nBins,'BinLimits',limBinHist);
         dataHist(iRow,iCol,:) = hist.Values;
    end
end

timeLapse = toc;
disp('elapsed time: ',timeLapse);

