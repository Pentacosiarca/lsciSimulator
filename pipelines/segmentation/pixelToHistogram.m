function dataHist = pixelToHistogram(data)

limBinHist = [20,150];
nBins = 100;

tic;

dataHist = calculateHistogram(data,nBins,limBinHist);

% dataHist = zeros(size(data,1),size(data,2),nBins);
% for iRow = 1:size(data,1)
%     for iCol = 1:size(data,2)
%          hist = histogram(data(iRow,iCol,:),nBins,'BinLimits',limBinHist);
%          dataHist(iRow,iCol,:) = hist.Values;
%     end
% end

timeLapse = toc;
disp(['elapsed time: ', num2str(timeLapse)]);

