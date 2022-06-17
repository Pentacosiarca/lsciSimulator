function clusterHist = clusterBasedHist (normHist)

load('averageHistograms.mat')



clusterHist = zeros(size(normHist,1),size(normHist,2));

for iRows = 1:size(normHist,1)
    for iCols = 1:size(normHist,2)
        pixelHist = squeeze(normHist(iRows,iCols,:));
        [~,pValBF] = kstest2(averageHistBF,pixelHist);
        [~,pValPR] = kstest2(averageHistPR,pixelHist);

        whichPval = find([pValBF,pValPR]==min([pValBF,pValPR]));
        clusterHist(iRows,iCols) = whichPval(1);

    end
end

