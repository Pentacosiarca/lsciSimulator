function vectorizedIndex = vectorizedCentroids (centroids,nRows)

vectorizedIndex = nan(size(centroids,1),1);
for iCentroids = 1:length(vectorizedIndex)
    vectorizedIndex(iCentroids) = nRows * centroids(iCentroids+size(centroids,1)) + centroids(iCentroids);
end
