function centroids = relocateCentroids(distanceP2C,nCentroids)

for iNewCentr = 1:nCentroids
    find(distanceP2C == iNewCentr);
end