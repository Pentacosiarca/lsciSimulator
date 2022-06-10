function centroids = seedCentroids(sizeX,sizeY)

factorNPoints = 10;
pointsX = round(linspace(1,sizeX,round(sizeX/factorNPoints)));
pointsX = pointsX(2:end-1);
pointsY = round(linspace(1,sizeY,round(sizeY/factorNPoints)));
pointsY = pointsY(2:end-1);

repX = repmat(pointsX,[length(pointsY) 1]);
repXreshaped = repX(:);
repY = repmat(pointsY,[1 length(pointsX)]);

centroids = [repXreshaped, repY'];