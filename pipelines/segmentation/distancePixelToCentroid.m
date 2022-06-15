function distanceP2C = distancePixelToCentroid (dataHist,vectorizedIndex)

Nh = 1;
Ns = 1;

% for iRow = 1:size()

dataHistReshape = single(reshape(dataHist,[size(dataHist,1)*size(dataHist,2), size(dataHist,3)]));
centroidsHist = single(dataHistReshape(vectorizedIndex,:));

dataHistReshapeRepMat = repmat(dataHistReshape,[1, 1, length(vectorizedIndex)]);


















