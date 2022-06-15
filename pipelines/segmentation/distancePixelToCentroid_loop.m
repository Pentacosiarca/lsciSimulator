%% Distance Superpixels to Centroids
function distanceP2C = distancePixelToCentroid_loop (normBF,centroids)

    distanceP2C = zeros(size(normBF,1),size(normBF,2));
    Nst = 0.0001;
    Nsp = 100;
    
    for iRows = 1:size(normBF,1)
        for iCols = 1:size(normBF,2)
            normBF_pixelHist = squeeze(normBF(iRows,iCols,:));
            Dst = distanceStats (normBF_pixelHist,normBF,centroids);
            Dsp = distanceEucl (iRows,iCols,centroids);
            DistanceMetric = sqrt((1/Dst*Nst).^2 + (Dsp/Nsp).^2);

            minDist = find(DistanceMetric == min(DistanceMetric));
            distanceP2C(iRows,iCols) = minDist(1);

        end
    end

end


%% Statistical Distance
function Dst = distanceStats (normBF_pixelHist,normBF,centroids)
    Dst = zeros(size(centroids,1),1);
    for iCentroids = 1:size(centroids,1)
        centroidHist = squeeze(normBF(centroids(iCentroids,2),centroids(iCentroids,1),:));
%         keyboard
%         h = kstest2(centroidHist,normBF_pixelHist);
%         figure,hold on,
%         plot(normBF_pixelHist)
%         plot(centroidHist)
%         hold off
        [~,Dst(iCentroids)] = kstest2(centroidHist,normBF_pixelHist);
%         Dst(iCentroids) = kstest2(centroidHist,normBF_pixelHist); %sum(abs(centroidHist - normBF_pixelHist));
    end
end

%% Euclidian Distance
function Dsp = distanceEucl (iRows,iCols,centroids)
        abSum = centroids - [iCols, iRows];
        abSquare = abSum.^2;
        aSumb = sum(abSquare,2);
        Dsp = sqrt(aSumb);
end





