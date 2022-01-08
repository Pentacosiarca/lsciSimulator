function averageMean = calculateAverageMean (bfiSP, getLocMinima)
    averageMean = mean(bfiSP(getLocMinima(1,1):getLocMinima(2,1)));
    for iSection = 2:1:size(getLocMinima,2)
        averageMean = (averageMean + ...
                      mean(bfiSP(getLocMinima(1,1):getLocMinima(2,1)))) ...
                      /2; 
    end
end