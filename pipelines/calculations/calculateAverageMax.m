function averageMax = calculateAverageMax (bfiSP, getLocMinima)
    averageMax = max (bfiSP (getLocMinima(1,1):getLocMinima(2,1)));
    for iSection = 2:1:size(getLocMinima,2)
        averageMax = (averageMax + ...
                      max(bfiSP (getLocMinima(1,iSection):getLocMinima(2,iSection)))) ...
                      /2; 
    end
end