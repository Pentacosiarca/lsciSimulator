function averageMin = calculateAverageMin (bfiSP, getLocMinima)
    averageMin = mean([bfiSP(getLocMinima(1,1)), bfiSP(getLocMinima(2,1))]);
    for iSection = 2:1:size(getLocMinima,2)
        averageMin = (averageMin + ...
                      mean([bfiSP(getLocMinima(1,iSection)), bfiSP(getLocMinima(2,iSection))])) ...
                      /2; 
    end
end


