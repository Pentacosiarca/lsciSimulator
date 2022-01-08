function resistivityIndex = calculateResistivityIndex (bfiSP, getLocMinima)
    
    averageMax = calculateAverageMax (bfiSP, getLocMinima);
    averageMin = calculateAverageMin (bfiSP, getLocMinima);
    
    resistivityIndex = (averageMax - averageMin)/averageMin;

end