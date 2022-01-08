function pulsatilityIndex = calculatePulsatilityIndex (bfiSP, getLocMinima)

    averageMax = calculateAverageMax (bfiSP, getLocMinima);
    averageMin = calculateAverageMin (bfiSP, getLocMinima);
    averageMean = calculateAverageMean (bfiSP, getLocMinima);

    pulsatilityIndex = (averageMax - averageMin)/averageMean;

end