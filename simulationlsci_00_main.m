%% main simulation script
clear all,clc,close all

nPeriods = 2;

lastFramePeriodic = 0;
lastFrameBrownian = 0;
nMoves = 0;

% load('dataSetPoint.mat')

for iPeriod = 1:1:nPeriods

    disp(['period number: ',num2str(iPeriod)])

kindOfMotion = 'periodic'; 

[periodicBlockOfParticles, nMoves] = simulationlsci_01_particlesPositions (kindOfMotion, lastFramePeriodic,nMoves);

kindOfMotion = 'brownian';

[brownianBlockOfParticles, nMoves] = simulationlsci_01_particlesPositions (kindOfMotion, lastFrameBrownian,nMoves);


[blockComposition, lastFramePeriodic, lastFrameBrownian] = simulationlsci_02_largerBlockComposition (periodicBlockOfParticles, brownianBlockOfParticles);


simulationlsci_03_specklePattern (blockComposition);

end

blocksFolder = './../data/simulation/';

simulationlsci_04_concatenatingBlocks (blocksFolder);






