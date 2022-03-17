%% main simulation script

nPeriods = 50;

lastFramePeriodic = 0;
lastFrameBrownian = 0;

for iPeriod = 1:1:nPeriods

kindOfMotion = 'periodic'; 

periodicBlockOfParticles = simulationlsci_01_particlesPositions (kindOfMotion, lastFramePeriodic);

kindOfMotion = 'brownian';

brownianBlockOfParticles = simulationlsci_01_particlesPositions (kindOfMotion, lastFrameBrownian);

[blockComposition, lastFramePeriodic, lastFrameBrownian] = simulationlsci_02_largerBlockComposition (periodicBlockOfParticles, brownianBlockOfParticles);

simulationlsci_03_specklePattern (blockComposition);

end






