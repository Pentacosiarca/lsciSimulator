function blockOfParticles = simulationlsci_01_particlesPositions (kindOfMotion, lastFrameParticles)

if lastFrameParticles == 0

% Modelling of particles in a volume
%% setup parameters:
% generating randomly distributed particles in a volume
nParticles = single(1000); % reference value: 100
particlesVolumeSizeXyz = single([10; 10; 0.01]); % reference value: [10; 10; 0.01]

% defining sensor
sensorDistanceToZ = single(1000); % reference value: 1000
sensorXY = 2; % [0.1,0.2,0.3,0.6,1,1.5,2.1,2.8.]
sensorSizeXy = single([sensorXY; sensorXY]); % reference value: [0.5; 0.5]
sensorResolutionXy = single([10; 10]); % reference value: [100; 100]

% if strcmp(kindOfMotion,'periodic')
    isAddNoise = 0;
    if isAddNoise
    end
    mainHeartFreq = 10; % reference: 10
    maxSpeedParticle = 0.004; % reference value: 0.001
    minSpeedParticle = 0.001;
    
    exposureTimeDivision = 0.001; % seconds, reference value = 0.005
    Fs = 200; % Sampling frequency (samples per second) 
    dt = (1/Fs) * exposureTimeDivision; % seconds per sample 
    nPeriods = 1; 
    tRecording = nPeriods / mainHeartFreq; % recording time in seconds
    t = (0:dt:tRecording)'; % seconds 
    signal = sawtooth (2*pi*mainHeartFreq*t,1/4);
    signal = abs(min(signal)) + signal;
    signal = signal / max(signal);

    signalSharp = sawtooth (2*pi*mainHeartFreq*t,0.5);
    signalSharp = abs(min(signalSharp)) + signalSharp;
    signalSharp = signalSharp / max(signalSharp);

    signalMask = square (2*pi*mainHeartFreq*t,60);
    signalMask = abs(min(signalMask)) + signalMask;
    signalMask = signalMask / max(signalMask);

    signalSharp = signalSharp .* signalMask;

    samplesToConsider = 20000;

    findThresholdSharpEdge = signalSharp(1:samplesToConsider);
    signalSharpFindZeros = find(findThresholdSharpEdge==0);
    thresholdSharpEdge = signalSharp(signalSharpFindZeros(2)-1);
    
    signalSharp = signalSharp - thresholdSharpEdge;
    signalSharp(signalSharp<0) = 0;
    
% figure,plot(signalSharp)
% figure,plot(signalSharp(1:samplesToConsider)), hold on
% plot(signalMask(1:samplesToConsider))

    signalDephased = zeros(1,length(signalSharp));

    signalDephasedAngle = 20; %angle in degrees
    signalDephasedAmplitude = 1.5; % bump on top of sawtooth signal

    samplesPerPeriod = length(signalSharp)/nPeriods;
    samplesDephased = floor((samplesPerPeriod/360)*signalDephasedAngle);
    signalDephased(samplesDephased+1:end) = signalDephasedAmplitude * ...
                                          signalSharp(1:end-samplesDephased);
    signalTwo = signal + signalDephased';
    
    filterPercentage = 10;
    medfiltCoeff = round(samplesPerPeriod/filterPercentage);

    factorFilter = 1;
    filterWindow = ones(1,medfiltCoeff)/factorFilter;
    
    signalTwo = filter(filterWindow, 1, signalTwo);


    
    periodicMovement = minSpeedParticle + ((maxSpeedParticle - minSpeedParticle) * signalTwo); % setting offset and scale
    nMoves = length(t);
% end

% particle movement
% nMoves = 20; % reference value: 1000
maxSpeedParticle = 0.001; % reference value: 0.001
directionParticleDegreesConstant = 0; % reference value: 0
% kindOfMotion = 'periodic'; % options: 'periodic' 'ordered' 'brownian' 'random'

% figure,plot(signal)

% %%%%%%%% PLOTS / SAVING%%%%%%%%%%
% save the signal
isSaveSignal = uint8(0);
saveFolder = './../data/simulation/';

% plot the figures
isPlot = uint8(0);
%     isSavePlot = uint8(0); % to implement


% create gif
isGif = uint8(0);
%     isSaveGif = uint8(0); % to implement

% Calculating Time Intensity Autocorrelation Function (tiaf)
isCalculateTiaf = uint8(0);
%     isSaveCalculatedTiaf = uint8(0); % to implement
        maxTau = single(50);
        
% To calibrate the camera position and or particles volume, etc use the
% following visualizer.
%   visualizing particles with respect to sensors:
isVisualizeParticlesSensors = uint8(0);


%% generating randomly distributed particles in a volume
% nParticles = 100;
% particlesVolumeSizeXyz = [0.5; 0.5; 0.1];
particlesPositionXyz = [rand(nParticles, 1)*particlesVolumeSizeXyz(1),...
                        rand(nParticles, 1)*particlesVolumeSizeXyz(2),...
                        -rand(nParticles, 1)*particlesVolumeSizeXyz(3)];

else
    particlesPositionXyz = lastFrameParticles;
end


%% defining sensor
% sensorDistanceToZ = 10;
% sensorSizeXy = [10;10];
% sensorResolutionXy = [20;20];

sensorCentreXyz = [particlesVolumeSizeXyz(1)/2;...
                   particlesVolumeSizeXyz(2)/2;...
                   sensorDistanceToZ];

sensorPixelsCoordinatesX = linspace(sensorCentreXyz(1)-sensorSizeXy(1),sensorCentreXyz(1)+sensorSizeXy(1),sensorResolutionXy(1));
sensorPixelsCoordinatesY = linspace(sensorCentreXyz(2)-sensorSizeXy(2),sensorCentreXyz(2)+sensorSizeXy(2),sensorResolutionXy(2));

sensorPixelsCoordinates_X_RepeatedYTimes = repmat(sensorPixelsCoordinatesX,[sensorResolutionXy(2) 1]);
sensorPixelsCoordinates_X_RepeatedYTimesReshaped = reshape(sensorPixelsCoordinates_X_RepeatedYTimes,[1 sensorResolutionXy(1)*sensorResolutionXy(2)]);
clear sensorPixelsCoordinates_X_RepeatedYTimes

sensorPixelsCoordinates_Y_RepeatedXTimes = repmat(sensorPixelsCoordinatesY,[1 sensorResolutionXy(1)]);
clear sensorPixelsCoordinatesX sensorPixelsCoordinatesY

sensorPixelsCoordinatesXyz = [sensorPixelsCoordinates_X_RepeatedYTimesReshaped;...
                              sensorPixelsCoordinates_Y_RepeatedXTimes;...
                              repmat(sensorCentreXyz(3),[1,sensorResolutionXy(1)*sensorResolutionXy(2)])]';
clear sensorPixelsCoordinates_X_RepeatedYTimesReshaped sensorPixelsCoordinates_Y_RepeatedXTimes


%% visualizing particles with respect to sensors
if isVisualizeParticlesSensors
    scatter3(sensorPixelsCoordinatesXyz(:,1),...
             sensorPixelsCoordinatesXyz(:,2),...
             sensorPixelsCoordinatesXyz(:,3));
    hold on;
    scatter3(particlesPositionXyz(:,1),...
             particlesPositionXyz(:,2),...
             particlesPositionXyz(:,3),'r');

    hold off;
    pause;
end


%% particle moving
% nMoves = 10000;
% 
% maxSpeedParticle = 0.01;
% directionParticleDegreesConstant = 0;
% isOrderedMotion = 1;

% the corners are (X,Y):
%   bottom left     = 0, 0
%   bottom right    = particlesVolumeSizeXyz(1), 0
%   upper left      = 0, particlesVolumeSizeXyz(2)
%   upper right     = particlesVolumeSizeXyz(1), particlesVolumeSizeXyz(2)
particlesFrameCorners = [0,0;...
                         particlesVolumeSizeXyz(1), 0;...
                         0, particlesVolumeSizeXyz(2);...
                         particlesVolumeSizeXyz(1), particlesVolumeSizeXyz(2)];
if strcmp(kindOfMotion,'periodic')
    nSamples = tRecording * Fs;
    exposureSamples = floor(nMoves / nSamples);
    exposureFrames = nan(sensorResolutionXy(1),sensorResolutionXy(2),exposureSamples);
    exposureCounter = 0;
    timeIntensityAutocorrelationFunction = nan(sensorResolutionXy(1),sensorResolutionXy(2),nSamples);
    timeIntensityCounter = 0;
else
    timeIntensityAutocorrelationFunction = nan(sensorResolutionXy(1),sensorResolutionXy(2),nMoves);
end

% euclideanDistanceParticlesToPixelsSumOfSquares = nan(nParticles,sensorResolutionXy(1)*sensorResolutionXy(2),nMoves,'single');
blockOfParticles = nan(nParticles,3,nMoves,'single');

wb = waitbar(0,'Please wait...');
for iter = 1:nMoves
%update waitbar
waitbar(single(iter)/single(nMoves),wb,['Processing simulation... iteration: ',num2str(single(iter)),' of ',num2str(single(nMoves))]);
    switch kindOfMotion
        case 'periodic'
            speedParticle = periodicMovement(iter);
            directionParticleDegrees = directionParticleDegreesConstant;
        case 'ordered'
            speedParticle = maxSpeedParticle;
            directionParticleDegrees = directionParticleDegreesConstant;
        case 'brownian'
            speedParticle = rand(1,nParticles)*maxSpeedParticle;
            directionParticleDegrees = randi([0 360],1,nParticles);
        case 'random'
            particlesPositionXyz = [rand(nParticles, 1)*particlesVolumeSizeXyz(1),...
                                    rand(nParticles, 1)*particlesVolumeSizeXyz(2),...
                                    -rand(nParticles, 1)*particlesVolumeSizeXyz(3)];
            speedParticle = 0;
            directionParticleDegrees = 0;
    end
        
    particlesPositionXyz(:,1) = particlesPositionXyz(:,1)' + ...
                                speedParticle .* cosd(directionParticleDegrees);
    particlesPositionXyz(:,2) = particlesPositionXyz(:,2)' + ...
                                speedParticle .* sind(directionParticleDegrees);
                                
    %% if particle is out of scope bring it to the begining of the frame
    %       in the same trajectory and direction
    
    % when particles are out of the X axis:
    % - from the side of the max volume side in the X axis
    
    if sum(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(1)) > 0
        particlesPositionXyz(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(1), 2) = ...
                 rand(1, sum(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(1))) * particlesVolumeSizeXyz(1);
        particlesPositionXyz(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(1), 1) = 0;
%                              particlesPositionXyz(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(1), 1) - ...
%                              cosd(directionParticleDegrees) * particlesVolumeSizeXyz(1);
    end
    % - below 0
    if sum(particlesPositionXyz(:,1) < 0) > 0
        particlesPositionXyz(particlesPositionXyz(:,1) < 0, 2) = ...
                 rand(1, sum(particlesPositionXyz(:,1) < 0)) * particlesVolumeSizeXyz(1);
        particlesPositionXyz(particlesPositionXyz(:,1) < 0, 1) = particlesVolumeSizeXyz(2); 
%                              particlesPositionXyz(particlesPositionXyz(:,1) < 0, 1) - ...
%                              cosd(directionParticleDegrees) * particlesVolumeSizeXyz(1);
    end
    
    % - from the side of the max volume side in the Y axis
    if sum(particlesPositionXyz(:,2) > particlesVolumeSizeXyz(2)) > 0
        particlesPositionXyz(particlesPositionXyz(:,2) > particlesVolumeSizeXyz(2), 2) = ...
                 rand(1, sum(particlesPositionXyz(:,2) > particlesVolumeSizeXyz(2))) * particlesVolumeSizeXyz(2);
        particlesPositionXyz(particlesPositionXyz(:,2) > particlesVolumeSizeXyz(2), 1) = 0; 
%                              particlesPositionXyz(particlesPositionXyz(:,1) > particlesVolumeSizeXyz(2), 2) - ...
%                              cosd(directionParticleDegrees) * particlesVolumeSizeXyz(2);
    end
    % - below 0
    if sum(particlesPositionXyz(:,2) < 0) > 0
        particlesPositionXyz(particlesPositionXyz(:,2) < 0, 2) = ...
                 rand(1, sum(particlesPositionXyz(:,2) < 0)) * particlesVolumeSizeXyz(2);
        particlesPositionXyz(particlesPositionXyz(:,2) < 0, 1) = particlesVolumeSizeXyz(1); 
%                              particlesPositionXyz(particlesPositionXyz(:,2) < 0, 1) - ...
%                              cosd(directionParticleDegrees) * particlesVolumeSizeXyz(2);
    end

    blockOfParticles(:,:,iter) = particlesPositionXyz;
                                
%     %% distance from particle to pixel
%     particlesAugmented = repmat(particlesPositionXyz,[1 1 sensorResolutionXy(1)*sensorResolutionXy(2)]);
%     % clear particlesPositionXyz
%     sensorPixelsCoordinatesXyzAugmented = repmat(sensorPixelsCoordinatesXyz',[1 1 nParticles]);
% %     clear sensorPixelsCoordinatesXyz 
% 
%     euclideanDistanceParticlesToPixelsDifference = bsxfun(@minus,...
%                                                 particlesAugmented,...
%                                                 permute(sensorPixelsCoordinatesXyzAugmented,[3 1 2]));
% %     clear particlesAugmented sensorPixelsCoordinatesXyzAugmented 
%     euclideanDistanceParticlesToPixelsSumOfSquares(:,:,iter) = squeeze(sum(euclideanDistanceParticlesToPixelsDifference.^2,2));
% %     clear euclideanDistanceParticlesToPixelsDifference


%     %% light scattered from particles for every pixel
%     % E=exp(i*k*r -1i*k*c*t)./r
%     %
%     % i - imaginary unit, r - distance from particle to pixel, 
%     % c - speed of light, k - wavenumber (2pi/lambda), wavelength lambda=0.785 um, time t=0
%     c = single(299792458000000); % um/s
%     lambda = single(0.785); % microns
%     k = 2*pi/lambda;
%     t = 0; %single(1);
%     r = euclideanDistanceParticlesToPixelsSumOfSquares;
% %     clear euclideanDistanceParticlesToPixelsSumOfSquares
%     
%     E=exp(1i*k*r -1i*k*c*t)./r;
% %     clear r
% 
%     superpositionOfScatteredLight = sum(E,1);
% %     clear E
%     pixelsIntensity = superpositionOfScatteredLight .* conj(superpositionOfScatteredLight);
% %     clear superpositionOfScatteredLight
%     pixelsIntensity = pixelsIntensity/ mean(pixelsIntensity);


%     %% Reshaping to make sensor image
%     sensorImage = reshape(pixelsIntensity,[sensorResolutionXy(1) sensorResolutionXy(2)]);
% %     clear pixelsIntensity

%     %% Save sensorImage
%     % 
%     if strcmp(kindOfMotion,'periodic')
%         exposureCounter = exposureCounter + 1;
%         exposureFrames(:,:,exposureCounter) = sensorImage;
%         modIter = ceil(iter/ exposureSamples) - floor(iter/ exposureSamples);
%         if modIter == 0
%             exposureCounter = 0;
%             timeIntensityCounter = timeIntensityCounter + 1;
%             timeIntensityAutocorrelationFunction(:,:,timeIntensityCounter) = mean(exposureFrames,3);
%         end
%     else
%         timeIntensityAutocorrelationFunction(:,:,iter) = sensorImage;
%     end
    
    
%     %% plot the figures
% %     isPlot = 0;
%     if isPlot
%         h = figure;
%         screenSize = get(0, 'ScreenSize');
%         set(gcf, 'Position',  [floor(screenSize(3)/9), floor(screenSize(4)/5), ...
%             floor(screenSize(3)/9)*7, floor(screenSize(4)/5)*3])
%         subplot(1,2,1),
%         scatter(particlesPositionXyz(:,1),particlesPositionXyz(:,2)),
%         axisXmin = 0;
%         axisXmax = particlesVolumeSizeXyz(1);
%         axisYmin = 0;
%         axisYmax = particlesVolumeSizeXyz(2);
%         axis([axisXmin axisXmax axisYmin axisYmax]),
%         subplot(1,2,2),
%         imagesc(sensorImage)
%         % caxis([prctile(sensorImage(:),0.1),prctile(sensorImage(:),99.9)])
%         pause;
%         close all;
%     end

%     %% calculate autocorrelation
%     isAutocorrelation = 0;
%     if isAutocorrelation
% 
%         variableContrastKernelAnalysis;
% 
% %         keyboard;
%     end
%     %% create gif
%     if isGif
%           % Capture the plot as an image 
% 
%         fileName = [savePath,kindOfMotion,'_lscigif.gif'];
%         h = figure('Visible','off');
% 
%         imagesc(sensorImage)
%         frame = getframe(h); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         if iter == 1 
%           imwrite(imind,cm,fileName,'gif', 'Loopcount',inf); 
%         else 
%           imwrite(imind,cm,fileName,'gif','WriteMode','append'); 
%         end 
%           
%     end
    
%     clear sensorImage 
    
%     if isCalculateTiaf || isGif
%         display(['iteration: ',num2str(iter),' of ',num2str(nMoves)])
%     end
end
close all
close(wb);
% keyboard;

%% Saving the signal 
if isSaveSignal

    fileName = [kindOfMotion,'_OnlyParticlesSignal_',num2str(sensorResolutionXy(1)),'_by_',num2str(sensorResolutionXy(2)),'_fs_',num2str(Fs),'_per_',num2str(nPeriods)];
    dirFile = dir([saveFolder,fileName,'*']);
    if ~isempty(dirFile)
        if length(num2str(length(dirFile)+1))<2
            numberWithZero = ['0',num2str(length(dirFile)+1)];
        else
            numberWithZero = num2str(length(dirFile)+1);
        end
        fileName = [fileName,'_',numberWithZero];
    end
    fileName = [fileName,'.mat'];
    save([saveFolder,fileName],'euclideanDistanceParticlesToPixelsSumOfSquares','exposureTimeDivision','Fs','nPeriods','mainHeartFreq','periodicMovement','-mat','-v7.3');
end

    
     

%% Calculating Time Intensity Autocorrelation Function (tiaf)
% isCalculateTiaf = 1;

if isCalculateTiaf

%     clear sensorPixelsCoordinatesXyz
%     maxTau = 100;
    tiafIntensityTimesIntensityDelayed = nan(sensorResolutionXy(1),...
                                             sensorResolutionXy(2),...
                                             nMoves-maxTau);
    tiaf = nan(1,maxTau);
    tiafMeanSquaredTime = mean(timeIntensityAutocorrelationFunction,3).^2;
    for iTau = 1:maxTau
        tiafIntensityTimesIntensityDelayed = mean(...
            timeIntensityAutocorrelationFunction(:,:,1:end-maxTau) .* ...
            timeIntensityAutocorrelationFunction(:,:,iTau:end-maxTau+iTau-1)...
            ,3);
        
        tiafDivision = tiafIntensityTimesIntensityDelayed ./ tiafMeanSquaredTime;
        tiaf(iTau) = mean(mean(tiafDivision,1),2);
    end

    decorrelationLag = linspace(1,maxTau,maxTau);
    figure,plot(decorrelationLag, tiaf)

    nameOrder = kindOfMotion;
    
    nameTitle = [nameOrder,' motion'];
    title(nameTitle)
    xlabel('decorrelation time \tau_{c}')
    ylabel('Time Intensity Autocorrelation')
    
    if isSaveCalculatedTiaf

        
    end
end




