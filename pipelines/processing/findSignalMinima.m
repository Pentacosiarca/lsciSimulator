%findSignalMinima - find the local minima of the data after calculating
%blood flow
%
% Syntax:  output1 = function_name(input1,input2)
%
% Inputs:
%    data           - blood flow signal
%    Fs             - sample frequency of the data
%
% Outputs:
%    locMinPrctile  - position of the local minima after discarding
%                     abnormal pulses
%
% Other m-files required: none 
% Subfunctions: none
% MAT-files required: none

% Authors: Alberto Gonzalez Olmos, DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 17-November-2021

%------------- BEGIN CODE --------------

function locMinMedian = findSignalMinima (data, Fs)

%% parameters             
% T = 1/Fs;             % Sampling period       
L = length(data);     % Length of signal
% t = (0:L-1)*T;        % Time vector

X = data - mean(data);
Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

heartRateRangeOfInterest = [5 14]; % in Hz
[~,findLow] = min(abs(f - heartRateRangeOfInterest(1)));
[~,findHigh] = min(abs(f - heartRateRangeOfInterest(2)));
freqOfInterest = P1 (findLow:findHigh);
freqGaussian = fit((1:1:length(freqOfInterest))',freqOfInterest,'gauss1');
freqPeak = round(freqGaussian.b1);
freqHalfWidth = round(freqGaussian.c1/2);
freqMainEnd   = freqPeak + freqHalfWidth;
freqMainMax = f (findLow + freqMainEnd);
freqNumWindow = ceil((1/freqMainMax)*Fs);
freqWindow = gausswin(freqNumWindow);

%% rough estimate of max peak
convValue = zeros(1, L - freqNumWindow);
iConv = 1;
while iConv<(L - freqNumWindow)
   sectionSignal = data(iConv:iConv+freqNumWindow-1); 
   sectionSignalNorm = (sectionSignal - min(sectionSignal)) ./ max(sectionSignal - min(sectionSignal));
   if max(sectionSignal) > max(sectionSignal(1),sectionSignal(end))

       [~,idxMax] = max(sectionSignal);
       convValue(iConv+idxMax-1) = sum(sectionSignalNorm.*freqWindow);
       iConv = iConv+ freqNumWindow ; %idxMax;
   else
       iConv = iConv + 1;
   end
end

convValueIdx = find(convValue>0);

%% find local minima
iPeaks = 1;
locMin = zeros (1, length(convValueIdx)-1);
countLocMin = 1;
while iPeaks < (length(convValueIdx)-1)
    sectionSignal = data(convValueIdx(iPeaks):convValueIdx(iPeaks+1));
    if (length(sectionSignal) > (0.5 * freqNumWindow)) && ...
       (length(sectionSignal) < (1.5 * freqNumWindow))
        
        [~,idxMin] = min(sectionSignal);
        locMin(countLocMin) = convValueIdx(iPeaks) + idxMin - 1;
        countLocMin = countLocMin + 1;
        iPeaks = iPeaks+1;
    else
        iPeaks = iPeaks+1;
    end
end
locMin(locMin==0) = [];

% figure,hold on,
% plot(data),
% scatter(locMin(1,:),data(locMin(1,:))),
% title('First estimate of minima based on window'),
% hold off,

%% Score pulses based on max / min values
pulsesScoreMaxMin = zeros(1,length(locMin)-1);
pulsesScoreMinMin = zeros(1,length(locMin)-1);
for iCheckPulses = 1:1:length(locMin)-1
    interval = data(locMin(iCheckPulses):locMin(iCheckPulses+1));
    pulsesScoreMaxMin(iCheckPulses) = max(interval) - min(interval);
    pulsesScoreMinMin(iCheckPulses) = abs(interval(1)-interval(end));
end


% figure,hold on,
% plot(data),
% scatter(locMin(1,:),data(locMin(1,:))),
% scatter(locMin(1,1:end-1),pulsesScoreMaxMin,'.'),
% scatter(locMin(1,1:end-1),pulsesScoreMinMin,'.'),
% legend('signal','local Minima','difference Max Min','difference Min Min')
% title('Score of Minima based on difference of Max, Min and Min, Min'),
% hold off,

%% Find pulses with scores within 2 times the median

medianMaxMin = median(pulsesScoreMaxMin);
medianMinMin = median(pulsesScoreMinMin);

highMedianMaxMin = medianMaxMin * 1.5;
highMedianMinMin = medianMinMin * 2;

normalPulses = find((pulsesScoreMaxMin < highMedianMaxMin) & ...
                    (pulsesScoreMinMin < highMedianMinMin));


locMinMedian = zeros(2,length(normalPulses)+1);

iPrctiles = 1;
while iPrctiles<length(normalPulses)
    if ((locMin(normalPulses(iPrctiles)+1) - locMin(normalPulses(iPrctiles))) < 1.5 * freqNumWindow) && ...
      ~(locMinMedian(iPrctiles) == locMin(normalPulses(iPrctiles)))

        locMinMedian(1,iPrctiles)    = locMin(normalPulses(iPrctiles));
        locMinMedian(2,iPrctiles)  = locMin(normalPulses(iPrctiles)+1);
    end
    iPrctiles = iPrctiles+1;
end

idxZero = find(locMinMedian(1,:)==0);
locMinMedian(:,idxZero) = [];

figure,hold on,
plot(data),
scatter(locMinMedian(1,:),data(locMinMedian(1,:)),'or'),
scatter(locMinMedian(2,:),data(locMinMedian(2,:)),'.b'),
legend('signal','Start Pulse, local minima','End Pulse, local minima')
title('Score of Minima based on difference of Max, Min and Min, Min'),
hold off,









