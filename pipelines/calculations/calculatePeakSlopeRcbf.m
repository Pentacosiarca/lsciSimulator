function [framePeak, slopeAngle, rCBF] = calculatePeakSlopeRcbf (bfiTS, sLsciBaseline, sLsciStartStim)

medWindow = 2001; % median window
polyfitDegree = 5; % degree of the polynomial
medfiltBfiTs = medfilt1(bfiTS,medWindow); % median filtering of bfi signal
curvePoly = polyfit((1:1:length(medfiltBfiTs))',medfiltBfiTs,polyfitDegree);
curveFitting = polyval(curvePoly,(1:1:length(medfiltBfiTs))');


%% finding the peak of the curve
framePeak = find(curveFitting == max(curveFitting));
framePeak = framePeak(1);


%% calculating the curve's slope's angle
curvePeak = curveFitting(1:framePeak); % curve before the peak
curveDeriv = curvePeak(2:end) - curvePeak(1:end-1); % derivative of the curve
minCurveDeriv = find(curveDeriv == max(curveDeriv));
medianMinCurveDeriv = median(minCurveDeriv);
sizeWindowGradient = 200;
windowGradient = curvePeak(medianMinCurveDeriv - sizeWindowGradient : medianMinCurveDeriv + sizeWindowGradient - 1);
A = [sizeWindowGradient,windowGradient(end)-windowGradient(1)];
B = [sizeWindowGradient,0];
slopeAngle = acosd(cos((A.*B)/(abs(A).*abs(B))));


%% check fitting results
isCheckFit = 0;
if isCheckFit
    figure,
    hold on,
    plot(bfiTS,'Color',[0.8,0.8,0.8])
    plot(medfiltBfiTs,'Color',[0.2,0.2,0.2],'LineWidth',3)
    plot(curveFitting,'LineWidth',2)
    title(['polyfit degree: ',num2str(polyfitDegree)])
    legend('blood flow index signal from spatial contrast',...
            'median averaged signal',...
            ['curve fitting with polyfit, degree: ',num2str(polyfitDegree)])
    hold off,
end


%% Computing rCBF
% get most stable sequence of frames in the baseline

bfiTSbase=squeeze(mean(sLsciBaseline,[1,2]));
bfiTSbase=1./(bfiTSbase.^2);
baselineBfi = bfiTSbase; 
baselineMedianDiff = abs(baselineBfi - median(baselineBfi));
baselineFindMedianPoints = find(baselineMedianDiff < median(baselineMedianDiff));
baselineDiffFindMedianPoints = baselineFindMedianPoints(2:end) - baselineFindMedianPoints(1:end-1);
baselineDFMPmovMedian = medfilt1(baselineDiffFindMedianPoints,3);
baselineDFMPMMZeros = zeros(1,length(baselineDFMPmovMedian));
baselineDFMPMMZeros(baselineDFMPmovMedian < 3) = 1;
baselineDiffZeros=reshape(find(diff([0;baselineDFMPMMZeros';0])~=0),2,[]);
[lgtmax,jmax]=max(diff(baselineDiffZeros));
istart=baselineDiffZeros(1,jmax);
baselineFramesIdx = baselineFindMedianPoints(istart:istart+lgtmax);


%% plot which frames considered in the baseline
isShowFramesBaseline = 0;
if isShowFramesBaseline
figure,plot(baselineBfi), hold on,
scatter(baselineFindMedianPoints(istart:istart+lgtmax),baselineBfi(baselineFindMedianPoints(istart:istart+lgtmax))),hold off
hold off,
end


% get most stable sequence of frames around the peak
peakFramesIdx = framePeak - round(length(baselineFramesIdx)/2) : framePeak + round(length(baselineFramesIdx)/2) - 1;

isRawFilesTemporalContrast = 0;
if isRawFilesTemporalContrast
    [rawDataBaselineStableFrames,~,~,~] = readRLS(fileName,baselineFramesIdx(1),baselineFramesIdx(end));
    peakFramesFromBaseline = peakFramesIdx + baselineFramesIdx(end);
    [rawDataPeakFrames,~,~,~] = readRLS(fileName,peakFramesFromBaseline(1),peakFramesFromBaseline(end));
    tLSCIbaseline=getTLSCI(rawDataBaselineStableFrames,25,'gpu');
    tLSCIpeak=getTLSCI(rawDataPeakFrames,25,'gpu');
    bfiTlsciBaseline = mean(tLSCIbaseline,3);
    bfiTlsciBaseline = squeeze(1./(bfiTlsciBaseline.^2));
    bfiTlsciPeak = mean(tLSCIpeak,3);
    bfiTlsciPeak = squeeze(1./(bfiTlsciPeak.^2));
    rCBF =  bfiTlsciPeak ./ bfiTlsciBaseline;
%     figure,imagesc(bfiTlsciPeak),title('tLSCI peak')
%     figure,imagesc(bfiTlsciBaseline),title('tLSCI baseline')
%     figure,imagesc(rCBF);
else
    baselineAverageFrames = mean(sLsciBaseline(:,:,baselineFramesIdx),3);
    peakAverageFrames = mean(sLsciStartStim(:,:,peakFramesIdx),3);
    rCBF = peakAverageFrames ./ baselineAverageFrames;
end









