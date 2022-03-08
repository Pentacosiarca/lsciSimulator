
% clear all, clc, close all
% load('sensorImage.mat')

kernelSizes = 1:20;

% count = 0;
for iKernels = 1:length(kernelSizes)
kernel = kernelSizes(iKernels);
filteredImage = boxFiltDown2D (sensorImage, kernel, 'average');

xcorrSensorImage = xcorr(filteredImage(50,:),'unbiased');
centeredSignal = xcorrSensorImage(round(length(xcorrSensorImage)/2)-10:round(length(xcorrSensorImage)/2)+10);
[pks,~,widths,~] = findpeaks(centeredSignal,1:length(centeredSignal));
maxPeakIdx = find(pks == max(pks));

speckleSize = widths(maxPeakIdx(1))/2;
K = std(filteredImage,1,'all') / mean(filteredImage,[1 2]);

% figure,plot(centeredSignal)
% figure,imagesc(filteredImage)

% 
% keyboard;
disp(['K: ',num2str(K),', Speckle Size: ',num2str(speckleSize)])
if exist("speckleVScontrastK",'var')
    count = count + 1;
else
    if exist("speckleVScontrastK.mat",'file') == 2
        load ("speckleVScontrastK")
        count = count + 1;
    else
        count = 1;
    end
end

speckleVScontrastK (:,count) = [K, speckleSize];

end

save('speckleVScontrastK','speckleVScontrastK','count')

figure,scatter(speckleVScontrastK(2,:),speckleVScontrastK(1,:)),
xlabel('Speckle Size')
ylabel('K')

% close(wb)

