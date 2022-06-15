Fs = 194;

%See read RLS file for detailed description of arguments
batchSize = 100; 
procType = 'fastcpu';
skipFrames = 0;
movMeanKernel = 1;
downSampling = 1; % set to 1 for no downsampling
ROI = [];

data = runBatchLSCI(fileName, 'temporal', 25, batchSize, procType, ...
                    skipFrames, movMeanKernel, downSampling, ROI);



[data,~,~,~]=readRLS(fileName);
data = single(data(:,:,1:1000));
dataSTD = std(data,[],3);
dataMEAN = mean(data,3);
dataLSCI_temporalContrast = dataSTD./dataMEAN;

[bloodFlow,parenchyma]=selectAreasToCompare(dataLSCI_temporalContrast);

[dataBloodFlow,~,~,~]=readRLS(fileName,0,[],bloodFlow);
[dataParenchyma,~,~,~]=readRLS(fileName,0,[],parenchyma);



dataSLSCI = getSLSCI(dataBloodFlow,7,procType);

figure,imagesc(dataSLSCI(:,:,1))

bfi = double(squeeze(mean(dataSLSCI,[1 2])));
bfi = 1./bfi.^2;
figure,plot(bfi)

    figure
    subplot(2,1,1)
    title('histogram inside vessel')
    histogram(dataBloodFlow(1,1,:),'BinLimits',[40,100])
    subplot(2,1,2)
    title('histogram parenchyma')
    histogram(dataParenchyma(1,1,:),'BinLimits',[40,100])

figure,histogram(dataBloodFlow,'r')
hold on
histogram(dataParenchyma,'b')
hold off


loopHist = 0;
if loopHist
    h = figure;
    for iter = 1:size(dataParenchyma,3)
        figure(h)
        subplot(2,1,1)
        histogram(dataBloodFlow(:,:,iter),'BinLimits',[40,150])
        subplot(2,1,2)
        histogram(dataParenchyma(:,:,iter),'BinLimits',[40,150])
    end
end