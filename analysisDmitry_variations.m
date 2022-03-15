clear all, clc, close all
addpath(genpath('./../'))

loadPath = 'C:\Users\Alberto\Desktop\github\data\simulation\';
load([loadPath,'periodic_signal_10_by_10_fs_200_per_50'])

% load([loadPath,'periodic_signal_10_by_10_fs_200_per_20_02.mat'])

dataLsci = timeIntensityAutocorrelationFunction;
clear timeIntensityAutocorrelationFunction

% dataLsci = dataLsci(:,:,1:3:end);
downsampleRatio = floor(length(periodicMovement)/size(dataLsci,3));
periodicDownSampled = downsample(periodicMovement,downsampleRatio);
periodicDownSampled(size(dataLsci,3)+1:end) = [];

x = 1:length(periodicDownSampled);
[~,loc] = findpeaks(-periodicDownSampled);

% locMinMedian = findSignalMinima (periodicDownSampled, 200);

locMinima = islocalmax(-periodicDownSampled);
locMinima(1) = 1;
locMinima(end) = 1;

isException = 1;
if isException
deleteLocMin = 11:20:loc(end);
locMinima(deleteLocMin) = 0;
end

figure,plot(x,periodicDownSampled,x(locMinima),periodicDownSampled(locMinima),'r*')

positionLocMin = find(locMinima>0);

sizeMedianPulse = median(positionLocMin(2:end)-positionLocMin(1:end-1));

lsLSCI=[];

kernel = length(positionLocMin)-1;


% for i=1:1:20
for i=1:1:sizeMedianPulse

% subdata=dataLsci(:,:,i:20:end);
% lsLSCI(:,:,i)=squeeze(std(subdata,0,3)./mean(subdata,3));

subdata=getTLSCI(dataLsci(:,:,i:20:end),kernel,'gpu');
lsLSCI(:,:,i)=mean(subdata,3);

end

tLSCI25=getTLSCI(dataLsci,5,'gpu');

sLSCI25=getSLSCI(dataLsci,5,'gpu');

tLSCI=zeros(10,10,20);

sLSCI=zeros(10,10,20);

for i=1:1:20

tLSCI(:,:,i)=mean(tLSCI25(:,:,i:20:end),3);

sLSCI(:,:,i)=mean(sLSCI25(:,:,i:20:end),3);

end

ts1=1./(squeeze(mean(lsLSCI,[1,2])).^2);

ts2=1./(squeeze(mean(tLSCI,[1,2])).^2);

ts3=1./(squeeze(mean(sLSCI,[1,2])).^2);

ts1=(ts1-min(ts1))./(max(ts1)-min(ts1));

ts2=(ts2-min(ts2))./(max(ts2)-min(ts2));

ts3=(ts3-min(ts3))./(max(ts3)-min(ts3));

origSignal = periodicDownSampled(1:20);
origSignal=(origSignal-min(origSignal))./(max(origSignal)-min(origSignal));

figure,
hold on


plot(origSignal,'b')

plot(ts1,'r')
plot(ts2,'g')
plot(ts3,'m')

hold off
legend('original','lossless','temporal contrast','spatial contrast')



