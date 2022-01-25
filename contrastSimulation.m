clear all
%close all
tic
%Computation/visualisation settings
rngKernel=110;
rng(rngKernel);
observationTime=100000; %mus
gpuBlockFrames=min(1000,observationTime); %frames per single GPU block. Reduce if getting out of memory error
%reset(gpuDevice(1));
storeXYZ=0; %do not use for large particlesN, or long observations
storeRawI=1; % only use when validating, not recommended for long simulations

animationDt=10;

%light source
lambda = 0.785; %mum
k=2*pi/lambda; %wavenumber
c=3*10^8; % mum/mus

%dynamics
dT=1; %time between instanteneous frames, mus
T=5000; %exposure time of the output frame, mus
tauc=50; % mus
motionType='ordered'; % 'unordered' or 'ordered'
switch motionType
    case 'unordered'
        D= 1/((tauc/dT)*k^2); % Diff coef
        d = sqrt(2*D); % random walk step length
    case 'ordered'
        d=1/((tauc/dT)/(lambda.^2));
end
d=ones(1,floor(observationTime/dT))*d;

%introduce signal to d here
%e.g. d=d+d*0.1*sin(1:1:length(d));


%volume
sizeX = 100; %mum
sizeY = 100; %mum
sizeZ = 10; %mum

%particles
particlesN = 50;
brN=floor(particlesN*1); % N of "brownian" particles
particleX=gpuArray(zeros(1,1,particlesN,'single'));
particleY=gpuArray(zeros(1,1,particlesN,'single'));
particleZ=gpuArray(zeros(1,1,particlesN,'single'));

%initial random positions
particleX(1,1,:) = rand(particlesN,1)*sizeX;
particleY(1,1,:) = rand(particlesN,1)*sizeY;
particleZ(1,1,:) = rand(particlesN,1)*sizeZ;

if storeXYZ==1
    storedXYZ=zeros(3,particlesN,observationTime);
end

%sensor
pixelsNx = 100;
pixelsNy = 100;
pixelSize = 0.5; %mum
sensorX = sizeX/2;
sensorY = sizeY/2;
sensorZ = sizeZ*10;
pixelPosX =gpuArray(single( sensorX + ones(pixelsNx,1)*((1:pixelsNx)-0.5*pixelsNx).*pixelSize));
pixelPosY =gpuArray(single( sensorY + ((1:pixelsNy)'-0.5*pixelsNy)*ones(1,pixelsNy).*pixelSize));

%pre-allocate memory for intensity values
expI=zeros(pixelsNx,pixelsNy,floor(observationTime/T),'single');
gpuExpI=gpuArray(zeros(pixelsNx,pixelsNy,'single'));
if storeRawI==1
    rawI=zeros(pixelsNx,pixelsNy,observationTime/dT,'single');
    subI=gpuArray(zeros(pixelsNx,pixelsNy,gpuBlockFrames,'single'));
end

%calculate field
tic
for i=1:1:floor(observationTime/dT)
    %Store particles position for visualisation
    if storeXYZ==1
        storedXYZ(1,:,i)=gather(particleX);
        storedXYZ(2,:,i)=gather(particleY);
        storedXYZ(3,:,i)=gather(particleZ);
    end
    
    r = sqrt((particleX - pixelPosX).^2 + (particleY - pixelPosY).^2 + (particleZ - sensorZ).^2);
    tk=exp(-1i.*k.*c.*((i-1)*dT)); %set to 1 if precision is not high enough
    E = exp(1i.*k.*r).*tk./r;
    E=sum(E,3);
    
    switch motionType
        case 'unordered'
            particleX(1:brN) = particleX(1:brN)+d(i)*randn(1,1,brN);
            particleY(1:brN) = particleY(1:brN)+d(i)*randn(1,1,brN);
            particleZ(1:brN) = particleZ(1:brN)+d(i)*randn(1,1,brN);
            particleX(1:brN)=mod(particleX(1:brN),sizeX);
            particleY(1:brN)=mod(particleY(1:brN),sizeY);
            particleZ(1:brN)=mod(particleZ(1:brN),sizeZ);
        case 'ordered'
            particleX(1:brN) = particleX(1:brN)+d(i);
            particleY(particleX>sizeX) = rand(sum(particleX>sizeX),1)*sizeY;
            particleZ(particleX>sizeX) =rand(sum(particleX>sizeX),1)*sizeZ;
            particleX(1:brN)=mod(particleX(1:brN),sizeX);
    end
    
   
    gpuExpI=gpuExpI+single(E.*conj(E));
    if mod(i,floor(T/dT))==0
        expI(:,:,floor(i/floor(T/dT)))=gather(gpuExpI);
        gpuExpI(:)=0;
    end
    
    if storeRawI==1
        subI(:,:,mod(i-1,gpuBlockFrames)+1) = single(E.*conj(E));  %intensity
        if mod(i-1,gpuBlockFrames)==gpuBlockFrames-1
            rawI(:,:,i-gpuBlockFrames+1:i)=gather(subI);
            toc
            tic
        end
    end
end

if storeRawI==1
    rawI = rawI/mean(rawI(:));
end
expI=expI/mean(expI(:)); %this is the final output of images with exposure time T

%% Validate decorrelation time by fitting the models to data (requires storing the raw intensity)
g2=getG2(rawI,500,size(rawI,3),[1,1]);
g2=g2./(mean(rawI,3).^2);
g2=squeeze(mean(g2,[1,2]));
t=squeeze((0:length(g2)-1).*dT);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0,0],...
               'Upper',[1,1,1000],...
               'StartPoint',[0.9 0,10]);
ftOrdered = fittype('1+beta*((exp(-(x/tauc)^2))^2)+c','options',fo);
[ftOrdered,ftOrdGof] = fit(t(2:end)',double(g2(2:end)),ftOrdered);
ftUnordered = fittype('1+beta*((exp(-(x/tauc)))^2)+c','options',fo);
[ftUnordered,ftUnordGof] = fit(t(2:end)',double(g2(2:end)),ftUnordered);


figure
plot(t,g2)
hold on
plot(t,1+ftOrdered.beta*((exp(-(t/ftOrdered.tauc).^2)).^2)+ftOrdered.c)
plot(t,1+ftUnordered.beta*((exp(-(t/ftUnordered.tauc))).^2)+ftUnordered.c)
hold off
legend({'Data',['Ordered R^2=',num2str(ftOrdGof.rsquare)],['Unordered R^2=',num2str(ftUnordGof.rsquare)]})
xlabel('Time lag')
ylabel('g2')
set(gcf,'color','w');
title(['Ordered tauc=',num2str(ftOrdered.tauc),' Unordered tauc=',num2str(ftUnordered.tauc)])
toc


%% make an intensity gif
visStep=1;
visDT=T;
dataToPlot=expI; %rawI;
h=figure
tmp=dataToPlot(:,:,1:10);
percI01=prctile(tmp(:),1);
percI99=prctile(tmp(:),99);
for i=1:visStep:size(dataToPlot,3)
    imagesc(dataToPlot(:,:,i))
    caxis([percI01,percI99])
    xlabel('pixels')
    ylabel('pixels')
    axis image
    set(gcf,'color','w');
    title([num2str(i*visDT),' mus'])
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,cm,'Intensity.gif','gif', 'Loopcount',inf, 'BackgroundColor', 0,'DelayTime',0.1);
    else
        imwrite(imind,cm,'Intensity.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end