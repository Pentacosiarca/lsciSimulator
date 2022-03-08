function G2=getG2(data,lagMax,sampleL,blocksN)

sampleL=sampleL-lagMax-1;
sizeY=size(data,1);
sizeX=size(data,2);
blockSizeY=floor(sizeY/blocksN(1));
blockSizeX=floor(sizeX/blocksN(2));
X=[1:blockSizeX:(sizeX-blockSizeX+1),sizeX];
Y=[1:blockSizeY:(sizeY-blockSizeY+1),sizeY];

G2=gpuArray(zeros(size(data,1),size(data,2),lagMax+1,'single'));
for i=1:1:length(X)-1
    for j=1:1:length(Y)-1
        subdata=gpuArray(single(data(Y(j):Y(j+1),X(i):X(i+1),:)));
        for lag=0:1:lagMax
            G2(Y(j):Y(j+1),X(i):X(i+1),lag+1)=mean(subdata(:,:,1:sampleL).*subdata(:,:,1+lag:sampleL+lag),3);
        end
    end
end
G2=gather(G2);
end
