% runBatchLSCI reads .rls data and calculates Spatial and/or Temporal
% contrast in batches. This function is used when the amount of available
% RAM memory is insufficient to store the whole LSCI sequence.
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Required Inputs:
%   fileName        - path to the .rls file
%   contrastType    - use 'spatial' or 'temporal' to calculate each kind of contrast.
%   contrastKernel  - size of the kernel to use for the contrast type.
%   batchSize       - size of the batch of data to be loaded for analysis.
%   procType        - choose the processor type, use: 
%                       'cpu', 'gpu', 'fastcpu', 'fastgpu'
%
% Optional Inputs:
%   skipFrames    - number of frames to skip at the begining of the recording. 
%                       Default is 0. Use 0 when you don't want to skip frames
%                       but need to enter other optional arguments.
%   movMeanKernel   - Kernel size for the move mean smoothing. Use 1 when 
%                       you don't want to use the algorithm but need to enter 
%                       other optional arguments.
%   downSampling    - Discard frames in between the downsampling value. Use
%                       0 when you don't want to downsample your data but
%                       want to select a ROI.
%   ROI             - Region of Interest:
%                       [firstRow,      lastRow ;
%                       firstColumn,   lastColumn]
%                       Default is the whole frame. 
%
% Outputs:
%    
%
% Example: 
%    
%
% Other m-files required: readRLS, getSLSCI, getTLSCI
% Subfunctions: readRLS, getSLSCI, getTLSCI
% MAT-files required: none
%
% See also: getTLSCI.m

% Authors: Alberto Gonzalez Olmos, Dmitry D Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 28-September-2021

%------------- BEGIN CODE --------------
function [dataOut,sampling,timeStampsOut,sizeT] = runBatchLSCI(fileName, contrastType, contrastKernel, batchSize, procType, varargin)
                                          

skipFrames = 0;
movMeanKernel = 1;
downSampling = 1;
ROI = [];
dataSize = 1;
maxFrames = 'all';

if ~isempty(varargin)
    skipFrames = varargin{1};
    for iVar = 2:1:length(varargin)
        if iVar==2, movMeanKernel = varargin{2}; end
        if iVar==3, downSampling = varargin{3}; end
        if iVar==4, ROI = varargin{4}; end
        if iVar==5, maxFrames = varargin{5}; end
    end
end

%% Reading .rls metadata
fileReadId = fopen(fileName, 'r');
fseek(fileReadId,0*1024,-1 );
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');
sizeT=fread(fileReadId,1,'*uint64');
sampling=fread(fileReadId,1,'*uint64');
version=fread(fileReadId,4,'*ubit8')';
if strcmp(version,'Ver.')
    nVer = fread(fileReadId,1,'*uint64');
    if nVer>1
        dataSize=fread(fileReadId,1,'*uint64');
    end
end


%% correct default values based on meta data
if isempty(ROI), ROI = [1,sizeY;1,sizeX]; end
switch dataSize
    case 1
        dataType='uint8';
    case 2
        dataType='uint16';
    otherwise
        error('Unindentified data type')
end
if ~strcmp(maxFrames,'all'), sizeT = maxFrames; end
numBatches = ceil((sizeT-skipFrames)/batchSize);


%% load batches of data after skipping frames
%pre-allocate memory for arrays
sizeTafterDownsampling = length(1:downSampling:sizeT-skipFrames);
batchSizeTafterDownsampling = length(1:downSampling:batchSize);
timeStampsOut=zeros(sizeTafterDownsampling,1,'int64');
dataOut=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),sizeTafterDownsampling,'single');


%move to the first timeStamp/frame location
firstByte=30*1024+sizeX*sizeY*skipFrames*dataSize+8*skipFrames;
fseek(fileReadId,firstByte,-1 );
temporalKernelIdx = 0;
if strcmp(contrastType,'temporal')
    temporalKernelIdx = contrastKernel;
end

wb = waitbar(0,'Please wait...');
counterFrames = 1;
tic
for iBatch=1:1:numBatches
    %update waitbar
    waitbar(single(iBatch)/single(numBatches),wb,['Processing ',contrastType,' contrast.']);
    %sizing batch to read in .lsc file
    firstIdx = (iBatch-1) * batchSize+1;
    lastIdx  = min(iBatch*batchSize, sizeT - skipFrames - temporalKernelIdx);
    batch = firstIdx:1:lastIdx + temporalKernelIdx;
    %memory allocation to read the batch
    data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),length(batch),'single');
    timeStamps = zeros(length(batch),1,'int64');
    %reading the batch from file
    for t=1:1:length(batch)
        timeStamps(t,1)=fread(fileReadId,1,'*uint64');
        frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
        data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
        counterFrames = counterFrames +1;
    end
    
    %rewind temporal kernel frames if calculating Temporal Contrast
    if strcmp(contrastType,'temporal')
        counterFrames = counterFrames - temporalKernelIdx;
        frewind(fileReadId)
        rewindFrames=firstByte+sizeX*sizeY*dataSize*(counterFrames-1)+8*(counterFrames-1);
        fseek(fileReadId,rewindFrames,-1 );
    end

    
    %% select contrast
    switch contrastType
        %spatial contrast
        case 'spatial'
            % use getSLSCI if there are changes in that function after the
            % last revision of this file
            usingGetSlsci = 0;
            if usingGetSlsci
            data = getSLSCI(data,contrastKernel,procType,size(data,3));
            else
                % processor type
                switch procType
                    case 'cpu'
                        for i=1:1:size(data,3)
                            frame=single(data(:,:,i));
                            frameMean=imfilter(frame,fspecial('average',[contrastKernel contrastKernel]));
                            frameSTD=stdfilt(frame,ones(contrastKernel));
                            data(:,:,i)=frameSTD./frameMean;
                        end
                    case 'gpu'
                        for i=1:1:size(data,3)
                            frame=gpuArray(single(data(:,:,i)));
                            frameMean=imfilter(frame,fspecial('average',[contrastKernel contrastKernel]));
                            frameSTD=stdfilt(frame,ones(contrastKernel));
                            data(:,:,i)=gather(frameSTD./frameMean);
                        end
                    case 'fastcpu'
                        kernel=ones(contrastKernel,contrastKernel,1);
                        data = stdfilt (data,kernel)...
                            ./ convn(data,kernel,'same') * length(kernel(:));
                    case 'fastgpu'
                        kernel=gpuArray(ones(contrastKernel,contrastKernel,1,'single'));
                        data = gather(stdfilt (gpuArray(data),kernel)...
                                   ./ convn(gpuArray(data),kernel,'same') * length(kernel(:)));
                end %procType
            end %usingGetSlsci
            
        %temporal contrast
        case 'temporal'
            % use getTLSCI if there are changes in that function after the
            % last revision of this file
            usingGetTlsci = 0;
            if usingGetTlsci
                data = getTLSCI(data,contrastKernel,procType,size(data,3));
            else
                switch procType
                    case 'cpu'
                        for i=1:1:size(data,3)-contrastKernel+1
                            frames=data(:,:,i:i+contrastKernel-1);
                            frameMean=squeeze(mean(frames,3));
                            frameSTD=squeeze(std(frames,0,3));
                            data(:,:,i)=frameSTD./frameMean;
                        end
                    case'gpu'
                        for i=1:1:size(data,3)-contrastKernel+1
                            frames=gpuArray(data(:,:,i:i+contrastKernel-1));
                            frameMean=squeeze(mean(frames,3));
                            frameSTD=squeeze(std(frames,0,3));
                            data(:,:,i)=frameSTD./frameMean;
                        end
                    case 'fastcpu'
                        data = movstd (data,[0,contrastKernel-1],0,3,'Endpoints','discard')...
                            ./ movmean(data,[0,contrastKernel-1],3,'Endpoints','discard');
                    case 'fastgpu'
                        data = gather (movstd (data,[0,contrastKernel-1],0,3,'Endpoints','discard')...
                            ./ movmean(data,[0,contrastKernel-1],3,'Endpoints','discard'));
                end %switch procType
            end %ifUsingTlsci
            
    end %contrastType
    
    %sizing the downsampled batch
    firstIdxDownsampled = (iBatch-1) * batchSizeTafterDownsampling+1;
    lastIdxDownsampled  = min(iBatch*batchSizeTafterDownsampling, sizeT - skipFrames - contrastKernel);
    batchDownsampled = firstIdxDownsampled:1:lastIdxDownsampled;
    dataDownsampled = round(linspace(1,size(data,3),length(batchDownsampled)));
    
    %calculating Moving Mean algorithm if the kernel is different from 1
    if movMeanKernel~=1, data = movmean(data,movMeanKernel); end
    
    %downsampling data and timeStamps
    dataOut (:,:,batchDownsampled) = data(:,:,dataDownsampled);
    timeStampsOut (batchDownsampled) = timeStamps(dataDownsampled);
end %for loop over the number of batches

fclose(fileReadId);
elapsed_time = toc;
close(wb);
disp(['Took ',num2str(elapsed_time),' seconds to compute the ',contrastType,' contrast with ',num2str(numBatches),' batches. The size of the data is: ',num2str(size(data,1)),'x',num2str(size(data,2)),'x',num2str(sizeT)]);

end
%------------- END OF CODE --------------
% Comments: note that using movstd leads to a rounding error up to e-06
% compared to convential std calculation. Switching to double precision
% decreases the error to the e-11 order.
%
% Large input data can lead to a memory overflow, particularly
% when uint8 input data is provided. This can be controlled by an additional
% outer loop and/or by converting sLSCI data to scaled integers or by
% allowing downsampling with a specific kernel size. Futhermore, memory overflow
% can be controlled when using 'fastcpu' or 'fastgpu' by specifying a
% smaller batch size.
%
% Finally, note that the 'gpu' and the 'fastgpu' processing types require
% a MATLAB supported GPU and their drivers to be installed.