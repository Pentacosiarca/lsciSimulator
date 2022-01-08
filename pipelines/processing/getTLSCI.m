%getTLSCI - calculates temporal Laser Speckle Contrast Images
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Inputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    kernelSize - number of pixels in a side of the kernel
%    procType   - choose the processor type: use 'cpu', 'gpu', 'fastcpu' or
%                 'fastgpu'. 'fastcpu' should be used by default.
% Optional Input:
%   batchSize   - size of the batch of data to be vectorized for analysis.
%
% Outputs:
%    tLSCI      - processed data as [y,x,t] 3d matrix
%
% Example:
%    tLSCI=getSLSCI(data,25,'gpu','none')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getSLSCI.m

% Authors: DD Postnov, Alberto Gonzalez Olmos
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 23-September-2021

%------------- BEGIN CODE --------------

function data=getTLSCI(data,kernelSize,procType,varargin)
data=single(data);
numBatches = 1;
batchSize=size(data,3);

if ~isempty(varargin)
    batchSize = varargin{1};
    numBatches = ceil(size(data,3)/batchSize);
end

switch procType
    case 'cpu'
        for i=1:1:size(data,3)-kernelSize+1
            frames=data(:,:,i:i+kernelSize-1);
            frameMean=squeeze(mean(frames,3));
            frameSTD=squeeze(std(frames,0,3));
            data(:,:,i)=frameSTD./frameMean;
        end
    case'gpu'
        for i=1:1:size(data,3)-kernelSize+1
            frames=gpuArray(data(:,:,i:i+kernelSize-1));
            frameMean=squeeze(mean(frames,3));
            frameSTD=squeeze(std(frames,0,3));
            data(:,:,i)=frameSTD./frameMean;
        end
    case 'fastcpu'
        for i = 1:1:numBatches
            firstIdx = (i-1) * batchSize+1;
            lastIdx  = min(i*batchSize, size(data,3) - kernelSize + 1);
            batch = firstIdx:1:lastIdx;
            batchKernel = firstIdx:1:(lastIdx+kernelSize-1);
            data(:,:,batch) = movstd (data(:,:,batchKernel),[0,kernelSize-1],0,3,'Endpoints','discard')...
                ./ movmean(data(:,:,batchKernel),[0,kernelSize-1],3,'Endpoints','discard');
        end
    case 'fastgpu'
        for i = 1:1:numBatches
            firstIdx = (i-1) * batchSize+1;
            lastIdx  = min(i*batchSize, size(data,3) - kernelSize + 1);
            batch = firstIdx:1:lastIdx;
            batchKernel = firstIdx:1:(lastIdx+kernelSize-1);
            dataBatch=gpuArray(data(:,:,batchKernel));
            data(:,:,batch) = movstd (dataBatch,[0,kernelSize-1],0,3,'Endpoints','discard')...
                ./ movmean(dataBatch,[0,kernelSize-1],3,'Endpoints','discard');
        end
end

data(:,:,size(data,3)-kernelSize+1+1:end)=[];

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






