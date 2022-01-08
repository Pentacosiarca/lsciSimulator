%getSLSCI - calculates spatial Laser Speckle Contrast Images using a square
%kernel
%
% Syntax:  output1 = function_name(input1,input2,input3,input4)
%
% Inputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    kernelSize - number of pixels in a side of the kernel
%    procType   - choose the processor type, use: 
%                   'cpu', 'gpu', 'fastcpu', 'fastgpu'
%
% Optional Input:
%   batchSize   - size of the batch of data to be vectorized for analysis.
%
% Outputs:
%    sLSCI      - processed data as [y,x,t] 3d matrix
%
% Example: 
%    sLSCI=getSLSCI(data,7,'gpu','none')
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getTLSCI.m

% Authors: DD Postnov, Alberto Gonzalez Olmos
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 23-September-2021

%------------- BEGIN CODE --------------

function data=getSLSCI(data,kernelSize,procType,varargin)
data=single(data);
numBatches = 1;
batchSize=size(data,3);

if ~isempty(varargin)
    batchSize = varargin{1};
    numBatches = ceil(size(data,3)/batchSize);
end

switch procType
    case 'cpu'
        for i=1:1:size(data,3)
            frame=single(data(:,:,i));
            frameMean=imfilter(frame,fspecial('average',[kernelSize kernelSize]));
            frameSTD=stdfilt(frame,ones(kernelSize));
            data(:,:,i)=frameSTD./frameMean;
        end
    case 'gpu'
        for i=1:1:size(data,3)
            frame=gpuArray(single(data(:,:,i)));
            frameMean=imfilter(frame,fspecial('average',[kernelSize kernelSize]));
            frameSTD=stdfilt(frame,ones(kernelSize));
            data(:,:,i)=gather(frameSTD./frameMean);
        end
    case 'fastcpu'
        kernel=ones(kernelSize,kernelSize,1);
        for i = 1:1:numBatches
            firstIdx = (i-1) * batchSize+1;
            lastIdx  = min(i*batchSize, size(data,3));
            batch = firstIdx:1:lastIdx;
            data(:,:,batch) = stdfilt (data(:,:,batch),kernel)...
                ./ convn(data(:,:,batch),kernel,'same') * length(kernel(:));
        end
    case 'fastgpu'
        kernel=gpuArray(ones(kernelSize,kernelSize,1,'single'));
        for i = 1:1:numBatches
            firstIdx = (i-1) * batchSize+1;
            lastIdx  = min(i*batchSize, size(data,3));
            batch = firstIdx:1:lastIdx;
            dataBatch=gpuArray(data(:,:,batch));
            data(:,:,batch) = stdfilt (dataBatch,kernel)...
                ./ convn(dataBatch,kernel,'same') * length(kernel(:));
        end
end
end

%------------- END OF CODE --------------
% Comments: note that there might be up to e-07 order rounding error  
% processing types.
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