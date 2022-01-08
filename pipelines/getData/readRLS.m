%readRLS - reads raw laser speckle (.rls) data files
%
% Syntax:  output1 = function_name(requiredInput1,optionalInput1,optionalInput2,optionalInput3)
%
% Required inputs:
%   fileName      - path to the .rls file
% Optional Inputs:
%   framesToSkip  - number of frames to skip at the begining of the recording. 
%                   Default is 0. Use 0 when you don't want to skip frames
%                   but need to enter other optional arguments
%   framesToRead  - number of frames to read. Use [] when you want to read 
%                   all frames and use ROI. By default reads all frames
%                   except the skipped ones.
%   ROI           - Region of Interest:
%                   [firstRow,      lastRow ;
%                    firstColumn,   lastColumn]
%                   Default is the whole frame.
% Outputs:
%    data       - raw laser speckle data as 3d [y,x,t] matrix
%    sampling   - sampling time in ms
%    timeStamps - time stamps in microseconds
%    sizeT      - total number of frames in the file
%
% Example 1:
%    [data,sampling,timeStamps,sizeT]=readRLS('D:\stroke_20201001.rls',10,100,[10,100;20,200]);
%
% Example 2:
%    [data,sampling,timeStamps,sizeT]=readRLS('D:\stroke_20201001.rls');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getTLSCI.m, getSLSCI.m

% Authors: Dmitry D Postnov, Alberto Gonzalez Olmos
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 23-September-2021

%------------- BEGIN CODE --------------
function [data,sampling,timeStamps,sizeT]=readRLS(fileName,varargin)

%check if the fileName ends with .rls, add .rls otherwise
C = strsplit(fileName,'.');
if ~strcmp(C{end},'rls')
    fileName=[fileName,'.rls'];
end

%initializing default values for optional and hidden parameters
framesToSkip=0;
framesToRead=[];
ROI=[];
dataSize=1; %initiate data size as default for 'uint8'

%checking that inputs are in the correct format:
for iVar = 1:length(varargin)
    if iVar==1
        framesToSkip = varargin{1}; 
        if ~(round(framesToSkip)==framesToSkip)
            framesToSkip=round(framesToSkip);
            warning('Rounding number of frames to nearest integer')
        end
    end
    if iVar==2 && ~isempty(varargin{2})
        framesToRead = varargin{2};
        if ~(round(framesToRead)==framesToRead)
            framesToRead=round(framesToRead);
            warning('Rounding number of frames to nearest integer')
        end
    end
    if iVar==3 
        ROI = varargin{3};
        if ~all(size(ROI)==[2,2])
            error(['ROI format is incorrect, please check that the matrix dimensions are: ',char(13),...
                    '[firstRow    --> , <--   lastRow   --> ; <--',char(13),...
                    ' firstColumn --> , <--   lastColumn      ]'])
        end
    end
end

%read the meta data
fileReadId = fopen(fileName, 'r');
fseek(fileReadId,0*1024,-1 );
sizeX=fread(fileReadId,1,'*uint64');
sizeY=fread(fileReadId,1,'*uint64');
sizeT=single(fread(fileReadId,1,'*uint64'));
sampling=single(fread(fileReadId,1,'*uint64'));
version=fread(fileReadId,4,'*ubit8')';
if strcmp(version,'Ver.')
    nVer = fread(fileReadId,1,'*uint64');
    if nVer>1
        dataSize=fread(fileReadId,1,'*uint64');
    end
end

%correct default values based on meta data
if isempty(ROI), ROI = [1,sizeY;1,sizeX]; end
if isempty(framesToRead), framesToRead = sizeT; end
switch dataSize
    case 1
        dataType='uint8';
    case 2
        dataType='uint16';
    otherwise
        error('Unindentified data type')
end

%pre-allocate memory for arrays
timeStamps=zeros(framesToRead,1,'int64');
data=zeros(length(ROI(1,1):1:ROI(1,2)),length(ROI(2,1):1:ROI(2,2)),framesToRead,dataType);

%move to the first timeStamp/frame location
firstByte=30*1024+sizeX*sizeY*uint64(framesToSkip)*dataSize+8*uint64(framesToSkip);
fseek(fileReadId,firstByte,-1 );

%read data and close the file
for t=1:1:framesToRead
    timeStamps(t)=fread(fileReadId,1,'*uint64');
    frame=fread(fileReadId,[sizeY,sizeX],['*',dataType]);
    data(:,:,t)=frame(ROI(1,1):1:ROI(1,2),ROI(2,1):1:ROI(2,2));
end
fclose(fileReadId);

end
%------------- END OF CODE --------------
% Comments: this code can be used only with files that are several times
% smaller than the available RAM memory. Use the batch version of readRLS
% (batchReadRLS) to process the data using less RAM memory.