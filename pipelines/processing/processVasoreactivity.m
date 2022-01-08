%calculateVasoreactivity - calculates several features from the vasoreactivity 
% from LSCI data
%
% Syntax:  [output1, output2, output3, output4] = function_name (input1,input2)
%
% Inputs:
%    bfiTimeSecondsStartStim    - time stamps for each frame from the start
%                                   of the stimulus
%    sLsciStartStim             - spatial contrast of LSCI data 
%                                   from the start of the stimulus
%
% Outputs:
%   ttp               - time to peak from the start of the stimulation
%   auc_beforePeak    - Area Under the Curve from time of stimulation 
%                           to the peak
%   slopeAngle        - the slope's angle of the "sigmoid" curve
%   rCBF              - relative cerebral blood flow map between an average 
%                           frame at the start of the simulation with an
%                           average frame at the peak.
%
% Other m-files required: none 
% Subfunctions: none
% MAT-files required: none

% Authors: Alberto Gonzalez Olmos, DD Postnov
% CFIN, Aarhus University
% email address: dpostnov@cfin.au.dk
% Last revision: 15-December-2021

%------------- BEGIN CODE --------------


function [ttp, auc_beforePeak, slopeAngle, rCBF] = processVasoreactivity (bfiTimeSecondsStartStim, sLsciStartStim)


% blood flow index signal from spatial contrast
% find infinite values and discard them
bfiTS=squeeze(mean(sLsciStartStim,[1,2]));
bfiTS=1./(bfiTS.^2);
bfiTsInf = find(bfiTS == inf);
sLsciStartStim(:,:,bfiTsInf) = [];
bfiTimeSecondsStartStim(bfiTsInf) = [];

% recalculate bfi without infinite values
bfiTS=squeeze(mean(sLsciStartStim,[1,2]));
bfiTS=1./(bfiTS.^2);


%% output parameters:
[framePeak, slopeAngle, rCBF] = calculatePeakSlopeRcbf (bfiTS, sLsciBaseline, sLsciStartStim);
ttp = bfiTimeSecondsStartStim(framePeak) - bfiTimeSecondsStartStim(1);
auc_beforePeak = sum(bfiTS(1:framePeak));

