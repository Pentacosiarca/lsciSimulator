function averagePulse = calculateAveragePulse (dataLsci,positionLocMin,cellTypeOfAveraging)

sizeMedianPulse = median(positionLocMin(2:end)-positionLocMin(1:end-1));
averagePulse = nan(size(dataLsci,1),size(dataLsci,2),sizeMedianPulse);

typeOfAveraging = cell2mat(cellTypeOfAveraging(1));

switch typeOfAveraging
    case 'meanContrast'
        typeContrast = cell2mat(cellTypeOfAveraging(2));
        if strcmp(typeContrast,'temporal')
            contrastKernel = 5; %length(positionLocMin)-1;
        else
            contrastKernel = 5;
        end
        batchSize = 100; 
        procType = 'fastgpu';
        switch typeContrast
            case 'temporal'
                dataLsci = getTLSCI(dataLsci,contrastKernel,procType);%,batchSize);
            case 'spatial'
                dataLsci = getSLSCI(dataLsci,contrastKernel,procType);%,batchSize);
        end

        averagePulse = ...
            imresize3(dataLsci(:,:,1:positionLocMin(2)-1),[size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');
        countSums = 0;

%         for iPulses = 2:length(positionLocMin)-1
        for iPulses = 1:sizeMedianPulse-1
            endPosition = min(positionLocMin(iPulses+1)-1,size(dataLsci,3));
            if length(positionLocMin(iPulses):endPosition)==sizeMedianPulse
%                 pulse = dataLsci(:,:,positionLocMin(iPulses):endPosition);

%                 pulse = imresize3(dataLsci(:,:,positionLocMin(iPulses):endPosition),...
%                                           [size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');

                pulse = dataLsci(:,:,positionLocMin(1:end-1)+iPulses-1);

                countSums = countSums + 1;
                averagePulse(:,:,iPulses) = mean(pulse,3);
            end
        end

%         averagePulse = averagePulse ./ (countSums+1);

    case 'temporalContrastMean'

%         resizedPulse = nan(size(dataLsci,1),size(dataLsci,2),sizeMedianPulse*(length(positionLocMin)-1));
%         for iPulses = 1:length(positionLocMin)-1
%             batchMedianPulse = (iPulses-1)*sizeMedianPulse+1:iPulses*sizeMedianPulse;
%             resizedPulse(:,:,batchMedianPulse) = ...
%                         imresize3(dataLsci(:,:,positionLocMin(iPulses):positionLocMin(iPulses+1)-1),...
%                        [size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');
%         end
        resizedPulse = dataLsci;

        contrastKernel = length(positionLocMin)-1;
        batchSize = 100; 
        procType = 'fastcpu';
        for iOnePulse = 1:sizeMedianPulse
            disp(num2str(positionLocMin(1:end-1)+iOnePulse-1));
%             pulse = ...
%                 getTLSCI(resizedPulse(:,:,positionLocMin(1:end-1)+iOnePulse-1),contrastKernel,procType);
%            averagePulse(:,:,iOnePulse) = mean(pulse,3);
        end
end



end