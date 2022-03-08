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
                dataLsci = getTLSCI(dataLsci,contrastKernel,procType,batchSize);
            case 'spatial'
                dataLsci = getSLSCI(dataLsci,contrastKernel,procType,batchSize);
        end

        averagePulse = ...
            imresize3(dataLsci(:,:,1:positionLocMin(2)),[size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');
        countSums = 0;

        for iPulses = 2:length(positionLocMin)-1


            endPosition = min(positionLocMin(iPulses+1),size(dataLsci,3));
            if endPosition>positionLocMin(iPulses)
                pulse = imresize3(dataLsci(:,:,positionLocMin(iPulses):endPosition),...
                                          [size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');
                countSums = countSums + 1;
                averagePulse = bsxfun(@plus,pulse,averagePulse);
            end

               
        end

        averagePulse = averagePulse ./ (countSums+1);

    case 'temporalContrastMean'

        resizedPulse = nan(size(dataLsci,1),size(dataLsci,2),sizeMedianPulse*(length(positionLocMin)-1));
        for iPulses = 1:length(positionLocMin)-1
            resizedPulse(:,:,positionLocMin(iPulses):positionLocMin(iPulses)+sizeMedianPulse-1) = ...
                        imresize3(dataLsci(:,:,positionLocMin(iPulses):positionLocMin(iPulses+1)),...
                       [size(dataLsci,1),size(dataLsci,2),sizeMedianPulse],'linear');
        end

        contrastKernel = length(positionLocMin)-1;
        batchSize = 100; 
        procType = 'fastgpu';
        for iOnePulse = 1:sizeMedianPulse
%             pulse = getTLSCI(resizedPulse(:,:,positionLocMin(1:end-1)+iOnePulse-1),contrastKernel,procType,batchSize);
%             averagePulse(:,:,iOnePulse) = mean(pulse,3);

            pulse = ...
                getTLSCI(resizedPulse(:,:,positionLocMin(1:end-1)+iOnePulse-1),contrastKernel,procType,batchSize);
           averagePulse(:,:,iOnePulse) = mean(pulse,3);
        end
end



end