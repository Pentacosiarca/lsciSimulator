function registrationCoordinates = registerMotion(dataLsci,referenceFrame)

% coordinates in Row, Columns to translate
registrationCoordinates = nan(2,size(dataLsci,3)-1); 
padding = 10;

batchSize = 100;
referenceFrame = getTLSCI(...
        dataLsci(:,:,100),...
        25,'fastgpu',batchSize);
referenceFrame = mean(referenceFrame,3);
referenceFrame = 1./(referenceFrame.^2);
% figure,imagesc(referenceFrame)

referenceFramePadded = zeros(size(dataLsci,1)+padding*2,size(dataLsci,2)+padding*2);
referenceFramePadded(padding+1:end-padding,padding+1:end-padding) = referenceFrame;
figure,imagesc(referenceFramePadded)
parfor iRegister = 1:1:size(dataLsci,3)-1

    disp([num2str(iRegister),' of ',num2str(size(dataLsci,3)-1)])
    frame = 1./(dataLsci(:,:,iRegister).^2);
    c = normxcorr2(frame,referenceFramePadded);

%     cropData = [200, 350,...
%                 100, 200];
%     data2 = dataLsci(cropData(1):cropData(2),cropData(3):cropData(4),iRegister);
% 
%     c = normxcorr2(data2,dataLsci(:,:,iRegister));

    % offset found by correlation
    [~,imax] = max(abs(c(:)));
    [rowPeak,colPeak] = ind2sub(size(c),imax(1));

    rowOffSet = rowPeak-size(dataLsci,1)+1;
    colOffSet = colPeak-size(dataLsci,2)+1;

    allowedOffset = 2;
    if (padding-allowedOffset < colOffSet) && (padding+allowedOffset > colOffSet)...
            && (padding-allowedOffset < rowOffSet) && (padding+allowedOffset > rowOffSet)
        registrationCoordinates(:,iRegister) = [0,0];
    else
        registrationCoordinates(:,iRegister) = [rowOffSet,colOffSet];
    end


end
