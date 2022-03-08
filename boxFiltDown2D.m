%% boxfiltering + downsamling

function filteredImage = boxFiltDown2D (im, kernel, typeOfFilter)

filteredImage = nan(floor(size(im,1)/kernel),floor(size(im,2)/kernel));

for iRow = 1:1:size(filteredImage,1)
    for iCol = 1:1:size(filteredImage,2)
        boxRows = 1+(iRow-1)*kernel : iRow*kernel;
        boxCols = 1+(iCol-1)*kernel : iCol*kernel;
        boxFrame = im(boxRows,boxCols);
        switch typeOfFilter
            case 'average'
                pixelVal = mean(boxFrame,[1 2]);
        end
        filteredImage(iRow,iCol) = pixelVal;
    end
end

% figure,imagesc(filteredImage)
% figure,imagesc(im)