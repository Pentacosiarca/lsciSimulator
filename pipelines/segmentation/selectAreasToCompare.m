function [bloodFlow,parenchyma] = selectAreasToCompare(dataLSCI_temporalContrast)

figure;
title('select blood flow region')
imagesc(dataLSCI_temporalContrast);
h_rect = drawrectangle();



% Rectangle position is given as [xmin, ymin, width, height]
pos_rect = h_rect.Position;
% Round off so the coordinates can be used as indices
pos_rect = round(pos_rect);
% Select part of the image
img_cropped = dataLSCI_temporalContrast(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));

bloodFlow = [pos_rect(2), pos_rect(2) + pos_rect(4); pos_rect(1), pos_rect(1) + pos_rect(3)];

figure,imagesc(img_cropped)
pause;

figure;
title('select parenchyma region')
imagesc(dataLSCI_temporalContrast);
h_rect = drawrectangle();



% Rectangle position is given as [xmin, ymin, width, height]
pos_rect = h_rect.Position;
% Round off so the coordinates can be used as indices
pos_rect = round(pos_rect);
% Select part of the image
img_cropped = dataLSCI_temporalContrast(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));
% Coordinates
parenchyma = [pos_rect(2), pos_rect(2) + pos_rect(4); pos_rect(1), pos_rect(1) + pos_rect(3)];

figure,imagesc(img_cropped)
pause;
close all