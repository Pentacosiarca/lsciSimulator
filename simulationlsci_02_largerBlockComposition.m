function [blockComposition, lastFramePeriodic, lastFrameBrownian] = simulationlsci_02_largerBlockComposition (periodicBlockOfParticles, brownianBlockOfParticles)

lastFramePeriodic = periodicBlockOfParticles(:,:,end);
lastFrameBrownian = brownianBlockOfParticles(:,:,end);

blockComposition = nan(size(brownianBlockOfParticles,1)+size(periodicBlockOfParticles,1),...
                     3,size(brownianBlockOfParticles,3));

stripeSizes = [0,1 ; 1,3 ; 3,6 ; 6,10];

sizePeriod = size(periodicBlockOfParticles,3);
wb = waitbar(0,'Please wait...');
for iter = 1:sizePeriod
%update waitbar
waitbar(single(iter)/single(sizePeriod),wb,['Processing simulation - large block composition -... iteration: ',num2str(single(iter)),' of ',num2str(single(sizePeriod))]);
brownianParticlesSlice = brownianBlockOfParticles(:,2,iter);
periodicParticlesSlice = periodicBlockOfParticles(:,2,iter);

cumulativeIdxBlockParticles = 0;
    for iStripes = 1:size(stripeSizes,1)

        % brownian particles

        brownianStripeIdx = find(brownianParticlesSlice >= stripeSizes(iStripes,1) &...
                                 brownianParticlesSlice <=  stripeSizes(iStripes,2));
        brownianStripe = double(brownianParticlesSlice(brownianStripeIdx));
        brownianStripe = abs(brownianStripe - min(brownianStripe));
        % normalize stripe
        brownianStripeNorm = brownianStripe./max(brownianStripe);
        % resize 
        brownianStripeResize = brownianStripeNorm * 0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1));
        brownianStripeResize = brownianStripeResize + stripeSizes(iStripes,1);

%         disp(['brownian max: ',num2str(0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1)))]);
%         disp(['starting at: ',num2str(stripeSizes(iStripes,1))]);

        blockComposition(cumulativeIdxBlockParticles+1 : cumulativeIdxBlockParticles+length(brownianStripeIdx),:,iter) = ...
            [ brownianBlockOfParticles(brownianStripeIdx,1,iter)' ;...
              brownianStripeResize' ;...
              brownianBlockOfParticles(brownianStripeIdx,3,iter)' ]';

        cumulativeIdxBlockParticles = cumulativeIdxBlockParticles + length(brownianStripeIdx);

        % periodic particles

        periodicStripeIdx = find(periodicParticlesSlice >= stripeSizes(iStripes,1) &...
                                 periodicParticlesSlice <=  stripeSizes(iStripes,2));
        periodicStripe = double(periodicParticlesSlice(periodicStripeIdx));
        periodicStripe = abs(periodicStripe - min(periodicStripe));
        % normalize stripe
        periodicStripeNorm = periodicStripe./max(periodicStripe);
        % resize 
        periodicStripeResize = periodicStripeNorm * 0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1));
        periodicStripeResize = periodicStripeResize + stripeSizes(iStripes,1) + 0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1));

%         disp(['periodic max: ',num2str(0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1)))]);
%         disp(['starting at: ',num2str(stripeSizes(iStripes,1) + 0.5 *(stripeSizes(iStripes,2) - stripeSizes(iStripes,1)))]);

        blockComposition(cumulativeIdxBlockParticles+1 : cumulativeIdxBlockParticles+length(periodicStripeIdx),:,iter) = ...
            [ periodicBlockOfParticles(periodicStripeIdx,1,iter)' ;...
              periodicStripeResize' ;...
              periodicBlockOfParticles(periodicStripeIdx,3,iter)' ]';

        cumulativeIdxBlockParticles = cumulativeIdxBlockParticles + length(periodicStripeIdx);

    end
end
close all
close(wb)


%% visualizations
isVisual = 0;
if isVisual
    % static representation
    figure,hold on,
    scatter(periodicBlockOfParticles(:,1,1),periodicBlockOfParticles(:,2,1),'b')
    scatter(brownianBlockOfParticles(:,1,1),brownianBlockOfParticles(:,2,1),'r')
    hold off,
    
    % animation of particles
    data = blockComposition;
    
    figure,
    for i=1:10:size(data,3)
        scatter(data(:,1,i),data(:,2,i))
        pause(0.01)
    end
end


