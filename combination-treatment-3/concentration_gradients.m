function [cgradX, cgradY, cgradZ] = concentration_gradients(concentration, coh_map)
% Input:
% 1. concentration - concentration field in 3d
% 2. coh_map - coherency image

% Output: (Note: gradient is normalized)
% 1. cgradX - x component of the gradient
% 2. cgradY - y component of the gradient
% 3. cgradZ - z component of the gradient

    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    concentration( coh_map==0 ) = 0; 
    
    % Centered differences in each direction 
    cgradX = zeros(size(X));
    cgradX(2:end-1,:,:) = concentration(3:end,:,:) - concentration(1:end-2,:,:);

    cgradY = zeros(size(X));
    cgradY(:,2:end-1,:) = concentration(:,3:end,:) - concentration(:,1:end-2,:);

    cgradZ = zeros(size(X));
    cgradZ(:,:,2:end-1) = concentration(:,:,3:end) - concentration(:,:,1:end-2);
    
    % Find largest gradient vector magnitude
    L = (cgradX.^2 + cgradY.^2 + cgradZ.^2).^(.5);
    maxL = max(max(max(L)));

    % Normalize the gradient 
    cgradX = cgradX./maxL;
    cgradY = cgradY./maxL;
    cgradZ = cgradZ./maxL;

end
