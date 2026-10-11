function [cgradX2, cgradY2, cgradZ2] = concentration_gradients2(concentration2, coh_map)

    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    concentration2( coh_map==0 ) = 0; 
    
    cgradX2 = zeros(size(X));
    cgradX2(2:end-1,:,:) = concentration2(3:end,:,:) - concentration2(1:end-2,:,:);

    cgradY2 = zeros(size(X));
    cgradY2(:,2:end-1,:) = concentration2(:,3:end,:) - concentration2(:,1:end-2,:);

    cgradZ2 = zeros(size(X));
    cgradZ2(:,:,2:end-1) = concentration2(:,:,3:end) - concentration2(:,:,1:end-2);
    
    % find largest vector magnitude
    L = (cgradX2.^2 + cgradY2.^2 + cgradZ2.^2).^(.5);
    maxL = max(max(max(L)));

    % normalize gradient
    cgradX2 = cgradX2./maxL;
    cgradY2 = cgradY2./maxL;
    cgradZ2 = cgradZ2./maxL;

end