function [num_cancers, Boundary_Ind] = find_boundary(coh_map, cancer_center)
% find_boundary.m creates a fixed interior target region for NSCs to move towards
% when they reach the boundary 

% Inputs:
% 1. coh_map - coherency image
% 2. cancer_center: contains coordinates of tumor center

% Outputs:
% 1. num_cancers - number of tumor sites
% 2. Boundary_Ind - linear coh_map indices inside target region.
%                   run_simulation_GBM samples these when directing NSC in out of 
%                   bounds region to the interior

    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    % Computes distance from every grid point to center of target region
    Boundary_Dist = sqrt(((X - 100)).^2+((Y - 140)).^2+((Z - 64)).^2);
    % Store matching voxels as linear indices into coh_map
    Boundary_Ind = find(Boundary_Dist < 11.42);

    % If multiple cancer sites exist, we assume they are all the same size
    num_cancers = size(cancer_center,1);

end

