function [EV, FA] = load_3D_run_tumor(EV_file, FA_file)
% load_3d_run_tumor opens and loads the eigenvector and fractional anisotropy files
% Inputs:
% 1. EV_file - data file containing principal diffusion vector at each voxel
% 2. FA_file - data file containing coherency value at each voxel

% Outputs:
% 1. EV - a four dimensional array of size nx x ny x nz x 3, where last dimension is eigenvector components
% 2. FA - a three dimansional array of size nx x ny x nz containing the coherency values 

    fid = fopen(EV_file,'r');
    EV = fread(fid,'float32','ieee-le');
    EV = permute(reshape(EV,3,200,280,128),[2,3,4,1]); % EV shape 200x280x128x3...this file contains principal eigenvector at each pixel
    fclose(fid);

    % fractional anisotropy 
    fid = fopen(FA_file,'r');
    FA = fread(fid,'float32','ieee-le');
    FA = reshape(FA,200,280,128);                      % FA contains 3D FA val at each pixel
    fclose(fid);

    % INTERPOLATE TO THE APPROPRIATE SIZE
    %limits = [21,180,21,160,36,105];    scale = 2;
    limits = [1,200,1,280,1,128];    scale = 0;
    EV = EV(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6),:);
    EV_x = interp3(EV(:,:,:,1),scale);  EV_y = interp3(EV(:,:,:,2),scale);  EV_z = interp3(EV(:,:,:,3),scale);
    clear EV; 
    EV(:,:,:,1) = EV_x; EV(:,:,:,2) = EV_y; EV(:,:,:,3) = EV_z;    
    clear EV_x EV_y EV_z

    FA = interp3(FA(limits(1):limits(2),limits(3):limits(4),limits(5):limits(6)),scale);
end 

