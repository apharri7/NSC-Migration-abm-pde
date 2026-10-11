function noiseGaussian = generate_3D_noise_gaussian(coh_map ) 
% generate_3D_noise_gaussian generates spatially correlated 3D noise

% Input:
% 1. coh_map - coherency image. Used to specify dimensions of noise field.

% Output:
% 1. noiseGaussian - spatially correlated noise field with same dimensions as coh_map

ngrid = size( coh_map ); 

% Define coordinate grid passed to get_KL_coeff
dt = .01; 
timegrid1 = [1:ngrid(1)]*dt; 
timegrid2 = [1:ngrid(2)]*dt; 
timegrid3 = [1:ngrid(3)]*dt; 

% Correlation length in each axis 
lc = [0.1, 0.1, 0.1]; 

% Generate Karhunen Loeve expansion of stochastic process 
[nphi1, ~, ~, rDim(1)] = get_KL_coeff(lc(1), timegrid1);
[nphi2, ~, ~, rDim(2)] = get_KL_coeff(lc(2), timegrid2);
[nphi3, ~, ~, rDim(3)] = get_KL_coeff(lc(3), timegrid3);

% Sets the noise amplitude
sig = 2.0245e-04;
noiseGaussian = zeros( length(timegrid1), length(timegrid2), length(timegrid3) ); 

for n = 1:rDim(1)
    for m = 1:rDim(2) 
        for k = 1:rDim(3) 
            rsample = randn(1); 
            noise2D = nphi1(:,n) * nphi2(:,m)'; 
            for i = 1:size(noiseGaussian,3) 
              noiseGaussian(:,:,i) = noiseGaussian(:,:,i) + sig* noise2D * nphi3(i,k) * rsample; 
            end
        end 
    end 
end 

