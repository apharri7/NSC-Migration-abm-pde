function op = setup_pde_operator(EV_file, FA_file, tumor_center_full, ...
    coh_limit, epsilon, crop_half, pad, use_full_domain, params)
% setup_pde_operator is a one-time assembly of the fixed discrete diffusion
% and drift operators, source field, etc.

% Input:
% 1. and 2.  EV_file, FA_file - eigenvector and coherency images
% 3. tumor_center_full - tumor center in original coh_map coordinate frame
% 4 coh_limit - coherency threshold 
% 5. epsilon - DDM interface thickness used in phi
% 6. crop_half - cropping half-width (ignored if use_full_domain=true) 
% 7. pad - number of zero-valued voxels added to every side of selected domain
% 8.  use_full_domain - true uses whole brain image; false uses a cropped version
% 9. params - struct containing .ic_type ('gaussian' or 'constant'),
%                              .chi_tumor - tumor-source strength
%                              .chi_bg - background logistic production strength
%                              .r_tumor - tumor source radius, .omega_tumor - source window transition width
%
% OUTPUT
% 1. op  - structure with the fixed PDE domain and reduced operators:
% a.        .L_diff_diag_r - diagonal diffusion operator
% b.        .L_corr_r - cross diffusion and drift correction operator
% c.        .dof - active voxel indices in a vectorized order
% d.        .perm - sparse-matrix ordering used prior to ILU/BiCGSTAB
% e.        .active - active brain mask
% f.        .source_fied_r - tumor-source fiede on active dof
% g.        coordinate maps, dimensions, params structure


[eigen_map, coh_map] = load_3D_run_tumor(EV_file, FA_file);
[nx_full, ny_full, nz_full, ~] = size(eigen_map);

% tumor center in original domain
cx0 = tumor_center_full(1); cy0 = tumor_center_full(2); cz0 = tumor_center_full(3);

% Can use the full domain or a cropped version centered around the tumor
if use_full_domain
    ix_full = 1:nx_full; iy_full = 1:ny_full; iz_full = 1:nz_full;
else
    ix_full = max(1, cx0-crop_half) : min(nx_full, cx0+crop_half);
    iy_full = max(1, cy0-crop_half) : min(ny_full, cy0+crop_half);
    iz_full = max(1, cz0-crop_half) : min(nz_full, cz0+crop_half);
end

% Get the selected region in eigen_map and coh_map
eigen_crop = double(eigen_map(ix_full, iy_full, iz_full, :));
coh_crop   = double(coh_map(ix_full, iy_full, iz_full));
clear eigen_map coh_map

% Extend our original domain (active brain region) to domain in which there are zeros 
% padded onto the original to create the computational domain
[ncx, ncy, ncz, ~] = size(eigen_crop);
nx = ncx + 2*pad; ny = ncy + 2*pad; nz = ncz + 2*pad;

eigen = zeros(nx, ny, nz, 3);
coh   = zeros(nx, ny, nz);
eigen(pad+(1:ncx), pad+(1:ncy), pad+(1:ncz), :) = eigen_crop;
coh(  pad+(1:ncx), pad+(1:ncy), pad+(1:ncz))    = coh_crop;
clear eigen_crop coh_crop

% Maps the tumor center in original domain to the padded domain
cx = (cx0 - ix_full(1) + 1) + pad;
cy = (cy0 - iy_full(1) + 1) + pad;
cz = (cz0 - iz_full(1) + 1) + pad;

[Y, X, Z] = meshgrid(1:ny, 1:nx, 1:nz);

% -- Construct anisotropic diffusion tensor
% Condition that separates what we consider white or grey matter
mask = coh > coh_limit;
vx = eigen(:,:,:,1); vy = eigen(:,:,:,2); vz = eigen(:,:,:,3);

% Normalize the DTI direction
vn = sqrt(vx.^2 + vy.^2 + vz.^2);
vn(vn == 0) = 1;
vx = vx./vn; vy = vy./vn; vz = vz./vn;
clear eigen

% Convert diffusion coefficients to voxels/day
scale_aniso = (6*10^3 * 24) / 13.5^2;
scale_iso   = (10^3 * 24) / 13.5^2;

% Constructs the diffusion tensor according to the white/grey matter regions
D11 = zeros(nx,ny,nz); D12 = D11; D13 = D11; D22 = D11; D23 = D11; D33 = D11;
D11(mask) = scale_aniso * vx(mask).^2;
D12(mask) = scale_aniso * vx(mask).*vy(mask);
D13(mask) = scale_aniso * vx(mask).*vz(mask);
D22(mask) = scale_aniso * vy(mask).^2;
D23(mask) = scale_aniso * vy(mask).*vz(mask);
D33(mask) = scale_aniso * vz(mask).^2;
iso = ~mask;
D11(iso) = scale_iso; D22(iso) = scale_iso; D33(iso) = scale_iso;
D21 = D12; D31 = D13; D32 = D23;        % diffusion tensor is symmetric
clear vx vy vz vn mask iso

% Constructs the phase field function (phi) and computational mask

% Nonzero coherency voxels define the brain
brain_mask = coh > 0;
% Computes the signed distance function used in phase field function
dist_outside = bwdist(brain_mask);
dist_inside  = bwdist(~brain_mask);
signed_dist  = double(dist_outside - dist_inside);
% Phase field function. 1 in brain, transitions to 0 outside the brain. 
phi_raw = 0.5*(1 - tanh((3*signed_dist)/epsilon));
% Stops division by 0 with parameter phi_floor when computing C=q/phi
phi = max(phi_raw, params.phi_floor);
% Solver domain includes brain and non-negligible diffuse interface decided
% by phi_cutoff
solve_mask = brain_mask | (phi_raw >= params.phi_cutoff);


% Construct tumor-sourced production field - used for Gaussian initial profile only

if isfield(params, 'ic_type') && strcmpi(params.ic_type, 'constant')
	source_field = zeros(nx, ny, nz);
else
    dist_to_tumor = sqrt((X-cx).^2 + (Y-cy).^2 + (Z-cz).^2);
    tumor_window = 0.5 * ...
        (1 - tanh((dist_to_tumor - params.r_tumor) / params.omega_tumor));

    source_field = params.chi_tumor * tumor_window;
    source_field(~brain_mask) = 0;

end

fprintf('source_field: min=%.4e, max=%.4e, mean=%.4e\n', ...
    min(source_field(:)), max(source_field(:)), mean(source_field(:)));


% Construct sparse shift and centered difference matrices

N = nx*ny*nz;
Ix = speye(nx); Iy = speye(ny); Iz = speye(nz);
Sx1 = circshift(Ix,[-1,0]); Sy1 = circshift(Iy,[-1,0]); Sz1 = circshift(Iz,[-1,0]);
Sx = kron(Iz, kron(Iy, Sx1));
Sy = kron(Iz, kron(Sy1, Ix));
Sz = kron(Sz1, kron(Iy, Ix));

% Centered first difference operator (dx=dy=dz=1)
Gx = 0.5*(Sx - Sx'); Gy = 0.5*(Sy - Sy'); Gz = 0.5*(Sz - Sz');

% Computes tensor drift, b = div(D)

bx = Gx*D11(:) + Gy*D21(:) + Gz*D31(:);
by = Gx*D12(:) + Gy*D22(:) + Gz*D32(:);
bz = Gx*D13(:) + Gy*D23(:) + Gz*D33(:);

% Construct phase field and diffusion tensor component matrices
Phi    = spdiags(phi(:), 0, N, N);
PhiInv = spdiags(1./max(phi(:),eps), 0, N, N);
D12m = spdiags(D12(:),0,N,N); D13m = spdiags(D13(:),0,N,N);
D21m = spdiags(D21(:),0,N,N); D23m = spdiags(D23(:),0,N,N);
D31m = spdiags(D31(:),0,N,N); D32m = spdiags(D32(:),0,N,N);

I = speye(N);

% Compute diagonal diffusion term fluxes 
D11x_face = 0.5*(D11(:) + Sx*D11(:));
phix_face = 0.5*(phi(:) + Sx*phi(:));
CoefX = spdiags(phix_face .* D11x_face, 0, N, N);
L_diff_xx = (I - Sx') * CoefX * (Sx - I) * PhiInv;

D22y_face = 0.5*(D22(:) + Sy*D22(:));
phiy_face = 0.5*(phi(:) + Sy*phi(:));
CoefY = spdiags(phiy_face .* D22y_face, 0, N, N);
L_diff_yy = (I - Sy') * CoefY * (Sy - I) * PhiInv;

D33z_face = 0.5*(D33(:) + Sz*D33(:));
phiz_face = 0.5*(phi(:) + Sz*phi(:));
CoefZ = spdiags(phiz_face .* D33z_face, 0, N, N);
L_diff_zz = (I - Sz') * CoefZ * (Sz - I) * PhiInv;

% Diagonal tensor portion of discretized diffusion term. Forms low order operator
L_diff_diag = L_diff_xx + L_diff_yy + L_diff_zz;

% Full discretized diffusion terms with cross terms added 
L_diff = L_diff_diag ...
       + Gx*Phi*(D12m*Gy + D13m*Gz)*PhiInv ...
       + Gy*Phi*(D21m*Gx + D23m*Gz)*PhiInv ...
       + Gz*Phi*(D31m*Gx + D32m*Gy)*PhiInv;

% Computes upwind discretization of drift term
bx_face = 0.5*(bx + Sx*bx);
Pplus_x  = spdiags(double(bx_face >= 0), 0, N, N);
Q_upwind_x = Pplus_x + (I-Pplus_x)*Sx;
Bx_face = spdiags(bx_face, 0, N, N);
L_drift_x = (I - Sx') * (Bx_face * Q_upwind_x);

by_face = 0.5*(by + Sy*by);
Pplus_y  = spdiags(double(by_face >= 0), 0, N, N);
Q_upwind_y = Pplus_y + (I-Pplus_y)*Sy;
By_face = spdiags(by_face, 0, N, N);
L_drift_y = (I - Sy') * (By_face * Q_upwind_y);

bz_face = 0.5*(bz + Sz*bz);
Pplus_z  = spdiags(double(bz_face >= 0), 0, N, N);
Q_upwind_z = Pplus_z + (I-Pplus_z)*Sz;
Bz_face = spdiags(bz_face, 0, N, N);
L_drift_z = (I - Sz') * (Bz_face * Q_upwind_z);

L_drift = L_drift_x + L_drift_y + L_drift_z;

clear Sx Sy Sz Gx Gy Gz Phi PhiInv D11x_face D22y_face D33z_face ...
      phix_face phiy_face phiz_face CoefX CoefY CoefZ D12m D13m D21m D23m D31m D32m ...
      bx_face by_face bz_face Pplus_x Pplus_y Pplus_z Q_upwind_x Q_upwind_y Q_upwind_z ...
      Bx_face By_face Bz_face L_drift_x L_drift_y L_drift_z I


% Reduce operators to brain and non-negligible interface
dof          = find(solve_mask(:));
inactive_idx = find(~solve_mask(:));
Ndof = numel(dof);

L_diff_diag_r = L_diff_diag(dof, dof);
L_diff_r = L_diff(dof, dof);
L_drift_r = L_drift(dof, dof);

% Return to boundary diagonals the contributions to discarded exterior voxels 
L_diff_diag_r = close_truncated_flux(L_diff_diag_r);
L_diff_r = close_truncated_flux(L_diff_r);
L_drift_r = close_truncated_flux(L_drift_r);

% symamd computes fill reducing permutation from the sparsity pattern of low order 
% diffusion matrix. 
perm = symamd(L_diff_diag_r);
% Contains terms exluded from low order matrix, which includes cross diffusion and drift term.
% This part is applied through a positivity correction
L_corr_r      = (L_diff_r  - L_diff_diag_r) + L_drift_r;

source_field_r = source_field(dof);

% Incorporate the linear tumor-source term into the implicit matrix
if strcmpi(params.ic_type, 'gaussian')
    L_diff_diag_r = L_diff_diag_r - ...
        spdiags(source_field_r, 0, Ndof, Ndof);
end

% Stores fixed operators, geometry, coordinate maps in op
op.nx = nx; op.ny = ny; op.nz = nz; op.N = N;
op.nx_full = nx_full; op.ny_full = ny_full; op.nz_full = nz_full;
op.X = X; op.Y = Y; op.Z = Z;
op.brain_mask = brain_mask;
op.phi = phi; op.solve_mask = solve_mask;
op.dof = dof; op.inactive_idx = inactive_idx; op.perm = perm;
op.L_diff_diag_r = L_diff_diag_r;
op.L_corr_r = L_corr_r;
op.cx = cx; op.cy = cy; op.cz = cz;
op.ix_full = ix_full; op.iy_full = iy_full; op.iz_full = iz_full;
op.pad = pad;
op.mid_z = round(nz/2);
op.source_field_r = source_field_r;
op.params = params;

fprintf('setup_pde_operator done: %d active DOF out of %d total voxels (%.1f%%)\n', ...
    numel(dof), N, 100*numel(dof)/N);
end


function L_r = close_truncated_flux(L_r)

column_leak = full(sum(L_r,1)).';

L_r = L_r - spdiags(column_leak, 0, size(L_r,1), size(L_r,2));

end
