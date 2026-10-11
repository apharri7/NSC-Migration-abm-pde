function [state, C_full] = step_pde_concentration(op, state, nsc_positions_full, dt)
% step_pde_concentration advances the pde by one timestep

% Inputs:
% 1. op - structure from setup_pde_operator.m. Contains fixed operators, functions, etc.
%           not changing in time
% 2. state - structure from init_pde_state.m. Stores q on active voxels and the initial guess for solver
% 3. nsc_positions_full - contains NSCs coordinates in the original domain
% 4. dt - pde timestep

% Outputs: 
% 1. state - updated pde state
% 2. C_full - concentration on the original grid with padding removed. 

% NSC internalization parameters
R_nsc       = op.params.R_nsc;
lambda_sink = op.params.lambda_sink;

% Settings for BiCGSTAB
tol   = 1e-8;
maxit = 200;

% Local grid dimensions, brian tissue voxels, and pde solver voxels
nx = op.nx; ny = op.ny; nz = op.nz;
dof = op.dof;
brain_mask = op.brain_mask;
solve_mask = op.solve_mask;
Ndof = numel(dof);

% Map NSC positions from the full coh_map frame into the local, padded frame
local_ijk = pde_full_to_local(op, nsc_positions_full);
local_ijk(:,1) = min(max(local_ijk(:,1), 1), nx);
local_ijk(:,2) = min(max(local_ijk(:,2), 1), ny);
local_ijk(:,3) = min(max(local_ijk(:,3), 1), nz);

% NSC sink field with truncated Gaussian distribution for ease of computation
N_NSC = size(local_ijk, 1);
sink_field = zeros(nx, ny, nz);
normfac = lambda_sink / ((2*pi*R_nsc^2)^(3/2));

% Evaluate the Gaussian only withing 6 stds of center. Contributions outside this 
% window are considered  negligible for computation
half_w = ceil(6*R_nsc);

for n = 1:N_NSC
    ic = local_ijk(n,1); jc = local_ijk(n,2); kc = local_ijk(n,3);
    ix = max(1, floor(ic-half_w)) : min(nx, ceil(ic+half_w));
    iy = max(1, floor(jc-half_w)) : min(ny, ceil(jc+half_w));
    iz = max(1, floor(kc-half_w)) : min(nz, ceil(kc+half_w));
    [xs, ys, zs] = ndgrid(ix, iy, iz);
    d2 = (xs-ic).^2 + (ys-jc).^2 + (zs-kc).^2;
    sink_field(ix,iy,iz) = sink_field(ix,iy,iz) + normfac*exp(-d2/(2*R_nsc^2));
end
sink_field(~brain_mask) = 0;

% Construct reduced low order backward Euler system

Sink_r = spdiags(sink_field(dof), 0, Ndof, Ndof);

% L_lo_r is montontone low order operator used in the implicit baseline solve
L_lo_r = op.L_diff_diag_r - Sink_r;

% Backward euler matrix
A_lo_r = speye(Ndof) - dt*L_lo_r;

% Applies fill-reducing permutation from setup_pde_operator.m. Same permutation is
% applied to rhs and initial guess to give an algebraically equivalent matrix eqn
p = op.perm;
A_lo_r_p = A_lo_r(p,p);

% Compute incomplete LU preconditioner for reordered system
ilu_setup.type = 'nofill';
[L1, U1] = ilu(A_lo_r_p, ilu_setup);


% Get the current concentration
C_r = state.q_r ./ max(state.phi_r, eps);

% Compute the logistic production term
logistic_prod_r = op.params.chi_bg .* state.phi_r .* C_r .* (1 - C_r);
logistic_forcing_r = dt .* logistic_prod_r;

% Solve low-order implicit system

% rhs involves the tumor source term (state.forcing_r) and the logistic production term (logistic_forcing_r)
rhs_r = state.q_r + state.forcing_r + logistic_forcing_r;

% Reorder the rhs and previous low order solution which is used as initial guess
rhs_r_p = rhs_r(p);
q_lo_guess_p = state.q_lo_guess(p);

% Solve A_lo_r_p * q_lo_p = rhs_r_p with no-fill ILU factors as preconditioners
[q_lo_p, flag, relres, iter] = bicgstab(A_lo_r_p, rhs_r_p, tol, maxit, L1, U1, q_lo_guess_p);
if flag ~= 0
    fprintf('  step_pde_concentration: bicgstab flag=%d (relres=%.2e, iters=%d)\n', flag, relres, iter);
end

% Undo permutation to put low order soluton in original order according to active DOF
q_lo = zeros(Ndof,1);
q_lo(p) = q_lo_p;

% Apply correction term
corr = dt * (op.L_corr_r * q_lo);
neg = corr < 0;

% Compute global limiter theta that prevents q from going negative
if any(neg)
    theta = min(1, 0.999 * min(q_lo(neg) ./ max(-corr(neg), eps)));
else
    theta = 1.0;
end
theta = max(0, theta);

% Add correction to low-order solution
q_r = q_lo + theta*corr;

% Store updated q and new low order guess
state.q_r = q_r;
state.q_lo_guess = q_lo;
state.step_count = state.step_count + 1;
state.last_mass_err = (sum(q_r) - state.mass0) / max(abs(state.mass0), eps);

% Compute C from q and reshape to 3d
q_full = zeros(op.N, 1);
q_full(dof) = q_r;
Q3d = reshape(q_full, nx, ny, nz);
C3d = Q3d ./ max(op.phi, eps);
C3d(~solve_mask) = 0;

% ABM only gets concentration in the brain 
C3d_output = C3d;
C3d_output(~brain_mask) = 0;

% Removes padding or cropping to original data dimensions
C_full = extract_full_field(op, C3d_output);
state.C_full = C_full;
end


function local_ijk = pde_full_to_local(op, full_ijk)
% Convert original coordinates to local coordinates
% to account for padding and cropping
local_ijk = zeros(size(full_ijk));
local_ijk(:,1) = full_ijk(:,1) - op.ix_full(1) + 1 + op.pad;
local_ijk(:,2) = full_ijk(:,2) - op.iy_full(1) + 1 + op.pad;
local_ijk(:,3) = full_ijk(:,3) - op.iz_full(1) + 1 + op.pad;
end


function C_full = extract_full_field(op, C3d)
% Go back to orginal grid from the padded and cropped grid
ncx = numel(op.ix_full); ncy = numel(op.iy_full); ncz = numel(op.iz_full);
C_crop = C3d(op.pad+(1:ncx), op.pad+(1:ncy), op.pad+(1:ncz));

C_full = zeros(op.nx_full, op.ny_full, op.nz_full);
C_full(op.ix_full, op.iy_full, op.iz_full) = C_crop;
end

