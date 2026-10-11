function state = init_pde_state(op, dt)
% init_pde_state: sets the initial concentration and terms needed for PDE solver

% Inputs:
% 1. op  - struct returned by setup_pde_operator.m
% 2. dt  - the PDE timestep 

% Output
% 1. state - structure containing: 
%           a. q_r - q=phi*C at active voxels
%           b. C_full - concentration in original domain
%           c. phi_r - phi at active voxels
%           d. forcing_r - fixed Gaussian IC forcing, zero for constant IC
%           e. q_lo_guess - initial BiCGSTAB guess vector 
%           f. mass0 - diagnostic (sum of q_r)
%           g. dt - pde timestep 
%           h. step_count - completed pde steps
%           f. last_mass_err - most recent mass diagnostic

% phi at active voxels
phi_r = op.phi(op.dof);


% Sets initial condition based on chosen constant or Gaussian
if isfield(op.params, 'ic_type') && strcmpi(op.params.ic_type, 'constant')
    c0 = op.params.ic_const * ones(op.nx, op.ny, op.nz);
else
    amp = op.params.ic_amp;
    sd  = op.params.ic_sd;
    c0 = amp .* exp( -((op.X-op.cx).^2 + (op.Y-op.cy).^2 + (op.Z-op.cz).^2) ./ (2*sd^2) );
end

% Add initial-condition noise.
if isfield(op.params, 'ic_noise_full') && ~isempty(op.params.ic_noise_full)
    if ~isfield(op.params, 'ic_noise_scale') || isempty(op.params.ic_noise_scale)
        error('init_pde_state:missingNoiseScale', ...
            'op.params.ic_noise_scale is required when op.params.ic_noise_full is provided.');
    end

    noise_local = initial_noise_to_local_frame(op, op.params.ic_noise_full);
    c0 = c0 + noise_local ./ op.params.ic_noise_scale;
end

% keep initial concentration only inside solve_mask
c0(~op.solve_mask) = 0;

% PDE solver advances q=phi*c
q0 = op.phi .* c0;
q0(~op.solve_mask) = 0;
q0_full = q0(:);

% Retain entries for active voxels only
q_r = q0_full(op.dof);

% Store initial reduced state 
state.q_r = q_r;
state.phi_r = phi_r;
state.mass0 = sum(q_r);

% Adds the fixed part of the tumor-sourced production term in Gaussian case
if isfield(op.params, 'ic_type') && strcmpi(op.params.ic_type, 'constant')
    state.forcing_r = zeros(size(q_r));
else
    state.forcing_r = dt * (op.source_field_r .* phi_r);
end

% Uses initial state as first guess, but this is updated in later pde steps with
% the low order solution
state.q_lo_guess = q_r;
state.dt = dt;
state.step_count = 0;
state.last_mass_err = 0;

% Converts q to C in the origianl domain
state.C_full = pde_state_to_full_concentration(op, q_r);
end


function C_full = pde_state_to_full_concentration(op, q_r)
% pde_state_to_full_concentration recovers C from q

% Insert active voxel q values. Inactive voxels still zero
q_full = zeros(op.N, 1);
q_full(op.dof) = q_r;

% Reshape into original 3d size
Q3d = reshape(q_full, op.nx, op.ny, op.nz);

% Undo the q=phi*C 
C3d = Q3d ./ max(op.phi, eps);
C3d(~op.solve_mask) = 0;
% ABM gets just the concentration in the brain tissue
C3d_output = C3d;
C3d_output(~op.brain_mask) = 0;

% Determines original domain size and removes padding from C
ncx = numel(op.ix_full); ncy = numel(op.iy_full); ncz = numel(op.iz_full);
C_crop = C3d_output(op.pad+(1:ncx), op.pad+(1:ncy), op.pad+(1:ncz));

% Puts the cropped C into array with dimensions mapping coh_map (original domain size)
C_full = zeros(op.nx_full, op.ny_full, op.nz_full);
C_full(op.ix_full, op.iy_full, op.iz_full) = C_crop;
end

function noise_local = initial_noise_to_local_frame(op, noise_full)
% initial_noise_to_local_frame - maps noise into padded, local pde frame

noise_full = double(noise_full);

local_size = [op.nx op.ny op.nz];
full_size = [op.nx_full op.ny_full op.nz_full];

if isequal(size(noise_full), local_size)
    % noise already represented in local padded pde frame
    noise_local = noise_full;
elseif isequal(size(noise_full), full_size)
    % Use the same original-frame crop
    ncx = numel(op.ix_full);
    ncy = numel(op.iy_full);
    ncz = numel(op.iz_full);

    noise_crop = noise_full(op.ix_full, op.iy_full, op.iz_full);
    
    % Put cropped noise into zero padded local array
    noise_local = zeros(op.nx, op.ny, op.nz);
    noise_local(op.pad+(1:ncx), op.pad+(1:ncy), op.pad+(1:ncz)) = noise_crop;
else
    % Print an error for mismatched dimensions
    error('init_pde_state:noiseSizeMismatch', ...
        ['Initial noise size is %s, but expected either local PDE size %s ' ...
         'or full coh_map size %s.'], mat2str(size(noise_full)), ...
        mat2str(local_size), mat2str(full_size));
end

% Noise outside solver mask  does not modify pde state
noise_local(~op.solve_mask) = 0;
end
