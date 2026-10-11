function [p, total_dead, n_seeds_current, daily_arrivals, pde_state] = run_simulation_GBM(alpha4chmtx, beta4dist, Boundary_Ind, cancer_center, cancer_size, ... 
op, pde_state, dt_pde, ...
chemo_sensitivity, chmtx_limit, coh_limit, coh_map, CONVERT2MICRON, d_g, d_w, d_noise_value, dmax, Doses, ...
dtime, eigen_map, Finaltimestep, has_cancer, inj_center, lbound, modelNum, ...
modelType, nbuffer, n_seeds, num_cancers, p, seed_ind, seed_sd, ubound)
% run_simulation_GBM advances the the coupled ABM and chemoattractant PDE

% Important inputs:
% 1. op - fixed PDE operator, 
% 2. pde_state - PDE state evolved
% 3. dt_pde - PDE timestep
% 4. coh_map - coherency 
% 5. eigen_map - principal DTI vectors
% 6. d_g - grey matter reference step lengths
% 7. d_w - white matter reference step length
% 8. alpha4chmtx, chemo_sensitivity, chmtx_limit - chemotactic bias parameters
% 11. Finaltimestep - final ABM time for each dose interval
% 12. n_seeds - cumulative NSC count after each injection
% 13. lbound, ubound - lower and upper brain z-surfaces
% 15. Boundary_Ind  - candidate interior target region used only when an NSC
%                   cannot be returned in bounds by the first retry
%
% Outputs:
% 1. p  - NSC trajectories structure
% 2. total_dead - zero because cell death is not modeled here
% 3. n_seeds_current - total number of cells at final dose
% 4. daily_arrivals - newly observed tumor occupying cells at the dose endpoints
% 5. pde_state - final PDE state

% For each seed
Tstart = tic;
disp('Running simulation...')

% No cells die in this ABM version; all injected cells remain
n_seeds_current = n_seeds(1,Doses);
total_dead = 0;

% Save z-slice heatmap at each simulated day
tumor_center = cancer_center(1,:);
run_tag = sprintf('base_cancer_%.0f_%.0f_%.0f', tumor_center(1), tumor_center(2), tumor_center(3));
save_concentration_slices = true;
concentration_slice_times = unique([1, 200:200:Finaltimestep(1,Doses)]);
concentration_slice_dir = sprintf('concentration_slices_%s', run_tag);


[Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));

    %%%% Spawn seed within injection site
i = 1;
for seed_loop = 1:n_seeds(1,1)
    seed = seed_ind(seed_loop);
    [seed(1,1),seed(1,2),seed(1,3)] = ind2sub(size(coh_map),seed);
    if seed(3) > ubound(seed(1),seed(2)) - nbuffer
        seed(3) = ubound(seed(1),seed(2)) - nbuffer;
    elseif seed(3) < lbound(seed(1),seed(2)) + nbuffer
        seed(3) = lbound(seed(1),seed(2)) + nbuffer;
    end
    p(seed_loop).coord(:,i) = seed'; % Coordinate of this starting point
    %%%% Initialize particular reference step length
    dmax(1,seed_loop) = betainv(rand(1), 1, beta4dist);
end

% i = 1 is the PDE initial condition from init_pde_state
concentration = pde_state.C_full;
save_concentration_slice_if_needed(i, concentration, tumor_center, ...
    save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);


% Advance concentration to i = 2 using the initial NSC positions at i = 1
seedposition = zeros(3, n_seeds(1,1));
for seed_loop = 1:n_seeds(1,1)
    seedx = round( p(seed_loop).coord(1,1) );
    seedy = round( p(seed_loop).coord(2,1) );
    seedz = round( p(seed_loop).coord(3,1) );
    seedposition(:, seed_loop) = [seedx, seedy, seedz]';
end
[pde_state, concentration] = step_pde_concentration(op, pde_state, seedposition', dt_pde);
concentration(coh_map == 0) = 0;
save_concentration_slice_if_needed(2, concentration, tumor_center, ...
    save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);

% Factor of 5000 scales to original VEGF data assumed maximum before gradient computed
[cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);

    %%%% First step
i = 2;
for seed_loop = 1:n_seeds(1,1)
    eigen_vect_noise = [rand(1)*2-1, rand(1), rand(1)]';
    eigen_vect = eigen_vect_noise/norm(eigen_vect_noise);
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(d_g*eigen_vect); % Find next coordinate at step d along EV direction
    seed = round( p(seed_loop).coord(:,i) );
    if seed(3) >= ubound(seed(1),seed(2)) - nbuffer || seed(3) <= lbound(seed(1),seed(2)) + nbuffer ...
        || isnan(ubound(seed(1),seed(2))) || isnan(lbound(seed(1),seed(2)))   
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
    end
end

% Matrix to store the position of the n_seeds for PDE's internalization term
seedposition = zeros(3, n_seeds(1,1));

 %%%% Advances the first dose of cells
for i = 3:Finaltimestep(1,1) 

    for seed_loop = 1:n_seeds(1,1)
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    [pde_state, concentration] = step_pde_concentration(op, pde_state, seedposition', dt_pde);
    concentration(coh_map == 0) = 0;
    % ------------------------------
    save_concentration_slice_if_needed(i, concentration, tumor_center, ...
        save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);
    % ------------------------------
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);

    % Move the seeds
    for seed_loop = 1:n_seeds(1,1)
        %%% Find index (Finds location of each of the seeds, store this)
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
                % --- Add this bounds check ---
        if seedx < 1 || seedx > size(coh_map,1) || ...
           seedy < 1 || seedy > size(coh_map,2) || ...
           seedz < 1 || seedz > size(coh_map,3) || ...
           isnan(seedx) || isnan(seedy) || isnan(seedz)

            % Put the cell back to the last valid location and skip
           p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
           continue;
        end

        ind = sub2ind(size(coh_map), seedx, seedy, seedz);
        
        %%% Calculate eigen vector, step length, and prefered direction +/-
        if coh_map(ind) > coh_limit % Get Eigen vector at the previous coordinate point
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) );
            d = dmax(1,seed_loop) * d_w;
            pref_mvmt = sign(dot( p(seed_loop).coord(:,i-1) - p(seed_loop).coord(:,i-2), eigen_vect) + eps ); % Move either in the +/- direction
            d_noise = d_noise_value;
        else
            eigen_vect = 0;
            d = dmax(1, seed_loop) * d_g; % Move with less magnitude if no WM present 
            pref_mvmt = 1;
            d_noise = 1;
        end
        eigen_vect_noise = (rand(3,1)*2-1);
        eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
        
        dstep = d;

        %%% Calculate chemotaxis vector and gradient strength
        chmtx_vect = [cgradX(ind); cgradY(ind); cgradZ(ind)];
        chmtx_bias = 0;
        switch modelNum
            case 1  % Chemotaxis
                % Move along chemotaxis (deterministic)
                dstep = dstep/dmax(1, seed_loop); % Keeps dstep deterministic (dmax is stochastic)
                if has_cancer
                    chmtx_bias = chemo_sensitivity;
                end
            case 2  % Chemotaxis, Stochasticity
                %  Move along chemotaxis (stochasticity)
                if has_cancer
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
            case 3  % Chemotaxis, Chemotax threshold, Stochasticity
                %  Move along chemotaxis once near cancer (stochasticity)
                if has_cancer && (concentration(ind) >= chmtx_limit)
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
        end
        %%% Take step: (curr coord) = (prev coord) + (coherency influence) + (cancer influence)
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);
        
        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);
            seednow = round(p(seed_loop).coord(:,i));
            countloop = 1;
            while ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                Center_Point = Boundary_Ind(ceil(rand*length(Boundary_Ind)));
                [Center_Point(1,1),Center_Point(1,2),Center_Point(1,3)] = ind2sub(size(coh_map),Center_Point);
                Center_Vector = Center_Point - p(seed_loop).coord(:,i)';
                CenterVector = Center_Vector./norm(Center_Vector);
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*CenterVector';
                countloop = countloop + 1;
                seednow = round(p(seed_loop).coord(:,i));
                if( countloop >= 50 );   break;   end 
            end
        end
    end
end

%% Run loop for doses
for j = 1:Doses-1
    i = Finaltimestep(1,j) + 1;

    seedposition = zeros(3, n_seeds(1,j));
    for seed_loop = 1:n_seeds(1,j)
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    [pde_state, concentration] = step_pde_concentration(op, pde_state, seedposition', dt_pde);
    concentration(coh_map == 0) = 0;
    % -----------------------------
    save_concentration_slice_if_needed(i, concentration, tumor_center, ...
        save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);
    % -----------------------------
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);

    for seed_loop = 1:n_seeds(1,j)
        %%% Find index
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
               
                % --- Add this bounds check ---
        if seedx < 1 || seedx > size(coh_map,1) || ...
           seedy < 1 || seedy > size(coh_map,2) || ...
           seedz < 1 || seedz > size(coh_map,3) || ...
           isnan(seedx) || isnan(seedy) || isnan(seedz)

            % Put the cell back to the last valid location and skip
           p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
           continue;
        end

        ind = sub2ind(size(coh_map), seedx, seedy, seedz);

        %%% Calculate eigen vector, step length, and prefered direction +/-
        if coh_map(ind) > coh_limit % Get Eigen vector at the previous coordinate point
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) ); 
            d = dmax(1, seed_loop) * d_w;
            pref_mvmt = sign(dot( p(seed_loop).coord(:,i-1) - p(seed_loop).coord(:,i-2), eigen_vect) + eps ); % Move either in the +/- direction
            d_noise = d_noise_value;
        else
            eigen_vect = 0;
            d = dmax(1, seed_loop) * d_g; % Move with less magnitude if no WM present 
            pref_mvmt = 1;
            d_noise = 1;
        end
        eigen_vect_noise = (rand(3,1)*2-1);
        eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);

        %%% Keep or reduce step length by some amount
        %%%%%% uniform[0, 1] 
%         dstep = rand(1)*dmax ; 
        %%%%%% beta[1, beta], if beta==1, uniform 
%         dstep = betainv( rand(1), 1, beta4dist )*d ; 
        %%%%%% deterministic
        dstep = d;

        %%% Calculate chemotaxis vector and gradient strength
        chmtx_vect = [cgradX(ind); cgradY(ind); cgradZ(ind)];
        chmtx_bias = 0;
        switch modelNum
            case 1  % Chemotaxis
                % Move along chemotaxis (deterministic)
                dstep = dstep/dmax(1,seed_loop); % Keeps dstep deterministic (dmax is stochastic)
                if has_cancer
                    chmtx_bias = chemo_sensitivity;
                end
            case 2  % Chemotaxis, Stochasticity 
                %  Move along chemotaxis (stochasticity)
                if has_cancer
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
            case 3  % Chemotaxis, Chemotax threshold, Stochasticity
                %  Move along chemotaxis once near cancer (stochasticity)
                if has_cancer && (concentration(ind) >= chmtx_limit)
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
        end
        %%% Take step: (curr coord) = (prev coord) + (coherency influence) + (cancer influence)
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);

        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);
            seednow = round(p(seed_loop).coord(:,i));
            countloop = 1;
            while ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                Center_Point = Boundary_Ind(ceil(rand*length(Boundary_Ind)));
                [Center_Point(1,1),Center_Point(1,2),Center_Point(1,3)] = ind2sub(size(coh_map),Center_Point);
                Center_Vector = Center_Point - p(seed_loop).coord(:,i)';
                CenterVector = Center_Vector./norm(Center_Vector);
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*CenterVector';
                countloop = countloop + 1;
                seednow = round(p(seed_loop).coord(:,i));
                if( countloop >= 50 );   break;   end 
            end
        end
    end
 % Injection
    for seed_loop = n_seeds(1,j)+1:n_seeds(1,j+1)
        seed = seed_ind(ceil(rand*length(seed_ind))); % Choose a random point index to start the simulation
        [seed(1,1),seed(1,2),seed(1,3)] = ind2sub(size(coh_map),seed);
        if seed(3) > ubound(seed(1),seed(2))
            seed(3) = ubound(seed(1),seed(2));
        elseif seed(3) < lbound(seed(1),seed(2))
            seed(3) = lbound(seed(1),seed(2));
        end
        p(seed_loop).coord(:,i) = seed'; % Coordinate of this starting point
        dmax(1,seed_loop) = betainv(rand(1), 1, beta4dist); 
    end
%%%% First step with last injection
    i = Finaltimestep(1,j) + 2;

    seedposition = zeros(3, n_seeds(1,j+1));
    for seed_loop = 1:n_seeds(1,j+1)
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    [pde_state, concentration] = step_pde_concentration(op, pde_state, seedposition', dt_pde);
    concentration(coh_map == 0) = 0;
    % ------------------------------
    save_concentration_slice_if_needed(i, concentration, tumor_center, ...
        save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);
    % ------------------------------
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);

    for seed_loop = 1:n_seeds(1,j)
        %%% Find index
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
               
        if seedx < 1 || seedx > size(coh_map,1) || ...
           seedy < 1 || seedy > size(coh_map,2) || ...
           seedz < 1 || seedz > size(coh_map,3) || ...
           isnan(seedx) || isnan(seedy) || isnan(seedz)

            % Put the cell back to the last valid location and skip
           p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
           continue;
        end

        ind = sub2ind(size(coh_map), seedx, seedy, seedz);

        %%% Calculate eigen vector, step length, and prefered direction +/-
        if coh_map(ind) > coh_limit % Get Eigen vector at the previous coordinate point
            eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) );
            d = dmax(1,seed_loop) * d_w;
            pref_mvmt = sign(dot( p(seed_loop).coord(:,i-1) - p(seed_loop).coord(:,i-2), eigen_vect) + eps );   % Move either in the +/- direction
            d_noise = d_noise_value;
        else
            eigen_vect = 0;
            d = dmax(1, seed_loop) * d_g; % Move with less magnitude if no WM present 
            pref_mvmt = 1;
            d_noise = 1;
        end
        eigen_vect_noise = (rand(3,1)*2-1);
        eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);

        %%% Keep or reduce step length by some amount
        %%%%%% uniform[0, 1] 
%         dstep = rand(1)*dmax ; 
        %%%%%% beta[1, beta], if beta==1, uniform 
%         dstep = betainv( rand(1), 1, beta4dist )*d ; 
        %%%%%% deterministic
        dstep = d;

        %%% Calculate chemotaxis vector and gradient strength
        chmtx_vect = [cgradX(ind); cgradY(ind); cgradZ(ind)];
        chmtx_bias = 0;
        switch modelNum
            case 1  % Chemotaxis
                % Move along chemotaxis (deterministic)
                dstep = dstep/dmax(1,seed_loop); % Keeps dstep deterministic (dmax is stochastic)
                if has_cancer
                    chmtx_bias = chemo_sensitivity;
                end
            case 2  % Chemotaxis, Stochasticity
                %  Move along chemotaxis (stochasticity)
                if has_cancer
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
            case 3  % Chemotaxis, Chemotax threshold, Stochasticity
                %  Move along chemotaxis once near cancer (stochasticity)
                if has_cancer && (concentration(ind) >= chmtx_limit)
                    chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                end
        end

        %%% Take step: (curr coord) = (prev coord) + (coherency influence) + (cancer influence)
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);

        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);
            seednow = round(p(seed_loop).coord(:,i));
            countloop = 1;
            while ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                Center_Point = Boundary_Ind(ceil(rand*length(Boundary_Ind)));
                [Center_Point(1,1),Center_Point(1,2),Center_Point(1,3)] = ind2sub(size(coh_map),Center_Point);
                Center_Vector = Center_Point - p(seed_loop).coord(:,i)';
                CenterVector = Center_Vector./norm(Center_Vector);
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*CenterVector';
                countloop = countloop + 1;
                seednow = round(p(seed_loop).coord(:,i));
                if( countloop >= 50 );   break;   end 
            end
        end
    end
  
%Run injected cells time 2
    for seed_loop = n_seeds(1,j)+1:n_seeds(1,j+1)
        eigen_vect_noise = [rand(1)*2-1, rand(1), rand(1)]';
        eigen_vect = eigen_vect_noise/norm(eigen_vect_noise);
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(d_g*eigen_vect); % Find next coordinate at step d along EV direction
        seed = round( p(seed_loop).coord(:,i) );
        if seed(3) >= ubound(seed(1),seed(2)) - nbuffer || seed(3) <= lbound(seed(1),seed(2)) + nbuffer ...
            || isnan(ubound(seed(1),seed(2))) || isnan(lbound(seed(1),seed(2)))   
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
        end
    end

%%%% Algorithm: while NSC stays in bounds

    % This loop should also have a PDE update for time
    for i = Finaltimestep(1,j)+3:Finaltimestep(1,j+1) % Run the loop a large number of steps 
        seedposition = zeros(3, n_seeds(1,j+1));
        for seed_loop = 1:n_seeds(1,j+1)
            seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
            seedy = round( p(seed_loop).coord(2,i-1) );
            seedz = round( p(seed_loop).coord(3,i-1) );
            seedposition(:, seed_loop) = [seedx, seedy, seedz]';
        end

        [pde_state, concentration] = step_pde_concentration(op, pde_state, seedposition', dt_pde);
        concentration(coh_map == 0) = 0;
        % ------------------------------
        save_concentration_slice_if_needed(i, concentration, tumor_center, ...
            save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag);
        % ------------------------------
        [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);


        for seed_loop = 1:n_seeds(1,j+1)
            %%% Find index
            seedx = round( p(seed_loop).coord(1,i-1) );  % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
            seedy = round(p(seed_loop).coord(2,i-1));
            seedz = round( p(seed_loop).coord(3,i-1) );
                    % --- Add this bounds check ---
            if seedx < 1 || seedx > size(coh_map,1) || ...
                seedy < 1 || seedy > size(coh_map,2) || ...
                seedz < 1 || seedz > size(coh_map,3) || ...
                isnan(seedx) || isnan(seedy) || isnan(seedz)

                % Put the cell back to the last valid location and skip
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
                continue;
            end

            ind = sub2ind(size(coh_map), seedx, seedy, seedz);

            %%% Calculate eigen vector, step length, and prefered direction +/-
            if coh_map(ind) > coh_limit % Get Eigen vector at the previous coordinate point
                eigen_vect = squeeze( eigen_map( seedx, seedy, seedz, : ) );
                d = dmax(1,seed_loop) * d_w;
                pref_mvmt = sign(dot( p(seed_loop).coord(:,i-1) - p(seed_loop).coord(:,i-2), eigen_vect) + eps ); % Move either in the +/- direction
                d_noise = d_noise_value;
            else  
                eigen_vect = 0;
                d = dmax(1, seed_loop) * d_g; % Move with less magnitude if no WM present 
                pref_mvmt = 1;
                d_noise = 1;
            end
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);

            %%% Keep or reduce step length by some amount
            %%%%%% uniform[0, 1] 
    %         dstep = rand(1)*dmax ; 
            %%%%%% beta[1, beta], if beta==1, uniform 
    %         dstep = betainv( rand(1), 1, beta4dist )*d ; 
            %%%%%% deterministic
            dstep = d;

          %%% Calculate chemotaxis vector and gradient strength
            chmtx_vect = [cgradX(ind); cgradY(ind); cgradZ(ind)];
            chmtx_bias = 0;
            switch modelNum
                case 1  % Chemotaxis
                    % Move along chemotaxis (deterministic)
                    dstep = dstep/dmax(1,seed_loop); % Keeps dstep deterministic (dmax is stochastic)
                    if has_cancer
                        chmtx_bias = chemo_sensitivity;
                    end
                case 2  % Chemotaxis, Stochasticity 
                    %  Move along chemotaxis (stochasticity)
                    if has_cancer
                        chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                    end
                case 3  % Chemotaxis, Chemotax threshold, Stochasticity
                %  Move along chemotaxis once near cancer (stochasticity)
                    if has_cancer && (concentration(ind) >= chmtx_limit)
                        chmtx_bias = betainv(rand(1),alpha4chmtx,1) * chemo_sensitivity;
                    end
            end

            %%% Take step: (curr coord) = (prev coord) + (coherency influence) + (cancer influence)
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);

            %%% If cell is out of bounds, put it back to previous coordinate
            seednow = round(p(seed_loop).coord(:,i));
            if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                eigen_vect_noise = (rand(3,1)*2-1);
                eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + dstep*d_noise*eigen_vect_noise);
                seednow = round(p(seed_loop).coord(:,i));
                countloop = 1;
                while ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                    ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                    ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                    Center_Point = Boundary_Ind(ceil(rand*length(Boundary_Ind)));
                    [Center_Point(1,1),Center_Point(1,2),Center_Point(1,3)] = ind2sub(size(coh_map),Center_Point);
                    Center_Vector = Center_Point - p(seed_loop).coord(:,i)';
                    CenterVector = Center_Vector./norm(Center_Vector);
                    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + dstep*CenterVector';
                    countloop = countloop + 1;
                    seednow = round(p(seed_loop).coord(:,i));
                    if( countloop >= 50 );   break;   end 
                end
            end
        end
    end
end

daily_arrivals = zeros(Doses, 1);
already_arrived = false(1, n_seeds_current);
for dose_loop = 1:Doses
    active_count = n_seeds(1,dose_loop);
    time_idx = Finaltimestep(1,dose_loop);
    arrived_now = false(1, n_seeds_current);
    for seed_loop = 1:active_count
        if size(p(seed_loop).coord, 2) >= time_idx
            arrived_now(seed_loop) = any(all(abs(cancer_center - p(seed_loop).coord(:,time_idx)') < cancer_size, 2));
        end
    end
    daily_arrivals(dose_loop) = sum(arrived_now & ~already_arrived);
    already_arrived = already_arrived | arrived_now;
end

toc( Tstart );
end

% Functions to save concentration slices each day
function save_concentration_slice_if_needed(i, concentration, tumor_center, ...
        save_concentration_slices, concentration_slice_times, concentration_slice_dir, run_tag)
if ~save_concentration_slices || ~any(i == concentration_slice_times)
    return;
end
save_concentration_slice(i, concentration, tumor_center, concentration_slice_dir, run_tag);
end

% just save images:
function save_concentration_slice(i, concentration, tumor_center, concentration_slice_dir, run_tag)
% save_concentration_slice saves an x-y concentration heat map through tumor center z-coordinate
if ~exist(concentration_slice_dir, 'dir')
    mkdir(concentration_slice_dir);
end

safe_tag = regexprep(run_tag, '[^A-Za-z0-9_.-]', '_');
z_slice = round(tumor_center(3));
z_slice = max(1, min(size(concentration, 3), z_slice));

concentration_slice = squeeze(concentration(:,:,z_slice));
filename = fullfile(concentration_slice_dir, ...
    sprintf('concentration_slice_z%d_i%04d_%s.png', z_slice, i, safe_tag));

fig = figure('Visible', 'off');
imagesc(concentration_slice');
set(gca, 'YDir', 'normal');
axis image;
colorbar;
colormap(parula);
hold on;
plot(tumor_center(1), tumor_center(2), 'wo', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title(sprintf('Concentration, z = %d, i = %d', z_slice, i));
saveas(fig, filename);
close(fig);
end