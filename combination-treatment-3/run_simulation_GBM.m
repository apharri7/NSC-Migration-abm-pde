function [p, total_dead, n_seeds_current, daily_arrivals, pde_state_vegf, pde_state_mcp1] = run_simulation_GBM(alpha4chmtx, beta4dist, Boundary_Ind, cancer_center, cancer_size, ... 
    op_vegf, pde_state_vegf, op_mcp1, pde_state_mcp1, dt_pde, ...
    chemo_sensitivity, chmtx_limit, coh_limit, coh_map, CONVERT2MICRON, d_g, d_w, d_noise_value, dmax, Doses, ...
    dtime, eigen_map, Finaltimestep, has_cancer, inj_center, lbound, modelNum, ...
    modelType, nbuffer, n_seeds, num_cancers, p, seed_ind, seed_sd, ubound)

% File is set up for 3 radiation fractions. Comments are added for
% additions for 5 fractions or 1 fraction
                        
%% For each seed
Tstart = tic;
disp('Running simulation...')

% LQ Elimination parameters
tumor_center = cancer_center(1,:)';
tumor_radius = 20;
alpha = 0.1; beta = 0.03;
dose = 2;
survival_prob = exp(-alpha * dose - beta * dose^2); 
total_dead = zeros(1,3);
total_alive = 0;
n_seeds_current = 0; % stores the current number of seeds

[Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));

    %%%% Spawn seed within injection site
i = 1;
for seed_loop = 1:n_seeds(1,1)
    %seed = seed_ind(ceil(rand*length(seed_ind))); % Choose a random point index to start the simulation
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

n_seeds_current = numel(p);

% Initial VEGF and MCP-1 fields come from init_pde_state.m.
concentration = pde_state_vegf.C_full;
concentration2 = pde_state_mcp1.C_full;

% FIRST STEP: SHOULD DETERMINE WHICH CELLS SURVIVE HERE AND UPDATE p and dmax accordingly nad num_current

    %%%% First step
i = 2;

for seed_loop = 1:n_seeds(1,1) % n_seeds(1,1) SHOULD BE CHANGED TO NUM_CURRENT
    eigen_vect_noise = [rand(1)*2-1, rand(1), rand(1)]';
    eigen_vect = eigen_vect_noise/norm(eigen_vect_noise);
    p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(d_g*eigen_vect); % Find next coordinate at step d along EV direction
    seed = round( p(seed_loop).coord(:,i) );
    if seed(3) >= ubound(seed(1),seed(2)) - nbuffer || seed(3) <= lbound(seed(1),seed(2)) + nbuffer ...
        || isnan(ubound(seed(1),seed(2))) || isnan(lbound(seed(1),seed(2)))   
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1);
    end
    
end



% SECOND STEP: SHOULD DETERMINE WHICH CELLS SURVIVE HERE AND UPDATE p and dmax accordingly nad num_current

% Matrix to store the position of the n_seeds for Internalization term
%seedposition = zeros(3, n_seeds_current); % MAY NEED TO MVE THIS INSIDE TIME LOOP

% NSC ARRIVAL TRACKING
has_arrived = false(1, n_seeds_current);    % For tracking NSC arrival, needs reset when n_seeds_current changees
alive_flags = true(1, n_seeds_current);
daily_arrivals = zeros(1, 17);               % 17 days of simulation, no change needed
day_index = 1;                               % Track which day we're on, shifts in outer time step loop

 %%%% Algorithm: while NSC stays in bounds
for i = 3:Finaltimestep(1,1) % Finaltimestep(1,1) runs to 600, so at last step update surv prob to n=3 fractions of radiation after this loop (runs to 599)

    % Time loop is where we should solve pde
    seedposition = zeros(3, n_seeds_current);
    
    
    
    
    % Find seeds current position (time i-1)
    for seed_loop = 1:n_seeds_current
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    % Advance VEGF and MCP-1 with the same implicit PDE solver.
    [pde_state_vegf, concentration] = step_pde_concentration(op_vegf, pde_state_vegf, seedposition', dt_pde);
    [pde_state_mcp1, concentration2] = step_pde_concentration(op_mcp1, pde_state_mcp1, seedposition', dt_pde);

    % Apply radiotherapy-driven MCP-1 boost after the PDE step and keep the PDE state synchronized.
    [concentration2, pde_state_mcp1] = apply_mcp1_boost_if_needed(i, concentration2, tumor_center, tumor_radius, op_mcp1, pde_state_mcp1);

    % Compute concentration gradients using physical concentration scales.
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);
    [cgradX2, cgradY2, cgradZ2] = concentration_gradients2(110*concentration2, coh_map);

    
    % Move the seeds
    for seed_loop = 1:n_seeds(1,1)   % (NEED TO CHANGE n_seeds(1,1) TO HOW MANY SEEDS ARE LEFT after ELIM)
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
        % Second chemotaxis vector
        chmtx_vect2 = [cgradX2(ind); cgradY2(ind); cgradZ2(ind)];
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
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);
        
        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);
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
        
        % Check for unique seed arrival at each step, sum them over a full day, and put the sum in the index according to day (arrivals matrix 1xnumdays)
        
        dist_to_tumor = norm(p(seed_loop).coord(:,i) - tumor_center); 
        
        if ~has_arrived(seed_loop) && dist_to_tumor <= tumor_radius
            has_arrived(seed_loop) = true;
            daily_arrivals(day_index) = daily_arrivals(day_index) + 1;
        end
        
    end
    
   % Update days
    if mod(i, 200) == 0  % if new day starts (e.g. i = 201, 401, etc.)
        day_index = day_index + 1;
    end 
    
end



n_seeds_current = numel(p);



%% run loop for doses
for j = 1:Doses-1
    seedposition = zeros(3, n_seeds_current); % HERE WE NEED TO CHANGE n_seeds(1,j) TO NUM CURRENT CELLS
    i = Finaltimestep(1,j) + 1;
    % Update concentration based on pde, update based  on j ------------
    

    
    % Find seeds current position (time i-1)
    for seed_loop = 1:n_seeds_current
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    % Advance VEGF and MCP-1 with the same implicit PDE solver.
    [pde_state_vegf, concentration] = step_pde_concentration(op_vegf, pde_state_vegf, seedposition', dt_pde);
    [pde_state_mcp1, concentration2] = step_pde_concentration(op_mcp1, pde_state_mcp1, seedposition', dt_pde);

    % Apply radiotherapy-driven MCP-1 boost after the PDE step and keep the PDE state synchronized.
    [concentration2, pde_state_mcp1] = apply_mcp1_boost_if_needed(i, concentration2, tumor_center, tumor_radius, op_mcp1, pde_state_mcp1);

    % Compute concentration gradients using physical concentration scales.
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);
    [cgradX2, cgradY2, cgradZ2] = concentration_gradients2(110*concentration2, coh_map);

    for seed_loop = 1:n_seeds_current  % AGAIN, CHANGE n_seeds(1,j) to NUM CURRENT CELL
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
        % Second chemotaxis vector
        chmtx_vect2 = [cgradX2(ind); cgradY2(ind); cgradZ2(ind)];
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
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);

        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect + chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);
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
        
        % Check for unique seed arrival at each step, sum them over a full day, and put the sum in the index according to day (arrivals matrix 1xnumdays)
        
        dist_to_tumor = norm(p(seed_loop).coord(:,i) - tumor_center); 
        
        if ~has_arrived(seed_loop) && dist_to_tumor <= tumor_radius
            has_arrived(seed_loop) = true;
            daily_arrivals(day_index) = daily_arrivals(day_index) + 1;
        end  
         
    end
    
    % HERE WE NEED TO UPDATE P AND D_MAX AND NUM CURRENT TO ELIM USING LQ
    % Update days
    if mod(i, 200) == 0  % if new day starts (e.g. i = 201, 401, etc.)
        day_index = day_index + 1;
    end

    
    
 % Injection, THIS JUST ADDS NEW CELL LOCATIONS
    for seed_loop = n_seeds_current+1:n_seeds(1,j+1) % HERE WE NEED TO CHANGE n_seeds(1,j) to NUM CURRENT
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
    
    % update size of has_arrived to match n_seeds current
    num_new_seeds = n_seeds(1,j+1) - n_seeds_current;
    has_arrived = [has_arrived, false(1, num_new_seeds)];
    
    % this should update the number of seeds from injection
    n_seeds_current = numel(p);
    
%%%% First step with last injection
    i = Finaltimestep(1,j) + 2;
    % Update concentration ---------------------------
    
    % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
    if i == 801 || i == 1001 || i == 1201
        alive_flags = true(1, n_seeds_current);
    end
    
    seedposition = zeros(3, n_seeds_current);
    
    % Find seeds current position (time i-1)
    for seed_loop = 1:n_seeds_current
        seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
        seedy = round( p(seed_loop).coord(2,i-1) );
        seedz = round( p(seed_loop).coord(3,i-1) );
        seedposition(:, seed_loop) = [seedx, seedy, seedz]';
    end

    % Advance VEGF and MCP-1 with the same implicit PDE solver.
    [pde_state_vegf, concentration] = step_pde_concentration(op_vegf, pde_state_vegf, seedposition', dt_pde);
    [pde_state_mcp1, concentration2] = step_pde_concentration(op_mcp1, pde_state_mcp1, seedposition', dt_pde);

    % Apply radiotherapy-driven MCP-1 boost after the PDE step and keep the PDE state synchronized.
    [concentration2, pde_state_mcp1] = apply_mcp1_boost_if_needed(i, concentration2, tumor_center, tumor_radius, op_mcp1, pde_state_mcp1);

    % Compute concentration gradients using physical concentration scales.
    [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);
    [cgradX2, cgradY2, cgradZ2] = concentration_gradients2(110*concentration2, coh_map);

    for seed_loop = 1:n_seeds_current % AND AGAIN n_seeds(1,j) SHOULD BE NUM CURRENT CELLS
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
        % Second chemotaxis vector
        chmtx_vect2 = [cgradX2(ind); cgradY2(ind); cgradZ2(ind)];
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
        p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);

        %%% If cell is out of bounds, put it back to previous coordinate
        seednow = round(p(seed_loop).coord(:,i));
        if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
            ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
            ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
            eigen_vect_noise = (rand(3,1)*2-1);
            eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);
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
        
         % CHECK FOR ARRIVAL HERE
        
        dist_to_tumor = norm(p(seed_loop).coord(:,i) - tumor_center); 
        
        if ~has_arrived(seed_loop) && dist_to_tumor <= tumor_radius
            has_arrived(seed_loop) = true;
            daily_arrivals(day_index) = daily_arrivals(day_index) + 1;
        end
        
       % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
        if i == 801 || i == 1001 || i == 1201 
            %%%%% LQ Addition! Check for cell survival
            cell_position = round(p(seed_loop).coord(:,i));
            dist_to_tumor = norm(cell_position - tumor_center);
    
            if dist_to_tumor <= tumor_radius && rand > survival_prob
                alive_flags(seed_loop) = false;
                continue;
            end
        end
        
        
        
        
    
    end
    
    % Update day for tracking arrivals per day

    if mod(i, 200) == 0  % if new day starts (e.g. i = 201, 401, etc.)
        day_index = day_index + 1;
    end

    % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
    if i == 801
        % update p, d_max, n_seeds_current
        total_dead(1,1) = total_dead(1,1) + sum(~alive_flags);
        p = p(alive_flags);   % remove dead seeds
        dmax = dmax(:, alive_flags); % update dmax if it's per-seed
        n_seeds_current = numel(p);
        % remove dead indices from has_arrived
        has_arrived = has_arrived(alive_flags);
    end
    
    if i == 1001
        % update p, d_max, n_seeds_current
        total_dead(1,2) = total_dead(1,2) + sum(~alive_flags);
        p = p(alive_flags);   % remove dead seeds
        dmax = dmax(:, alive_flags); % update dmax if it's per-seed
        n_seeds_current = numel(p);
        % remove dead indices from has_arrived
        has_arrived = has_arrived(alive_flags);
    end
    
    if i == 1201
        % update p, d_max, n_seeds_current
        total_dead(1,3) = total_dead(1,3) + sum(~alive_flags);
        p = p(alive_flags);   % remove dead seeds
        dmax = dmax(:, alive_flags); % update dmax if it's per-seed
        n_seeds_current = numel(p);
        % remove dead indices from has_arrived
        has_arrived = has_arrived(alive_flags);
    end
    
    
    
    
 
    
%Run injected cells time 2
    for seed_loop = 1:n_seeds_current % THIS SHOULD JUST BE NUM_CURRENT!
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
        % Update concentration in this loop -------------------------------
        
        % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
        if i == 801 || i == 1001 || i == 1201 
            alive_flags = true(1, n_seeds_current);
        end
       
        seedposition = zeros(3, n_seeds_current);

        % Find seeds current position (time i-1)
        for seed_loop = 1:n_seeds_current
            seedx = round( p(seed_loop).coord(1,i-1) ); % Initialize (seedx,seedy,seedz) as (x,y,z) coords of prev step
            seedy = round( p(seed_loop).coord(2,i-1) );
            seedz = round( p(seed_loop).coord(3,i-1) );
            seedposition(:, seed_loop) = [seedx, seedy, seedz]';
        end

        % Advance VEGF and MCP-1 with the same implicit PDE solver.
        [pde_state_vegf, concentration] = step_pde_concentration(op_vegf, pde_state_vegf, seedposition', dt_pde);
        [pde_state_mcp1, concentration2] = step_pde_concentration(op_mcp1, pde_state_mcp1, seedposition', dt_pde);

        % Apply radiotherapy-driven MCP-1 boost after the PDE step and keep the PDE state synchronized.
        [concentration2, pde_state_mcp1] = apply_mcp1_boost_if_needed(i, concentration2, tumor_center, tumor_radius, op_mcp1, pde_state_mcp1);

        % Compute concentration gradients using physical concentration scales.
        [cgradX, cgradY, cgradZ] = concentration_gradients(5000*concentration, coh_map);
        [cgradX2, cgradY2, cgradZ2] = concentration_gradients2(110*concentration2, coh_map);
        
        

        for seed_loop = 1:n_seeds_current % HERE n_seeds(1,j+1) SHOULD BE NUM CURRENT
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
            % Second chemotaxis vector
            chmtx_vect2 = [cgradX2(ind); cgradY2(ind); cgradZ2(ind)];
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
            p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);

            %%% If cell is out of bounds, put it back to previous coordinate
            seednow = round(p(seed_loop).coord(:,i));
            if ~( all(seednow' - size(coh_map)<=-1) && all( seednow'>1 ) ) || ...                                           %Seed reaches edge of coh_map
                ( seednow(3) > ubound(seednow(1),seednow(2)) || seednow(3) < lbound(seednow(1),seednow(2)) ) || ...         %Seed is outside of the lower or upper bounds (z)
                ( isnan(ubound(seednow(1),seednow(2))) || isnan(lbound(seednow(1),seednow(2))) )                            %Seed is outside of the side bounds of brain (x,y)
                eigen_vect_noise = (rand(3,1)*2-1);
                eigen_vect_noise = eigen_vect_noise/norm(eigen_vect_noise);
                p(seed_loop).coord(:,i) = p(seed_loop).coord(:,i-1) + 1/dtime*(-pref_mvmt*dstep*eigen_vect + chmtx_bias*chmtx_vect+ chmtx_bias*chmtx_vect2 + dstep*d_noise*eigen_vect_noise);
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
            
            % TRACK ARRIVAL (non-repeated) before DEATH
            
            % CHECK FOR ARRIVAL HERE
        
            dist_to_tumor = norm(p(seed_loop).coord(:,i) - tumor_center); 
        
            if ~has_arrived(seed_loop) && dist_to_tumor <= tumor_radius
                has_arrived(seed_loop) = true;
                daily_arrivals(day_index) = daily_arrivals(day_index) + 1;
            end
          
          
            % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
            if i == 801 || i == 1001 || i == 1201 
            %%%%% LQ Addition! Check for cell survival
                cell_position = round(p(seed_loop).coord(:,i));
                dist_to_tumor = norm(cell_position - tumor_center);
    
                if dist_to_tumor <= tumor_radius && rand > survival_prob
                    alive_flags(seed_loop) = false;
                    continue;
                end
            end 
           
        end
        
        % Update day for tracking arrivals per day

        if mod(i, 200) == 0  % if new day starts (e.g. i = 201, 401, etc.)
            day_index = day_index + 1;
        end
        
        % add i==1601 for 5 fractions or remove 1001, 1201 for 1 fraction
        if i == 801
            % update p, d_max, n_seeds_current
            total_dead(1,1) = total_dead(1,1) + sum(~alive_flags);
            p = p(alive_flags);   % remove dead seeds
            dmax = dmax(:, alive_flags); % update dmax if it's per-seed
            n_seeds_current = numel(p);
            % remove dead indices from has_arrived
            has_arrived = has_arrived(alive_flags);
        end
    
        if i == 1001
            % update p, d_max, n_seeds_current
            total_dead(1,2) = total_dead(1,2) + sum(~alive_flags);
            p = p(alive_flags);   % remove dead seeds
            dmax = dmax(:, alive_flags); % update dmax if it's per-seed
            n_seeds_current = numel(p);
            % remove dead indices from has_arrived
            has_arrived = has_arrived(alive_flags);
        end
    
        if i == 1201
            % update p, d_max, n_seeds_current
            total_dead(1,3) = total_dead(1,3) + sum(~alive_flags);
            p = p(alive_flags);   % remove dead seeds
            dmax = dmax(:, alive_flags); % update dmax if it's per-seed
            n_seeds_current = numel(p);
            % remove dead indices from has_arrived
            has_arrived = has_arrived(alive_flags);
        end
    
        
    end
end

total_alive = total_alive + sum(alive_flags);

toc( Tstart );
end

function [concentration2, pde_state_mcp1] = apply_mcp1_boost_if_needed(i, concentration2, tumor_center, tumor_radius, op_mcp1, pde_state_mcp1)
boost_factor = 1;

if i == 801
    boost_factor = 1.699427737;
    % remove 1201 for 1 fraction
elseif i == 1201
    boost_factor = 1.033365813;
 % add boost at 1601 for 5 fractions
% elseif i == 1601
    % boost_factor = 1.307860307;
end

if boost_factor ~= 1
    [Xgrid, Ygrid, Zgrid] = ndgrid(1:size(concentration2,1), ...
                                   1:size(concentration2,2), ...
                                   1:size(concentration2,3));
    distance_from_tumor = sqrt((Xgrid - tumor_center(1)).^2 + ...
                               (Ygrid - tumor_center(2)).^2 + ...
                               (Zgrid - tumor_center(3)).^2);
    tumor_mask = (distance_from_tumor <= tumor_radius);

    concentration2(tumor_mask) = boost_factor * concentration2(tumor_mask);
    pde_state_mcp1 = apply_mcp1_boost_to_state(op_mcp1, pde_state_mcp1, tumor_mask, boost_factor);

    pde_state_mcp1.C_full = concentration2;
end
end

function state = apply_mcp1_boost_to_state(op, state, tumor_mask_full, boost_factor)
ncx = numel(op.ix_full);
ncy = numel(op.iy_full);
ncz = numel(op.iz_full);

% Crop the full-grid tumor mask in the same way as the PDE domain.
tumor_mask_crop = tumor_mask_full( ...
    op.ix_full, op.iy_full, op.iz_full);

% Insert the cropped mask into the padded PDE grid.
tumor_mask_local = false(op.nx, op.ny, op.nz);
tumor_mask_local( ...
    op.pad+(1:ncx), ...
    op.pad+(1:ncy), ...
    op.pad+(1:ncz)) = tumor_mask_crop;

% Apply the boost only within the brain
tumor_mask_local = tumor_mask_local & ...
                   op.brain_mask & ...
                   op.solve_mask;

% 3-D mask to the reduced PDE ordering
tumor_mask_r = tumor_mask_local(op.dof);

% Multiplying q by the boost factor multiplies C by the same factor
% because q = phi*C and phi does not change
state.q_r(tumor_mask_r) = ...
    boost_factor .* state.q_r(tumor_mask_r);

% Update the next solver guess
state.q_lo_guess = state.q_r;

state.mass0 = sum(state.q_r);
end
