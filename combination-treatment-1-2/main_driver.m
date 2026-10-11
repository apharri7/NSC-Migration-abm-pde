%% NSC Migration in Response to GBM

% Loads the DTI data. Defines the ABM and PDE parameters.
% Initializes PDE solver, runs the coupled model simulation, 
% and saves summary statistics and plots in the current directory.

%% Load data


disp('Gathering data...')
% EV_file = ... write path to DTI eigenvector image
% FA_file = ... write path to fractional anisotropy image
[eigen_map, coh_map] = load_3D_run_tumor(EV_file, FA_file); 


%% Set ABM parameters 

% 300 NSCs injected at start of simulation. Additional doses of 300 cells 
% injected on days 3,5, 7, 10, 12.   17 days simulated.
n_seeds = [300 600 900 1200 1500 1800];     % NSC counts 
Finaltimestep = [3 5 7 10 12 17];           % time (Days)

timeinterval = 200;                         % Number of ABM/PDE updates per simulated day
Doses = numel(n_seeds);
Finaltimestep = Finaltimestep*timeinterval; 
Finaltimestep(1:end-1) = Finaltimestep(1:end-1) - 1;  
CONVERT2MICRON =    54;                     % Average micrometers per one voxel
coh_limit =         0.4;                    % White matter coherency threshold
chmtx_limit =       0.1;                    % Chemotaxis threshold
d_w =               2;                      % Reference step size on white matter
d_g =             1;                        % Reference step size on grey matter
cancer_size =       [20, 20, 20];           % Half-widths of tumor arrival box
chemo_sensitivity = 1;                    % chemoattractant sensitivity level
%%% Injection Center 
modelType =         'Intranasal';           % Determines injection center
% modelType =         'Intracerebral'; 
%%% Cancer Location
% cancer_center =     [120, 400, 150];      % Right, back putamen tumor center 
cancer_center = [79 84 78];                 % Right, front putamen tumor center
%%% Conditional
has_cancer =        1;
modelNum =          2;                      % Chooses the stochastic chemotaxis model
alpha4chmtx = 1;                            % Stochastic chemosensitivity parameter

beta4dist =         1;                      % Movement scale parameter
p = struct('coord',{});                     % Stores each NSC trajectory as 3xtimepts array
dmax = zeros(1, n_seeds(1,Doses));          % per-NSC stochastic reference step multiplier
dtime = 1;                                  % ABM movement divisor - 1 leaves unchanges
nbuffer = 5;                                % Keeps injected NSCs away from upper and lower z bounds
d_noise_value = .1/d_w;                     % White matter directional noise coefficient


% Generates injection-site indices and upper and lower brain bounds
[seed_ind,seed_sd,ubound,lbound,inj_center] = set_initial_GBM(coh_map, modelType, Doses, n_seeds); 

% Identify boundary region to redirect NSC outside the brain 
[num_cancers, Boundary_Ind] = find_boundary(coh_map, cancer_center);

% One-time setup of the PDE solver for chemoattractant VEGF

% dt_pde. One PDE step happens per ABM update (200 per day)
dt_pde = 1/timeinterval;

pde_params = struct();                  

% Production term parameters
pde_params.chi_tumor = 24;              % tumor-localized production rate
pde_params.chi_bg_gaussian = 0.5;       % weak background logistic for Gaussian i.c.
pde_params.chi_bg_constant = 24;        % full background logistic for constant i.c.

% Tumor-source and NSC-sink parameters
pde_params.r_tumor     = 20;          % Tumor source radius (voxels)
pde_params.omega_tumor = 1;           % Transition width of tumor window
pde_params.R_nsc       = 7.5/13.5;      % NSC sink radius (voxels)
pde_params.lambda_sink = 100*24;        % NSC internalization rate

% Diffuse domain reduction parameters
pde_params.phi_floor = 1e-12;
pde_params.phi_cutoff = 1e-8;

% Initial-condition selection
vegf_ic_type = 'gaussian';              % constant or gaussian
add_vegf_noise = true;                  

pde_params.ic_type = vegf_ic_type;
pde_params.ic_const     = 0.6094;   % Normalized constant initial profile   
pde_params.ic_amp      = 0.6294;      % Gaussian initial profile peak concentration (normalized)
pde_params.ic_sd       = 67.9025;     % Gaussian initial profile standard deviation (voxels)

% Select background logistic production based on initial condition type
if strcmpi(pde_params.ic_type, 'constant')
    pde_params.chi_bg = pde_params.chi_bg_constant;        
else
    pde_params.chi_bg = pde_params.chi_bg_gaussian;         
end

% Adds spatially correlated initial profile noise in full brain
if add_vegf_noise
    if strcmpi(vegf_ic_type, 'constant')
        pde_params.ic_noise_full = generate_3D_noise_constant(coh_map);
    else
        pde_params.ic_noise_full = generate_3D_noise_gaussian(coh_map);
    end
    pde_params.ic_noise_scale = 5000;                                   % VEGF concentration scaling 
end

% Diffuse domain method parameters
pde_coh_limit = coh_limit;                                        
pde_epsilon   = 2;                              % Diffuse interface thickness (voxels), 2 used for simulations
pde_crop_half = 5;                              % Unused while pde_use_full_domain = true
pde_pad       = 15;                              % Zero padding around brain to create computational domain
pde_use_full_domain = false;                    % Set to 'false' for faster cropped tests, or 'true' for the full run

% Setup of all fixed PDE operators one time and initializes q = phi*C
op = setup_pde_operator(EV_file, FA_file, cancer_center, ...
    pde_coh_limit, pde_epsilon, pde_crop_half, pde_pad, pde_use_full_domain, pde_params);          
pde_state = init_pde_state(op, dt_pde);                                                            

% Run coupled ABM-PDE simulation
[p, total_dead, n_seeds_current, daily_arrivals, pde_state] = run_simulation_GBM(alpha4chmtx, beta4dist, Boundary_Ind, cancer_center, cancer_size, ... 
    op, pde_state, dt_pde, ...
    chemo_sensitivity, chmtx_limit, coh_limit, coh_map, CONVERT2MICRON, d_g, d_w, d_noise_value, dmax, Doses, ...
    dtime, eigen_map, Finaltimestep, has_cancer, inj_center, lbound, modelNum, ...
    modelType, nbuffer, n_seeds, num_cancers, p, seed_ind, seed_sd, ubound);                       

  
% Save simulation summaries.                       
% Save total dead cells (0 unless radiotherapy treatment applied)                  
filename = 'total_dead.txt';
fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file %s for writing.', filename);
end
fprintf(fid, 'Total dead seeds: ');
fprintf(fid, '%d ', total_dead);
fprintf(fid, '\n');
fclose(fid);  

% Save daily arrivals to file
filename = 'daily_arrivals.txt';
fid = fopen(filename, 'w');
if fid == -1
    error('Failed to open file %s for writing.', filename);
end
fprintf(fid, '%d\n', daily_arrivals);
fclose(fid);


%% Downsample trajectories for plotting

disp('Plotting graphs...')


% Get the NSC positions at every 100 time steps
acc = 100;                                                                                                 
p_original = p;
p_new = struct('coord',{});
for i=1:n_seeds_current
    p_new(i).coord = p(i).coord(:,1:acc:Finaltimestep(1,Doses)); 
end
p = p_new; clear p_new;
xvals = (1:acc:Finaltimestep(1, Doses));

% NSC injected after first dose have leading zero coordinates before their injection time.
% Replace placeholders with downsamples copy so we do not get artificial paths from origin 
% to injection site
for j = 1:n_seeds_current
    Injection_Loc = find(p(j).coord(1,:)== 0,1,'Last') + 1;
    for k = 1:Injection_Loc 
        p(j).coord(:,k) = p(j).coord(:,Injection_Loc);
    end    
end    


%% Plot downsampled path of seeds
figure;
plotenvironment_run_tumor(.7,coh_map, ubound, lbound, seed_sd, inj_center, cancer_size, cancer_center, has_cancer, num_cancers);            

for i=1: n_seeds_current
    plot3(p(i).coord(1,:),p(i).coord(2,:),p(i).coord(3,:),'linewidth',3);                                       %Plot path
    plot3(p(i).coord(1, end),p(i).coord(2, end),p(i).coord(3, end),'ok','Markersize', 6,'markerfacecolor','w'); %Plot final coords
end

filename = 'path.png';
saveas(gcf, filename);

%% Boxplots of NSC's distance from center of injection site
clear distInit;
distInit = zeros(n_seeds_current,Finaltimestep(1, Doses)/acc);
for n = 1:n_seeds_current
    distInit(n,:) = sqrt( sum ( (p(n).coord - inj_center'*ones(1,size(p(n).coord,2))).^2, 1 ) );
end

num_plots = 10;
data = distInit(:, (1:num_plots)*floor(Finaltimestep(1, Doses)/(num_plots*acc)))*CONVERT2MICRON;
labels = {num_plots};
for i=1:num_plots
    labels{i} = num2str( floor(Finaltimestep(1, Doses)/num_plots)*i /100 );
end
figure;     boxplot( data, labels );
hold on;    plot( (1:acc:Finaltimestep(1, Doses))*(num_plots/Finaltimestep(1, Doses)), median( distInit )*CONVERT2MICRON );
% hold on;    plot( (1:acc:Finaltimestep(1, Doses))*(num_plots/Finaltimestep(1, Doses)), mean( distInit )*CONVERT2MICRON );

ylim([0 max(max(data))+500]);
xlim([0 num_plots+1]);
xlabel( 'Time Intervals (days)' ); ylabel( 'Distance from injection site (\mu)' );

% Save with a unique name
filename = 'distance.png';
saveas(gcf, filename);

%% Determine whether NSC made it within a certain radius to cancer center
if has_cancer
    for k=1: num_cancers
        percAtCancerGraph = zeros(Finaltimestep(1, Doses)/acc,1);
        for i = 1:Finaltimestep(1, Doses)/acc
            numAtCancer = 0;
            
            for j = 1:n_seeds_current
                if all( abs(cancer_center(k,:) - p(j).coord(:,i)') < cancer_size )
                    numAtCancer = numAtCancer + 1;
                end
            end
            
            percAtCancerGraph(i,1) = numAtCancer;        
        end
        percAtCancerGraph = 100*percAtCancerGraph/n_seeds_current;
        
        % compute max perc at cancer and save to .txt
        max_percent = max(percAtCancerGraph);
        
        % compute average percent over sim
        average_percent = mean(percAtCancerGraph);
        
        % compute the percentage of unique NSCs ever arrivad
        ever_reached = false(n_seeds_current,1);
        n_full_timesteps = size(p_original(1).coord,2);
        
        for i=1:n_full_timesteps
            for j=1:n_seeds_current
                if all( ...
                        abs(cancer_center(k,:) - ...
                        p_original(j).coord(:,i)') < cancer_size)
                    ever_reached(j) = true;
                end
            end
        end
        
        % convert number of unique NSCs arrived into %
        percent_ever_arrived = 100*sum(ever_reached)/n_seeds_current;
        
        % save max percent
        filename = sprintf('max_percent_on_cancer_site%d.txt', k);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Maximum percent on cancer: %.2f%%\n', max_percent);
        fclose(fid);
        
        % save avg percent
        filename = sprintf('average_percent_on_cancer_site%d.txt', k);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Average percent on cancer: %.2f%%\n', average_percent);
        fclose(fid);
        
        % savae number of unique NSCs ever arrived
        filename = sprintf('percent_ever_arrived_site%d.txt', k);
        fid = fopen(filename, 'w');
        fprintf(fid, 'Percent of unique NSCs that ever arrived: %.2f%%\n', percent_ever_arrived);
        fclose(fid);

        indAtCancer = []; 
        for i = 1:n_seeds_current
            if all( abs(cancer_center - p(i).coord(:,end)') < cancer_size )
                indAtCancer = [indAtCancer,i]; 
            end
        end

        figure; plot(xvals,percAtCancerGraph,'r');
        xlabel('Time Intervals');   ylabel('Percent of NSC that Reach Cancer Site');

        curr_site_num=k;

        % Save with a unique name
        filename = 'oncancer.png';
        saveas(gcf, filename);
    end
end


% Plot for drift
x_boundary_value = 100;     % x-location of boundary line
final_idx = size(p(1).coord, 2);  % final time index
z_slice = 78;               % z-slice to project onto

figure;
plotenvironment_run_tumor(0.7, coh_map, ubound, lbound, seed_sd, inj_center, cancer_size, cancer_center, has_cancer, num_cancers);
hold on;

x_positions = zeros(1, n_seeds_current);
y_positions = zeros(1, n_seeds_current);
z_positions = zeros(1, n_seeds_current);

for i = 1:n_seeds_current
    x = p(i).coord(1, final_idx);
    y = p(i).coord(2, final_idx);
    z = p(i).coord(3, final_idx);

    x_positions(i) = x;
    y_positions(i) = y;
    z_positions(i) = z;

    % Plot only if near choice z-slice
    if abs(z - z_slice) < 1
        plot3(x, y, z_slice, 'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w');
    end
end

y_range = [1, size(coh_map,2)];  % full y-axis range
plot3([x_boundary_value x_boundary_value], y_range, [z_slice z_slice], 'r-', 'LineWidth', 2);

left_count = sum(x_positions < x_boundary_value);
right_count = sum(x_positions >= x_boundary_value);
disp(['Seeds left of boundary (x < ', num2str(x_boundary_value), '): ', num2str(left_count)]);
disp(['Seeds right of boundary (x ≥ ', num2str(x_boundary_value), '): ', num2str(right_count)]);

xlabel('x'); ylabel('y'); zlabel('z');
title(['Final Seed Positions (z = ', num2str(z_slice), ') with Boundary at x = ', num2str(x_boundary_value)]);
axis tight; grid on;

% Save with a unique name
filename = 'drift.png';
saveas(gcf, filename);

eps = 1e-6;
ratio = left_count / max(right_count, eps);  % Avoid divide by zero

filename = 'seed_distribution_ratio.mat';
save(filename, 'left_count', 'right_count', 'ratio');




% Plot for each day
n_days = 17;
step_indices = 2:2:34;  % These are the saved time steps in data (1-based indices)

figure('Units','normalized','Position',[0.05 0.05 0.9 0.85]);  % Large window
tiledlayout(4, 5, 'Padding', 'compact', 'TileSpacing', 'compact');  % 4x5 grid

for day_idx = 1:n_days
    step_idx = step_indices(day_idx);

    nexttile;

    plotenvironment_run_tumor(0.7, coh_map, ubound, lbound, seed_sd, inj_center, ...
                               cancer_size, cancer_center, has_cancer, num_cancers);
    hold on;

    for i = 1:n_seeds_current
        if size(p(i).coord, 2) >= step_idx
            plot3(p(i).coord(1, step_idx), ...
                  p(i).coord(2, step_idx), ...
                  p(i).coord(3, step_idx), ...
                  'ok', 'MarkerSize', 4, 'MarkerFaceColor', 'w');  % white circle with black center
        end
    end

    title(['Day ', num2str(day_idx), ' (step ', num2str(step_idx), ')']);
    axis tight;
    daspect([1 1 1]);  % equal axis scaling
end

sgtitle('Final Seed Positions at End of Each Day');

% Save with a unique name
filename = 'position17.png';
saveas(gcf, filename);
