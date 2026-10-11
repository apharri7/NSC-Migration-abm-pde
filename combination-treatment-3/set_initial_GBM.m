function [seed_ind,seed_sd,indup,inddown,inj_center] = set_initial_GBM(coh_map, modelType, Doses, n_seeds) 

% set_initial_gbm.m selects NSC injection locations and z-boundaries for brain

% Input:
%   1. coh_map - coherency image. Voxels where coh_map ~= 0 are treated as part of brain domain.
%   2. modelType - injection location. These simulations use only 'Intranasal'
%   3. Doses - number of NSC injection groups
%   4. n_seeds - cumulative number of NSC after each new dose

% Output: 
%   1. seed_ind - linear coh_map indices selected as initial NSC positions
%   2. seed_sd - radius of spherical injection region
%   3. indup - upper occupied z-index for every (x,y) column of coh_map
%   4. inddown - lower occupied z-index for every (x,y) column of coh_map
%   5. inj_center - center of injection region

    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    % Select injection center and spherical injection region
    switch modelType
        case 'Intranasal'
            inj_center =    [100.5, 31, 60];
        case 'Intracerebral'
            inj_center =    [133.25, 88.25, 63.25];
    end
    seed_sd = 16; 
    seed_ROI = sqrt((X-inj_center(1)).^2 + (Y-inj_center(2)).^2 + (Z-inj_center(3)).^2)<seed_sd;   
    seed_ind = find(seed_ROI); 
 
   % Set injection positions 
   p = sobolset(1,'Skip',1e3,'Leap',1e2); 
   ind = ceil( net(p, n_seeds(Doses)) *length(seed_ind) ); 
   seed_ind = seed_ind(ind); 

   % Create better upper and lower bounds (z-axis) for mouse brain
    for n = 1:size(coh_map,1)
        for m = 1:size(coh_map,2)
            ind = find( squeeze(coh_map(n,m,:)) ~=0 );
            if( isempty(ind) )
                indup(n,m) = NaN;
                inddown(n,m) = NaN;  
            else
                indup(n,m) = max(ind);
                inddown(n,m) = min(ind);            
            end
        end
    end

    if size(coh_map,1) == 200 
        indup(:,1:5) = NaN; inddown(:,1:5) = NaN;
        indup(:,261:end) = NaN; inddown(:,261:end) = NaN;
        [indx,indy] = find( indup > 120); 
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        while( ~isempty( indx ) ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( indup( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) < 128*4 ) 
                    jj = find( nbd~=128 ); njj = setdiff( 1:4, jj ); 
                    indup( indx(n), indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) );   
        end

        [indx,indy] = find( inddown == 1 ); 
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        nrepeat = 0; 
        while( ~isempty( indx ) && nrepeat < 10 ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( inddown( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) > 5 ) 
                    jj = intersect( find( nbd > 1 ), find( ~isnan(nbd) ) ); njj = setdiff( 1:4, jj ); 
                    inddown( indx(n) , indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) ); 
            nrepeat = nrepeat + 1; 
        end
        
        [indx,indy] = find( inddown <= 2 );
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        nrepeat = 0; 
        while( ~isempty( indx ) && nrepeat < 10 ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( inddown( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) > 10 ) 
                    jj = intersect( find( nbd > 2 ), find( ~isnan(nbd) ) ); njj = setdiff( 1:4, jj ); 
                    inddown( indx(n) , indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) ); 
            nrepeat = nrepeat + 1; 
        end
        
        [indx,indy] = find( inddown <= 3 );
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        nrepeat = 0; 
        while( ~isempty( indx ) && nrepeat < 10 ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( inddown( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) > 12 ) 
                    jj = intersect( find( nbd > 3 ), find( ~isnan(nbd) ) ); njj = setdiff( 1:4, jj ); 
                    inddown( indx(n) , indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) ); 
            nrepeat = nrepeat + 1; 
        end
        
        [indx,indy] = find( inddown <= 4 );
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        nrepeat = 0; 
        while( ~isempty( indx ) && nrepeat < 10 ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( inddown( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) > 16 ) 
                    jj = intersect( find( nbd > 4 ), find( ~isnan(nbd) ) ); njj = setdiff( 1:4, jj ); 
                    inddown( indx(n) , indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) ); 
            nrepeat = nrepeat + 1; 
        end        

        [indx,indy] = find( inddown(:,1:200) <= 8 );
        ii = [1,1; -1, 1; 1,-1; -1,-1]; 
        nrepeat = 0; 
        while( ~isempty( indx ) && nrepeat < 10 ) 
            nn = []; 
            for n = 1:length(indx) 
                nbd = diag( inddown( indx(n)+ii(:,1), indy(n)+ii(:,2) ) ); 
                if( sum(nbd) > 32 ) 
                    jj = intersect( find( nbd > 8 ), find( ~isnan(nbd) ) ); njj = setdiff( 1:4, jj ); 
                    inddown( indx(n) , indy(n) ) = mean( nbd(jj) ); 
                    nn = [nn, n]; 
                end
            end
            indx = indx( setdiff( 1:length(indx), nn ) ); 
            indy = indy( setdiff( 1:length(indy), nn ) ); 
            nrepeat = nrepeat + 1; 
        end              
    end
end
