function plotenvironment_run_tumor(coherency_value,coh_map, ubound, lbound, seed_sd, inj_center, cancer_size, cancer_center, has_cancer, num_cancers)
% plotenvironment_run_tumor plots the brian, injection region, and tumor
% Inputs: 
% 1. coherency_value - threshold for distinguishing between white and grey matter
% 2. coh_map - coherency image
% 3. ubound, lbound - upper and lower z-coordinate for each (x,y) location
% 4. seed_sd - radius of displayed injection region
% 5. inj_center - center of the injection region
% 6. cancer_size, cancer_center - size has three components displaying half-axis width in x,y,z directions, center of each tumor
% 7. has_cancer - logical check to see if simulation includes tumor region
% 8. num_cancers - number of cancer sites considered in simulation

% Output: plot

    [Y,X,Z] = meshgrid(1:size(coh_map,2),1:size(coh_map,1),1:size(coh_map,3));
    
    hold on; surf( ubound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot upper bound of mouse brain
    hold on; surf( lbound', 'linestyle', 'none' , 'FaceColor', [0 0.4470 0.7410] , 'facealpha', 0.3 );      %Plot lower bound of mouse brain

    % Plot the injection region 
    
    [x,y,z] = sphere;
    x = x*seed_sd + inj_center(1);
    y = y*seed_sd + inj_center(2);
    z = z*seed_sd + inj_center(3);
    h = surf(x, y, z);                                                                                      %Plot injection site
    set(h,'FaceColor',[1 0 1],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');

    if has_cancer
        for i=1:num_cancers
            [x,y,z] = sphere;
            x = x*cancer_size(1) + cancer_center(i,1);
            y = y*cancer_size(2) + cancer_center(i,2);
            z = z*cancer_size(3) + cancer_center(i,3);
            h = surf(x, y, z);                                                                              %Plot cancer site
            set(h,'FaceColor',[1 .7 0],'FaceAlpha',0.4,'FaceLighting','gouraud','EdgeColor','none');
        end
    end
    
    ii = find( coh_map > coherency_value ); 
    hold on; plot3( X(ii(1:100:end)),Y(ii(1:100:end)),Z(ii(1:100:end)),'b.');                               %Plot WM track
    
    xlabel('x'); ylabel('y'); zlabel('z'); grid on; axis equal; daspect([1 1 1]); camlight;
end
