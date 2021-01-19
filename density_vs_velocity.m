% Created by: Alexandros A. Fragkopoulos, Johannes Frey, Flora-Maud le Menn
% Max Planck institute for dynamics and self-organization
%
% Last modified on: 12/01/2021
%
% It uses the results from the processing_raw.m to calculate the local
% densities and velocities using voronoi tessellation.
% 
% The inputs to the program are:
% pathname: a string containing the directory of the processed results
%
% filename: the name of the file of the proccessed results
%
% saut: the number of frames over which to calculate the average v and rho.
% Typically this is 33, which corresponds to 1 second.
%
% n_bins: this is the number of radial bins. The radius is binned such that 
% the area is perserved, meaning that (r*dr) is constant. This is important
% for small radii, which in the case of constant dr would correspond to 
% small areas, thus particles and therefore low statistics. The constant 
% (r*dr) allows for all bins to have high statistics. Typically the value 
% of n_bins is 50.

function density_vs_velocity(pathname,filename,saut,n_bins)
%% Parameters
fps = 33;
conversion_64=500/751;
conversion_10=500/1187;
starting_image = 1;
%% Load Data
load(strcat(pathname,filename),'all_info','n_images');
%%
min_r = min(all_info(:,6)); % find the minimum radius in the data
max_r = max(all_info(:,6)); % find the maximum radius in the data
u0 = (max_r^2-min_r^2)./(n_bins); % calculate the associated binning width
r_bin = sqrt((0:1:n_bins).*u0+min_r^2)';% An array with the upper side of all the bins
%%
mean_v = zeros([n_bins,floor((n_images-starting_image)/saut)]);
mean_std_v = zeros([n_bins,floor((n_images-starting_image)/saut)]);
mean_rho = zeros([n_bins,floor((n_images-starting_image)/saut)]);
mean_std_rho = zeros([n_bins,floor((n_images-starting_image)/saut)]);
%%
for j = 1:floor((n_images-starting_image)/saut) %look at the frame j, with a saut between each frame we look at.
    disp(strcat('Velocity vs Density: ', num2str(100*j/floor((n_images-starting_image)/saut)),'% Done'))
    for i=1:n_bins
        % Find all cells in the i-th bin
        ind = find(all_info(:,6)>=(r_bin(i)) & all_info(:,6)<(r_bin(i+1)) & all_info(:,10)>=(j-1)*saut+starting_image & all_info(:,10)<j*saut+starting_image);
        % the all_info(:,10)==j tells to look only in the frame j, if the
        % step is different than 1, we need to look in all the frames
        % between two "saut"
        if ~isempty(ind)
        % Calculate values associated to the i-th bin    
            mean_v(i,j) = sqrt(mean(all_info(ind,5).^2)); % average velocity
            mean_std_v(i,j) = std(all_info(ind,5).^2)/(2.*mean_v(i,j).*sqrt(length(ind))); %mean standart error = division by the total number of particles taken into account
            mean_rho(i,j) = mean(all_info(ind,8)); % average density
            mean_std_rho(i,j) = std(all_info(ind,8))/sqrt(length(ind)); %mean standart deviation of density
        end
    end
end
%% Save everything
save(strcat(pathname,'density_vs_velocity.mat'),'mean_v','mean_std_v','mean_rho','mean_std_rho','n_bins','r_bin','saut','-v7.3')
end