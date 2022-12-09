% Function that outputs the extent of spatial overlap between an
% ROI in a study and a Yeo network, as well as the average sign of the
% correlation between the region of overlap and the seed in seed-to-voxel
% analyses.

% Inputs:
% - datatype: A study's method for measuring functional connectivity. This
%   is a string: 'other', 'seed_to_voxel', or 'seed_to_voxel_CIFTI.'
%   The last of these three options refers to studies that produced spatial
%   maps of seed-to-voxel connectivity in CIFTI format.
% - filename: The filename containing the spatial map of the ROI. The ROI
%   is either a 10mm-radius sphere around a set of coordinates or a mask of
%   a seed, unless the datatype is 'seed_to_voxel_CIFTI'; in this case, the
%   file is a matrix of x,y,z coordinates (in MNI space) of voxels that 
%   have strong connectivity with the seed of interest (columns 1-3), as 
%   well as the sign of the connectivity of each voxel.
% - atlas: The atlas of the networks that the ROIs are desired to overlap
%   with. This is a string: 'yeo', 'yeo+tian', 'yeo+cerebellum', or
%   'yeo+tian+cerebellum.'

% Outputs:
% - signs: The sign of the connectivity between the seed and the set of 
%   voxels that overlap with each  network in seed-to-voxel analyses.
%   This is an n x 1 vector if the datatype is 'seed_to_voxel', where n is
%   the number of networks that contain at least one voxel in the 
%   spatial map of seed-to-voxel connectivity, or it is an empty array if 
%   the datatype is 'other'.
% - overlap_yeo_p: The overlaps of the ROI with the Yeo networks according
%   to the partial paradigm. This is a 1 x n_networks array in which each element is
%   the proportion of voxels in the ROI that overlap with the corresponding
%   Yeo network.
% - overlap_yeo_w: The overlaps of the ROI with the Yeo networks according
%   to the winner-takes-all (WTA) paradigm. This is a scalar representing
%   the single Yeo network that contains the most voxels in the ROI.

% Kenneth Shinozuka 2022

function [signs, overlaps_p, overlaps_w] = calculate_overlap(datatype, filename, atlas)

    root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/';
    
    % Load the combined Yeo and Tian networks in 2mm resolution. Values 1-7
    % correspond to the 7 cortical Yeo networks, and values 8-23 correspond 
    % to the 16 subcortical Tian structures.
    if strcmp(atlas, 'yeo+tian')
        atlas_vol = spm_vol(fullfile(root_path, 'Atlases', 'yeo_tian_atlas.nii'));
        num_networks = 15;
    elseif strcmp(atlas, 'yeo+tian+cerebellum')
        atlas_vol = spm_vol(fullfile(root_path, 'Atlases', 'yeo_tian_cerebellum_atlas.nii'));
        num_networks = 16;
    elseif strcmp(atlas, 'yeo+cerebellum')
        atlas_vol = spm_vol(fullfile(root_path, 'Atlases', 'yeo_cerebellum_atlas.nii'));
        num_networks = 8;
    elseif strcmp(atlas, 'yeo')
        atlas_vol = spm_vol(fullfile(root_path, 'Atlases', 'yeo', 'MNI', 'Yeo2011_7Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz'));
        num_networks = 7;
    else
        error('Atlas must be one of the following: yeo+tian, yeo+tian+cerebellum, yeo+cerebellum, yeo');
    end
    atlas_img = spm_read_vols(atlas_vol);
    
    % yeo_vals contains the Yeo network that each voxel in the ROI belongs 
    % to. yeo_vals is a 1 x nvoxels vector, 
    % except for seed-to-voxel analyses; in this case, yeo_vals is a 2 x 
    % nvoxels vector. Here, the second row contains the t-value at each 
    % voxel. We then average these values across each region that overlaps 
    % with each Yeo network and take the sign of each average. These signs 
    % are stored in signs_w and inputted into the scoring algorithms. For 
    % studies that do not measure seed-to-voxel connectivity, we use the 
    % sign that was reported  in the study, so signs_p and signs_w are 
    % empty arrays in this case. 
    atlas_vals = [];
    signs = [];
   
    if strcmp(datatype, 'other') | strcmp(datatype, 'seed_to_voxel')      
        % load the ROI
        nii = spm_vol(filename);
        img = spm_read_vols(nii);
        [d1, d2, d3] = ind2sub(size(img), find(img~=0)); % locations of voxels inside the ROI       
    elseif strcmp(datatype, 'seed_to_voxel_CIFTI')
        coords = readmatrix(filename);
        vox = [coords(:,1) coords(:,2) coords(:,3) ones(size(coords,1),1)]*(inv(atlas_vol.mat))'; % convert MNI coordinates to voxels in 2mm resolution
        vox = round(vox(:,1:3));
        d1 = vox(:,1);
        d2 = vox(:,2);
        d3 = vox(:,3);
    else
        error('Datatype must be other, seed_to_voxel, or seed_to_voxel_CIFTI.');
    end
    
    for i = 1:size(d1,1)
        atlas_val = atlas_img(d1(i), d2(i), d3(i)); % the Yeo/Tian network corresponding to each voxel inside the ROI
        if atlas_val > 0
            % The index of a left-hemispheric portion of a subcortical 
            % structure is always 8 greater than the index of the
            % right-hemispheric portion. Subtracting 8 ensures that we treat
            % overlaps with the right-hemispheric and left-hemispheric
            % portions as equivalent overlaps.
            if strcmp(atlas, 'yeo+tian') | strcmp(atlas, 'yeo+tian+cerebellum')
                if atlas_val > num_networks
                    atlas_val = atlas_val - 8;
                end
            end
            % In the Yeo+cerebellum atlas and Yeo+cerebellum+Tian atlas, 
            % some cerebellar regions overlap with other networks. 
            % Here, we set all such overlapping regions to be part of 
            % the cerebellum (the cerebellum has the highest index).
            if strcmp(atlas, 'yeo+cerebellum') | strcmp(atlas, 'yeo+tian+cerebellum')
                cerebellum_index = num_networks;
                if atlas_val > cerebellum_index
                    atlas_val = cerebellum_index; 
                end
            end
            if strcmp(datatype, 'seed_to_voxel')
                s = sign(img(d1(i), d2(i), d3(i))); % sign of the t-value of the voxel at hand
                atlas_vals = [atlas_vals [atlas_val; s]];
            elseif strcmp(datatype, 'seed_to_voxel_CIFTI')
                s = coords(i,4); % sign of the t-value of the voxel at hand
                atlas_vals = [atlas_vals [atlas_val; s]];
            else
                atlas_vals = [atlas_vals atlas_val];
            end
        end
    end
    
    
    % In each region that overlaps with a Yeo network, compute the sign of
    % the average of the t-values of the voxels in that region, but only
    % for seed-to-voxel analyses, as discussed above.
    if strcmp(datatype, 'seed_to_voxel') | strcmp(datatype, 'seed_to_voxel_CIFTI')
        unique_atlas_vals = unique(atlas_vals(1,:)); % all the Yeo networks that the spatial map overlaps with
        for r = 1:size(unique_atlas_vals,2)
            u = unique_atlas_vals(r);
            t = [];
            for i = 1:size(atlas_vals(1,:),2) % sign of the t-value of each voxel
                if atlas_vals(1,i) == u 
                    t = [t atlas_vals(2,i)]; % signs of the t-values of all the voxels that are contained within a particular Yeo newtork
                end
            end
            % Assign a negative sign to a Yeo network if the mean of
            % t-values is less than zero, but a positive sign otherwise.
            if mean(t) < 0
                append = [u; -1];
            else
                append = [u; 1];
            end
            signs = [signs append];
        end
    end

    if isempty(atlas_vals)
        overlaps_p = zeros(1,num_networks);
        overlaps_w = 0;
    else
        counts_network = hist(atlas_vals(1,:), [1:num_networks]); % histogram of the number of voxels overlapping with each network
        overlaps_p = counts_network/sum(counts_network); % partial paradigm computes the proportion of voxels overlapping with each Yeo network
        overlaps_w = mode(atlas_vals(1,:)); % WTA paradigm computes the single Yeo network that overlaps with the most voxels
    end
        
end
