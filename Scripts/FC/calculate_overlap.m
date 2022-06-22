% Function that outputs the extent of spatial overlap between an
% ROI in a study and a Yeo network, as well as the average sign of the
% correlation between the region of overlap and the seed in seed-to-voxel
% analyses.

% Inputs:
% - datatype: A study's method for measuring functional connectivity. This
%   is either 'other' or 'seed-to-voxel.'
% - filename: The filename containing the spatial map of the ROI. The ROI
%   is either a 10mm-radius sphere around a set of coordinates or a mask of
%   a seed.

% Outputs:
% - signs: The sign of the connectivity between the seed and the set of 
%   voxels that overlap with each Yeo network in seed-to-voxel analyses.
%   This is an nx1 vector if the datatype is 'seed-to-voxel', where n is
%   the number of Yeo networks that contain at least one voxel in the 
%   spatial map of seed-to-voxel connectivity, or it is an empty array if 
%   the datatype is 'other'.
% - overlap_yeo_p: The overlaps of the ROI with the Yeo networks according
%   to the partial paradigm. This is a 1x7 array in which each element is
%   the proportion of voxels in the ROI that overlap with the corresponding
%   Yeo network.
% - overlap_yeo_w: The overlaps of the ROI with the Yeo networks according
%   to the winner-takes-all (WTA) paradigm. This is a scalar representing
%   the single Yeo network that contains the most voxels in the ROI.

function [signs, overlap_yeo_p, overlap_yeo_w] = calculate_overlap(datatype, filename)

    root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';
    
    % load the Yeo networks in 2mm resolution
    yeo = spm_vol(fullfile(root_path, 'Yeo2011_7Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz'));
    yeo_img = spm_read_vols(yeo);
    
    % yeo_vals contains the Yeo network that each voxel in the ROI belongs to 
    % yeo_vals is a 1 x nvoxels vector, except for seed-to-voxel analyses;
    % in this case, yeo_vals is a 2 x nvoxels vector. Here, the second row
    % contains the t-value at each voxel. We then average these values
    % across each region that overlaps with each Yeo network and take the
    % sign of each average. These signs are stored in signs_w
    % and inputted into the scoring algorithms. For studies that do not 
    % measure seed-to-voxel connectivity, we use the sign that was reported 
    % in the study, so signs_p and signs_w are empty arrays in this case. 
    yeo_vals = [];
    signs = [];
   
    % load the ROI
    nii = spm_vol(filename);
    img = spm_read_vols(nii);
    
    % locations of voxels inside the ROI
    [d1, d2, d3] = ind2sub(size(img), find(img~=0));
    
    for i = 1:size(d1,1)
        yeo_val = yeo_img(d1(i), d2(i), d3(i)); % the Yeo network corresponding to each voxel inside the ROI
        if yeo_val > 0
            if strcmp(datatype, 'seed_to_voxel')
                voxel_val = img(d1(i), d2(i), d3(i)); % t-value of the voxel at hand
                yeo_vals = [yeo_vals [yeo_val; voxel_val]];
            else
                yeo_vals = [yeo_vals yeo_val];
            end
        end
    end
    
    % In each region that overlaps with a Yeo network, compute the sign of
    % the average of the t-values of the voxels in that region, but only
    % for seed-to-voxel analyses, as discussed above.
    if strcmp(datatype, 'seed_to_voxel')
        unique_yeos = unique(yeo_vals(1,:)); % all the Yeo networks that the spatial map overlaps with
        for y = 1:size(unique_yeos,2)
            u = unique_yeos(y);
            t = [];
            for i = 1:size(yeo_vals(1,:),2) % the t-value of each voxel
                if yeo_vals(1,i) == u 
                    t = [t yeo_vals(2,i)]; % t-values of all the voxels that are contained within a particular Yeo newtork
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

    if isempty(yeo_vals)
        overlap_yeo_p = [0 0 0 0 0 0 0];
        overlap_yeo_w = 0;
    else
        counts_yeo = hist(yeo_vals(1,:), [1:7]); % histogram of the number of voxels overlapping with each Yeo network
        overlap_yeo_p = counts_yeo/sum(counts_yeo); % partial paradigm computes the proportion of voxels overlapping with each Yeo network
        overlap_yeo_w = mode(yeo_vals(1,:)); % WTA paradigm computes the single Yeo network that overlaps with the most voxels
    end
        
end