% This code computes a scoring matrix on different parcellations of the
% same resting-state functional connectivity data from the NYU Test-Retest
% dataset. The goal is to determine whether the scoring matrices will be
% consistent across these different parcellations. Functional connectivity
% is computed via CONN.

% NOTE: There is one section ("Extract each ROI in the parcellations") that
% requires you to an execute a script in Terminal. More instructions in
% the comments of that section.

root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';
cd(fullfile(root_path, 'Data/FC/Parcellation Consistency'))
addpath(fullfile(root_path, 'Scripts/FC'))

%% Threshold each matrix at the top 5% of values

% We do this because the scoring algorithm updates the scoring matrix for
% all pairs of ROIs that are considered to be connected. If almost every
% pair of ROIs has a nonzero value of connectivity, as is the case for all
% of the correlation matrices, then all pairs of Yeo networks will receive
% the same score according to the algorithm. Therefore, we need to set a
% threshold for connectivity. When scoring the psychedelic dataset, this
% threshold is statistical significance (as defined by each study), whereas
% here, we have (arbitrarily) decided that the threshold is the top 5% of
% correlation values of each connectivity matrix.

files = dir('*matrix.csv'); % the FC matrices for all the parcellations
for file = files'
    FC_matrix = readmatrix(file.name);
    FC_matrix(isnan(FC_matrix)) = min(min(FC_matrix)); % set all the NaN values of the matrix to the minimum value
    topten_val = prctile(FC_matrix, 95);
    FC_matrix = FC_matrix .* (FC_matrix > topten_val);
    writematrix(FC_matrix, fullfile(append(erase(file.name,'.csv'),'_thr.csv'))); % e.g. schaefer_FC_matrix_thr.csv
end

%% Get indices of ROIs that have nonzero FC values, and the sign of their FC

% (This data will be inputted into the scoring matrix.)

thr_files = dir('*_thr.csv'); % thresholded FC matrices for all the parcellations
for thr_file = thr_files'
    thr_matrix = readmatrix(thr_file.name);
    inds = [];
    signs = [];
    % loop through all entries of the thresholded FC matrices
    for r = 1:size(thr_matrix,1)
        for c = 1:size(thr_matrix,2)
            if thr_matrix(r,c) ~= 0
                inds = [inds; [r c]];
                if thr_matrix(r,c) > 0
                    signs = [signs; 1];
                else
                    signs = [signs; -1];
                end
            end
        end
    end
    indsigns = [inds signs]; % combine all the indices and signs into a single nx2 matrix
    writematrix(indsigns, fullfile(append(erase(thr_file.name,'_matrix_thr.csv'),'_indsigns.csv'))); % e.g. schaefer_FC_indsigns.csv
end

%% Extract each ROI in the parcellations

cd(fullfile(root_path, 'Scripts/FC/Parcellation Consistency'))

% NOTE: To run this section, you will need to execute extract_ROI_maps.sh
% in Scripts/FC/extract_ROI_maps.sh in Terminal. For some reason, while it is
% posible to run Bash scripts from MATLAB, I get an error:
% "fslmaths: command not found" when I try to do this.

% Extract a mask of the ROI if it is defined as a spatial map rather than a
% set of coordinates. This is the case for the Brainnetome, Craddock,
% Harvard-Oxford, Power, Schaefer, and Shen parcellations.

% Make sure your version of bash is newer than 4.0 in order to run the code
% below. Note we use bash because we use FSLmaths to extract the spatial
% maps.


%% Score the networks

% Note: the code below corresponds to the seed-to-seed segment of the
% "Score the networks" section of FC_meta_analysis.m

cd(fullfile(root_path, 'Data/FC/Parcellation Consistency'))

% Struct to store the scoring matrices computed on each parcellation.
p.parcellation = {'brainnetome', 'craddock', 'harvard_oxford', 'power', 'raichle', 'schaefer', 'shen'};
p.partial_scores = {zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7)};
p.wta_scores = {zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7), zeros(7,7)};

for par = 1:length(p.parcellation)
    atlas_root = p.parcellation{par}; % name of parcellation, e.g. "schaefer"
    FC_matrix = readmatrix(fullfile([atlas_root '_FC_indsigns.csv']));
    for r = 1:size(FC_matrix,1)
        mask1_ind = FC_matrix(r,1) - 1; % subtract 1 because FSL is 0-indexed
        mask2_ind = FC_matrix(r,2) - 1;
        mask1_filename = fullfile([atlas_root '_mask_' num2str(mask1_ind) '.nii.gz']);
        mask2_filename = fullfile([atlas_root '_mask_' num2str(mask2_ind) '.nii.gz']);
        signs = FC_matrix(r,3);
        [s, overlaps_p_mask1, overlaps_w_mask1] = calculate_overlap('other', mask1_filename);
        [s, overlaps_p_mask2, overlaps_w_mask2] = calculate_overlap('other', mask2_filename);   
        [p.partial_scores{par}, p.wta_scores{par}] = scoring('other', signs, overlaps_p_mask1, overlaps_p_mask2, overlaps_w_mask1, overlaps_w_mask2, p.partial_scores{par}, p.wta_scores{par}, size(FC_matrix,1)); 
    end
end 

%% Compare score matrices for different parcellations

pvals = ones(length(p.parcellation));

for i = 1:length(p.partial_scores)
    for j = (i+1):length(p.partial_scores)
        
        % Define the "true distance" between score matrices as the sum of 
        % the absolute value of the difference between each corresponding
        % element
        true_dist = 0;
        for r = 1:7
            for c = 1:7
                true_dist = true_dist + abs(p.partial_scores{i}(r,c) - p.partial_scores{j}(r,c));
            end
        end
        
        % Permute the elements of one of the scoring matrices and then
        % calculate the distance between the score matrices. Repeat this
        % process 1000 times. These "permuted" distances will form a null
        % distribution.
        dists = [];
        for perm = 1:1000
            perm_inds = randperm(numel(p.partial_scores{j}));
            shuffled_j = reshape(p.partial_scores{j}(perm_inds), size(p.partial_scores{j}));
            dist = 0;
            for r = 1:7
                for c = 1:7
                    dist = dist + abs(p.partial_scores{i}(r,c) - shuffled_j(r,c));
                end
            end
            dists = [dists dist];
        end
        
        % Calculate the p-value of the true distance as its percentile
        % within the null distribution
        nlessp = sum(dists < true_dist);
        nequalp = sum(dists == true_dist);
        pp = (nlessp + 0.5*nequalp) / length(dists);
        pvals(i,j) = pp;
        pvals(j,i) = pp;
        
    end
end

pvals