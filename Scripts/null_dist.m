% Function that computes a null distribution for assessing the statistical
% significance of the functional connectivity meta-analysis results. null_dist
% calls compute_surrogate_study; this computes correlations between
% randomly chosen windows of resting-state Human Connectome Project data
% parcellated with the Dessikan-Killiany (DK) atlas (either 62 ROIs for cortical
% connectivity or 80 ROIs for cortical and subcortical connectivity.)
% null_dist then runs compute_surrogate_study 22 times (22 is the number of
% studies in the psychedelic FC dataset), producing a "surrogate study
% set." null_dist calculates the FC score matrix (num_networks x num_networks
% where num_networks = 7 for only cortical FC and num_networks = 15 for
% both cortical and subcortical FC) for each surrogate study set. (Note that 
% it scores FC without taking the sign of connectivity into account, as the 
% sign in the scoring algorithm refers to whether or not connectivity 
% increases or decreases under psychedelics.) null_dist repeats this process 
% num_iter times, resulting in a num_networks  x num_networks x num_iter 
% tensor. This matrix forms a null distribution for each entry of the score 
% matrix. In FC_meta_analysis.m, if an entry of the score matrix is above 
% the 95th percentile or below the 5th percentile  of the score matrix, then 
% it is considered significant. 

% Inputs: 
% - num_iter: The number of surrogate study sets that null_dist computes.
%   Typically this is set to 1000. 
% - num_rois: The number of regions of interest in the DK atlas whose FC
%   is calculated by compute_surrogate_study.m. num_rois must be either 62
%   (for only cortical connectivity) or 80 (for cortical and subcortical
%   connectivity.)

% Outputs:
% - all_Cp: A num_networks x num_networks x num_iter tensor, in which each
%   num_networks x num_networks matrix is an FC score matrix computed with
%   the partial paradigm.
% - all_Cw: The same, but for winner-takes-all (WTA) paradigm.

% Kenneth Shinozuka 2022

function [all_Cp, all_Cw] = null_dist(num_iter,num_rois)

    root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/';
    cd(root_path)
    addpath(fullfile(root_path, 'Scripts', 'natsortfiles')) % Scripts for numerically sorting files in a directory

    % load HCP data
    load /Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/FC/Null/DK80_BOLD_HCP_100subs.mat

    % calculate the overlaps between the ROIs in the DK parcellation and
    % the Yeo networks
    OP = []; % partial overlaps (1x7 vector for each of the 80 ROIs)
    OW = []; % winner-takes-all (WTA) overlaps (scalar for each ROI)
    if num_rois == 62
        cd(fullfile(root_path, 'Atlases', 'dk62'))
        atlas = 'yeo';
        num_networks = 7;
    elseif num_rois == 80
        cd(fullfile(root_path, 'Atlases', 'dk80'))  
        atlas = 'yeo+tian';
        num_networks = 15;
    else
        error('num_rois must be either 62 (cortical) or 80 (subcortical and cortical).');
    end
    roi_files = natsortfiles(dir('*.nii.gz'));
    for roi_file = roi_files'
        [signs_p, overlaps_p, overlaps_w]  = calculate_overlap('other', roi_file.name, atlas);
        OP = [OP; overlaps_p];
        OW = [OW; overlaps_w];
    end
    
    if num_networks ~= 7 | num_networks ~= 15
        error('num_networks must be either 7 (cortical) or 15 (subcortical and cortical).');
    end
    
    num_studies = 22; % number of studies in the psychedelic dataset
    N = ceil(mean([15, 15, 19, 15, 10, 20, 20, 20, 12, 20, 9, 22, 22, 16, 15, 20, 15, 20, 18, 24, 15])); % average number of subjects in the psychedelic dataset
    TR = ceil(mean([210, 220, 220, 300, 300, 300, 180, 220, 300, 150, 213, 156, 156, 220, 300, 240, 240, 240])); % average number of TRs in the psychedelic dataset (sometimes the number of TRs was not reported)
    N_HCP = 100; % number of subjects in the HCP dataset
    T_HCP = 1189; % number of TRs in the HCP dataset

    all_Cp = zeros(num_networks,num_networks,num_iter); % tensor of surrogate study sets, each element of which is a num_networks x num_networks partial scoring matrix
    all_Cw = zeros(num_networks,num_networks,num_iter); % tensor of surrogate study sets, each element of which is a num_networks x num_networks winner-takes-all (WTA) scoring matirx

    for a = 1:num_iter

        C_p = zeros(num_networks); % partial scoring matrix
        C_w = zeros(num_networks); % WTA scoring matrix
        C_p_nosign = zeros(num_networks); % scored without taking the sign of connectivity into account
        C_w_nosign = zeros(num_networks); 

        for i = 1:num_studies
            surrogate_study = compute_surrogate_study(T_HCP, TR, ts, num_rois, N_HCP, N); % matrix consisting of indices of ROIs that are significantly connected in resting-state
            % Compute scoring matrix for the surrogate study set
            for row = 1:size(surrogate_study,1)
                ind1 = surrogate_study(row,1); % ROI #1
                ind2 = surrogate_study(row,2); % ROI #2
                overlaps_1_p = OP(ind1,:); % partial overlap between ROI #1 and Yeo+Tian networks
                overlaps_2_p = OP(ind2,:); % " ROI #2
                overlaps_1_w = OW(ind1); % WTA overlap between ROI #1 and Yeo+Tian networks
                overlaps_2_w = OW(ind2); % " ROI #2
                [~, ~, C_p_nosign, C_w_nosign] = scoring('other', 0, overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(surrogate_study,1), 0);
            end
            clear surrogate_study
        end

        all_Cp(:,:,a) = C_p_nosign;
        all_Cw(:,:,a) = C_w_nosign;

        clear C_p C_w C_p_nosign C_w_nosign

        a

    end

%     save(fullfile(root_path, 'FC', 'Null', 'all_Cp_null_DK62_061222_sign.mat'), 'all_Cp_null')
%     save(fullfile(root_path, 'FC', 'Null', 'all_Cw_null_DK62_061222_sign.mat'), 'all_Cw_null')
