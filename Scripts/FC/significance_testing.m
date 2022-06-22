% This code performs significance testing on the functional connectivity
% results (FC_meta_analysis.m). The goal is to compare the scoring matrices
% computed on the psychedelic dataset with scoring matrices computed on
% resting-state functional connectivity (null dataset), which is
% parcellated using the Schaefer100 atlas. As a first step to computing the
% null scoring matrices, we calculate the overlap between each Schaefer100
% ROI and the Yeo networks. We then determine the Schaefer100 ROIs that are
% functionally connected in the null dataset and refer to this as a
% "surrogate study." We iterate this process once for each study in the
% psychedelic dataset, resulting in a "surrogate study set" whose scoring 
% matrix we can compare with that of the psychedelic dataset. We compute 
% the null scoring matrix on the surrogate study set using the 
% pre-calculated overlaps with the Yeo networks. In total, we generate 500 
% surrogate study sets and their corresponding null scoring matrices. This 
% yields a histogram of 500 scores for each element of the 7x7 null scoring 
% matrix; note that this histogram serves as a null distribution. Based on
% this histogram, we then determine p-values for each element of the 
% scoring matrix for the psychedelic dataset.

root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';
cd(fullfile(root_path, 'Data', 'FC', 'Null'))

%% Calculate overlap between each Schaefer100 ROI with the 7 Yeo networks in advance

% First, define a 10mm-radius sphere around each of the Schaefer100 
% coordinates. Adapted from: 
% http://jpeelle.net/mri/misc/marsbar_roi.html

mkdir 'Spheres'
mkdir 'Spheres/Schaefer100'
outDir = 'Spheres/Schaefer100'; 
sphereRadius = 10;
coordinates_file = fullfile(root_path, 'Data', 'FC', 'Null', 'schaefer100_coordinates.xlsx');
coordinates = readmatrix(coordinates_file);
for row = 1:size(coordinates,1)
    x = coordinates(row,1);
    y = coordinates(row,2);
    z = coordinates(row,3);
    roiLabel = sprintf('%i_%i_%i', x, y, z);
    sphereROI = maroi_sphere(struct('centre', [x y z], 'radius', sphereRadius));
    outName = fullfile(root_path, outDir, sprintf('sphere_mask_%s', roiLabel));
    save_as_image(sphereROI, [outName '.nii']);
end

% Now calculate the overlap between each sphere and the Yeo networks.

OP = []; % partial overlaps (7x1 vector for each sphere)
OW = []; % winner-takes-all (WTA) overlaps (scalar for each sphere)
num_rois = 100;
cd(fullfile(root_path, 'Data', 'FC', 'Null', 'Spheres', 'Schaefer100'))
sphere_files = dir('*.xlsx');
for sphere_file = sphere_files'
    sphere = load(sphere_file.name);
    [signs_p, overlaps_yeo_p, overlaps_yeo_w]  = calculate_overlap('other', sphere);
    OP = [OP; overlaps_yeo_p];
    OW = [OW; overlaps_yeo_w];
end
        

%% Generate surrogate study set

% Load HCP data
load(fullfile(root_path, 'Data', 'FC', 'Null', 'Schaefer100_BOLD_HCP.mat'));
ts = BOLD_timeseries_HCP;

num_studies = 18; % number of studies in the psychedelic dataset
N = ceil(mean([15, 15, 19, 15, 10, 20, 20, 20, 12, 20, 9, 22, 16, 15, 20, 15, 20])); % average number of subjects in the psychedelic dataset
TR = ceil(mean([210, 220, 220, 300, 300, 300, 180, 220, 300, 150, 213, 156, 220, 300, 240])); % average TR of studies in the psychedelic dataset
N_HCP = 100; % number of subjects in the HCP dataset
T_HCP = 1189; % number of TRs in the HCP dataset

n_iter = 500; % number of surrogate study sets

all_Cp = zeros(7,7,n_iter); % tensor of surrogate study sets, each of which is a 7x7 partial scoring matrix
all_Cw = zeros(7,7,n_iter); % tensor of surrogate study sets, each of which is a 7x7 winner-takes-all (WTA) scoring matirx

% NOTE: You may have to run the for loop below several times and
% concatenate the outputs of all_Cp; for some reason, the for loop seems to
% get indefinitely stuck after 85 to 120 iterations.
for a = 1:n_iter
    
    surrogate_study_set = []; 

    C_p = zeros(7); % partial scoring matrix
    C_w = zeros(7); % WTA scoring matrix

    for i = 1:num_studies
        surrogate_study_set = compute_surrogate_study(T_HCP, TR, ts, num_rois, N_HCP, N, surrogate_study_set);
    end

    % Compute scoring matrix for the surrogate study
    for row = 1:size(surrogate_study_set,1)
        ind1 = surrogate_study_set(row,1); % ROI #1
        ind2 = surrogate_study_set(row,2); % ROI #2
        signs = surrogate_study_set(row,3); % sign of connectivity between ROIs
        overlaps_yeo1_p = OP(ind1,:); % partial overlap between ROI #1 and Yeo networks
        overlaps_yeo2_p = OP(ind2,:); % " ROI #2
        overlaps_yeo1_w = OW(ind1); % WTA overlap between ROI #1 and Yeo networks
        overlaps_yeo2_w = OW(ind2); % " ROI #2
        [C_p, C_w] = scoring('other', signs, overlaps_yeo1_p, overlaps_yeo2_p, overlaps_yeo1_w, overlaps_yeo2_w, C_p, C_w, zeros(num_studies,1));
    end
    
    clear surrogate_study_set
    
    all_Cp(:,:,a) = C_p;
    all_Cw(:,:,a) = C_w;
    
    clear C_p C_w

end

save('all_Cp_null.mat', 'all_Cp');
save('all_Cw_null.mat', 'all_Cw');

% figure;
% h1 = heatmap(pvalsp);
% caxis([0.05 1]);
% h1.XDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
% h1.YDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
% 
% figure;
% h2 = heatmap(pvalsw);
% caxis([0.05 1]);
% h2.XDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
% h2.YDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};