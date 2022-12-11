%% Introduction

% Broadly speaking, the studies define ROIs in three ways: (1) around
% individual sets of coordinates; (2) as networks of coordinates; (3) as 
% spatial maps of seeds; and (4) as whole-brain maps of voxels that are 
% functionally connected to seeds, i.e. in seed-to-voxel analyses. 
% (Note that the seeds could be defined anatomically, or in a data-driven 
% way such as ICA.) 

% In cases (1) and (2), we defined 10mm-radius spheres around the
% coordinates of interest, as performed in the first two sections of code. 
% Additionally, for each of these studies, we created a file ending in 
% '_coordinates.xlsx'that contains the (x,y,z) MNI coordinates of each ROI. 
% In cases (3) and (4),  we obtained spatial maps of the seeds/voxels and 
% saved them as files ending in '_seed.nii'.

% Note that some studies defined ROIs in multiple different ways. For
% instance, one study defined a portion of its ROIs as seeds and the
% remaining ROIs as networks of coordinates. Therefore, in total, there were
% six types of connectivity in the studies: coordinate-to-coordinate,
% network-to-network, seed-to-seed, seed-to-voxel, seed-to-network, and
% seed-to-coordinate.

% We created a file 'FC_indices.xlsx' that reports which ROIs are
% functionally connected in each study. The first row of the file
% contains the name of the study and the second row displays the type of
% connectivity in that study (there are six possibilities, as discussed
% above). The remaining rows either show the indices of the functionally
% connected coordinates/networks when the ROIs fall into categories (1) or (2) 
% (where the index represents the row number of the coordinate in 
% '*_coordinates.xlsx'), or they show the names of the functionally connected 
% seeds when the ROIs fall into categories (3) or (4). For each study, the
% third column also shows the sign of the change in connectivity between
% the corresponding ROIs; that is, if an entry in the column is 1, then
% the change in connectivity is positive, but if the entry is -1, then the
% change is negative.

% We then score the ROIs in two parts. Firstly, we calculate the
% overlap between the ROI and the Yeo networks (calculate_overlap.m or
% network_overlaps in the case of networks). Then, we update the scoring 
% matrix, representing the overall score of each pair of Yeo networks. 

% Kenneth Shinozuka 2022

% NOTE: You must update this to reflect the root path to the data on your
% device.
root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/';
cd(root_path)
addpath(fullfile(root_path, 'Scripts'))

%% Indicate whether a study defines ROIs based on coordinates or masks, and the number of secondary analyses associated with a dataset

studies.name = {'araujo_2012', 'barrett_2020_claustrum', 'bershad_2020', 'carhart-harris_2016_s2v', 'carhart-harris_2016_ICA', 'gaddis_2022', 'grimm_2018', 'kaelen_2016', 'lebedev_2015', 'luppi_2021_segregated', 'luppi_2021_time-averaged', 'madsen_2021', 'mason_2020', 'mason_2021', 'mcculloch_2021', 'muller_2017_s2s', 'muller_2017_s2v', 'muller_2018_ICA', 'muller_2018_s2v', 'palhano-fontes_2015_s2v', 'pasquini_2020_s2v', 'pasquini_2020_s2s', 'preller_2018', 'roseman_2014', 'sampedro_2017', 'smigielski_2019_ICA1', 'smigielski_2019_ICA2', 'tagliazucchi_2014'};
studies.datatype = {'seeds', 'seeds', 'coordinates', 'seeds', 'seeds', 'seeds', 'coordinates', 'seeds', 'coordinates', 'coordinates', 'coordinates', 'coordinates', 'seeds', 'seeds', 'coordinates', 'coordinates', 'seeds', 'seeds', 'seeds', 'coordinates', 'seeds', 'seeds', 'neither', 'seeds', 'coordinates', 'coordinates', 'seeds', 'seeds'};
studies.secondaries = {0, 0, 0, 4, 4, 0, 0, 4, 2, 4, 4, 0, 1, 1, 0, 3, 3, 3, 3, 0, 1, 1, 0, 2, 0, 1, 1, 2};    

%% Define a sphere around the coordinates

% Adapted from: http://jpeelle.net/mri/misc/marsbar_roi.html

addpath(fullfile(root_path, 'Toolboxes', 'marsbar'))

outDir = fullfile(root_path, 'Spheres');
sphereRadius = 10;

for i = 1:length(studies.datatype) % each study
    if strcmp(studies.datatype{i}, 'coordinates')
        % get the file with the relevant coordinates for that study
        study_coordinates_file = fullfile(root_path, 'FC', [studies.name{i} '_' studies.datatype{i} '.xlsx']);
        study_coordinates = readmatrix(study_coordinates_file);
        for j = 1:size(study_coordinates,1)
            row = study_coordinates(j,:);
            if ~isnan(row(1))
                x = row(1);
                y = row(2);
                z = row(3);
            end
            roiLabel = sprintf('%d_%d_%d', x, y, z); 
            sphereROI = maroi_sphere(struct('centre', [x y z], 'radius', sphereRadius)); % construct the sphere
            outName = fullfile(outDir, sprintf('sphere_mask_%s', roiLabel));
            save_as_image(sphereROI, [outName '.nii']);
        end
    end
end

%% Score the networks

% which atlas to use to score the networks
atlas = 'yeo';

if strcmp(atlas, 'yeo+tian+cerebellum')
    num_networks = 16;
    network_labels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FPN', 'DMN', 'Hipp', 'Amyg', 'pThal', 'aThal', 'NuAc', 'GP', 'Put', 'Caud', 'Cerebellum'};
    atlas_file = spm_vol(fullfile(root_path, 'Atlases', 'yeo_tian_cerebellum_atlas.nii'));
    atlas_img = spm_read_vols(atlas_file);
elseif strcmp(atlas, 'yeo+tian')
    num_networks = 15;
    network_labels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FPN', 'DMN', 'Hipp', 'Amyg', 'pThal', 'aThal', 'NuAc', 'GP', 'Put', 'Caud'};
    atlas_file = spm_vol(fullfile(root_path, 'Atlases', 'yeo_tian_atlas.nii'));
    atlas_img = spm_read_vols(atlas_file);
elseif strcmp(atlas, 'yeo+cerebellum')
    num_networks = 8;
    network_labels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FPN', 'DMN', 'Cerebellum'};
    atlas_file = spm_vol(fullfile(root_path, 'Atlases', 'yeo_cerebellum_atlas.nii'));
    atlas_img = spm_read_vols(atlas_file);
else
    num_networks = 7;
    network_labels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FPN', 'DMN'};
    atlas_file = spm_vol(fullfile(root_path, 'Atlases', 'yeo_atlas.nii'));
    atlas_img = spm_read_vols(atlas_file);
end

% matrix of scores, where C_ij is the connectivity between Yeo network i
% and Yeo network j. C_p is the scoring matrix based on the partial
% paradigm; C_w is the scoring matrix based on the winners-take-all
% paradigm.
C_p = zeros(num_networks);
C_w = zeros(num_networks);
C_p_nosign = zeros(num_networks);
C_w_nosign = zeros(num_networks);

% load file indicating which pairs of seeds/coordinates are connected to each other
[num_pairs, txt_pairs, raw_pairs] = xlsread(fullfile(root_path, 'FC', 'FC_indices.xlsx'));

for j = 1:3:size(raw_pairs, 2) % cycle through the studies
    study = char(raw_pairs(1,j)); % name of the study 
    num_secondaries = studies.secondaries{find(contains(studies.name,study))};
    if strcmp(raw_pairs(2,j), 'network_to_network')
        col = cell2mat(raw_pairs(3:end,j)); % indices of networks
        col = col(~isnan(col)); 
        coordinate_file = fullfile(root_path, 'FC', [study '_' 'coordinates.xlsx']); % coordinates of each region in each network
        coordinates = readmatrix(coordinate_file);
        [networks_p, networks_w] = network_overlaps(coordinates, atlas); % note that there is a separate function for calcluating overlaps in networks, where the overlaps get averaged across the individual ROIs in the network
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC 
        signs = signs(~isnan(signs)); 
        % for each pair of functionally connected networks, get their
        % overlaps with the Yeo networks and then update the scoring matrix
        % accordingly
        for r = 1:size(col,1)
            overlaps_1_p = networks_p(cell2mat(raw_pairs(r+2,j)),:);  
            overlaps_2_p = networks_p(cell2mat(raw_pairs(r+2,j+1)),:);
            overlaps_1_w = networks_w(cell2mat(raw_pairs(r+2,j)),:);
            overlaps_2_w = networks_w(cell2mat(raw_pairs(r+2,j+1)),:);
            [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(r), overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(col,1), num_secondaries);
        end
    elseif strcmp(raw_pairs(2,j), 'coord_to_coord')
        coordinate_file = fullfile(root_path, 'FC', [study '_' 'coordinates.xlsx']);
        coordinates = readmatrix(coordinate_file);
        col = cell2mat(raw_pairs(3:end,j)); % indices of coordinates
        col = col(~isnan(col));  
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
        signs = signs(~isnan(signs));
        for r = 1:size(col,1)
            ind1 = cell2mat(raw_pairs(r+2,j)); % index of coordinate 1
            ind2 = cell2mat(raw_pairs(r+2,j+1)); % index of coordinate 2
            % x, y, z values of coordinate 1
            x1 = coordinates(ind1,1);
            y1 = coordinates(ind1,2);
            z1 = coordinates(ind1,3);
            % if x, y, and z are integers, then the filename of the sphere
            % around the coordinate will be expressed in integers; otherwise, 
            % the filename will be expressed in floating point decimals
            if floor(x1)==x1 & floor(y1)==y1 & floor(z1)==z1
                sphere_file1 = fullfile(root_path, 'Spheres', ['sphere_mask_' sprintf('%d', x1) '_' sprintf('%d', y1) '_' sprintf('%d', z1) '.nii']);
            else
                sphere_file1 = fullfile(root_path, 'Spheres', ['sphere_mask_' sprintf('%.4f', x1) '_' sprintf('%.4f', y1) '_' sprintf('%.4f', z1) '.nii']);
            end
            % x, y, z values of coordinate 2
            x2 = coordinates(ind2,1);
            y2 = coordinates(ind2,2);
            z2 = coordinates(ind2,3);
            if floor(x2)==x2 & floor(y2)==y2 & floor(z2)==z2
                sphere_file2 = fullfile(root_path, 'Spheres', ['sphere_mask_' sprintf('%d', x2) '_' sprintf('%d', y2) '_' sprintf('%d', z2) '.nii']);                
            else
                sphere_file2 = fullfile(root_path, 'Spheres', ['sphere_mask_' sprintf('%.4f', x2) '_' sprintf('%.4f', y2) '_' sprintf('%.4f', z2) '.nii']);
            end
            [signs_p, overlaps_1_p, overlaps_1_w]  = calculate_overlap('other', sphere_file1, atlas);
            [signs_p, overlaps_2_p, overlaps_2_w] = calculate_overlap('other', sphere_file2, atlas);
            [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(r), overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(col,1), num_secondaries);
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_voxel') | strcmp(raw_pairs(2,j), 'seed_to_voxel_CIFTI')
        seeds = txt_pairs(3:end,j); % name of the seed
        seeds = seeds(~cellfun('isempty',seeds));        
        for s = 1:size(seeds,1)
            % I was unable to obtain the V1 seed used by Carhart-Harris et
            % al. (2016), so instead I assumed that the seed overlaps
            % completely with Yeo network 1, which is the visual network.
            if strcmp(study, 'carhart-harris_2016_s2v') & strcmp(char(seeds(s)), 'V1')
                overlaps_1_p = zeros(1,num_networks);
                overlaps_1_p(1) = 1;
                overlaps_1_w = 1;
            else
                seed1_filename = fullfile(root_path, 'FC', [study '_' char(seeds(s)) '_seed.nii']);
                [signs_p, overlaps_1_p, overlaps_1_w] = calculate_overlap('other', seed1_filename, atlas);
            end
            % Note that the spatial map of the functionally connected
            % voxels indicates the sign of the connectivity at each voxel,
            % hence FC_indices.xlsx does not specify the sign.
            if strcmp(raw_pairs(2,j), 'seed_to_voxel')
                voxels_file = fullfile(root_path, 'FC', [study '_' char(seeds(s)) '.nii']); % map of voxels that are functionally connected to the seed
                [signs2_p, overlaps_2_p, overlaps_2_w] = calculate_overlap('seed_to_voxel', voxels_file, atlas);
                [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('seed_to_voxel', signs2_p, overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds,1), num_secondaries);
            else
                voxels_coordinate_file = fullfile(root_path, 'FC', [study '_' 'coordinates.xlsx']);
                [signs2_p, overlaps_2_p, overlaps_2_w] = calculate_overlap('seed_to_voxel_CIFTI', voxels_coordinate_file, atlas);               
                [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('seed_to_voxel', signs2_p, overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds,1), num_secondaries);
            end
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_network')
        coordinate_file = fullfile(root_path, 'FC', [study '_' 'coordinates.xlsx']);
        coordinates = readmatrix(coordinate_file); % coordinates of the regions within each network 
        [networks_p, networks_w] = network_overlaps(coordinates, atlas);
        seeds = txt_pairs(3:end,j); % names of the seeds
        seeds = seeds(~cellfun('isempty',seeds));  
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
        signs = signs(~isnan(signs));
        for s = 1:size(seeds,1)
            seed_filename = fullfile(root_path, 'FC', [study '_' char(seeds(s)) '_seed.nii']);
            [signs_p, overlaps_mask_p, overlaps_mask_w] = calculate_overlap('other', seed_filename, atlas);
            % overlaps of the network that is functionally connected to the seed at hand
            overlaps_network_p = networks_p(cell2mat(raw_pairs(s+2,j+1)),:);
            overlaps_network_w = networks_w(cell2mat(raw_pairs(s+2,j+1)),:);
            [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(s), overlaps_mask_p, overlaps_network_p, overlaps_mask_w, overlaps_network_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds,1), num_secondaries);
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_coordinate')
        coordinate_file = fullfile(root_path, 'FC', [study '_' 'coordinates.xlsx']);
        coordinates = readmatrix(coordinate_file);
        seeds = txt_pairs(3:end,j); % names of the seeds
        seeds = seeds(~cellfun('isempty',seeds));  
        % since the same seed is used repeatedly in studies that measure
        % seed-to-coordinate connectivity, we calculate the overlaps of
        % that seed with the Yeo networks a single time.
        unique_seeds = unique(seeds);
        array_overlaps_mask_p = [];
        array_overlaps_mask_w = [];
        % assign an index to each seed
        unique_mask_inds = 1:length(unique_seeds); 
        maskMappings = containers.Map(unique_seeds, unique_mask_inds); 
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
        signs = signs(~isnan(signs));
        for u = 1:size(unique_seeds,1)
            seed_filename = fullfile(root_path, 'FC', [study '_' char(unique_seeds(u)) '_seed.nii']);
            [signs_p, overlaps_mask_p, overlaps_mask_w] = calculate_overlap('other', seed_filename, atlas);
            array_overlaps_mask_p = vertcat(array_overlaps_mask_p, overlaps_mask_p);
            array_overlaps_mask_w = vertcat(array_overlaps_mask_w, overlaps_mask_w);
        end
        for s = 1:size(seeds,1)
            seed_filename = fullfile(root_path, 'FC', [study '_' char(seeds(s)) '_seed.nii']);
            % retrieve the overlaps that were already calculated above
            % between the seed at hand and the Yeo networks 
            overlaps_mask_p = array_overlaps_mask_p(maskMappings(char(seeds(s))),:);
            overlaps_mask_w = array_overlaps_mask_w(maskMappings(char(seeds(s))),:);
            coordinate_ind = cell2mat(raw_pairs(s+2,j+1)); % index of the coordinate that is connected to the seed
            % x, y, and z values of the coordinate
            x = coordinates(coordinate_ind,1);
            y = coordinates(coordinate_ind,2);
            z = coordinates(coordinate_ind,3);
            sphere_file = fullfile(root_path, 'Spheres', ['sphere_mask_' num2str(x) '_' num2str(y) '_' num2str(z) '.nii']); % 10mm sphere around the coordinate
            [signs_p, overlaps_coordinates_p, overlaps_coordinates_w] = calculate_overlap('other', sphere_file, atlas);
            [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(s), overlaps_mask_p, overlaps_coordinates_p, overlaps_mask_w, overlaps_coordinates_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds,1), num_secondaries);
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_seed')
        % determine whether the seed index is a number or a string of
        % characters
        if all(cellfun('isempty', txt_pairs(3:end,j))) % seed index is a number
            seeds1 = cell2mat(raw_pairs(3:end,j)); % index of seed 1
            seeds1 = seeds1(~isnan(seeds1));   
            seeds2 = cell2mat(raw_pairs(3:end,j+1)); % index of seed 2
            seeds2 = seeds2(~isnan(seeds2)); 
            signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
            signs = signs(~isnan(signs));
            for s = 1:size(seeds1,1)
                seed1_filename = fullfile(root_path, 'FC', [study '_' num2str(seeds1(s)) '_seed.nii']);
                seed2_filename = fullfile(root_path, 'FC', [study '_' num2str(seeds2(s)) '_seed.nii']);
                [signs_p, overlaps_1_p, overlaps_1_w] = calculate_overlap('other', seed1_filename, atlas);
                [signs_p, overlaps_2_p, overlaps_2_w] = calculate_overlap('other', seed2_filename, atlas);
                [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(s), overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds1,1), num_secondaries);
            end
        else % seed index is a string
            seeds1 = txt_pairs(3:end,j); % index of seed 1 
            seeds1 = seeds1(~cellfun('isempty',seeds1));   
            seeds2 = txt_pairs(3:end,j+1); % index of seed 2 
            seeds2 = seeds2(~cellfun('isempty',seeds2)); 
            signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
            signs = signs(~isnan(signs));
            for s = 1:size(seeds1,1)
                seed1_filename = fullfile(root_path, 'FC', [study '_' seeds1{s} '_seed.nii']);
                seed2_filename = fullfile(root_path, 'FC', [study '_' seeds2{s} '_seed.nii']);
                [signs_p, overlaps_1_p, overlaps_1_w] = calculate_overlap('other', seed1_filename, atlas);
                [signs_p, overlaps_2_p, overlaps_2_w] = calculate_overlap('other', seed2_filename, atlas);
                [C_p, C_w, C_p_nosign, C_w_nosign] = scoring('other', signs(s), overlaps_1_p, overlaps_2_p, overlaps_1_w, overlaps_2_w, C_p, C_w, C_p_nosign, C_w_nosign, size(seeds1,1), num_secondaries);
            end
         end
    end
end

%% Determine which elements of the scoring matrix are significant/insignificant

% Create null distribution
n_iter = 1000;
num_rois = 62;
[all_Cp, all_Cw] = null_dist(n_iter,num_rois);
all_Cp_null = all_Cp;
all_Cw_null = all_Cw;

n_iter = size(all_Cp_null, 3); % number of scoring matrices in the surrogate study set
pvalsp = zeros(num_networks);
pvalsw = zeros(num_networks);

for r = 1:num_networks
    for c = 1:num_networks
        % Extract an n_iter vector of all the values for each entry of the scoring
        % matrices in the surrogate study set. This vector forms the null 
        % distribution for that entry.
        null_Cp_vals = squeeze(all_Cp_null(r,c,:));
        null_Cw_vals = squeeze(all_Cw_null(r,c,:));
        % Calculate the percentile of each cell in the scoring matrix of
        % the meta-analysis relative to the null distribution. Code adapted 
        % from here: 
        % https://uk.mathworks.com/matlabcentral/answers/182131-percentile-of-a-value-based-on-array-of-data
        nlessp = sum(null_Cp_vals < C_p_nosign(r,c));
        nequalp = sum(null_Cp_vals == C_p_nosign(r,c));
        pp = (nlessp + 0.5*nequalp) / length(null_Cp_vals);
        pvalsp(r,c) = pp;
        nlessw = sum(null_Cw_vals < C_w_nosign(r,c));
        nequalw = sum(null_Cw_vals == C_w_nosign(r,c));
        pw = (nlessw + 0.5*nequalw) / length(null_Cw_vals);
        pvalsw(r,c) = pw;       
    end
end    

% Eliminate elements of the scoring matrix if they are below the 5th
% percentile or above the 95th percentile of the null distribution
C_p_sig = C_p;
C_w_sig = C_w;
for r = 1:num_networks
    for c = 1:num_networks
        if pvalsp(r,c) > 0.05 & pvalsp(r,c) < 0.95
            C_p_sig(r,c) = NaN;
        end
        if pvalsw(r,c) > 0.05 & pvalsw(r,c) < 0.95
            C_w_sig(r,c) = NaN;
        end
    end
end

%% Normalise the scoring matrix by the sizes of the Yeo/Tian networks

% Note: I am not sure whether to normalise the scoring matrix because it
% may overinflate the functional connectivity between the limbic network,
% which is relatively very small, and other networks.
    
% number of voxels in each Yeo network
atlas_sizes = [];
for i = 1:num_networks
    numvox = size(find(atlas_img==i),1);
    atlas_sizes = [atlas_sizes; numvox];        
end

% Account for the fact that subcortical structures are lateralised to only
% one hemisphere (in our use of the Tian parcellation) and are a lot
% smaller than cortical structures
if strcmp(atlas, 'yeo+tian') | strcmp(atlas, 'yeo+tian+cerebellum')
    cort_avg_size = mean(atlas_sizes(1:7));
    subcort_avg_size = mean(atlas_sizes(8:15)) * 2;
    ratio_c2sc = cort_avg_size/subcort_avg_size;
    atlas_sizes(8:15) = atlas_sizes(8:15)/(10000/ratio_c2sc);
end
atlas_sizes(1:7) = atlas_sizes(1:7)/10000;

norm_C_p = zeros(num_networks);
norm_C_w = zeros(num_networks);
norm_C_p_nosig = zeros(num_networks);
norm_C_w_nosig = zeros(num_networks);

% Normalise the score of the connectivity between Yeo networks y1 and y2 by
% the average of the sizes of y1 and y2
for y1 = 1:num_networks
    for y2 = 1:num_networks
        norm = (atlas_sizes(y1) + atlas_sizes(y2))/2;
        if ~isnan(C_p_sig(y1,y2))
            norm_C_p(y1,y2) = C_p_sig(y1,y2)/norm;
        else
            norm_C_p(y1,y2) = C_p_sig(y1,y2);
        end
        norm_C_p_nosig(y1,y2) = C_p(y1,y2)/norm;
        if ~isnan(C_w_sig(y1,y2))
            norm_C_w(y1,y2) = C_w_sig(y1,y2)/norm; 
        else
            norm_C_w(y1,y2) = C_w_sig(y1,y2);
        end
        norm_C_w_nosig(y1,y2) = C_w(y1,y2)/norm;
    end
end


%% Make a heatmap of the results

 
f = figure;
h = heatmap(norm_C_p_nosig);
h.XDisplayLabels = network_labels;
h.YDisplayLabels = network_labels;
h.MissingDataColor = [0 0 0];

% Make the heatmap square. Code taken from here: https://uk.mathworks.com/matlabcentral/answers/578140-huge-white-space-around-the-plot-after-saving
% Temporarily change figure units 
originalUnits.fig = f.Units;  % save original units (probaly normalized)
originalUnits.ax = h.Units; 
f.Units = 'centimeters';  % any unit that will result in squares
h.Units = 'Normalize'; 
% Get number of rows & columns
sz = size(h.ColorData); 
% make axes square (not the table cells, just the axes)
h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
f.InnerPosition(3:4) = min(f.InnerPosition(3:4))*[1,1]; 
% Change figure position;
if sz(1)>sz(2)
    % make the figure size more narrow
    f.InnerPosition(3) = f.InnerPosition(3)*(sz(2)/sz(1)); 
else
    % make the figure size shorter
    f.InnerPosition(4) = f.InnerPosition(4)*(sz(1)/sz(2)); 
end
% return original figure units
f.Units = originalUnits.fig; 
h.Units = originalUnits.ax; 

% Set background color to white
set(f, 'color', 'white');

f = figure;
h = heatmap(norm_C_w_nosig);
h.XDisplayLabels = network_labels;
h.YDisplayLabels = network_labels;
h.MissingDataColor = [0 0 0];

% Make the heatmap square - massive faff!

% Temporarily change figure units 
originalUnits.fig = f.Units;  % save original units (probaly normalized)
originalUnits.ax = h.Units; 
f.Units = 'centimeters';  % any unit that will result in squares
h.Units = 'Normalize'; 
% Get number of rows & columns
sz = size(h.ColorData); 
% make axes square (not the table cells, just the axes)
h.Position(3:4) = min(h.Position(3:4))*[1,1]; 
f.InnerPosition(3:4) = min(f.InnerPosition(3:4))*[1,1]; 
% Change figure position;
if sz(1)>sz(2)
    % make the figure size more narrow
    f.InnerPosition(3) = f.InnerPosition(3)*(sz(2)/sz(1)); 
else
    % make the figure size shorter
    f.InnerPosition(4) = f.InnerPosition(4)*(sz(1)/sz(2)); 
end
% return original figure units
f.Units = originalUnits.fig; 
h.Units = originalUnits.ax; 

% Set background color to white
set(f, 'color', 'white');
