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

% NOTE: You must update this to reflect the root path to the data on your
% device.
root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';
cd(root_path)

%% Indicate whether a study defines ROIs based on coordinates or seeds

studyKeys = {'araujo_2012', 'barrett_2020_claustrum', 'bershad_2020', 'carhart-harris_2016_s2v', 'carhart-harris_2016_ica', 'kaelen_2016', 'lebedev_2015', 'luppi_2021_segregated', 'luppi_2021_time-averaged', 'mason_2021', 'madsen_2021', 'mcculloch_2021', 'muller_2017_s2s', 'muller_2017_s2v', 'muller_2018', 'palhano-fontes_2015_s2v', 'pasquini_2020_intra', 'pasquini_2020_inter', 'roseman_2014', 'sampedro_2017', 'smigielski_2019_ICA1', 'smigielski_2019_ICA2'};
values = {'seeds', 'coordinates', 'coordinates', 'seeds', 'seeds', 'seeds', 'coordinates', 'coordinates', 'coordinates', 'seeds', 'coordinates', 'coordinates', 'coordinates', 'seeds', 'seeds', 'coordinates', 'coordinates', 'seeds', 'seeds', 'coordinates', 'coordinates', 'seeds'};
    
studyMappings = containers.Map(studyKeys, values);

%% Define a sphere around the coordinates

% Adapted from: http://jpeelle.net/mri/misc/marsbar_roi.html

outDir = fullfile(root_path, 'Data', 'FC', 'Spheres');
sphereRadius = 10;

for i = 1:length(studyKeys) % each study
    if strcmp(studyMappings(studyKeys{i}), 'coordinates')
        % get the file with the relevant coordinates for that study
        study_coordinates_file = fullfile(root_path, 'Data', 'FC', [studyKeys{i} '_' studyMappings(studyKeys{i}) '.xlsx']);
        study_coordinates = readmatrix(study_coordinates_file);
        for j = 1:size(study_coordinates,1)
            row = study_coordinates(j,:);
            if ~isnan(row(1))
                x = row(1);
                y = row(2);
                z = row(3);
            end
            roiLabel = sprintf('%.4f_%.4f_%.4f', x, y, z); 
            sphereROI = maroi_sphere(struct('centre', [x y z], 'radius', sphereRadius)); % construct the sphere
            outName = fullfile(outDir, sprintf('sphere_mask_%s', roiLabel));
            save_as_image(sphereROI, [outName '.nii']);
        end
    end
end

%% Score the networks

% matrix of scores, where C_ij is the connectivity between Yeo network i
% and Yeo network j. C_p is the scoring matrix based on the partial
% paradigm; C_w is the scoring matrix based on the winners-take-all
% paradigm.
C_p = zeros(7);
C_w = zeros(7);

% load file indicating which pairs of seeds/coordinates are connected to each other
[num_pairs, txt_pairs, raw_pairs] = xlsread('Data/FC/FC_indices.xlsx');

for j = 1:3:size(raw_pairs, 2) % cycle through the studies
    study = char(raw_pairs(1,j)); % name of the study 
    if strcmp(raw_pairs(2,j), 'network_to_network')
        col = cell2mat(raw_pairs(3:end,j)); % indices of networks
        col = col(~isnan(col)); 
        coordinate_file = fullfile('Data', 'FC', [study '_' 'coordinates.xlsx']); % coordinates of each region in each network
        coordinates = readmatrix(coordinate_file);
        [networks_yeo_p, networks_yeo_w] = network_overlaps(coordinates); % note that there is a separate function for calcluating overlaps in networks, where the overlaps get averaged across the individual ROIs in the network
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC 
        signs = signs(~isnan(signs)); 
        % for each pair of functionally connected networks, get their
        % overlaps with the Yeo networks and then update the scoring matrix
        % accordingly
        for r = 1:size(col,1)
            overlaps_yeo1_p = networks_yeo_p(cell2mat(raw_pairs(r+2,j)),:);  
            overlaps_yeo2_p = networks_yeo_p(cell2mat(raw_pairs(r+2,j+1)),:);
            overlaps_yeo1_w = networks_yeo_w(cell2mat(raw_pairs(r+2,j)),:);
            overlaps_yeo2_w = networks_yeo_w(cell2mat(raw_pairs(r+2,j+1)),:);
            [C_p, C_w] = scoring('other', signs(r), overlaps_yeo1_p, overlaps_yeo2_p, overlaps_yeo1_w, overlaps_yeo2_w, C_p, C_w, size(col,1));
        end
    elseif strcmp(raw_pairs(2,j), 'coord_to_coord')
        coordinate_file = fullfile('Data', 'FC', [study '_' 'coordinates.xlsx']);
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
                sphere_file1 = fullfile('Data', 'FC', 'Spheres', ['sphere_mask_' sprintf('%d', x1) '_' sprintf('%d', y1) '_' sprintf('%d', z1) '.nii']);
            else
                sphere_file1 = fullfile('Data', 'FC', 'Spheres', ['sphere_mask_' sprintf('%.4f', x1) '_' sprintf('%.4f', y1) '_' sprintf('%.4f', z1) '.nii']);
            end
            % x, y, z values of coordinate 2
            x2 = coordinates(ind2,1);
            y2 = coordinates(ind2,2);
            z2 = coordinates(ind2,3);
            if floor(x2)==x2 & floor(y2)==y2 & floor(z2)==z2
                sphere_file2 = fullfile('Data', 'FC', 'Spheres', ['sphere_mask_' sprintf('%d', x2) '_' sprintf('%d', y2) '_' sprintf('%d', z2) '.nii']);                
            else
                sphere_file2 = fullfile('Data', 'FC', 'Spheres', ['sphere_mask_' sprintf('%.4f', x2) '_' sprintf('%.4f', y2) '_' sprintf('%.4f', z2) '.nii']);
            end
            [signs_p, overlaps_yeo1_p, overlaps_yeo1_w]  = calculate_overlap('other', sphere_file1);
            [signs_p, overlaps_yeo2_p, overlaps_yeo2_w] = calculate_overlap('other', sphere_file2);
            [C_p, C_w] = scoring('other', signs(r), overlaps_yeo1_p, overlaps_yeo2_p, overlaps_yeo1_w, overlaps_yeo2_w, C_p, C_w, size(col,1));
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_voxel')
        seeds = txt_pairs(3:end,j); % name of the seed
        seeds = seeds(~cellfun('isempty',seeds));        
        for s = 1:size(seeds,1)
            % I was unable to obtain the V1 seed used by Carhart-Harris et
            % al. (2016), so instead I assumed that the seed overlaps
            % completely with Yeo network 1, which is the visual network.
            if strcmp(study, 'carhart-harris_2016_s2v') & strcmp(char(seeds(s)), 'V1')
                overlap_yeo1_p = [1 0 0 0 0 0 0];
                overlap_yeo1_w = 1;
            else
                seed1_filename = fullfile('Data', 'FC', [study '_' char(seeds(s)) '_seed.nii']);
                [signs_p, overlap_yeo1_p, overlap_yeo1_w] = calculate_overlap('other', seed1_filename);
            end
            % Note that the spatial map of the functionally connected
            % voxels indicates the sign of the connectivity at each voxel,
            % hence FC_indices.xlsx does not specify the sign.
            voxels_filename = fullfile('Data', 'FC', [study '_' char(seeds(s)) '.nii']); % map of voxels that are functionally connected to the seed
            [signs2_p, overlap_yeo2_p, overlap_yeo2_w] = calculate_overlap('seed_to_voxel', voxels_filename);
            [C_p, C_w] = scoring('seed_to_voxel', signs2_p, overlap_yeo1_p, overlap_yeo2_p, overlap_yeo1_w, overlap_yeo2_w, C_p, C_w, size(seeds,1));
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_network')
        coordinate_file = fullfile('Data', 'FC', [study '_' 'coordinates.xlsx']);
        coordinates = readmatrix(coordinate_file); % coordinates of the regions within each network 
        [networks_yeo_p, networks_yeo_w] = network_overlaps(coordinates);
        seeds = txt_pairs(3:end,j); % names of the seeds
        seeds = seeds(~cellfun('isempty',seeds));  
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
        signs = signs(~isnan(signs));
        for s = 1:size(seeds,1)
            seed_filename = fullfile('Data', 'FC', [study '_' char(seeds(s)) '_seed.nii']);
            [signs_p, overlaps_mask_yeo_p, overlaps_mask_yeo_w] = calculate_overlap('other', seed_filename);
            % overlaps of the network that is functionally connected to the seed at hand
            overlaps_network_yeo_p = networks_yeo_p(cell2mat(raw_pairs(s+2,j+1)),:);
            overlaps_network_yeo_w = networks_yeo_w(cell2mat(raw_pairs(s+2,j+1)),:);
            [C_p, C_w] = scoring('other', signs(s), overlaps_mask_yeo_p, overlaps_network_yeo_p, overlaps_mask_yeo_w, overlaps_network_yeo_w, C_p, C_w, size(seeds,1));
        end
    elseif strcmp(raw_pairs(2,j), 'seed_to_coordinate')
        coordinate_file = fullfile('Data', 'FC', [study '_' 'coordinates.xlsx']);
        coordinates = readmatrix(coordinate_file);
        seeds = txt_pairs(3:end,j); % names of the seeds
        seeds = seeds(~cellfun('isempty',seeds));  
        % since the same seed is used repeatedly in studies that measure
        % seed-to-coordinate connectivity, we calculate the overlaps of
        % that seed with the Yeo networks a single time.
        unique_seeds = unique(seeds);
        array_overlaps_mask_yeo_p = [];
        array_overlaps_mask_yeo_w = [];
        % assign an index to each seed
        unique_mask_inds = 1:length(unique_seeds); 
        maskMappings = containers.Map(unique_seeds, unique_mask_inds); 
        signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
        signs = signs(~isnan(signs));
        for u = 1:size(unique_seeds,1)
            seed_filename = fullfile('Data', 'FC', [study '_' char(unique_seeds(u)) '_seed.nii']);
            [signs_p, overlaps_mask_yeo_p, overlaps_mask_yeo_w] = calculate_overlap('other', seed_filename);
            array_overlaps_mask_yeo_p = vertcat(array_overlaps_mask_yeo_p, overlaps_mask_yeo_p);
            array_overlaps_mask_yeo_w = vertcat(array_overlaps_mask_yeo_w, overlaps_mask_yeo_w);
        end
        for s = 1:size(seeds,1)
            seed_filename = fullfile('Data', 'FC', [study '_' char(seeds(s)) '_seed.nii']);
            % retrieve the overlaps that were already calculated above
            % between the seed at hand and the Yeo networks 
            overlaps_mask_yeo_p = array_overlaps_mask_yeo_p(maskMappings(char(seeds(s))),:);
            overlaps_mask_yeo_w = array_overlaps_mask_yeo_w(maskMappings(char(seeds(s))),:);
            coordinate_ind = cell2mat(raw_pairs(s+2,j+1)); % index of the coordinate that is connected to the seed
            % x, y, and z values of the coordinate
            x = coordinates(coordinate_ind,1);
            y = coordinates(coordinate_ind,2);
            z = coordinates(coordinate_ind,3);
            sphere_file = fullfile('Data', 'FC', 'Spheres', ['sphere_mask_' num2str(x) '_' num2str(y) '_' num2str(z) '.nii']); % 10mm sphere around the coordinate
            [signs_p, overlaps_coordinates_yeo_p, overlaps_coordinates_yeo_w] = calculate_overlap('other', sphere_file);
            [C_p, C_w] = scoring('other', signs(s), overlaps_mask_yeo_p, overlaps_coordinates_yeo_p, overlaps_mask_yeo_w, overlaps_coordinates_yeo_w, C_p, C_w, size(seeds,1));
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
                seed1_filename = fullfile('Data', 'FC', [study '_' num2str(seeds1(s)) '_seed.nii']);
                seed2_filename = fullfile('Data', 'FC', [study '_' num2str(seeds2(s)) '_seed.nii']);
                [signs_p, overlaps_yeo1_p, overlaps_yeo1_w] = calculate_overlap('other', seed1_filename);
                [signs_p, overlaps_yeo2_p, overlaps_yeo2_w] = calculate_overlap('other', seed2_filename);
                [C_p, C_w] = scoring('other', signs(s), overlaps_yeo1_p, overlaps_yeo2_p, overlaps_yeo1_w, overlaps_yeo2_w, C_p, C_w, size(seeds1,1));
            end
        else % seed index is a string
            seeds1 = txt_pairs(3:end,j); % index of seed 1 
            seeds1 = seeds1(~cellfun('isempty',seeds1));   
            seeds2 = txt_pairs(3:end,j+1); % index of seed 2 
            seeds2 = seeds2(~cellfun('isempty',seeds2)); 
            signs = cell2mat(raw_pairs(3:end,j+2)); % sign of the change in FC
            signs = signs(~isnan(signs));
            for s = 1:size(seeds1,1)
                seed1_filename = fullfile('Data', 'FC', [study '_' seeds1{s} '_seed.nii']);
                seed2_filename = fullfile('Data', 'FC', [study '_' seeds2{s} '_seed.nii']);
                [signs_p, overlaps_yeo1_p, overlaps_yeo1_w] = calculate_overlap('other', seed1_filename);
                [signs_p, overlaps_yeo2_p, overlaps_yeo2_w] = calculate_overlap('other', seed2_filename);
                [C_p, C_w] = scoring('other', signs(s), overlaps_yeo1_p, overlaps_yeo2_p, overlaps_yeo1_w, overlaps_yeo2_w, C_p, C_w, size(seeds1,1));
            end
         end
    end
end

%% Determine which elements of the scoring matrix are significant/insignificant

% 'null_surrogate_set.m' calculates a surrogate study set, in which each
% element of the set is a scoring matrix computed on sober resting-state
% functional connectivity data. all_Cp contains all 500 of the partial
% scoring matrices in the surrogate study sets, and all_Cw contains all
% 500 of the WTA matrices. They are each 7x7x500 tensors.
load(fullfile('Data', 'FC', 'Null', 'all_Cp_null.mat'))
load(fullfile('Data', 'FC', 'Null', 'all_Cw_null.mat'))

n_iter = 500; % number of scoring matrices in the surrogate study set
pvalsp = zeros(7);
pvalsw = zeros(7);

for r = 1:7
    for c = 1:7
        % Produce a histogram of values for each cell of the scoring
        % matrices in the surrogate study set (Cp_vals and Cw_vals are 
        % 1x500 vectors). This histogram is the null distribution for that 
        % cell.
        null_Cp_vals = [];
        null_Cw_vals = []; 
        for m = 1:n_iter
            null_Cp_vals = [null_Cp_vals all_Cp(r,c,m)];
            null_Cw_vals = [null_Cw_vals all_Cw(r,c,m)];
        end
        % Calculate the percentile of each cell in the scoring matrix of
        % the meta-analysis relative to the null distribution. Code adapted 
        % from here: 
        % https://uk.mathworks.com/matlabcentral/answers/182131-percentile-of-a-value-based-on-array-of-data
        nlessp = sum(null_Cp_vals < C_p(r,c));
        nequalp = sum(null_Cp_vals == C_p(r,c));
        pp = (nlessp + 0.5*nequalp) / length(null_Cp_vals);
        pvalsp(r,c) = pp;
        nlessw = sum(null_Cw_vals < C_w(r,c));
        nequalw = sum(null_Cw_vals == C_w(r,c));
        pw = (nlessw + 0.5*nequalw) / length(null_Cw_vals);
        pvalsw(r,c) = pw;       
    end
end    

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

% eliminate elements of the scoring matrix if their corresponding p-value
% is greater than 0.05
C_w_sig = C_w;
C_p_sig = C_p;
for r = 1:7
    for c = 1:7
        if pvalsp(r,c) >= 0.05
            C_p_sig(r,c) = NaN;
        end
        if pvalsw(r,c) >= 0.05
            C_w_sig(r,c) = NaN;
        end
    end
end

%% Normalise the scoring matrix by the sizes of the Yeo networks

% Note: I am not sure whether to normalise the scoring matrix because it
% may overinflate the functional connectivity between the limbic network,
% which is relatively very small, and other networks.

yeo = spm_vol(fullfile(root_path, 'Yeo2011_7Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz'));
yeo_img = spm_read_vols(yeo);
    
% number of voxels in each Yeo network
yeo_sizes = [];
for i = 1:7
    numvox = size(find(yeo_img==i),1);
    yeo_sizes = [yeo_sizes numvox];
end
yeo_sizes = yeo_sizes/10000; % scale down the sizes since we will be dividing the values of C by the sizes

norm_C_p = zeros(7,7);
norm_C_w = zeros(7,7);

% divide the score of the connectivity between Yeo networks y1 and y2 by
% the average of the sizes of y1 and y2
for y1 = 1:size(C_p,1)
    for y2 = 1:size(C_p,2)
        norm = (yeo_sizes(y1)+yeo_sizes(y2))/2;
        if ~isnan(C_p_sig(y1,y2))
            norm_C_p(y1,y2) = C_p_sig(y1,y2)/norm;
        else
            norm_C_p(y1,y2) = C_p_sig(y1,y2);
        end
        if ~isnan(C_w_sig(y1,y2))
            norm_C_w(y1,y2) = C_w_sig(y1,y2)/norm; 
        else
            norm_C_w(y1,y2) = C_w_sig(y1,y2);
        end
    end
end


%% Make a heatmap of the results

figure;
h1 = heatmap(norm_C_p);
h1.XDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
h1.YDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
h1.MissingDataColor = [0 0 0];

figure;
h2 = heatmap(norm_C_w);
h2.XDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
h2.YDisplayLabels = {'Visual', 'SMN', 'VAN', 'DAN', 'Limbic', 'FP', 'DMN'};
h2.MissingDataColor = [0 0 0];

%% Visualise the results using OSL

figure;
spatial_basis_file = fullfile(root_path, 'Yeo2011_7Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz');
p = parcellation(spatial_basis_file);
p.labels = [{'Visual'}, {'SMN'}, {'VAN'}, {'DAN'}, {'Limbic'}, {'FP'}, {'DMN'}];
[h_patch,h_scatter] = p.plot_network(norm_C_w,0,NaN,true);
osl_spinning_brain('FC_brain_WTA_150322.mp4');

figure;
p = parcellation(spatial_basis_file);
p.labels = [{'Visual'}, {'SMN'}, {'VAN'}, {'DAN'}, {'Limbic'}, {'FP'}, {'DMN'}];
[h_patch,h_scatter] = p.plot_network(norm_C_p,0,NaN,true);
osl_spinning_brain('FC_brain_partial_150322.mp4');


