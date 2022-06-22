% This function calculates the mean expression level of each gene of 
% interest in each Yeo network. It is run in 'genes_to_MNI.m'.

% Inputs:
% - coords: MNI coordinates of the samples at which gene expression is
%   measured by a single probe (chosen in genes_to_MNI.m)
% - explevels: expression levels at those coordinates

% Outputs:
% - mean_yeo_explevels: mean expression level of the gene in each Yeo
%   network (7x1 vector)

function mean_yeo_explevels = calculate_overlap_genes(coords, explevels)

    root_path = '/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Literature/Meta-Analysis/Code Review/';

    yeo_explevels = {[], [], [], [], [], [], []};
     
    % find the Yeo network that a coordinate overlaps with
    yeo_atlas = spm_vol(fullfile(root_path, 'Yeo2011_7Networks_MNI152_FreeSurferConformed2mm_LiberalMask.nii.gz'));
    yeo_img = spm_read_vols(yeo_atlas);
    vox = [coords(:,1) coords(:,2) coords(:,3) ones(size(coords,1),1)]*(inv(yeo_atlas.mat))'; % convert MNI coordinates from mm to voxels
    vox = vox(:,1:3); % columns 1-3 of vox are the x,y,z coordinates in voxel-space
    for c = 1:size(coords,1)
        yeo_val = yeo_img(round(vox(c,1)),round(vox(c,2)),round(vox(c,3))); % the Yeo network associated with a coordinate
        if yeo_val ~= 0
            explevel = explevels(c); % the expression level of the gene at that coordinate
            yeo_explevels{yeo_val} = [yeo_explevels{yeo_val}; explevel]; % add that expression level to the vector of expression levels for the corresponding Yeo network
        end
    end
    
    % get the mean expression level of a gene across each Yeo network
    mean_yeo_explevels = [];
    for y = 1:7
        if isempty(yeo_explevels{y})
            mean_yeo_explevels = [mean_yeo_explevels; 0];
        else
            mean_yeo_explevels = [mean_yeo_explevels; mean(yeo_explevels{y})];
        end
    end

end